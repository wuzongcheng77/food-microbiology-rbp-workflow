"""Independent YRBP intake prototype. No original pipeline writes or accession bonuses.
FASTA -> sequence checks -> optional prepared-structure docking -> contacts -> report.
Missing structures yield explicit work requests, never invented docking scores.
"""
import argparse,csv,hashlib,json,re,subprocess,time
from pathlib import Path
AA=set('ACDEFGHIKLMNPQRSTVWY')
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def fasta(path):
    records=[]; name=None;seq=[]
    for line in Path(path).read_text(encoding='utf-8-sig').splitlines():
        if line.startswith('>'):
            if name is not None:records.append((name,''.join(seq).upper()))
            name=line[1:].split()[0];seq=[]
        elif line.strip():
            if name is None:raise ValueError('Sequence before FASTA header')
            seq.append(line.strip())
    if name is not None:records.append((name,''.join(seq).upper()))
    if not records:raise ValueError('Empty FASTA')
    if len({n for n,s in records})!=len(records):raise ValueError('Duplicate candidate IDs')
    for n,s in records:
        if not s or set(s)-AA:raise ValueError('Invalid/ambiguous amino acids: '+n)
    return records
def atoms(path,first_model=False):
    result=[]
    for l in Path(path).read_text().splitlines():
        if first_model and l.startswith('ENDMDL'):break
        if l.startswith(('ATOM  ','HETATM')):
            typ=l.split()[-1]
            if typ in ('H','HD','HS'):continue
            result.append((l[17:20].strip(),l[21:22].strip(),l[22:26].strip(),tuple(float(l[a:b]) for a,b in [(30,38),(38,46),(46,54)])))
    return result
def contacts(protein,pose):
    residues={};pairs=0
    for rn,ch,ri,p in atoms(protein):
        for _,_,_,q in atoms(pose,True):
            d=sum((a-b)**2 for a,b in zip(p,q))**.5
            if d<=4.0:
                pairs+=1;k=(rn,ch,ri);residues[k]=min(d,residues.get(k,999))
    return [{'residue_name':k[0],'chain':k[1],'residue_number':k[2],'minimum_distance_A':round(v,3)} for k,v in sorted(residues.items())],pairs
def dumpcsv(path,rows,fields=None):
    with Path(path).open('w',encoding='utf-8',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields or list(rows[0]));w.writeheader();w.writerows(rows)
def _fingerprint(args,protein,ligand,spec):
    # Exact-execution identity: same engine, receptor, ligand, grid, seed and exhaustiveness.
    payload={'command':args,'protein_sha256':sha(protein),'ligand_sha256':sha(ligand),
             'grid':{k:spec.get(k) for k in ['center_x','center_y','center_z','size_x','size_y','size_z','exhaustiveness']},
             'seed':'20260929','engine_sha256':sha(args[0]) if Path(args[0]).is_file() else 'unavailable'}
    return hashlib.sha256(json.dumps(payload,sort_keys=True).encode()).hexdigest()

def run(input_path,manifest_path,out,vina=None,resume=False):
    # resume=False (default): output dir must be new (never overwrite history).
    # resume=True: continue an existing run dir; candidates whose completed.json
    # fingerprint matches are reused as cache_reuse, everything else (re)executes.
    out=Path(out).resolve()
    if resume:
        out.mkdir(parents=True,exist_ok=True)
    else:
        out.mkdir(parents=True,exist_ok=False)
    manifest=json.loads(Path(manifest_path).read_text()) if manifest_path else {}
    base=Path(manifest_path).resolve().parent if manifest_path else Path.cwd()
    rows=[];evidence=[];seen={}
    for n,s in fasta(input_path):
        h=hashlib.sha256(s.encode()).hexdigest();dup=seen.get(h,'');seen[h]=n
        row={'candidate_id':n,'sequence_sha256':h,'length':len(s),'duplicate_sequence_of':dup,'hydrophobic_fraction':sum(s.count(a) for a in 'AVILMFWY')/len(s),'cysteine_fraction':s.count('C')/len(s),'literature_status':'not_searched','docking_status':'structure_required','execution_mode':'not_run_no_structure','affinity_kcal_mol':'','contact_residues':'','contact_atom_pairs':''}
        spec=manifest.get(n)
        if spec:
            if spec['sequence_sha256']!=h:raise ValueError('Manifest sequence mismatch: '+n)
            protein=(base/spec['protein_pdbqt']).resolve();ligand=(base/spec['ligand_pdbqt']).resolve()
            for p,k in [(protein,'protein_sha256'),(ligand,'ligand_sha256')]:
                if sha(p)!=spec[k]:raise ValueError('Structure checksum mismatch: '+n)
            row['docking_status']='prepared_structure_available'
            row['execution_mode']='not_run_no_engine' if not vina else 'pending'
            if vina:
                folder=out/h[:16];folder.mkdir(exist_ok=True);pose=folder/'pose.pdbqt'
                args=[str(vina),'--receptor',str(protein),'--ligand',str(ligand),'--out',str(pose),'--seed','20260929','--cpu','2','--exhaustiveness',str(spec.get('exhaustiveness',8)),'--num_modes','3']
                for key in ['center_x','center_y','center_z','size_x','size_y','size_z']:args.extend(['--'+key,str(spec[key])])
                fp=_fingerprint(args,protein,ligand,spec)
                done=folder/'completed.json'
                reused=False
                if resume and done.exists():
                    try:
                        rec=json.loads(done.read_text(encoding='utf-8'))
                        if (rec.get('fingerprint')==fp and pose.is_file() and sha(pose)==rec.get('pose_sha256')
                            and (folder/'contacts.csv').is_file() and sha(folder/'contacts.csv')==rec.get('contacts_sha256')):
                            row.update(docking_status='completed_cache_reuse',execution_mode='cache_reuse',
                                       affinity_kcal_mol=rec['affinity_kcal_mol'],contact_residues=rec['contact_residues'],
                                       contact_atom_pairs=rec['contact_atom_pairs'])
                            evidence.append({'candidate_id':n,'sequence_sha256':h,'stage':'cache_reuse','fingerprint':fp,'pose_sha256':rec['pose_sha256'],'source':'prior_completed_run'})
                            reused=True
                    except Exception:
                        reused=False
                if not reused:
                    start=time.monotonic()
                    try:
                        proc=subprocess.run(args,capture_output=True,timeout=600)
                        (folder/'vina.stdout.bin').write_bytes(proc.stdout)
                        (folder/'vina.stderr.bin').write_bytes(proc.stderr)
                        (folder/'vina.log').write_text((proc.stdout+b'\n'+proc.stderr).decode('utf-8',errors='replace'),encoding='utf-8')
                        if proc.returncode or not pose.exists():raise RuntimeError('Vina failed; inspect vina.log')
                        m=re.search(r'REMARK VINA RESULT:\s*([-\d.]+)',pose.read_text())
                        if not m:raise RuntimeError('No Vina result')
                        ct,pairs=contacts(protein,pose)
                        dumpcsv(folder/'contacts.csv',ct,['residue_name','chain','residue_number','minimum_distance_A'])
                        affinity=float(m.group(1))
                        row.update(docking_status='completed_new_run',execution_mode='real_docking',affinity_kcal_mol=affinity,contact_residues=len(ct),contact_atom_pairs=pairs)
                        rec={'contacts_sha256':sha(folder/'contacts.csv'),'fingerprint':fp,'affinity_kcal_mol':affinity,'contact_residues':len(ct),'contact_atom_pairs':pairs,'pose_sha256':sha(pose),'pose_sha256_path':str(pose),'seconds':time.monotonic()-start}
                        done.write_text(json.dumps(rec,indent=2),encoding='utf-8')
                        evidence.append({'candidate_id':n,'sequence_sha256':h,'stage':'real_docking','protein_sha256':sha(protein),'ligand_sha256':sha(ligand),'pose_sha256':sha(pose),'seconds':time.monotonic()-start,'command':args,'structure_origin':spec.get('origin','user_supplied'),'scope':'Prepared structure supplied; folding not performed in this intake run; contacts are geometric proximity, not validated bonds.'})
                    except (subprocess.TimeoutExpired,RuntimeError) as e:
                        row['docking_status']='failed';row['execution_mode']='failed'
                        evidence.append({'candidate_id':n,'sequence_sha256':h,'stage':'failed','error':str(e)})
        rows.append(row)
        # Incremental flush: a crash after candidate N never loses candidates 1..N.
        dumpcsv(out/'candidate_results.csv',rows)
    modes={}
    for r in rows:modes[r['execution_mode']]=modes.get(r['execution_mode'],0)+1
    (out/'run_manifest.json').write_text(json.dumps({'input_sha256':sha(input_path),'candidate_count':len(rows),'evidence':evidence,'execution_modes':modes,'scope':'Execution demonstration, not an external predictive benchmark. Vina affinities across different protein structures are not calibrated probabilities or a sole candidate-ranking criterion.'},indent=2))
    (out/'structure_requests.json').write_text(json.dumps([{'candidate_id':r['candidate_id'],'sequence_sha256':r['sequence_sha256'],'needed':'sequence-matched protein structure, prepared ligand and justified grid'} for r in rows if r['docking_status']=='structure_required'],indent=2))
    return rows
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--fasta',required=True);p.add_argument('--manifest');p.add_argument('--out',required=True);p.add_argument('--vina');a=p.parse_args()
    print(json.dumps(run(a.fasta,a.manifest,a.out,a.vina),indent=2))
