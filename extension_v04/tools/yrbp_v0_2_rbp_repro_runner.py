"""Versioned RBP reproduction orchestrator (v0.2) wired into the original YRBP tree.

This is the called entry that connects the versioned intake/evidence modules to the
frozen historical RBP execution path. It is deliberately narrow and honest:

  * Protein-FASTA intake (never genome ORF calling in this entry).
  * Read-only lookup of the frozen v5.4.8z2 827-candidate registry for membership,
    duplicate/development-set tracking and historical component scores.
  * Optional real docking on caller-supplied prepared structures (Vina), with
    cache-reuse and raw-byte logging.
  * Optional evidence-bounded ranking under a fixed, pre-declared policy.

It NEVER writes to the frozen rank table, the release pointer, other branches or
frozen model artifacts. Stages it does not perform (ORF calling, folding, full
E6/E7 recomputation) are recorded as not_done, never silently treated as zero.
"""
import argparse, csv, hashlib, json
from pathlib import Path

from yrbp_candidate_intake_support_v0_2 import run, fasta
from yrbp_candidate_intake_v0_2 import audit
from yrbp_sequence_registry_v0_2 import load_registry, membership
from yrbp_input_contract_v0_3 import sequence_features, receptor_contract
from yrbp_evidence_bridge_v0_4 import e6_annotations, structure_sequence_audit, evidence_request

VERSION = '0.4.0'
REPO = Path(__file__).resolve().parent.parent
FROZEN_RANK_CSV = REPO / 'data/derived/YRBP_branch_runs/v5_4_8z2_G_pre1_target_engineering_composite_rank/YRBP_v5_4_8z2_G_pre1_target_engineering_composite_full_rank_table.csv'
PROTECTED = {
    'frozen_rank_script': REPO / 'tools/yrbp_v5_4_8z2_g_pre1_target_engineering_composite_rank.py',
    'release_pointer': REPO / 'config/yamps_current_release_v1.json',
    'frozen_rank_table': FROZEN_RANK_CSV,
}


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def load_frozen_registry(rank_path=FROZEN_RANK_CSV):
    with Path(rank_path).open(encoding='utf-8') as f:
        rows = list(csv.DictReader(f))
    by_acc = {r['source_accession']: r for r in rows}
    return rows, by_acc


def protected_hashes():
    return {name: sha(p) for name, p in PROTECTED.items() if p.exists()}


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--fasta', required=True)
    p.add_argument('--manifest', help='prepared-structure manifest (json); omit for sequence-only intake')
    p.add_argument('--vina', help='path to vina executable; omit to skip docking')
    p.add_argument('--out', required=True, help='fresh run directory')
    p.add_argument('--policy')
    p.add_argument('--evidence')
    p.add_argument('--resume', action='store_true')
    p.add_argument('--receptor-evidence', help='Independent accession-keyed JSON evidence, usable without structures')
    p.add_argument('--rank-registry', help='Explicit frozen rank CSV for portable reproduction')
    p.add_argument('--sequence-registry', help='Corresponding B0 CSV for portable reproduction')
    args = p.parse_args()

    if args.evidence and not args.policy:p.error('--evidence requires --policy')
    if bool(args.rank_registry) != bool(args.sequence_registry):
        p.error('--rank-registry and --sequence-registry must be provided together')
    registry_paths = {'rank_path': args.rank_registry, 'b0_path': args.sequence_registry} if args.rank_registry else {}
    protected_before = protected_hashes()
    supplied_registry_hashes = {k: sha(v) for k,v in registry_paths.items()}
    sequence_registry = load_registry(**registry_paths)
    out = Path(args.out).resolve()
    # run() owns the run-dir lifecycle (refuse existing unless --resume).

    # ---- Stage 1: intake + optional docking (real execution / cache / not_done) ----
    records = fasta(args.fasta)
    original_manifest = json.loads(Path(args.manifest).read_text(encoding='utf-8-sig')) if args.manifest else {}
    declarations = json.loads(Path(args.receptor_evidence).read_text(encoding='utf-8-sig')) if args.receptor_evidence else {}
    unknown_ids = set(declarations) - {n for n,seq in records}
    if unknown_ids: raise ValueError('Receptor evidence has IDs absent from FASTA: ' + ','.join(sorted(unknown_ids)))
    combined = {}
    for n,seq in records:
        spec = dict(original_manifest.get(n, {}))
        if n in declarations:
            if 'receptor_evidence' in spec and spec['receptor_evidence'] != declarations[n]:
                raise ValueError('Conflicting receptor declarations: ' + n)
            spec['receptor_evidence'] = declarations[n]
        combined[n] = spec
    gates = {n: receptor_contract(sequence_features(seq)['sequence_sha256'], combined[n]) for n, seq in records}
    structure_base = Path(args.manifest).resolve().parent if args.manifest else Path.cwd()
    structures = {n: structure_sequence_audit(seq, original_manifest.get(n), structure_base) for n,seq in records}
    for n,gate in gates.items():
        if args.vina and not structures[n]['automatic_acceptance']:
            gate['status'] = 'requires_review'
            gate['reasons'].append(structures[n]['status'])
    # Gate all requested docking as one batch: no partial hidden execution.
    if args.vina and any(g['status'] != 'declaration_consistent' for g in gates.values()):
        if out.exists() and not args.resume:
            raise FileExistsError(out)
        out.mkdir(parents=True, exist_ok=True)
        (out / 'receptor_preflight_blocked.json').write_text(json.dumps(gates, indent=2), encoding='utf-8')
        raise ValueError('Docking not started: receptor declarations need review; see receptor_preflight_blocked.json')
    intake_rows = run(args.fasta, args.manifest, out, args.vina, resume=args.resume)
    annotations = [e6_annotations(n,seq) for n,seq in records]
    requests = [evidence_request(n, gates[n], structures[n]) for n,seq in records]
    for filename,rows in [('fresh_e6_annotations.csv', annotations), ('evidence_requests.csv', requests)]:
        with (out / filename).open('w',encoding='utf-8',newline='') as f:
            w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    (out / 'structure_sequence_audit.json').write_text(json.dumps(structures,indent=2),encoding='utf-8')
    fresh = [{'candidate_id': n, **sequence_features(seq)} for n, seq in records]
    with (out / 'fresh_sequence_features.csv').open('w', encoding='utf-8', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(fresh[0])); w.writeheader(); w.writerows(fresh)
    (out / 'receptor_contract_audit.json').write_text(json.dumps(gates, indent=2), encoding='utf-8')
    for row in intake_rows:
        row['receptor_contract_status'] = gates[row['candidate_id']]['status']
    feature_manifest = {'schema': 'raw_composition_v1', 'records': len(fresh),
                        'unique_sequences': len({r['sequence_sha256'] for r in fresh}),
                        'input_sha256': sha(args.fasta), 'module_sha256': sha(Path(__file__).with_name('yrbp_input_contract_v0_3.py')),
                        'scope': 'Fresh descriptive features only. Not substitutes for normalized E6/E7 predictors; not used to retune frozen ranking.'}
    (out / 'fresh_feature_manifest.json').write_text(json.dumps(feature_manifest, indent=2), encoding='utf-8')

    # ---- Stage 2: read-only historical registry membership ----
    frozen_rows, frozen_by_acc = load_frozen_registry(args.rank_registry or FROZEN_RANK_CSV)
    registry_count = len(frozen_rows)

    source_register = []
    for r in intake_rows:
        acc = r['candidate_id']
        identity = membership(acc, r['sequence_sha256'], sequence_registry)
        r.update(identity)
        hist = frozen_by_acc.get(acc) if identity['historical_identity_status']=='accession_and_sequence_match' else None
        r['in_historical_827'] = identity['sequence_in_historical_827']
        r['historical_rank'] = hist['G_pre1_target_engineering_composite_rank'] if hist else ''
        r['historical_score'] = hist['G_pre1_target_engineering_composite_score'] if hist else ''
        r['development_set_participation'] = 'same_sequence_as_frozen_827' if identity['sequence_in_historical_827'] else 'exact_sequence_absent_independence_unverified'
        source_register.append({
            'accession': acc,
            **identity,
            'sequence_sha256': r['sequence_sha256'],
            'in_historical_827': r['in_historical_827'],
            'historical_rank': r['historical_rank'],
            'duplicate_sequence_of': r['duplicate_sequence_of'],
            'docking_status': r['docking_status'],
            'execution_mode': r['execution_mode'],
            'publication_source': 'frozen v5.4.8z2 candidate set (read-only); UniProt accession used as identifier; retrieval version/date not re-verified in this run',
        })

    # rewrite candidate_results.csv with membership columns
    fields = list(intake_rows[0].keys())
    with (out / 'candidate_results.csv').open('w', encoding='utf-8', newline='') as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(intake_rows)
    with (out / 'source_register.csv').open('w', encoding='utf-8', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(source_register[0].keys()))
        w.writeheader()
        w.writerows(source_register)

    # ---- Stage 3: evidence audit (computed) or not_done ----
    audited = None
    evidence_stage = 'not_done_no_policy'
    if args.policy and args.evidence:
        policy = json.loads(Path(args.policy).read_text(encoding='utf-8-sig'))
        evidence = json.loads(Path(args.evidence).read_text(encoding='utf-8-sig'))
        audited = audit(intake_rows, policy, evidence)
        (out / 'evidence_audit.json').write_text(json.dumps(audited, ensure_ascii=False, indent=2), encoding='utf-8')
        evidence_stage = 'computed_fixed_policy'

    # ---- Stage ledger: real / cache / mock / not_done ----
    modes = {}
    for r in intake_rows:
        modes[r['execution_mode']] = modes.get(r['execution_mode'], 0) + 1
    ledger = [
        {'stage': 'E6_raw_annotation_export', 'execution_mode': 'real_execution', 'count': len(annotations), 'note': 'Input-derived length/C/K connected to E6 annotation schema; scoring_ready=false; missing predictors not filled'},
        {'stage': 'structure_sequence_correspondence', 'execution_mode': 'real_execution', 'note': 'First-model ATOM residues; partial or mismatching chains require review'},
        {'stage': 'fresh_sequence_descriptors', 'execution_mode': 'real_execution', 'count': len(fresh), 'note': 'AA20 composition and entropy; not a binding predictor'},
        {'stage': 'receptor_declaration_check', 'execution_mode': 'real_execution', 'note': 'Missing or inconsistent declarations block requested docking; not biological validation'},
        {'stage': 'input_intake_protein_fasta', 'execution_mode': 'real_execution', 'count': len(intake_rows), 'note': 'Protein FASTA only; genome/ORF calling not wired'},
        {'stage': 'historical_827_registry_lookup', 'execution_mode': 'read_only_cache', 'count': registry_count, 'note': 'Frozen v5.4.8z2 table read; not modified'},
        {'stage': 'structure_resolution', 'execution_mode': 'real_execution' if args.manifest else 'not_done_no_manifest', 'count': sum(1 for r in intake_rows if r['docking_status'] != 'structure_required'), 'note': 'Prepared structures validated by declared sha256'},
        {'stage': 'docking_vina', 'execution_mode': 'mixed' if modes else 'not_done', 'detail': modes, 'count': modes.get('real_docking', 0), 'note': 'real_docking=%d cache_reuse=%d failed=%d no_structure=%d' % (modes.get('real_docking', 0), modes.get('cache_reuse', 0), modes.get('failed', 0), modes.get('not_run_no_structure', 0))},
        {'stage': 'contact_qc_4A_proximal', 'execution_mode': 'real_execution' if modes.get('real_docking') or modes.get('cache_reuse') else 'not_done', 'note': 'Geometric proximity in first pose; not validated bonds'},
        {'stage': 'evidence_audit_rank_bounds', 'execution_mode': evidence_stage, 'note': 'Unknown kept as interval; never coerced to 0; fixed weights declared in policy'},
        {'stage': 'genome_orf_calling', 'execution_mode': 'not_done_not_implemented', 'note': 'No end-to-end genome claim'},
        {'stage': 'automatic_structure_folding', 'execution_mode': 'not_done_not_implemented', 'note': 'Only caller-supplied prepared structures docked'},
        {'stage': 'full_E6_E7_rerank', 'execution_mode': 'not_done_not_implemented', 'note': 'Historical component scores used read-only; no re-weighting to preserve Gp17 rank'},
    ]
    (out / 'stage_ledger.json').write_text(json.dumps(ledger, ensure_ascii=False, indent=2), encoding='utf-8')

    # ---- Provenance / protection record ----
    protected_after = protected_hashes()
    if supplied_registry_hashes != {k:sha(v) for k,v in registry_paths.items()}:
        raise RuntimeError('Supplied registry changed during execution')
    provenance = {
        'version': VERSION,
        'tool_sha256': sha(__file__),
        'inputs': {'fasta': str(Path(args.fasta).resolve()), 'fasta_sha256': sha(args.fasta),
                   'manifest': args.manifest, 'vina': args.vina,
                   'receptor_evidence_sha256': sha(args.receptor_evidence) if args.receptor_evidence else None,
                   'supplied_registry_hashes': supplied_registry_hashes,
                   'policy_sha256': sha(args.policy) if args.policy else None,
                   'evidence_sha256': sha(args.evidence) if args.evidence else None,
                   'manifest_sha256': sha(args.manifest) if args.manifest else None},
        'protected_files_before': protected_before,
        'protected_files_after': protected_after,
        'protected_files_unchanged': protected_before == protected_after,
        'frozen_registry_rows': registry_count,
        'candidates_in_827': sum(1 for r in intake_rows if r['in_historical_827']),
        'candidates_outside_827': sum(1 for r in intake_rows if not r['in_historical_827']),
        'docking_modes': modes,
        'evidence_stage': evidence_stage,
        'scope': 'Computational execution and evidence bookkeeping. No accuracy or biological-hit claim without independent labels. Rank bounds are conditional mathematical intervals, not confidence intervals.',
    }
    (out / 'repro_manifest.json').write_text(json.dumps(provenance, ensure_ascii=False, indent=2), encoding='utf-8')
    if protected_before != protected_after:raise RuntimeError('Protected file changed during execution; inspect repro_manifest.json')
    print(json.dumps({'out': str(out), 'candidates': len(intake_rows), 'in_827': provenance['candidates_in_827'], 'outside_827': provenance['candidates_outside_827'], 'docking_modes': modes, 'version': VERSION}))


if __name__ == '__main__':
    main()
