import csv
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys
import numpy as np
import pytest
from PIL import Image

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from scoring import load_numeric,components,scores,rounded,rank,frozen_ties,WEIGHTS,analysis_ties
from feedback import reproduce_feedback,register_verified

def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):
    with p.open(encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f))

def test_input_hashes():
    assert sha(ROOT/'inputs/frozen_top50.csv')=='8d7e43b5fea103a3cbc38f62bda845d8ebe3d3e57ade6c591090f41aad3b45c2'
    assert sha(ROOT/'inputs/frozen_full.csv')=='85c3c3b9b8e5db3827fbdc7945426666226beb46b8461dd8605e187e50c94886'
    for item in json.loads((ROOT/'manifests/input_provenance.json').read_text()):
        assert sha(ROOT/item['destination'])==item['sha256'],item['destination']

@pytest.mark.parametrize('n',[50,827])
def test_frozen_scores_and_ranks(n):
    a=load_numeric(ROOT/'inputs/pre_numeric.csv')[:n]
    rows=read(ROOT/'inputs/pre_labels.csv')[:n]
    ids=[r['source_accession'] for r in rows]; lex={s:i for i,s in enumerate(sorted(ids))}
    s=rounded(scores(components(a),WEIGHTS)); _,r=rank(s,frozen_ties(a,[lex[i] for i in ids]))
    assert list(s)==[float(x['G_pre1_target_engineering_composite_score']) for x in rows]
    assert list(r)==[int(x['G_pre1_target_engineering_composite_rank']) for x in rows]
    i=ids.index('Q9T0Z9');assert s[i]==.856763 and r[i]==1

def test_rank_integrity_and_trial_accounting():
    for scope,n in [('reference_full827',827),('reference_top50',50)]:
        out=ROOT/'outputs'/scope
        summary=read(out/'Supplementary_Table_S4.csv')
        for row in summary:
            r=read(out/(('baseline' if row['analysis']=='full' else row['analysis'])+'.csv'))
            assert len({x['accession'] for x in r})==n
            assert sorted(int(x['rank']) for x in r)==list(range(1,n+1))
        freq=read(out/'candidate_inclusion_frequency.csv')
        assert sum(float(r['top5_frequency']) for r in freq)==pytest.approx(5)
        assert len(read(out/'weight_perturbation_10000.csv'))==10000

def test_no_accession_forcing(tmp_path):
    for p in (ROOT/'src').glob('*.py'):
        text=p.read_text();assert 'Q9T0Z9' not in text and 'Gp17' not in text
    before,_=reproduce_feedback(ROOT)
    ids=[r['source_accession'] for r in read(ROOT/'inputs/frozen_top50.csv')]
    mapping={s:f'ANON_{i:03d}' for i,s in enumerate(ids)}
    # Relabel all source fields, including the registry and anchor. Scores must follow data.
    for p in (ROOT/'inputs').rglob('*.csv'):
        dest=tmp_path/p.relative_to(ROOT);dest.parent.mkdir(parents=True,exist_ok=True)
        text=p.read_text(encoding='utf-8-sig')
        for old,new in sorted(mapping.items(),key=lambda kv:-len(kv[0])):text=text.replace(old,new)
        dest.write_text(text,encoding='utf-8')
    after,_=reproduce_feedback(tmp_path)
    for b,a in zip(before,after):
        assert a['accession']==mapping[b['accession']]
        for key in ['G_pro1_rank','G_pro1_score','C_verified_score','wet_lab_status']:assert a[key]==b[key]
    values=np.array([.2,.9,.1]);order,_=rank(values,analysis_ties(3,20260906))
    assert order[0]==1

@pytest.mark.parametrize('column',['C_verified','wet_lab','ELISA','food_matrix','enrichment','G_pro1'])
def test_no_wet_lab_leakage(tmp_path,column):
    p=tmp_path/'tainted.csv';lines=(ROOT/'inputs/pre_numeric.csv').read_text().splitlines()
    p.write_text(lines[0]+','+column+'\n'+lines[1]+',1\n',encoding='utf-8')
    with pytest.raises(ValueError,match='exact frozen numeric schema'):load_numeric(p)
    # Pure pre scorer has no file I/O except its one schema-checked numeric loader.
    text=(ROOT/'src/scoring.py').read_text()
    assert 'feedback' not in text.replace('pre-feedback','').replace('outcome inputs','')

def test_not_tested_semantics():
    rows,_=reproduce_feedback(ROOT)
    for r in rows:
        if r['wet_lab_status']=='not_tested': assert r['C_verified_status_score']=='0.500000'
    assert sum(r['wet_lab_status']=='not_tested' for r in rows)==49
    assert register_verified([])=={}
    with pytest.raises(ValueError):
        register_verified([dict(evidence_layer='C_verified',accession='x',evidence_status='unknown')])

def test_gpro1_reproduction():
    rows,comps=reproduce_feedback(ROOT)
    frozen=read(ROOT/'inputs/feedback/expected_top50.csv')
    for a,b in zip(rows,frozen):
        for k in ['accession','G_pro1_rank','G_pro1_score','deployment_readiness_score',
                  'C_verified_status_score','C_verified_anchor_similarity_score','C_verified_score','wet_lab_status']:
            assert a[k]==b[k],(a['accession'],k)

def test_same_seed_byte_identical_outputs(tmp_path):
    for name in ['a','b']:
        subprocess.run([sys.executable,'-B',str(ROOT/'scripts/run_analysis.py'),'--out',str(tmp_path/name)],check=True,capture_output=True)
    expected=json.loads((ROOT/'manifests/expected_output_hashes.json').read_text())
    for p in (tmp_path/'a').iterdir():
        assert p.read_bytes()==(tmp_path/'b'/p.name).read_bytes(),p.name
        assert sha(p)==expected['outputs/reference_full827/'+p.name]

def test_expected_output_hashes():
    for path,h in json.loads((ROOT/'manifests/expected_output_hashes.json').read_text()).items():
        assert sha(ROOT/path)==h,path

def test_figure_export():
    d=ROOT/'figures/Figure6'
    with Image.open(d/'Figure6.png') as im:
        assert im.size==(7200,5400);assert im.info['dpi'][0]==pytest.approx(600,abs=.01)
    with Image.open(d/'Figure6.tiff') as im:
        assert im.size==(7200,5400);assert im.tag_v2[259]==5
    assert '<text' in (d/'Figure6.svg').read_text()

def test_panel_c_label_layout_full827(tmp_path,monkeypatch):
    sys.path.insert(0,str(ROOT/'scripts'))
    import figure6
    from matplotlib.figure import Figure
    captured=[]
    monkeypatch.setattr(Figure,'savefig',lambda fig,*args,**kwargs: captured.append(fig))
    figure6.make(ROOT/'outputs/reference_full827',tmp_path/'figure')
    fig=captured[0];fig.canvas.draw()
    ax=fig.axes[2];renderer=fig.canvas.get_renderer()
    frequencies=read(ROOT/'outputs/reference_full827/candidate_inclusion_frequency.csv')
    assert len(frequencies)==827
    values=sorted([100*float(r['top5_frequency']) for r in frequencies if float(r['top5_frequency'])>0],reverse=True)
    assert len(ax.texts)==len(values)==9
    for bar,label,value in zip(ax.patches,ax.texts,values):
        bbox=label.get_window_extent(renderer);barbox=bar.get_window_extent(renderer)
        assert label.get_text()==f'{value:.2f}%'
        assert ax.bbox.contains(bbox.x0,bbox.y0) and ax.bbox.contains(bbox.x1,bbox.y1)
        if value>=95:
            assert label.get_ha()=='right'
            assert barbox.contains(bbox.x0,bbox.y0) and barbox.contains(bbox.x1,bbox.y1)
        else:
            assert label.get_ha()=='left' and bbox.x0>barbox.x1

def test_panel_c_label_threshold():
    sys.path.insert(0,str(ROOT/'scripts'))
    import figure6
    fig,ax=figure6.plt.subplots(figsize=(6,3))
    values=[0.02,94.99,95.0,100.0]
    bars=ax.barh(range(len(values)),values);ax.set_xlim(0,120)
    texts=figure6.label_panel_c(ax,bars,values)
    fig.canvas.draw();renderer=fig.canvas.get_renderer()
    for value,bar,text in zip(values,bars,texts):
        bounds=text.get_window_extent(renderer);edge=bar.get_window_extent(renderer).x1
        assert bounds.x0>=ax.bbox.x0 and bounds.x1<=ax.bbox.x1
        assert (bounds.x1<edge and text.get_ha()=='right') if value>=95 else (bounds.x0>edge and text.get_ha()=='left')
    figure6.plt.close(fig)
