import test_fixture_support
import sys,unittest,copy,tempfile,subprocess
from pathlib import Path
from unittest.mock import patch
sys.path.insert(0,str(test_fixture_support.ROOT/'tools'))
from yrbp_candidate_intake_v0_2 import audit
class EvidenceTests(unittest.TestCase):
    def setUp(self):
        self.rows=[{'candidate_id':'a','sequence_sha256':'abc'},{'candidate_id':'b','sequence_sha256':'def'}]
        self.policy={'version':'synthetic-contract-test','normalization_reference':'test-only normalized [0,1] values','weights':{'x':1}}
        self.evidence={'a':{'sequence_sha256':'abc','features':{'x':{'status':'observed','value':0,'source':'synthetic test fixture'}}}}
    def test_unknown_not_zero(self):
        r=audit(self.rows,self.policy,self.evidence)
        self.assertEqual((r[0]['lower'],r[0]['upper']),(0,0));self.assertEqual((r[1]['lower'],r[1]['upper']),(0,1))
    def test_order_invariance(self):
        a={x['candidate_id']:x for x in audit(self.rows,self.policy,self.evidence)}
        b={x['candidate_id']:x for x in audit(list(reversed(self.rows)),self.policy,self.evidence)}
        self.assertEqual(a,b)
    def test_sequence_mismatch_rejected(self):
        self.evidence['a']['sequence_sha256']='wrong'
        with self.assertRaises(ValueError):audit(self.rows,self.policy,self.evidence)
    def test_source_required(self):
        del self.evidence['a']['features']['x']['source']
        with self.assertRaises(ValueError):audit(self.rows,self.policy,self.evidence)
    def test_unknown_numeric_rejected(self):
        self.evidence['a']['features']['x']['status']='unknown'
        with self.assertRaises(ValueError):audit(self.rows,self.policy,self.evidence)
    def test_nans_and_invalid_weights(self):
        for w in [0,-1,float('nan'),float('inf'),True]:
            self.policy['weights']['x']=w
            with self.assertRaises(ValueError):audit(self.rows,self.policy,self.evidence)
    def test_foreign_candidate_rejected(self):
        self.evidence['unexpected']={}
        with self.assertRaises(ValueError):audit(self.rows,self.policy,self.evidence)
    def test_rename_invariance(self):
        a=audit(self.rows,self.policy,self.evidence)
        rows=copy.deepcopy(self.rows);rows[0]['candidate_id']='Q9T0Z9'
        b=audit(rows,self.policy,{'Q9T0Z9':self.evidence['a']})
        for key in ['lower','upper','coverage','best_possible_rank','worst_possible_rank']:self.assertEqual(a[0][key],b[0][key])
    def test_windows_log_bytes_preserved(self):
        from yrbp_candidate_intake_support_v0_2 import run
        demo=Path(__file__).resolve().parents[1]/'02_workflow/demo'
        def fake(args,**kwargs):
            self.assertNotIn('text',kwargs)
            Path(args[args.index('--out')+1]).write_text('REMARK VINA RESULT: -1.0\n')
            return subprocess.CompletedProcess(args,0,b'Windows locale: \xb2\xe2\xca\xd4',b'')
        with tempfile.TemporaryDirectory() as d,patch('yrbp_candidate_intake_support_v0_2.subprocess.run',side_effect=fake):
            rows=run(demo/'candidates.fasta',demo/'structures.json',Path(d)/'run','test-only-fake-vina')
            self.assertEqual(len(rows),2)
            files=list((Path(d)/'run').rglob('vina.stdout.bin'))
            self.assertEqual(len(files),2)
            self.assertTrue(all(p.read_bytes()==b'Windows locale: \xb2\xe2\xca\xd4' for p in files))
if __name__=='__main__':unittest.main(verbosity=2)
