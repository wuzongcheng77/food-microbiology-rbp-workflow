import test_fixture_support
import sys,unittest,tempfile,hashlib,json,subprocess
from pathlib import Path
from unittest.mock import patch
sys.path.insert(0,str(test_fixture_support.ROOT/'tools'))
from yrbp_candidate_intake_support_v0_2 import fasta,run
from yrbp_sequence_registry_v0_2 import load_registry,membership
class IdentityTests(unittest.TestCase):
    def test_exact_pool_counts(self):
        rows=load_registry();self.assertEqual(len(rows),827)
        self.assertEqual(len({r['sequence_sha256'] for r in rows}),744)
    def test_renamed_sequence_not_external(self):
        rows=load_registry();r=rows[0]
        result=membership('renamed_test_probe',r['sequence_sha256'],rows)
        self.assertTrue(result['sequence_in_historical_827'])
        self.assertEqual(result['historical_identity_status'],'sequence_match_under_other_id')
    def test_title_BZ_allowed_sequence_BZ_rejected(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'fixture.fasta'
            p.write_text('>BZ_identifier\nACDEFGHIKLMNPQRSTVWY\n');self.assertEqual(len(fasta(p)),1)
            for aa in ['B','Z']:
                p.write_text('>invalid_fixture\nMALK'+aa+'\n')
                with self.assertRaises(ValueError):fasta(p)
    def test_tampered_pose_forces_reexecution(self):
        demo=Path(__file__).resolve().parents[1]/'02_workflow/demo';calls=[]
        def fake(args,**kw):
            calls.append(args);Path(args[args.index('--out')+1]).write_text('REMARK VINA RESULT: -3.25\n')
            return subprocess.CompletedProcess(args,0,b'fake-test-only',b'')
        with tempfile.TemporaryDirectory() as d,patch('yrbp_candidate_intake_support_v0_2.subprocess.run',side_effect=fake):
            out=Path(d)/'run';run(demo/'candidates.fasta',demo/'structures.json',out,'fake-test-only')
            n=len(calls);next(out.rglob('pose.pdbqt')).write_text('REMARK VINA RESULT: -99.0\n')
            rows=run(demo/'candidates.fasta',demo/'structures.json',out,'fake-test-only',resume=True)
            self.assertEqual(len(calls),n+1);self.assertNotIn(-99,[r['affinity_kcal_mol'] for r in rows])
if __name__=='__main__':unittest.main(verbosity=2)
