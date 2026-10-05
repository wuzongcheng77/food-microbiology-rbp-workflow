import test_fixture_support
import sys, unittest, tempfile, json
from pathlib import Path
from unittest.mock import patch
sys.path.insert(0,str(test_fixture_support.ROOT/'tools'))
from yrbp_input_contract_v0_3 import sequence_features, receptor_contract
import yrbp_v0_2_rbp_repro_runner as runner

class ContractTests(unittest.TestCase):
    def test_raw_features_no_score_substitution(self):
        f=sequence_features('ACAC')
        self.assertEqual(f['count_A'],2)
        self.assertEqual(f['fraction_A'],0.5)
        self.assertEqual(f['composition_entropy_bits'],1)
        self.assertNotIn('docking_support_score',f)
        for aa in ['B','Z']:
            with self.assertRaises(ValueError):sequence_features('AC'+aa)
    def test_missing_not_negative_evidence(self):
        self.assertEqual(receptor_contract('hash',{})['status'],'requires_review')
    def test_consistency_is_not_validation(self):
        c=dict.fromkeys(['target_taxon','target_context','receptor_identity','evidence_source','reviewer','reviewed_on','ligand_identity','grid_rationale_source','structure_mapping_source'],'synthetic-test-only')
        c.update(sequence_sha256='hash',receptor_class='glycan',ligand_class='glycan',method='vina_prepared_ligand',evidence_scope='candidate_specific',review_status='reviewed')
        self.assertEqual(receptor_contract('hash',{'receptor_evidence':c})['status'],'declaration_consistent')
        c['ligand_class']='protein'
        self.assertIn('receptor_ligand_class_mismatch',receptor_contract('hash',{'receptor_evidence':c})['reasons'])
        c['ligand_class']='glycan';c['evidence_scope']='phage_level'
        self.assertEqual(receptor_contract('hash',{'receptor_evidence':c})['status'],'requires_review')
    def test_runner_stops_before_docking(self):
        with tempfile.TemporaryDirectory() as d:
            d=Path(d);f=d/'input.fasta';f.write_text('>testBZ\nACDEFGHIKLMNPQRSTVWY\n')
            with patch.object(sys,'argv',['runner','--fasta',str(f),'--out',str(d/'run'),'--vina','must-not-execute']),patch.object(runner,'run') as call:
                with self.assertRaises(ValueError):runner.main()
                call.assert_not_called()
            self.assertTrue((d/'run/receptor_preflight_blocked.json').exists())

if __name__=='__main__':unittest.main()
