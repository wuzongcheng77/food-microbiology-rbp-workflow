import test_fixture_support
import sys, unittest, tempfile, hashlib, json, csv
from pathlib import Path
from unittest.mock import patch
sys.path.insert(0, str(test_fixture_support.ROOT/'tools'))
from yrbp_evidence_bridge_v0_4 import e6_annotations, structure_sequence_audit
import yrbp_v0_2_rbp_repro_runner as runner

class BridgeTests(unittest.TestCase):
    def test_counts_do_not_impute_predictors(self):
        row = e6_annotations('probe', 'ACCKKK')
        self.assertEqual((row['sequence_length_aa'],row['cysteine_count'],row['lysine_count']), (6,2,3))
        self.assertEqual(row['cysteine_handle_count'], '')
        self.assertEqual(row['predicted_solubility_label'], '')
        self.assertEqual(row['engineering_scoring_ready'], 'false')

    def test_structure_residue_correspondence_and_hash(self):
        with tempfile.TemporaryDirectory() as tmp:
            path=Path(tmp)/'probe.pdbqt'
            path.write_text('ATOM      1  CA  ALA A   1       0.000   0.000   0.000\nATOM      2  CA  CYS A   2       1.000   1.000   1.000\n')
            spec={'protein_pdbqt':path.name,'protein_sha256':hashlib.sha256(path.read_bytes()).hexdigest()}
            self.assertTrue(structure_sequence_audit('AC',spec,tmp)['automatic_acceptance'])
            self.assertFalse(structure_sequence_audit('ACK',spec,tmp)['automatic_acceptance'])
            self.assertFalse(structure_sequence_audit('AD',spec,tmp)['automatic_acceptance'])
            spec['protein_sha256']='wrong'
            self.assertEqual(structure_sequence_audit('AC',spec,tmp)['status'],'structure_hash_mismatch')

    def test_receptor_evidence_without_structure(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);f=root/'input.fasta';f.write_text('>external_test\nACDEFGHIKLMNPQRSTVWY\n')
            ev=root/'evidence.json';ev.write_text(json.dumps({'external_test':{'receptor_class':'protein','evidence_scope':'phage_level'}}))
            with patch.object(sys,'argv',['runner','--fasta',str(f),'--receptor-evidence',str(ev),'--out',str(root/'out')]):
                runner.main()
            audit=json.loads((root/'out/receptor_contract_audit.json').read_text())
            self.assertIn('receptor_class_outside_this_adapter',audit['external_test']['reasons'])
            self.assertTrue((root/'out/fresh_e6_annotations.csv').exists())

    def test_unknown_evidence_id_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);f=root/'input.fasta';f.write_text('>one\nACDE\n')
            ev=root/'ev.json';ev.write_text('{"two":{}}')
            with patch.object(sys,'argv',['runner','--fasta',str(f),'--receptor-evidence',str(ev),'--out',str(root/'out')]),patch.object(runner,'run') as called:
                with self.assertRaises(ValueError): runner.main()
                called.assert_not_called()

if __name__=='__main__': unittest.main()
