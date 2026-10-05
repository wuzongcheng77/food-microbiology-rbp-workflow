import test_fixture_support
"""Integration tests for the 2026-10-03 RBP wiring fixes.

Covers: resume/cache-reuse, incremental flush, execution_mode tracking,
orchestrator historical-827 membership and stage ledger. Synthetic fixture
sequences are software probes, never biological evidence.
"""
import sys, unittest, tempfile, subprocess, json, csv
from pathlib import Path
from unittest.mock import patch

TOOLS = str(test_fixture_support.ROOT/'tools')
sys.path.insert(0, TOOLS)
from yrbp_candidate_intake_support_v0_2 import run  # noqa: E402

DEMO = test_fixture_support.ROOT/'02_workflow/demo'


def fake_vina_factory(calls):
    def fake(args, **kw):
        calls.append(list(args))
        Path(args[args.index('--out') + 1]).write_text('REMARK VINA RESULT: -3.25\n')
        return subprocess.CompletedProcess(args, 0, b'utf8 ok', b'')
    return fake


class ResumeTests(unittest.TestCase):
    def test_resume_reuses_completed_pose(self):
        calls = []
        with tempfile.TemporaryDirectory() as d:
            out = Path(d) / 'run'
            with patch('yrbp_candidate_intake_support_v0_2.subprocess.run', side_effect=fake_vina_factory(calls)):
                run(DEMO / 'candidates.fasta', DEMO / 'structures.json', out, 'fake-vina', resume=False)
            n_first = len(calls)
            self.assertGreaterEqual(n_first, 2)
            # second pass with resume: nothing should re-execute
            with patch('yrbp_candidate_intake_support_v0_2.subprocess.run', side_effect=fake_vina_factory(calls)):
                rows = run(DEMO / 'candidates.fasta', DEMO / 'structures.json', out, 'fake-vina', resume=True)
            self.assertEqual(len(calls), n_first, 'resume must not re-invoke vina')
            modes = {r['candidate_id']: r['execution_mode'] for r in rows}
            self.assertTrue(all(m == 'cache_reuse' for m in modes.values()), modes)

    def test_incremental_results_always_written(self):
        with tempfile.TemporaryDirectory() as d:
            out = Path(d) / 'run'

            def boom(args, **kw):
                Path(args[args.index('--out') + 1]).write_text('REMARK VINA RESULT: -1.0\n')
                return subprocess.CompletedProcess(args, 0, b'', b'')

            with patch('yrbp_candidate_intake_support_v0_2.subprocess.run', side_effect=boom):
                rows = run(DEMO / 'candidates.fasta', DEMO / 'structures.json', out, 'fake-vina', resume=False)
            self.assertTrue((out / 'candidate_results.csv').exists())
            saved = list(csv.DictReader((out / 'candidate_results.csv').open(encoding='utf-8')))
            self.assertEqual(len(saved), len(rows))


class OrchestratorTests(unittest.TestCase):
    def test_membership_and_stage_ledger(self):
        sys.path.insert(0, TOOLS)
        import importlib
        mod = importlib.import_module('yrbp_v0_2_rbp_repro_runner')
        with tempfile.TemporaryDirectory() as d:
            fa = Path(d) / 'in.fasta'
            # Q9T0Z9 is a real frozen-827 accession name; ZZ9NOTIN is not.
            # Sequences are synthetic software probes only.
            fa.write_text('>Q9T0Z9\nMALKDEFGHIVLWQA\n>ZZ9NOTIN8\nMALKDEFGHIVLWQATYVAN\n')
            out = Path(d) / 'out'
            sys.argv = ['x', '--fasta', str(fa), '--out', str(out)]
            mod.main()
            rows = list(csv.DictReader((out / 'candidate_results.csv').open(encoding='utf-8')))
            by = {r['candidate_id']: r for r in rows}
            self.assertEqual(by['Q9T0Z9']['in_historical_827'], 'False')
            self.assertEqual(by['Q9T0Z9']['accession_in_historical_827'], 'True')
            self.assertEqual(by['Q9T0Z9']['historical_identity_status'], 'accession_sequence_mismatch')
            self.assertEqual(by['ZZ9NOTIN8']['in_historical_827'], 'False')
            ledger = json.loads((out / 'stage_ledger.json').read_text(encoding='utf-8'))
            stages = {s['stage']: s['execution_mode'] for s in ledger}
            self.assertEqual(stages['genome_orf_calling'], 'not_done_not_implemented')
            self.assertEqual(stages['automatic_structure_folding'], 'not_done_not_implemented')
            man = json.loads((out / 'repro_manifest.json').read_text(encoding='utf-8'))
            self.assertEqual(man['protected_files_before'], man['protected_files_after'])


if __name__ == '__main__':
    unittest.main(verbosity=2)
