from pathlib import Path
import sys
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'tools'))
import yrbp_sequence_registry_v0_2 as registry
registry.load_registry.__defaults__=(ROOT/'data/frozen_rank.csv',ROOT/'data/sequence_registry.csv')
import yrbp_v0_2_rbp_repro_runner as runner
runner.FROZEN_RANK_CSV=ROOT/'data/frozen_rank.csv'
runner.PROTECTED={'supplied_frozen_rank_table':runner.FROZEN_RANK_CSV}
