import sys
from pathlib import Path

# make tests/r_compare.py importable from the tests and the benchmarks
sys.path.insert(0, str(Path(__file__).parent))
