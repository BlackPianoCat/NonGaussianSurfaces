"""Run the full test suite without installing the package.

Usage:
    python src/run_tests.py            # all tests
    python src/run_tests.py -k plots   # any extra args go to pytest
Credits: Sebastian Korsak, with Claude.
"""
import sys
from pathlib import Path

import pytest

SRC = Path(__file__).resolve().parent


def main():
    sys.path.insert(0, str(SRC))                       # make `import nongauss` work from a clone
    return pytest.main([str(SRC / 'tests'), '-q', '-p', 'no:cacheprovider', *sys.argv[1:]])


if __name__ == '__main__':
    sys.exit(main())
