"""Module entry: ``python -m knottidy``."""

from __future__ import annotations

import sys

from knottidy.cli import main

if __name__ == "__main__":
    sys.exit(main())
