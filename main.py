"""Launch MendelPRIME from this source checkout."""

import sys
from pathlib import Path

if __name__ == "__main__":
    sys.path.insert(0, str(Path(__file__).resolve().parent / "src"))
    from mendelprime.runtime.cli import main

    raise SystemExit(main())
