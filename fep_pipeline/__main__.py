"""Allow `python -m fep_pipeline` to run the CLI."""
from .prepare_fep import main

if __name__ == "__main__":
    raise SystemExit(main())
