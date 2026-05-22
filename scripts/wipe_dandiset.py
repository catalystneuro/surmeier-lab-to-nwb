"""Delete all assets from dandiset 001538 (the embargoed draft) before re-uploading.

Per-mouse-day refactor changes the asset paths, so the new upload would otherwise
coexist with the old per-modality files. This script wipes the draft so the new
upload represents the full dandiset state.

Safety: requires explicit --confirm. Prints what it WOULD delete first.
"""

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

from dandi.dandiapi import DandiAPIClient
from dotenv import load_dotenv

REPO_ROOT = Path(__file__).resolve().parents[1]
load_dotenv(REPO_ROOT / ".env")

DANDISET_ID = "001538"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--version", default="draft", help="Dandiset version (default: draft)")
    parser.add_argument(
        "--confirm", action="store_true", help="Actually perform the deletion. Without this, dry-run only."
    )
    args = parser.parse_args()

    token = os.environ.get("DANDI_API_TOKEN")
    if not token:
        print("DANDI_API_TOKEN not set in env or .env", file=sys.stderr)
        return 1

    client = DandiAPIClient(token=token)
    client.authenticate(token=token)
    ds = client.get_dandiset(DANDISET_ID, args.version)
    assets = list(ds.get_assets())

    print(f"Dandiset {DANDISET_ID} (version={args.version}): {len(assets)} assets")
    if not assets:
        print("Nothing to delete.")
        return 0

    # Show breakdown by top-level subject directory
    from collections import Counter

    breakdown = Counter()
    for a in assets:
        first_dir = a.path.split("/")[0]
        breakdown[first_dir] += 1
    print()
    print("Breakdown by subject:")
    for sub, n in sorted(breakdown.items()):
        print(f"  {sub:60s} {n} assets")

    total_size_mb = sum(a.size for a in assets) / 1024 / 1024
    print()
    print(f"Total size: {total_size_mb:.1f} MB ({total_size_mb/1024:.2f} GB)")

    if not args.confirm:
        print()
        print("Dry-run only. Re-run with --confirm to actually delete.")
        return 0

    print()
    print(f"DELETING {len(assets)} assets from {DANDISET_ID}/{args.version}...")
    failed = 0
    for i, a in enumerate(assets):
        try:
            a.delete()
            if (i + 1) % 25 == 0:
                print(f"  deleted {i+1}/{len(assets)}")
        except Exception as e:
            failed += 1
            print(f"  FAIL {a.path[:80]}: {e}")

    print()
    print(f"Done. {len(assets) - failed} deleted, {failed} failed.")
    return 0 if failed == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
