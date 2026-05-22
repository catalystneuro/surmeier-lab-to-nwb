"""Compare locally-converted NWB files against the embargoed dandiset 001538.

Reads DANDI_API_TOKEN from the project's .env. Walks ./nwb_files/ recursively
for local NWB files, then lists all assets currently on the dandiset and
compares. Output is filename-only (not byte-for-byte content), plus a per-
category breakdown so we can confirm session counts match per figure.

Usage:
    uv run scripts/compare_local_vs_dandi.py
"""

from __future__ import annotations

import os
import sys
from collections import Counter
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
NWB_FILES_DIR = REPO_ROOT / "nwb_files"
DANDISET_ID = "001538"


def _load_env() -> None:
    env_path = REPO_ROOT / ".env"
    if not env_path.exists():
        return
    for line in env_path.read_text().splitlines():
        line = line.strip()
        if not line or line.startswith("#") or "=" not in line:
            continue
        k, v = line.split("=", 1)
        os.environ.setdefault(k.strip(), v.strip())


def _local_files() -> dict[str, Path]:
    """Map basename -> full path for all local .nwb files."""
    out: dict[str, Path] = {}
    if not NWB_FILES_DIR.exists():
        print(f"Local nwb_files dir not found: {NWB_FILES_DIR}", file=sys.stderr)
        return out
    for p in NWB_FILES_DIR.rglob("*.nwb"):
        out[p.name] = p
    return out


def _session_token_from_name(nwb_path: str) -> str:
    """Extract the ses-XXXXX token from a DANDI-style filename.

    Filenames look like:
        sub-<id>/sub-<id>_ses-<token>_<modalities>.nwb
    """
    fname = nwb_path.split("/")[-1]
    parts = fname.split("_")
    for p in parts:
        if p.startswith("ses-"):
            return p[len("ses-") :]
    return fname


def main() -> int:
    _load_env()
    from dandi.dandiapi import DandiAPIClient  # type: ignore

    print(f"Local: scanning {NWB_FILES_DIR} ...")
    local = _local_files()
    print(f"Local: found {len(local)} .nwb files")

    print(f"DANDI: fetching asset list for dandiset {DANDISET_ID} ...")
    token = os.environ.get("DANDI_API_TOKEN")
    if not token:
        print("WARNING: DANDI_API_TOKEN not set; dandiset 001538 is embargoed and will likely 403", file=sys.stderr)
    with DandiAPIClient(token=token) as client:
        ds = client.get_dandiset(DANDISET_ID, "draft")
        dandi_paths: dict[str, str] = {}
        for asset in ds.get_assets():
            fname = asset.path.split("/")[-1]
            dandi_paths[fname] = asset.path
    print(f"DANDI: found {len(dandi_paths)} assets")

    local_names = set(local.keys())
    dandi_names = set(dandi_paths.keys())

    both = local_names & dandi_names
    only_local = local_names - dandi_names
    only_dandi = dandi_names - local_names

    print()
    print(f"In both:       {len(both)}")
    print(f"Only local:    {len(only_local)}")
    print(f"Only on DANDI: {len(only_dandi)}")

    # Session-token-based comparison (more robust to subject-id changes).
    local_sessions = {_session_token_from_name(str(p)) for p in local.values()}
    dandi_sessions = {_session_token_from_name(p) for p in dandi_paths.values()}
    both_sessions = local_sessions & dandi_sessions
    only_local_sessions = local_sessions - dandi_sessions
    only_dandi_sessions = dandi_sessions - local_sessions

    print()
    print("By session token (more robust to subject-id changes):")
    print(f"  In both:       {len(both_sessions)}")
    print(f"  Only local:    {len(only_local_sessions)}")
    print(f"  Only on DANDI: {len(only_dandi_sessions)}")

    # Per-category counts (local).
    print()
    print("Local files by category/figure:")
    cat_counts: Counter[str] = Counter()
    for p in local.values():
        rel = p.relative_to(NWB_FILES_DIR)
        if len(rel.parts) >= 2:
            cat_counts[f"{rel.parts[0]}/{rel.parts[1]}"] += 1
        else:
            cat_counts[str(rel.parts[0])] += 1
    for cat, n in sorted(cat_counts.items(), key=lambda kv: -kv[1]):
        print(f"  {n:4d}  {cat}")

    # Sample of differences.
    if only_local:
        print()
        print(f"Sample of only-local filenames ({min(len(only_local), 20)} of {len(only_local)}):")
        for n in sorted(only_local)[:20]:
            print(f"  {n}")
    if only_dandi:
        print()
        print(f"Sample of only-DANDI filenames ({min(len(only_dandi), 20)} of {len(only_dandi)}):")
        for n in sorted(only_dandi)[:20]:
            print(f"  {n}")
    if only_local_sessions:
        print()
        print(f"Session tokens only in local (first 20 of {len(only_local_sessions)}):")
        for s in sorted(only_local_sessions)[:20]:
            print(f"  {s}")
    if only_dandi_sessions:
        print()
        print(f"Session tokens only on DANDI (first 20 of {len(only_dandi_sessions)}):")
        for s in sorted(only_dandi_sessions)[:20]:
            print(f"  {s}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
