#!/usr/bin/env python3
"""Summarize returned outputs for a prepared Z+jet HTCondor campaign."""

import argparse
import json
from pathlib import Path
from typing import Dict, List


REPOSITORY = Path(__file__).resolve().parents[1]


def campaign_path(value: str) -> Path:
    direct = Path(value).expanduser()
    if direct.is_dir():
        return direct.resolve()
    return (REPOSITORY/"condor"/"jobs"/value).resolve()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("campaign", help="campaign name or directory")
    parser.add_argument("--details", action="store_true",
                        help="list every job that is not ready")
    args = parser.parse_args()

    campaign_dir = campaign_path(args.campaign)
    metadata_path = campaign_dir/"campaign.json"
    if not metadata_path.is_file():
        raise FileNotFoundError(f"campaign metadata not found: {metadata_path}")
    metadata = json.loads(metadata_path.read_text(encoding="utf-8"))

    counts: Dict[str, Dict[str, int]] = {
        sample: {"expected": 0, "ready": 0, "empty": 0, "missing": 0,
                 "unvalidated": 0}
        for sample in ("mc","data")
    }
    incomplete: List[str] = []
    for job in metadata["jobs"]:
        sample = job["sample"]
        chunk_id = job["chunk_id"]
        result = Path(job["result_path"])
        log = campaign_dir/"logs"/f"{sample}_{chunk_id}.out"
        counts[sample]["expected"] += 1

        if not result.exists():
            state = "missing"
        elif result.stat().st_size == 0:
            state = "empty"
        elif not log.is_file() or "Job finished successfully" not in \
                log.read_text(encoding="utf-8",errors="replace"):
            state = "unvalidated"
        else:
            state = "ready"
        counts[sample][state] += 1
        if state != "ready":
            incomplete.append(f"{sample}_{chunk_id}: {state}")

    print(f"Campaign: {metadata['campaign']}")
    print(f"Input files: {metadata['mc_files']} MC, "
          f"{metadata['data_files']} data")
    for sample in ("mc","data"):
        item = counts[sample]
        print(f"{sample.upper():4s}: {item['ready']}/{item['expected']} ready; "
              f"{item['missing']} missing, {item['empty']} empty, "
              f"{item['unvalidated']} unvalidated")

    if incomplete:
        shown = incomplete if args.details else incomplete[:10]
        print("Not ready to merge:")
        for item in shown:
            print(f"  {item}")
        if not args.details and len(incomplete)>len(shown):
            print(f"  ... and {len(incomplete)-len(shown)} more "
                  f"(pass --details to list all)")
        raise SystemExit(1)

    print("READY TO MERGE")
    print(f"Run: python3 scripts/merge_condor.py {metadata['campaign']} --force")


if __name__ == "__main__":
    main()
