#!/usr/bin/env python3
"""Summarize returned outputs for a prepared Z+jet HTCondor campaign."""

import argparse
import json
from pathlib import Path
from typing import Dict, List

from condor_storage import is_remote, remote_file_sizes


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
    parser.add_argument(
        "--additional-log-dir", action="append", type=Path, default=[],
        help=("also search this directory for worker .out logs; repeat for "
              "campaigns recovered from an AFS quota hold"))
    args = parser.parse_args()

    campaign_dir = campaign_path(args.campaign)
    metadata_path = campaign_dir/"campaign.json"
    if not metadata_path.is_file():
        raise FileNotFoundError(f"campaign metadata not found: {metadata_path}")
    metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
    result_directory = metadata.get("storage", {}).get("result_directory", "")
    configured_log_directory = Path(metadata.get("storage", {}).get(
        "log_directory",campaign_dir/"logs"))
    log_directories = [configured_log_directory]
    for value in args.additional_log_dir:
        directory = value.expanduser().resolve()
        if directory not in log_directories:
            log_directories.append(directory)
    remote_results = is_remote(result_directory)
    remote_sizes = remote_file_sizes(result_directory) if remote_results else {}

    counts: Dict[str, Dict[str, int]] = {
        sample: {"expected": 0, "ready": 0, "empty": 0, "missing": 0,
                 "unvalidated": 0}
        for sample in ("mc","data")
    }
    incomplete: List[str] = []
    for job in metadata["jobs"]:
        sample = job["sample"]
        chunk_id = job["chunk_id"]
        result_value = job["result_path"]
        logs = [directory/f"{sample}_{chunk_id}.out"
                for directory in log_directories]
        counts[sample]["expected"] += 1

        if remote_results:
            exists = job["output_file"] in remote_sizes
            empty = exists and remote_sizes[job["output_file"]] == 0
        else:
            result = Path(result_value)
            exists = result.exists()
            empty = exists and result.stat().st_size == 0

        if not exists:
            state = "missing"
        elif empty:
            state = "empty"
        elif not any(log.is_file() and "Job finished successfully" in
                     log.read_text(encoding="utf-8",errors="replace")
                     for log in logs):
            state = "unvalidated"
        else:
            state = "ready"
        counts[sample][state] += 1
        if state != "ready":
            incomplete.append(f"{sample}_{chunk_id}: {state}")

    print(f"Campaign: {metadata['campaign']}")
    print(f"Input files: {metadata['mc_files']} MC, "
          f"{metadata['data_files']} data")
    print(f"Partial results: {result_directory}")
    print("Job logs: " + ", ".join(str(path) for path in log_directories))
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
