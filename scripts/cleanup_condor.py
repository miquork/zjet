#!/usr/bin/env python3
"""Safely inspect or remove one completed Z+jet HTCondor campaign."""

import argparse
import json
import shutil
from pathlib import Path

from condor_storage import is_remote, remove_remote_files


REPOSITORY = Path(__file__).resolve().parents[1]
CAMPAIGN_ROOT = (REPOSITORY/"condor"/"jobs").resolve()


def campaign_path(value: str) -> Path:
    direct = Path(value).expanduser()
    path = direct.resolve() if direct.is_dir() else (CAMPAIGN_ROOT/value).resolve()
    try:
        path.relative_to(CAMPAIGN_ROOT)
    except ValueError as error:
        raise ValueError(f"campaign must be below {CAMPAIGN_ROOT}") from error
    return path


def directory_size(path: Path) -> int:
    return sum(item.stat().st_size for item in path.rglob("*") if item.is_file())


def human_size(value: int) -> str:
    amount = float(value)
    for unit in ("B","KiB","MiB","GiB","TiB"):
        if amount < 1024.0 or unit == "TiB":
            return f"{amount:.1f} {unit}"
        amount /= 1024.0
    return f"{amount:.1f} TiB"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("campaign", help="campaign name or directory")
    parser.add_argument("--delete", action="store_true",
                        help="archive provenance and remove bulky AFS intermediates")
    parser.add_argument(
        "--delete-remote-results", action="store_true",
        help="also remove exactly the EOS partial outputs listed in campaign.json")
    parser.add_argument("--archive-dir", type=Path,
                        default=REPOSITORY/"rootfiles"/"campaigns")
    args = parser.parse_args()

    campaign_dir = campaign_path(args.campaign)
    metadata_path = campaign_dir/"campaign.json"
    if not metadata_path.is_file():
        raise FileNotFoundError(f"campaign metadata not found: {metadata_path}")
    metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
    if metadata.get("campaign") != campaign_dir.name:
        raise RuntimeError("campaign name does not match campaign.json")

    result_directory = metadata.get("storage", {}).get("result_directory", "")
    remote_results = is_remote(result_directory)
    provenance_name = f"zjet_{metadata['campaign']}_provenance.json"
    provenance = campaign_dir/provenance_name
    merge_log = campaign_dir/f"zjet_{metadata['campaign']}_merge.log"

    print(f"Campaign directory: {campaign_dir}")
    print(f"AFS space used: {human_size(directory_size(campaign_dir))}")
    print(f"Partial result storage: {result_directory}")
    print(f"Jobs: {len(metadata['jobs'])}")
    if not args.delete:
        print("Dry run only. Pass --delete after a successful merge.")
        if remote_results:
            print("Add --delete-remote-results to remove EOS partial ROOT files too.")
        return

    if not provenance.is_file() or not merge_log.is_file():
        raise RuntimeError(
            "merge provenance is missing; run scripts/merge_condor.py first")
    archive_dir = args.archive_dir.expanduser().resolve()
    archive_dir.mkdir(parents=True,exist_ok=True)
    shutil.copy2(provenance,archive_dir/provenance.name)
    shutil.copy2(merge_log,archive_dir/merge_log.name)

    if args.delete_remote_results:
        if not remote_results:
            raise ValueError("campaign does not use remote result storage")
        urls = [job["result_path"] for job in metadata["jobs"]]
        print(f"Removing {len(urls)} manifest-listed EOS partial outputs...")
        remove_remote_files(urls)

    removed = []
    for name in ("chunks","logs","results"):
        target = campaign_dir/name
        if target.is_dir():
            shutil.rmtree(target)
            removed.append(name)
    print(f"Removed campaign intermediates: {', '.join(removed) or 'none'}")
    print(f"Retained the small manifest in {campaign_dir}")
    print(f"Archived provenance in {archive_dir}")


if __name__ == "__main__":
    main()
