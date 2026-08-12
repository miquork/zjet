#!/usr/bin/env python3
"""Safely inspect or remove one completed Z+jet HTCondor campaign."""

import argparse
import getpass
import json
import os
import shutil
import subprocess
from pathlib import Path
from typing import Optional

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


def external_log_directory(metadata: dict, campaign: str,
                           campaign_dir: Path) -> Optional[Path]:
    value = metadata.get("storage", {}).get("log_directory", "")
    if not value:
        return None
    path = Path(value).expanduser().resolve()
    internal = (campaign_dir/"logs").resolve()
    if path == internal:
        return None
    user = getpass.getuser()
    safe_root = Path(
        f"/afs/cern.ch/work/{user[0]}/{user}/zjet-condor").resolve()
    if path.parent != safe_root or path.name != campaign:
        raise ValueError(
            f"refusing external log directory outside {safe_root}/CAMPAIGN: "
            f"{path}")
    return path


def require_cms_proxy(minimum_seconds: int = 600) -> Path:
    """Fail before any EOS deletion unless a usable CMS proxy is available."""
    if shutil.which("voms-proxy-info") is None:
        raise RuntimeError("voms-proxy-info is unavailable; run cleanup on lxplus")
    candidates = []
    configured = os.environ.get("X509_USER_PROXY")
    if configured:
        candidates.append(Path(configured).expanduser())
    candidates.append(
        Path.home()/"private"/f"x509up_u{os.getuid()}")
    default = subprocess.run(
        ["voms-proxy-info","-path"],check=False,text=True,
        stdout=subprocess.PIPE,stderr=subprocess.DEVNULL)
    if default.returncode == 0 and default.stdout.strip():
        candidates.append(Path(default.stdout.strip()))
    for candidate in candidates:
        if not candidate.is_file():
            continue
        lifetime = subprocess.run(
            ["voms-proxy-info","-file",str(candidate),"-timeleft"],
            check=False,text=True,stdout=subprocess.PIPE,
            stderr=subprocess.DEVNULL)
        fqan = subprocess.run(
            ["voms-proxy-info","-file",str(candidate),"-fqan"],
            check=False,text=True,stdout=subprocess.PIPE,
            stderr=subprocess.DEVNULL)
        try:
            seconds = int(lifetime.stdout.strip())
        except ValueError:
            seconds = 0
        if (lifetime.returncode == 0 and fqan.returncode == 0 and
                seconds >= minimum_seconds and "/cms/" in fqan.stdout):
            os.environ["X509_USER_PROXY"] = str(candidate.resolve())
            print(f"Using CMS proxy {candidate.resolve()}; "
                  f"lifetime {seconds/3600:.1f} h")
            return candidate.resolve()
    raise RuntimeError(
        "remote deletion requires a CMS proxy with at least 10 minutes "
        "remaining; run: voms-proxy-init --rfc --voms cms --valid 192:00")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("campaign", help="campaign name or directory")
    parser.add_argument("--delete", action="store_true",
                        help="archive provenance and remove bulky AFS intermediates")
    parser.add_argument(
        "--delete-remote-results", action="store_true",
        help="also remove exactly the EOS partial outputs listed in campaign.json")
    parser.add_argument(
        "--discard-unmerged", action="store_true",
        help="allow removal of an intentionally abandoned, unmerged test campaign")
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
    external_logs = external_log_directory(
        metadata,campaign_dir.name,campaign_dir)
    remote_results = is_remote(result_directory)
    provenance_name = f"zjet_{metadata['campaign']}_provenance.json"
    provenance = campaign_dir/provenance_name
    merge_log = campaign_dir/f"zjet_{metadata['campaign']}_merge.log"

    print(f"Campaign directory: {campaign_dir}")
    print(f"AFS space used: {human_size(directory_size(campaign_dir))}")
    print(f"Partial result storage: {result_directory}")
    if external_logs:
        print(f"External job logs: {external_logs}")
    print(f"Jobs: {len(metadata['jobs'])}")
    if not args.delete:
        print("Dry run only. Pass --delete after a successful merge.")
        if remote_results:
            print("Add --delete-remote-results to remove EOS partial ROOT files too.")
        return

    merged = provenance.is_file() and merge_log.is_file()
    if not merged and not args.discard_unmerged:
        raise RuntimeError(
            "merge provenance is missing; merge first or explicitly pass "
            "--discard-unmerged for an abandoned test campaign")
    if args.delete_remote_results:
        if not remote_results:
            raise ValueError("campaign does not use remote result storage")
        # Authenticate before creating archives or removing any local/remote
        # object. This prevents a half-started cleanup due to a missing proxy.
        require_cms_proxy()
    archive_dir = args.archive_dir.expanduser().resolve()
    archive_dir.mkdir(parents=True,exist_ok=True)
    if merged:
        shutil.copy2(provenance,archive_dir/provenance.name)
        shutil.copy2(merge_log,archive_dir/merge_log.name)

    if args.delete_remote_results:
        urls = [job["result_path"] for job in metadata["jobs"]]
        print(f"Removing {len(urls)} manifest-listed EOS partial outputs...")
        removed_count, missing_count = remove_remote_files(urls)
        print(f"EOS cleanup: removed {removed_count}; already absent "
              f"{missing_count}")

    removed = []
    for name in ("chunks","logs","results"):
        target = campaign_dir/name
        if target.is_dir():
            shutil.rmtree(target)
            removed.append(name)
    if external_logs and external_logs.is_dir():
        shutil.rmtree(external_logs)
        removed.append(str(external_logs))
    print(f"Removed campaign intermediates: {', '.join(removed) or 'none'}")
    print(f"Retained the small manifest in {campaign_dir}")
    if merged:
        print(f"Archived provenance in {archive_dir}")
    else:
        print("Discarded an unmerged campaign; no merge provenance was available.")


if __name__ == "__main__":
    try:
        main()
    except (RuntimeError,FileNotFoundError,ValueError) as error:
        print(f"ERROR: {error}")
        raise SystemExit(1)
