#!/usr/bin/env python3
"""Prepare a chunked CERN HTCondor campaign for the Z+jet analysis."""

import argparse
import json
import os
import re
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Optional, Tuple


REPOSITORY = Path(__file__).resolve().parents[1]
CAMPAIGN_PATTERN = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_.-]*$")
SAFE_FILENAME_PATTERN = re.compile(r"^[A-Za-z0-9_.+-]+$")


def read_file_list(path: Path, maximum: int) -> List[str]:
    values: List[str] = []
    for line in path.read_text(encoding="utf-8").splitlines():
        value = line.strip()
        if not value or value.startswith("#"):
            continue
        values.append(value)
        if maximum >= 0 and len(values) >= maximum:
            break
    if not values:
        raise ValueError(f"no input files found in {path}")
    return values


def chunks(values: List[str], size: int) -> List[List[str]]:
    return [values[index:index + size] for index in range(0, len(values), size)]


def optional_input(path_value: str, label: str) -> Tuple[Optional[Path], str]:
    if not path_value:
        return None, ""
    path = Path(path_value).expanduser().resolve()
    if not path.is_file():
        raise FileNotFoundError(f"{label} does not exist: {path}")
    if not SAFE_FILENAME_PATTERN.fullmatch(path.name):
        raise ValueError(f"{label} basename is not batch-safe: {path.name}")
    return path, path.name


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mc-list", required=True, type=Path)
    parser.add_argument("--data-list", required=True, type=Path)
    parser.add_argument("--campaign", required=True,
                        help="new campaign name (letters, digits, dot, dash, underscore)")
    parser.add_argument("--files-per-job", type=int, default=10)
    parser.add_argument("--max-mc-files", type=int, default=-1)
    parser.add_argument("--max-data-files", type=int, default=-1)
    parser.add_argument("--golden-json", default="")
    parser.add_argument("--lumi-pileup", default="")
    parser.add_argument("--pileup-weights", default="")
    parser.add_argument("--job-flavour", default="longlunch",
                        choices=("espresso", "microcentury", "longlunch",
                                 "workday", "tomorrow", "testmatch", "nextweek"))
    parser.add_argument("--retries", type=int, default=2)
    args = parser.parse_args()

    if not CAMPAIGN_PATTERN.fullmatch(args.campaign):
        raise ValueError("invalid campaign name")
    if args.files_per_job <= 0:
        raise ValueError("--files-per-job must be positive")
    if args.retries < 0:
        raise ValueError("--retries cannot be negative")

    proxy_value = os.environ.get("X509_USER_PROXY","")
    if not proxy_value:
        raise RuntimeError("X509_USER_PROXY is not set")
    proxy_path = Path(proxy_value).expanduser().resolve()
    if not proxy_path.is_file():
        raise FileNotFoundError(f"X509_USER_PROXY does not exist: {proxy_path}")
    if not str(proxy_path).startswith("/afs/"):
        raise ValueError(
            "X509_USER_PROXY must point to a protected AFS file for CERN "
            "HTCondor; worker schedds cannot read an lxplus-local /tmp proxy")
    try:
        proxy_path.relative_to(REPOSITORY)
    except ValueError:
        pass
    else:
        raise ValueError("do not store an X.509 proxy inside the public repository")

    mc_list = args.mc_list.expanduser().resolve()
    data_list = args.data_list.expanduser().resolve()
    mc_files = read_file_list(mc_list,args.max_mc_files)
    data_files = read_file_list(data_list,args.max_data_files)

    golden_path, golden_name = optional_input(args.golden_json,"golden JSON")
    lumi_path, lumi_name = optional_input(args.lumi_pileup,"lumi pileup file")
    weights_path, weights_name = optional_input(args.pileup_weights,
                                                 "pileup-weight file")
    optional_paths = [path for path in (golden_path,lumi_path,weights_path)
                      if path is not None]
    optional_names = [path.name for path in optional_paths]
    if len(optional_names) != len(set(optional_names)):
        raise ValueError("optional input files must have distinct basenames")
    reserved_names = {
        "zjet.C", "zjet.h", "ZJetLumi.h", "mk_compile.C",
        "run_zjet_job.C",
    }
    if reserved_names.intersection(optional_names):
        raise ValueError("an optional input basename conflicts with source files")

    campaign_dir = REPOSITORY / "condor" / "jobs" / args.campaign
    if campaign_dir.exists():
        raise FileExistsError(
            f"campaign directory already exists: {campaign_dir}; use a new name")
    chunk_dir = campaign_dir / "chunks"
    log_dir = campaign_dir / "logs"
    result_dir = campaign_dir / "results"
    for directory in (chunk_dir,log_dir,result_dir):
        directory.mkdir(parents=True,exist_ok=False)

    jobs: List[Dict[str, object]] = []
    for sample, values, output_tag in (
            ("mc",mc_files,"MC"),("data",data_files,"DATA")):
        for index, group in enumerate(chunks(values,args.files_per_job)):
            chunk_name = f"{sample}_{index:04d}.txt"
            chunk_path = chunk_dir / chunk_name
            chunk_path.write_text("".join(f"{value}\n" for value in group),
                                  encoding="utf-8")
            output_file = f"zjet_{output_tag}_{index:04d}.root"
            jobs.append({
                "sample": sample,
                "chunk_id": f"{index:04d}",
                "chunk_name": chunk_name,
                "chunk_path": str(chunk_path.relative_to(REPOSITORY)),
                "output_file": output_file,
                "result_path": str((result_dir/output_file).resolve()),
                "input_files": len(group),
            })

    manifest_path = campaign_dir / "jobs.tsv"
    manifest_path.write_text(
        "".join(
            "{sample}\t{chunk_id}\t{chunk_path}\t{chunk_name}\t"
            "{output_file}\t{result_path}\n".format(**job)
            for job in jobs),
        encoding="utf-8")

    metadata = {
        "campaign": args.campaign,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "files_per_job": args.files_per_job,
        "mc_files": len(mc_files),
        "data_files": len(data_files),
        "jobs": jobs,
    }
    (campaign_dir/"campaign.json").write_text(
        json.dumps(metadata,indent=2) + "\n",encoding="utf-8")

    common_inputs = [
        "zjet.C", "zjet.h", "ZJetLumi.h", "mk_compile.C",
        "run_zjet_job.C",
    ]
    transfer_inputs = common_inputs + [str(path) for path in optional_paths]
    transfer_inputs.append("$(chunk_path)")
    golden_argument = golden_name or "-"
    lumi_argument = lumi_name or "-"
    weights_argument = weights_name or "-"
    submit_text = f"""universe = vanilla
executable = condor/run_zjet_job.sh
initialdir = {REPOSITORY}
arguments = $(sample) $(chunk_name) $(output_file) {golden_argument} {lumi_argument} {weights_argument}

output = {log_dir}/$(sample)_$(chunk_id).out
error = {log_dir}/$(sample)_$(chunk_id).err
log = {log_dir}/condor.log

getenv = True
MY.WantOS = \"el9\"
request_cpus = 1
request_memory = 2GB
request_disk = 2GB
+JobFlavour = \"{args.job_flavour}\"

use_x509userproxy = True
x509userproxy = {proxy_path}

should_transfer_files = YES
when_to_transfer_output = ON_EXIT
transfer_input_files = {','.join(transfer_inputs)}
transfer_output_files = $(output_file)
transfer_output_remaps = \"$(output_file)=$(result_path)\"

on_exit_remove = (ExitBySignal == False) && (ExitCode == 0)
max_retries = {args.retries}
requirements = (Machine =!= split(LastRemoteHost, \"@\")[1])

queue sample,chunk_id,chunk_path,chunk_name,output_file,result_path from {manifest_path}
"""
    submit_path = campaign_dir / "zjet.sub"
    submit_path.write_text(submit_text,encoding="utf-8")

    mc_jobs = sum(job["sample"] == "mc" for job in jobs)
    data_jobs = sum(job["sample"] == "data" for job in jobs)
    print(f"Prepared campaign {args.campaign}: {len(mc_files)} MC files in "
          f"{mc_jobs} jobs and {len(data_files)} data files in {data_jobs} jobs.")
    print(f"Submit with: condor_submit {submit_path.relative_to(REPOSITORY)}")
    print(f"Results will return to: {result_dir.relative_to(REPOSITORY)}")


if __name__ == "__main__":
    main()
