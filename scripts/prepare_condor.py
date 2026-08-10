#!/usr/bin/env python3
"""Prepare a chunked CERN HTCondor campaign for the Z+jet analysis."""

import argparse
import hashlib
import json
import os
import re
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Optional, Tuple

from condor_storage import normalize_remote_directory


REPOSITORY = Path(__file__).resolve().parents[1]
CAMPAIGN_PATTERN = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_.-]*$")
SAFE_FILENAME_PATTERN = re.compile(r"^[A-Za-z0-9_.+-]+$")
DEFAULT_JEC_L2 = (REPOSITORY/"CondFormats/JetMETObjects/data"/
                  "RunIII2024Summer24_V2_MC_L2Relative_AK4PUPPI.txt")
DEFAULT_JEC_RESIDUAL = (REPOSITORY/"CondFormats/JetMETObjects/data"/
                        "Prompt24_Run2024I_nib1_V11M_DATA_"
                        "L2L3Residual_AK4PFPuppi.txt")


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


def sha256_bytes(value: bytes) -> str:
    return hashlib.sha256(value).hexdigest()


def file_description(path: Optional[Path]) -> Dict[str, object]:
    if path is None:
        return {"enabled": False}
    return {
        "enabled": True,
        "basename": path.name,
        "sha256": sha256_bytes(path.read_bytes()),
    }


def git_description() -> Dict[str, object]:
    try:
        commit = subprocess.run(
            ["git", "rev-parse", "HEAD"], cwd=REPOSITORY, check=True,
            text=True, stdout=subprocess.PIPE,
            stderr=subprocess.DEVNULL).stdout.strip()
        dirty = bool(subprocess.run(
            ["git", "status", "--porcelain", "--untracked-files=no"],
            cwd=REPOSITORY, check=True, text=True, stdout=subprocess.PIPE,
            stderr=subprocess.DEVNULL).stdout.strip())
        return {"commit": commit, "tracked_files_modified": dirty}
    except (FileNotFoundError, subprocess.CalledProcessError):
        return {"commit": None, "tracked_files_modified": None}


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
    parser.add_argument("--jec-l2", default=str(DEFAULT_JEC_L2),
                        help="L2 correction applied to raw MC and data jets")
    parser.add_argument(
        "--jec-residual", default=str(DEFAULT_JEC_RESIDUAL),
        help="additional residual correction applied to raw data jets")
    parser.add_argument("--job-flavour", default="longlunch",
                        choices=("espresso", "microcentury", "longlunch",
                                 "workday", "tomorrow", "testmatch", "nextweek"))
    parser.add_argument("--retries", type=int, default=2)
    parser.add_argument(
        "--eos-results", default="",
        help=("EOS root:// directory for partial ROOT files. This avoids "
              "returning large outputs to AFS; for example "
              "root://eosuser.cern.ch//eos/user/v/voutila/zjet/CAMPAIGN"))
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
    jec_l2_path, _ = optional_input(args.jec_l2,"L2 JEC file")
    jec_residual_path, _ = optional_input(args.jec_residual,
                                           "data residual JEC file")
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

    eos_results = (normalize_remote_directory(args.eos_results)
                   if args.eos_results else "")
    campaign_dir = REPOSITORY / "condor" / "jobs" / args.campaign
    if campaign_dir.exists():
        raise FileExistsError(
            f"campaign directory already exists: {campaign_dir}; use a new name")
    chunk_dir = campaign_dir / "chunks"
    log_dir = campaign_dir / "logs"
    result_dir = campaign_dir / "results"
    directories = (chunk_dir,log_dir) if eos_results else \
        (chunk_dir,log_dir,result_dir)
    for directory in directories:
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
            result_path = (eos_results + output_file if eos_results else
                           str((result_dir/output_file).resolve()))
            jobs.append({
                "sample": sample,
                "chunk_id": f"{index:04d}",
                "chunk_name": chunk_name,
                "chunk_path": str(chunk_path.relative_to(REPOSITORY)),
                "output_file": output_file,
                "result_path": result_path,
                "input_files": len(group),
            })

    # Keep the submit manifest limited to values referenced by the submit
    # description.  Descriptive values such as chunk_name and the per-job EOS
    # result URL remain available in campaign.json.  HTCondor warns about Queue
    # columns which are not expanded anywhere in the submit description.
    manifest_path = campaign_dir / "jobs.tsv"
    if eos_results:
        manifest_fields = ("sample", "chunk_id", "chunk_path", "output_file")
    else:
        manifest_fields = (
            "sample", "chunk_id", "chunk_path", "output_file", "result_path")
    manifest_path.write_text(
        "".join(
            "\t".join(str(job[field]) for field in manifest_fields) + "\n"
            for job in jobs),
        encoding="utf-8")

    selected_mc = "".join(f"{value}\n" for value in mc_files).encode("utf-8")
    selected_data = "".join(f"{value}\n" for value in data_files).encode("utf-8")
    metadata = {
        "campaign": args.campaign,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "files_per_job": args.files_per_job,
        "mc_files": len(mc_files),
        "data_files": len(data_files),
        "jobs": jobs,
        "storage": {
            "result_mode": "eos" if eos_results else "afs",
            "result_directory": eos_results or str(result_dir.resolve()),
        },
        "source": git_description(),
        "inputs": {
            "mc_list": {"basename": mc_list.name,
                        "selected_sha256": sha256_bytes(selected_mc)},
            "data_list": {"basename": data_list.name,
                          "selected_sha256": sha256_bytes(selected_data)},
            "golden_json": file_description(golden_path),
            "lumi_pileup": file_description(lumi_path),
            "pileup_weights": file_description(weights_path),
            "jec_l2": file_description(jec_l2_path),
            "jec_residual": file_description(jec_residual_path),
        },
        "analysis": {
            "method": ("all accepted Z-jet pairs; +90 and -90 degree "
                       "sidebands, each with weight 0.5"),
            "jet_pt": "recomputed from raw pT with the configured JEC chain",
            "jec_recomputed_from_raw_pt": True,
            "stored_residual_profile": "inverse data residual correction",
            "jer_smearing": {"enabled": False},
            "jet_veto_map": {"enabled": False},
        },
        "command": " ".join(sys.argv),
    }
    (campaign_dir/"campaign.json").write_text(
        json.dumps(metadata,indent=2) + "\n",encoding="utf-8")

    common_inputs = [
        "zjet.C", "zjet.h", "ZJetLumi.h", "mk_compile.C",
        "run_zjet_job.C",
        "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h",
        "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h",
        "CondFormats/JetMETObjects/interface/SimpleJetCorrector.h",
        "CondFormats/JetMETObjects/src/Utilities.cc",
        "CondFormats/JetMETObjects/src/JetCorrectorParameters.cc",
        "CondFormats/JetMETObjects/src/SimpleJetCorrector.cc",
        "CondFormats/JetMETObjects/src/FactorizedJetCorrector.cc",
    ]
    transfer_inputs = list(common_inputs)
    staged_optional_names = []
    for path in optional_paths:
        try:
            relative = path.relative_to(REPOSITORY)
        except ValueError:
            transfer_inputs.append(str(path))
            staged_optional_names.append(path.name)
        else:
            transfer_inputs.append(str(relative))
            staged_optional_names.append(str(relative))
    correction_arguments = []
    for path in (jec_l2_path,jec_residual_path):
        if path is None:
            correction_arguments.append("-")
            continue
        try:
            relative = path.relative_to(REPOSITORY)
        except ValueError:
            transfer_inputs.append(str(path))
            correction_arguments.append(path.name)
        else:
            transfer_inputs.append(str(relative))
            correction_arguments.append(str(relative))
    transfer_inputs.append("$(chunk_path)")
    staged_optional = iter(staged_optional_names)
    golden_argument = next(staged_optional) if golden_path else "-"
    lumi_argument = next(staged_optional) if lumi_path else "-"
    weights_argument = next(staged_optional) if weights_path else "-"
    jec_l2_argument, jec_residual_argument = correction_arguments
    output_directive = (f"output_destination = {eos_results}\n"
                        "MY.XRDCP_CREATE_DIR = True\n"
                        if eos_results else
                        'transfer_output_remaps = "$(output_file)=$(result_path)"\n')
    submit_text = f"""universe = vanilla
executable = condor/run_zjet_job.sh
initialdir = {REPOSITORY}
arguments = $(sample) $(chunk_path) $(output_file) {golden_argument} {lumi_argument} {weights_argument} {jec_l2_argument} {jec_residual_argument}

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
preserve_relative_paths = True
transfer_input_files = {','.join(transfer_inputs)}
transfer_output_files = $(output_file)
{output_directive}

on_exit_remove = (ExitBySignal == False) && (ExitCode == 0)
max_retries = {args.retries}
requirements = (Machine =!= split(LastRemoteHost, \"@\")[1])

queue {','.join(manifest_fields)} from {manifest_path}
"""
    submit_path = campaign_dir / "zjet.sub"
    submit_path.write_text(submit_text,encoding="utf-8")

    mc_jobs = sum(job["sample"] == "mc" for job in jobs)
    data_jobs = sum(job["sample"] == "data" for job in jobs)
    print(f"Prepared campaign {args.campaign}: {len(mc_files)} MC files in "
          f"{mc_jobs} jobs and {len(data_files)} data files in {data_jobs} jobs.")
    print(f"Submit with: condor_submit {submit_path.relative_to(REPOSITORY)}")
    print(f"Results will be written to: "
          f"{eos_results or result_dir.relative_to(REPOSITORY)}")
    print(f"Check readiness with: python3 scripts/status_condor.py "
          f"{args.campaign}")


if __name__ == "__main__":
    main()
