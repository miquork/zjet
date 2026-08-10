#!/usr/bin/env python3
"""Validate and merge all outputs from a prepared Z+jet HTCondor campaign."""

import argparse
import json
import shutil
import subprocess
import tempfile
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Union

from condor_storage import (ensure_remote_directory, is_remote,
                            normalize_remote_directory, remote_basenames,
                            remote_file_sizes)


REPOSITORY = Path(__file__).resolve().parents[1]


def campaign_path(value: str) -> Path:
    direct = Path(value).expanduser()
    if direct.is_dir():
        return direct.resolve()
    return (REPOSITORY/"condor"/"jobs"/value).resolve()


def root_macro_argument(value: Path) -> str:
    return str(value).replace("\\", "\\\\").replace('"', '\\"')


def public_provenance(metadata: Dict[str, object]) -> Dict[str, object]:
    """Drop operational absolute paths and per-file URLs from the archive."""
    return {
        "campaign": metadata["campaign"],
        "created_utc": metadata["created_utc"],
        "merged_utc": datetime.now(timezone.utc).isoformat(),
        "files_per_job": metadata["files_per_job"],
        "mc_files": metadata["mc_files"],
        "data_files": metadata["data_files"],
        "mc_jobs": sum(job["sample"] == "mc" for job in metadata["jobs"]),
        "data_jobs": sum(job["sample"] == "data" for job in metadata["jobs"]),
        "source": metadata.get("source", {}),
        "inputs": metadata.get("inputs", {}),
        "analysis": metadata.get("analysis", {}),
    }


def provenance_log(provenance: Dict[str, object]) -> str:
    inputs = provenance["inputs"]
    analysis = provenance["analysis"]
    source = provenance["source"]
    return "\n".join([
        f"Campaign: {provenance['campaign']}",
        f"Created UTC: {provenance['created_utc']}",
        f"Merged UTC: {provenance['merged_utc']}",
        f"Git commit: {source.get('commit')}",
        f"Tracked files modified at preparation: "
        f"{source.get('tracked_files_modified')}",
        f"MC list: {inputs.get('mc_list')}",
        f"Data list: {inputs.get('data_list')}",
        f"Golden JSON: {inputs.get('golden_json')}",
        f"Lumisection pileup: {inputs.get('lumi_pileup')}",
        f"Pileup weights: {inputs.get('pileup_weights')}",
        f"L2 JEC: {inputs.get('jec_l2')}",
        f"Data residual JEC: {inputs.get('jec_residual')}",
        f"Analysis settings: {analysis}",
        f"Inputs: {provenance['mc_files']} MC files in "
        f"{provenance['mc_jobs']} jobs; {provenance['data_files']} data "
        f"files in {provenance['data_jobs']} jobs",
        "",
    ])


def upload(local: Path, remote: str, force: bool) -> None:
    xrdcp = shutil.which("xrdcp")
    if not xrdcp:
        raise RuntimeError("xrdcp is unavailable; run this command on lxplus")
    command = [xrdcp]
    if force:
        command.append("-f")
    command.extend([str(local), remote])
    subprocess.run(command, check=True)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("campaign", help="campaign name or directory")
    parser.add_argument(
        "--output-dir", default=str(REPOSITORY/"rootfiles"),
        help=("local directory or EOS root:// directory for merged files; "
              "temporary hadd output is kept outside AFS"))
    parser.add_argument("--scratch-dir", type=Path, default=None,
                        help="local temporary directory (default: system TMPDIR)")
    parser.add_argument("--force", action="store_true",
                        help="replace existing zjet_MC.root and zjet_DATA.root")
    args = parser.parse_args()

    campaign_dir = campaign_path(args.campaign)
    metadata_path = campaign_dir/"campaign.json"
    if not metadata_path.is_file():
        raise FileNotFoundError(f"campaign metadata not found: {metadata_path}")
    metadata = json.loads(metadata_path.read_text(encoding="utf-8"))

    result_directory = metadata.get("storage", {}).get("result_directory", "")
    remote_results = is_remote(result_directory)
    remote_sizes = remote_file_sizes(result_directory) if remote_results else {}
    groups: Dict[str, List[str]] = {"mc": [], "data": []}
    missing: List[str] = []
    for job in metadata["jobs"]:
        value = job["result_path"]
        if remote_results:
            valid = remote_sizes.get(job["output_file"],0)>0
        else:
            path = Path(value)
            valid = path.is_file() and path.stat().st_size > 0
        if not valid:
            missing.append(value)
        else:
            groups[job["sample"]].append(value)
    if missing:
        examples = "\n".join(f"  {path}" for path in missing[:10])
        raise RuntimeError(
            f"{len(missing)} expected job outputs are missing or empty:\n{examples}")

    hadd = shutil.which("hadd")
    if not hadd:
        raise RuntimeError("hadd is not available in the current ROOT environment")
    destination_value = args.output_dir
    remote_destination = is_remote(destination_value)
    if remote_destination:
        destination: Union[str, Path] = normalize_remote_directory(destination_value)
        ensure_remote_directory(destination)
        existing_outputs = remote_basenames(destination)
    else:
        destination = Path(destination_value).expanduser().resolve()
        destination.mkdir(parents=True,exist_ok=True)
        existing_outputs = {path.name for path in destination.iterdir()}

    provenance = public_provenance(metadata)
    provenance_name = f"zjet_{metadata['campaign']}_provenance.json"
    campaign_provenance = campaign_dir/provenance_name
    campaign_provenance.write_text(
        json.dumps(provenance,indent=2) + "\n",encoding="utf-8")
    log_name = f"zjet_{metadata['campaign']}_merge.log"
    campaign_log = campaign_dir/log_name
    campaign_log.write_text(provenance_log(provenance),encoding="utf-8")

    scratch_parent = (args.scratch_dir.expanduser().resolve()
                      if args.scratch_dir else None)
    if scratch_parent:
        scratch_parent.mkdir(parents=True,exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="zjet_merge_",
                                     dir=scratch_parent) as temporary_name:
        temporary_dir = Path(temporary_name)
        for sample, output_name in (("mc","zjet_MC.root"),
                                    ("data","zjet_DATA.root")):
            inputs = groups[sample]
            if not inputs:
                continue
            if output_name in existing_outputs and not args.force:
                raise FileExistsError(
                    f"output exists at destination: {output_name}; "
                    "pass --force to replace it")
            temporary = temporary_dir/output_name
            subprocess.run([hadd,"-f",str(temporary),*inputs],check=True)
            macro = (f'embedCampaignMetadata.C('
                     f'"{root_macro_argument(temporary)}",'
                     f'"{root_macro_argument(campaign_provenance)}")')
            subprocess.run(["root","-l","-b","-q",macro],
                           cwd=REPOSITORY,check=True)
            if remote_destination:
                output = destination + output_name
                upload(temporary,output,args.force)
            else:
                output = destination/output_name
                shutil.copy2(temporary,output)
            print(f"Merged {len(inputs)} {sample} outputs into {output}")

    if remote_destination:
        upload(campaign_provenance,destination+provenance_name,args.force)
        upload(campaign_log,destination+log_name,args.force)
        print(f"Wrote provenance to {destination+provenance_name}")
    else:
        output_provenance = destination/provenance_name
        shutil.copy2(campaign_provenance,output_provenance)
        shutil.copy2(campaign_log,destination/log_name)
        print(f"Wrote provenance to {output_provenance}")


if __name__ == "__main__":
    main()
