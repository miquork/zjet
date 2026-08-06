#!/usr/bin/env python3
"""Validate and merge all outputs from a prepared Z+jet HTCondor campaign."""

import argparse
import json
import os
import shutil
import subprocess
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
    parser.add_argument("--output-dir", type=Path,
                        default=REPOSITORY/"rootfiles")
    parser.add_argument("--force", action="store_true",
                        help="replace existing zjet_MC.root and zjet_DATA.root")
    args = parser.parse_args()

    campaign_dir = campaign_path(args.campaign)
    metadata_path = campaign_dir/"campaign.json"
    if not metadata_path.is_file():
        raise FileNotFoundError(f"campaign metadata not found: {metadata_path}")
    metadata = json.loads(metadata_path.read_text(encoding="utf-8"))

    groups: Dict[str, List[Path]] = {"mc": [], "data": []}
    missing: List[Path] = []
    for job in metadata["jobs"]:
        path = Path(job["result_path"])
        if not path.is_file() or path.stat().st_size == 0:
            missing.append(path)
        else:
            groups[job["sample"]].append(path)
    if missing:
        examples = "\n".join(f"  {path}" for path in missing[:10])
        raise RuntimeError(
            f"{len(missing)} expected job outputs are missing or empty:\n{examples}")

    hadd = shutil.which("hadd")
    if not hadd:
        raise RuntimeError("hadd is not available in the current ROOT environment")
    output_dir = args.output_dir.expanduser().resolve()
    output_dir.mkdir(parents=True,exist_ok=True)

    for sample, output_name in (("mc","zjet_MC.root"),
                                ("data","zjet_DATA.root")):
        inputs = groups[sample]
        if not inputs:
            continue
        output = output_dir/output_name
        if output.exists() and not args.force:
            raise FileExistsError(f"output exists: {output}; pass --force to replace it")
        temporary = output.with_name(f".{output.name}.tmp.{os.getpid()}")
        try:
            subprocess.run([hadd,"-f",str(temporary),
                            *(str(path) for path in inputs)],check=True)
            temporary.replace(output)
        finally:
            if temporary.exists():
                temporary.unlink()
        print(f"Merged {len(inputs)} {sample} outputs into {output}")


if __name__ == "__main__":
    main()
