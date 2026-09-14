#!/usr/bin/env python3
"""Prepare (never submit) one missing/empty job after the UUID-audit repair.

Keep the original physics configuration, input chunk and result/log destinations.
Only ZJetInputCounters.h may differ among recorded worker sources. Original logs
are backed up, and hashes record the limited mixed-version campaign repair.
"""
import argparse
import hashlib
import json
import re
import shlex
import shutil
import subprocess
from datetime import datetime, timezone
from pathlib import Path

from condor_storage import is_remote, remote_file_sizes
from status_condor import campaign_path


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def build_submit(original, row, retry_dir):
    matches = list(re.finditer(r"(?mi)^queue\s+([^\n]+)$", original))
    if len(matches) != 1:
        raise ValueError("Expected exactly one generated queue statement")
    match = matches[0]
    fields = match.group(1).split(" from ")[0].strip().split(",")
    if fields not in (["sample", "chunk_id", "chunk_path", "output_file"],
                      ["sample", "chunk_id", "chunk_path", "output_file", "result_path"]):
        raise ValueError("Unrecognized queue schema; refusing to guess")
    values = [str(row[field]) for field in fields]
    if any(re.search(r"[\s\"()]", value) for value in values):
        raise ValueError("Queue values contain unsafe submit-file characters")
    queue = "queue " + ",".join(fields) + " from (\n" + "\t".join(values) + "\n)"
    updated = original[:match.start()] + queue + original[match.end():]
    # Separate event log allows condor_wait to wait specifically for this retry.
    updated, n = re.subn(r"(?mi)^log\s*=.*$", "log = " + str(retry_dir/"condor.log"), updated)
    if n != 1:
        raise ValueError("Expected exactly one event log setting")
    return updated


def setting(text, name):
    matches = re.findall(r"(?mi)^" + re.escape(name) + r"\s*=\s*(.+)$", text)
    if len(matches) != 1:
        raise ValueError("Missing or ambiguous submit setting: " + name)
    return matches[0].strip()


def prepare(campaign, job_name, cluster, label):
    directory = campaign_path(campaign)
    metadata = json.loads((directory/"campaign.json").read_text())
    jobs = [j for j in metadata["jobs"] if j["sample"]+"_"+j["chunk_id"] == job_name]
    if len(jobs) != 1:
        raise ValueError("Expected exactly one matching job")
    job = jobs[0]
    original = (directory/"zjet.sub").read_text()
    repo = Path(setting(original, "initialdir")).resolve()
    if repo != Path(__file__).resolve().parents[1]:
        raise ValueError("Run this from the original campaign checkout")
    if re.search(r"[\s\"()]", str(directory)):
        raise ValueError("Unsupported campaign path characters")
    active = subprocess.check_output(["condor_q", str(cluster), "-autoformat", "ClusterId", "ProcId"], text=True)
    if active.strip():
        raise ValueError("Original cluster still has queued jobs; do not race a retry")
    result = job["result_path"]
    if is_remote(result):
        size = remote_file_sizes(metadata["storage"]["result_directory"]).get(job["output_file"], 0)
    else:
        path = Path(result)
        size = path.stat().st_size if path.exists() else 0
    if size:
        raise ValueError("Refusing to overwrite a nonempty result: " + job["output_file"])
    proxy = Path(setting(original, "x509userproxy"))
    if not str(proxy).startswith("/afs/") or not proxy.is_file() or proxy.stat().st_mode & 0o077:
        raise ValueError("Original submit proxy must be a readable, protected AFS file")
    seconds = int(subprocess.check_output(["voms-proxy-info", "-file", str(proxy), "-timeleft"], text=True).strip())
    fqan = subprocess.check_output(["voms-proxy-info", "-file", str(proxy), "-fqan"], text=True)
    if seconds < 7200 or "/cms/" not in fqan:
        raise ValueError("Renew the CMS proxy at the original submit path (need at least 2 hours)")
    transferred = setting(original, "transfer_input_files").split(",")
    changed = {}
    hashes = {}
    for name, description in metadata["source_files"].items():
        if name not in transferred and name != setting(original, "executable"):
            continue
        value = digest(repo/name)
        hashes[name] = value
        if value != description.get("sha256"):
            changed[name] = {"original": description.get("sha256"), "retry": value}
    if set(changed) - {"ZJetInputCounters.h"}:
        raise ValueError("Other worker source changes would mix physics versions: " + ", ".join(changed))
    # Also check all calibration/JSON inputs for which the campaign saved hashes.
    for description in metadata.get("inputs", {}).values():
        if not isinstance(description, dict) or not description.get("sha256"):
            continue
        candidates = [repo/value for value in transferred
                      if Path(value).name == description.get("basename")]
        if len(candidates) != 1 or digest(candidates[0]) != description["sha256"]:
            raise ValueError("Calibration input changed or missing: " + str(description.get("basename")))
    chunk = repo/job["chunk_path"]
    lines = [line for line in chunk.read_text().splitlines() if line.strip() and not line.lstrip().startswith("#")]
    if len(lines) != job["input_files"]:
        raise ValueError("Original input chunk length changed")
    sample_lines = []
    for entry in metadata["jobs"]:
        if entry["sample"] == job["sample"]:
            sample_lines.extend(line.strip() for line in (repo/entry["chunk_path"]).read_text().splitlines()
                                if line.strip() and not line.lstrip().startswith("#"))
    expected = metadata["inputs"][job["sample"]+"_list"]["selected_sha256"]
    if hashlib.sha256(("\n".join(sample_lines)+"\n").encode()).hexdigest() != expected:
        raise ValueError("Campaign input chunks no longer match the original selected-list hash")
    retry = directory/"retries"/(job_name+"_"+label)
    submission = build_submit(original, job, retry)
    # Fail if this label was used before; never overwrite a previous retry.
    retry.mkdir(parents=True, exist_ok=False)
    for suffix in ("out", "err"):
        log = Path(metadata["storage"]["log_directory"])/(job_name+"."+suffix)
        if log.is_file():
            shutil.copy2(log, retry/("original."+suffix))
    record = {"kind": "prepared_uuid_retry", "job": job_name,
              "prepared_utc": datetime.now(timezone.utc).isoformat(),
              "original_cluster": cluster, "changed_worker_sources": changed,
              "worker_source_sha256": hashes, "chunk_sha256": digest(chunk),
              "original_submit_sha256": digest(directory/"zjet.sub"),
              "note": "Prepared only; worker .out log records actual execution hashes."}
    (retry/"repair.json").write_text(json.dumps(record, indent=2)+"\n")
    (retry/"zjet.sub").write_text(submission)
    print("Prepared exactly one job: " + job_name)
    print("Original physics settings/results retained; original worker logs backed up.")
    print("Proxy remaining: " + str(seconds) + " seconds")
    print("Submit explicitly: condor_submit -maxjobs 1 " + shlex.quote(str(retry/"zjet.sub")))
    print("Then wait: condor_wait " + shlex.quote(str(retry/"condor.log")))
    return retry


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("campaign")
    parser.add_argument("job")
    parser.add_argument("--cluster", type=int, required=True, help="original completed cluster")
    parser.add_argument("--label", default="uuidfix")
    args = parser.parse_args()
    if not re.fullmatch(r"[A-Za-z0-9_-]+", args.label):
        parser.error("Invalid retry label")
    try:
        prepare(args.campaign, args.job, args.cluster, args.label)
    except (OSError, ValueError, RuntimeError, subprocess.CalledProcessError) as error:
        parser.exit(1, "ERROR: " + str(error) + "\nNo jobs submitted.\n")


if __name__ == "__main__":
    main()
