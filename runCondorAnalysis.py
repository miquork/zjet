#!/usr/bin/env python3
"""Run a checkpointed and explicitly approved Z+jet HTCondor workflow."""

import argparse
import getpass
import json
import math
import os
import re
import shlex
import shutil
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Optional


REPOSITORY = Path(__file__).resolve().parent
WORKFLOW_DIR = REPOSITORY / "condor" / "jobs" / "_workflows"
STAGES = [
    "created", "preflight_complete", "smoke_prepared", "smoke_submitted",
    "smoke_ready", "full_prepared", "full_submitted", "full_ready",
    "merged", "compatibility_written",
]
PRESETS = {
    "run2024i": {
        "label": "Run 2024I data and Summer24 DY MC",
        "mc_list": "textfiles/generated/summer24_mc.txt",
        "data_list": "textfiles/generated/run2024i_data.txt",
        "golden_json": "data/Cert_Collisions2024_378981_386951_Golden.json",
        "jec_l2": (
            "CondFormats/JetMETObjects/data/"
            "RunIII2024Summer24_V2_MC_L2Relative_AK4PUPPI.txt"
        ),
        "jec_residual": (
            "CondFormats/JetMETObjects/data/"
            "Prompt24_Run2024I_nib1_V11M_DATA_"
            "L2L3Residual_AK4PFPuppi.txt"
        ),
        # Deliberately conservative, order-of-magnitude planning estimates.
        "mc_cpu_minutes_per_file": 2.0,
        "data_cpu_minutes_per_file": 1.0,
    },
}


def run(command: List[str], *, capture: bool = False, check: bool = True,
        environment: Optional[Dict[str, str]] = None
        ) -> subprocess.CompletedProcess:
    print("+ " + shlex.join(command), flush=True)
    return subprocess.run(
        command, cwd=REPOSITORY,
        env=environment if environment is not None else os.environ.copy(),
        check=check,
        text=True, stdout=subprocess.PIPE if capture else None,
        stderr=subprocess.STDOUT if capture else None,
    )


def confirm(question: str, default: bool = False) -> bool:
    suffix = " [Y/n] " if default else " [y/N] "
    while True:
        answer = input(question + suffix).strip().lower()
        if not answer:
            return default
        if answer in ("y", "yes"):
            return True
        if answer in ("n", "no"):
            return False
        print("Please answer y or n.")


def require_commands(names: List[str]) -> None:
    missing = [name for name in names if shutil.which(name) is None]
    if missing:
        raise RuntimeError("missing command(s): " + ", ".join(missing))


def nonempty_lines(path: Path) -> int:
    if not path.is_file():
        raise FileNotFoundError(f"required input is missing: {path}")
    count = sum(bool(line.strip() and not line.lstrip().startswith("#"))
                for line in path.read_text(encoding="utf-8").splitlines())
    if count == 0:
        raise ValueError(f"input list is empty: {path}")
    return count


def resolve_preset(name: Optional[str]) -> str:
    if name:
        return name
    if not sys.stdin.isatty():
        raise RuntimeError("pass --preset when stdin is not interactive")
    names = list(PRESETS)
    print("Supported analyses:")
    for index, key in enumerate(names, 1):
        print(f"  {index}. {key}: {PRESETS[key]['label']}")
    while True:
        answer = input("Select an analysis: ").strip()
        if answer.isdigit() and 1 <= int(answer) <= len(names):
            return names[int(answer) - 1]
        if answer in PRESETS:
            return answer
        print("Select a listed number or preset name.")


def default_storage(campaign: str) -> Dict[str, str]:
    user = getpass.getuser()
    initial = user[0]
    return {
        "log_root": f"/afs/cern.ch/work/{initial}/{user}/zjet-condor",
        "eos_root": (
            "root://eosuser.cern.ch//eos/user/"
            f"{initial}/{user}/zjet"
        ),
        "smoke_campaign": campaign + "_smoke",
        "full_campaign": campaign + "_full",
    }


def create_state(args: argparse.Namespace, preset_name: str) -> Dict[str, object]:
    campaign = args.campaign or (
        preset_name + "_" + datetime.now(timezone.utc).strftime("%Y%m%d_%H%M")
    )
    storage = default_storage(campaign)
    if args.log_root:
        storage["log_root"] = str(args.log_root.expanduser().resolve())
    if args.eos_root:
        storage["eos_root"] = args.eos_root.rstrip("/")
    return {
        "workflow": campaign,
        "preset": preset_name,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "stage": "created",
        "files_per_job": args.files_per_job,
        "poll_seconds": args.poll_seconds,
        "compatibility_output": str(args.compatibility_output),
        "flavor_placeholders": args.flavor_placeholders,
        **storage,
        "clusters": {},
    }


def state_path(name: str) -> Path:
    direct = Path(name).expanduser()
    if direct.is_file():
        return direct.resolve()
    return WORKFLOW_DIR / f"{name}.json"


def save_state(path: Path, state: Dict[str, object]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(".json.tmp")
    temporary.write_text(json.dumps(state, indent=2) + "\n", encoding="utf-8")
    temporary.replace(path)


def advance(path: Path, state: Dict[str, object], stage: str) -> None:
    state["stage"] = stage
    state["updated_utc"] = datetime.now(timezone.utc).isoformat()
    save_state(path, state)
    print(f"Checkpoint saved: {stage} ({path})")


def command_output(command: List[str]) -> str:
    result = run(command, capture=True)
    output = result.stdout.strip()
    if output:
        print(output)
    return output


def valid_proxy(path: Path, minimum_seconds: int) -> bool:
    if not path.is_file():
        return False
    lifetime = subprocess.run(
        ["voms-proxy-info", "-file", str(path), "-timeleft"],
        text=True, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL,
    )
    if lifetime.returncode != 0:
        return False
    try:
        seconds = int(lifetime.stdout.strip())
    except ValueError:
        return False
    fqan = subprocess.run(
        ["voms-proxy-info", "-file", str(path), "-fqan"],
        text=True, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL,
    )
    return (fqan.returncode == 0 and "/cms/" in fqan.stdout and
            seconds >= minimum_seconds)


def ensure_proxy(minimum_hours: float) -> Path:
    minimum_seconds = int(minimum_hours * 3600)
    candidates = []
    configured = os.environ.get("X509_USER_PROXY")
    if configured:
        candidates.append(Path(configured).expanduser())
    shared = Path.home() / "private" / f"x509up_u{os.getuid()}"
    candidates.append(shared)
    default = subprocess.run(
        ["voms-proxy-info", "-path"], text=True, stdout=subprocess.PIPE,
        stderr=subprocess.DEVNULL,
    )
    if default.returncode == 0 and default.stdout.strip():
        candidates.append(Path(default.stdout.strip()))
    proxy = next((path.resolve() for path in candidates
                  if valid_proxy(path, minimum_seconds)), None)
    if proxy is None:
        raise RuntimeError(
            f"no CMS VOMS proxy with at least {minimum_hours:g} h remaining; "
            "run: voms-proxy-init --rfc --voms cms --valid 192:00"
        )
    if not str(proxy).startswith("/afs/"):
        if not confirm(f"Copy the valid proxy {proxy} to protected AFS storage?",
                       default=True):
            raise RuntimeError("HTCondor requires an AFS-visible proxy")
        shared.parent.mkdir(mode=0o700, parents=True, exist_ok=True)
        shared.parent.chmod(0o700)
        if shutil.which("fs"):
            run(["fs", "setacl", "-dir", str(shared.parent),
                 "-acl", "system:anyuser", "none"])
            run(["fs", "setacl", "-dir", str(shared.parent),
                 "-acl", getpass.getuser(), "all"])
        shutil.copyfile(proxy, shared)
        shared.chmod(0o600)
        proxy = shared.resolve()
    if not valid_proxy(proxy, minimum_seconds):
        raise RuntimeError(f"copied proxy failed validation: {proxy}")
    os.environ["X509_USER_PROXY"] = str(proxy)
    lifetime = command_output(
        ["voms-proxy-info", "-file", str(proxy), "-timeleft"])
    print(f"Using AFS proxy {proxy}; lifetime {int(lifetime)/3600:.1f} h")
    return proxy


def report_space(path: Path) -> None:
    existing = path
    while not existing.exists() and existing != existing.parent:
        existing = existing.parent
    print(f"Storage check for {path} (nearest existing path: {existing})")
    if shutil.which("fs") and str(existing).startswith("/afs/"):
        run(["fs", "listquota", "-path", str(existing)], check=False)
    run(["df", "-h", str(existing)], check=False)


def local_compiler_environment() -> Dict[str, str]:
    """Keep preflight compiler temporaries away from the AFS home quota."""
    cache_root = REPOSITORY / ".cache" / "condor-preflight"
    ccache_temporary = cache_root / "ccache" / "tmp"
    temporary = cache_root / "tmp"
    ccache_temporary.mkdir(parents=True, exist_ok=True)
    temporary.mkdir(parents=True, exist_ok=True)
    environment = os.environ.copy()
    environment.update({
        "CCACHE_DISABLE": "1",
        "CCACHE_DIR": str(cache_root / "ccache"),
        "CCACHE_TEMPDIR": str(ccache_temporary),
        "XDG_CACHE_HOME": str(cache_root),
        "TMPDIR": str(temporary),
    })
    return environment


def preflight(state: Dict[str, object], preset: Dict[str, object],
              skip_pull: bool) -> None:
    require_commands([
        "git", "root", "condor_submit", "condor_q", "klist",
        "voms-proxy-info", "gfal-ls", "hadd", "xrdcp", "xrdfs",
    ])
    run(["klist", "-s"])
    tracked = command_output(
        ["git", "status", "--porcelain", "--untracked-files=no"])
    if tracked:
        raise RuntimeError("tracked files are modified; commit or stash before running")
    head_before_pull = command_output(["git", "rev-parse", "HEAD"])
    if not skip_pull and confirm("Run git pull --ff-only before preparing jobs?",
                                 default=True):
        run(["git", "pull", "--ff-only"])
        head_after_pull = command_output(["git", "rev-parse", "HEAD"])
        if head_after_pull != head_before_pull:
            print("The pull updated the workflow code. Restarting is required "
                  "so this Python process also uses the new version.")
            print(f"Resume with: {resume_command(state)}")
            raise SystemExit(0)
    report_space(Path(str(state["log_root"])))
    report_space(REPOSITORY)
    for key in ("mc_list", "data_list", "golden_json", "jec_l2",
                "jec_residual"):
        path = REPOSITORY / str(preset[key])
        if not path.is_file():
            raise FileNotFoundError(f"preset input is missing: {path}")
    run(["root", "-l", "-b", "-q", "mk_compile.C"],
        environment=local_compiler_environment())
    for key in ("mc_list", "data_list"):
        path = REPOSITORY / str(preset[key])
        first = next(line.strip() for line in path.read_text(
            encoding="utf-8").splitlines() if line.strip() and
            not line.lstrip().startswith("#"))
        run(["gfal-ls", first])


def prepare_command(state: Dict[str, object], preset: Dict[str, object],
                    smoke: bool) -> List[str]:
    campaign_key = "smoke_campaign" if smoke else "full_campaign"
    campaign = str(state[campaign_key])
    command = [
        sys.executable, "scripts/prepare_condor.py",
        "--mc-list", str(preset["mc_list"]),
        "--data-list", str(preset["data_list"]),
        "--campaign", campaign,
        "--files-per-job", "1" if smoke else str(state["files_per_job"]),
        "--job-flavour", "espresso" if smoke else "workday",
        "--golden-json", str(preset["golden_json"]),
        "--jec-l2", str(preset["jec_l2"]),
        "--jec-residual", str(preset["jec_residual"]),
        "--log-dir", f"{state['log_root']}/{campaign}",
        "--eos-results", f"{state['eos_root']}/{campaign}",
    ]
    if smoke:
        command.extend(["--max-mc-files", "1", "--max-data-files", "1"])
    return command


def prepare(state_path_value: Path, state: Dict[str, object],
            preset: Dict[str, object], smoke: bool) -> None:
    campaign_key = "smoke_campaign" if smoke else "full_campaign"
    campaign = str(state[campaign_key])
    metadata = REPOSITORY / "condor" / "jobs" / campaign / "campaign.json"
    if metadata.is_file():
        print(f"Reusing prepared campaign {campaign}")
    else:
        run(prepare_command(state, preset, smoke))
    advance(state_path_value, state,
            "smoke_prepared" if smoke else "full_prepared")


def submit(state_path_value: Path, state: Dict[str, object],
           smoke: bool) -> None:
    campaign_key = "smoke_campaign" if smoke else "full_campaign"
    campaign = str(state[campaign_key])
    kind = "two-job smoke test" if smoke else "full analysis campaign"
    if not confirm(f"Submit the {kind} {campaign} to HTCondor?", default=False):
        print(f"Stopped before submission. Resume with: {resume_command(state)}")
        raise SystemExit(0)
    submit_file = REPOSITORY / "condor" / "jobs" / campaign / "zjet.sub"
    output = command_output(["condor_submit", str(submit_file)])
    match = re.search(r"submitted to cluster\s+(\d+)", output)
    if not match:
        raise RuntimeError("could not parse the HTCondor cluster ID")
    state["clusters"]["smoke" if smoke else "full"] = int(match.group(1))
    advance(state_path_value, state,
            "smoke_submitted" if smoke else "full_submitted")


def queue_states(cluster: int) -> List[List[str]]:
    result = run(
        ["condor_q", str(cluster), "-autoformat", "ClusterId", "ProcId",
         "JobStatus", "HoldReason"], capture=True, check=False,
    )
    if result.returncode != 0:
        raise RuntimeError(result.stdout.strip() or "condor_q failed")
    return [line.split(maxsplit=3) for line in result.stdout.splitlines()
            if line.strip()]


def wait_for_campaign(state_path_value: Path, state: Dict[str, object],
                      smoke: bool) -> None:
    key = "smoke" if smoke else "full"
    campaign = str(state["smoke_campaign" if smoke else "full_campaign"])
    cluster = int(state["clusters"][key])
    if not confirm(f"Wait for and validate cluster {cluster} now?", default=True):
        print(f"The jobs keep running. Resume with: {resume_command(state)}")
        raise SystemExit(0)
    poll_seconds = max(15, int(state["poll_seconds"]))
    while True:
        jobs = queue_states(cluster)
        held = [job for job in jobs if len(job) >= 3 and job[2] == "5"]
        if held:
            print(f"Cluster {cluster} has {len(held)} held job(s):")
            for job in held[:10]:
                print("  " + " ".join(job))
            print(f"Inspect with: condor_q {cluster} -hold")
            print(f"Resume later with: {resume_command(state)}")
            raise SystemExit(2)
        if not jobs:
            break
        counts: Dict[str, int] = {}
        for job in jobs:
            status = job[2] if len(job) >= 3 else "unknown"
            counts[status] = counts.get(status, 0) + 1
        timestamp = datetime.now().astimezone().strftime("%Y-%m-%d %H:%M:%S")
        print(f"{timestamp}: cluster {cluster} still active; states {counts}")
        time.sleep(poll_seconds)
    run([sys.executable, "scripts/status_condor.py", campaign])
    advance(state_path_value, state, "smoke_ready" if smoke else "full_ready")


def resource_summary(state: Dict[str, object], preset: Dict[str, object]) -> None:
    mc_files = nonempty_lines(REPOSITORY / str(preset["mc_list"]))
    data_files = nonempty_lines(REPOSITORY / str(preset["data_list"]))
    chunk = int(state["files_per_job"])
    jobs = math.ceil(mc_files / chunk) + math.ceil(data_files / chunk)
    cpu_minutes = (mc_files * float(preset["mc_cpu_minutes_per_file"]) +
                   data_files * float(preset["data_cpu_minutes_per_file"]))
    print("Full-campaign planning estimate:")
    print(f"  {mc_files} MC files and {data_files} data files")
    print(f"  {jobs} HTCondor jobs at {chunk} files/job")
    print(f"  roughly {cpu_minutes/60:.1f} CPU hours in total")
    print("  expected wall time is tens of minutes with normal CERN slot availability")
    print("These are empirical order-of-magnitude estimates, not guarantees.")


def merge(state_path_value: Path, state: Dict[str, object]) -> None:
    campaign = str(state["full_campaign"])
    destination = f"{state['eos_root']}/{campaign}/merged"
    if not confirm(f"Merge validated outputs into {destination}?", default=False):
        print(f"Stopped before merge. Resume with: {resume_command(state)}")
        raise SystemExit(0)
    run([sys.executable, "scripts/merge_condor.py", campaign,
         "--output-dir", destination])
    state["merged_directory"] = destination
    advance(state_path_value, state, "merged")


def write_compatibility(state_path_value: Path,
                        state: Dict[str, object]) -> None:
    output = Path(str(state["compatibility_output"]))
    if not output.is_absolute():
        output = REPOSITORY / output
    placeholders = bool(state["flavor_placeholders"])
    warning = " with EMPTY Z+flavor placeholders" if placeholders else ""
    if not confirm(f"Write or replace {output}{warning}?", default=False):
        print(f"Stopped before compatibility output. Resume with: {resume_command(state)}")
        raise SystemExit(0)
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary_output = output.with_name(output.name + ".part")
    if temporary_output.exists():
        temporary_output.unlink()
    merged = str(state["merged_directory"]).rstrip("/")
    macro = (
        f'writeJecsys3.C("{merged}/zjet_DATA.root",'
        f'"{merged}/zjet_MC.root","{temporary_output}",'
        f'{str(placeholders).lower()})'
    )
    try:
        run(["root", "-l", "-b", "-q", macro])
        if not temporary_output.is_file() or temporary_output.stat().st_size == 0:
            raise RuntimeError(
                f"compatibility writer did not create {temporary_output}")
        temporary_output.replace(output)
    finally:
        if temporary_output.exists():
            temporary_output.unlink()
    advance(state_path_value, state, "compatibility_written")


def resume_command(state: Dict[str, object]) -> str:
    return f"python3 runCondorAnalysis.py --resume {state['workflow']}"


def print_plan(state: Dict[str, object], preset: Dict[str, object]) -> None:
    print(f"Workflow: {state['workflow']}")
    print(f"Preset: {state['preset']} ({preset['label']})")
    print(f"Smoke campaign: {state['smoke_campaign']}")
    print(f"Full campaign: {state['full_campaign']}")
    print(f"Logs: {state['log_root']}")
    print(f"EOS: {state['eos_root']}")
    resource_summary(state, preset)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--preset", choices=sorted(PRESETS))
    parser.add_argument("--list-presets", action="store_true")
    parser.add_argument("--campaign", help="immutable workflow/campaign prefix")
    parser.add_argument("--resume", help="workflow name or state JSON path")
    parser.add_argument("--plan", action="store_true",
                        help="print the plan without checks, writes or submissions")
    parser.add_argument("--files-per-job", type=int, default=10)
    parser.add_argument("--poll-seconds", type=int, default=60)
    parser.add_argument("--log-root", type=Path)
    parser.add_argument("--eos-root")
    parser.add_argument("--compatibility-output", type=Path,
                        default=Path("rootfiles/zjet_JMENANO_compat.root"))
    parser.add_argument("--flavor-placeholders", action="store_true",
                        help="write explicitly marked empty Z+flavor objects")
    parser.add_argument("--skip-pull", action="store_true")
    args = parser.parse_args()

    if args.list_presets:
        for name, preset in PRESETS.items():
            print(f"{name}: {preset['label']}")
        return
    if args.files_per_job <= 0:
        parser.error("--files-per-job must be positive")
    if args.poll_seconds <= 0:
        parser.error("--poll-seconds must be positive")

    if args.resume:
        path = state_path(args.resume)
        if not path.is_file():
            raise FileNotFoundError(f"workflow state not found: {path}")
        state = json.loads(path.read_text(encoding="utf-8"))
        preset_name = str(state["preset"])
    else:
        preset_name = resolve_preset(args.preset)
        state = create_state(args, preset_name)
        path = state_path(str(state["workflow"]))
        if path.exists():
            raise FileExistsError(f"workflow already exists; use --resume: {path}")
    preset = PRESETS[preset_name]

    print_plan(state, preset)
    if args.plan:
        return
    if not sys.stdin.isatty():
        raise RuntimeError("workflow execution requires an interactive terminal")
    if not args.resume:
        save_state(path, state)

    try:
        require_commands(["voms-proxy-info"])
        ensure_proxy(24.0)
        while state["stage"] != "compatibility_written":
            stage = str(state["stage"])
            if stage == "created":
                preflight(state, preset, args.skip_pull)
                advance(path, state, "preflight_complete")
            elif stage == "preflight_complete":
                prepare(path, state, preset, True)
            elif stage == "smoke_prepared":
                submit(path, state, True)
            elif stage == "smoke_submitted":
                wait_for_campaign(path, state, True)
            elif stage == "smoke_ready":
                if not confirm("Prepare the full campaign now?", default=True):
                    print(f"Resume with: {resume_command(state)}")
                    return
                prepare(path, state, preset, False)
            elif stage == "full_prepared":
                resource_summary(state, preset)
                submit(path, state, False)
            elif stage == "full_submitted":
                wait_for_campaign(path, state, False)
            elif stage == "full_ready":
                merge(path, state)
            elif stage == "merged":
                write_compatibility(path, state)
            else:
                raise RuntimeError(f"unknown workflow stage: {stage}")
    except KeyboardInterrupt:
        print(f"\nInterrupted safely at checkpoint {state['stage']}.")
        print(f"Resume with: {resume_command(state)}")
        raise SystemExit(130)
    except (subprocess.CalledProcessError, RuntimeError, FileNotFoundError,
            ValueError) as error:
        print(f"\nERROR at checkpoint {state['stage']}: {error}",
              file=sys.stderr)
        print("No checkpoint was advanced; existing successful outputs were "
              "left in place.",file=sys.stderr)
        print(f"Resume with: {resume_command(state)}",file=sys.stderr)
        raise SystemExit(1)

    print("Workflow complete.")
    print(f"Compatibility file: {state['compatibility_output']}")


if __name__ == "__main__":
    main()
