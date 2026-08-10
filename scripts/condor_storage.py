#!/usr/bin/env python3
"""Shared storage helpers for Z+jet HTCondor campaign scripts."""

import shutil
import subprocess
import re
from pathlib import Path
from typing import Dict, List, Set, Tuple
from urllib.parse import urlsplit


def is_remote(value: str) -> bool:
    return value.startswith("root://")


def normalize_remote_directory(value: str) -> str:
    if not is_remote(value):
        raise ValueError("EOS destination must be a root:// URL")
    parts = urlsplit(value)
    if not parts.netloc or not parts.path.startswith("//eos/"):
        raise ValueError(
            "EOS destination must look like "
            "root://eosuser.cern.ch//eos/user/<initial>/<user>/...")
    return value.rstrip("/") + "/"


def split_remote(value: str) -> Tuple[str, str]:
    parts = urlsplit(value)
    endpoint = f"{parts.scheme}://{parts.netloc}"
    path = parts.path[1:] if parts.path.startswith("//") else parts.path
    return endpoint, path


def remote_basenames(directory: str) -> Set[str]:
    """List one CERN EOS directory with Kerberos-capable XRootD tools."""
    normalized = normalize_remote_directory(directory)
    endpoint, path = split_remote(normalized)
    if shutil.which("xrdfs"):
        command = ["xrdfs", endpoint, "ls", path]
    else:
        command = ["gfal-ls", normalized]
    try:
        result = subprocess.run(command, check=True, text=True,
                                stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE)
    except FileNotFoundError as error:
        raise RuntimeError(
            "neither xrdfs nor gfal-ls is available; run this on lxplus") from error
    except subprocess.CalledProcessError as error:
        detail = error.stderr.strip() or error.stdout.strip()
        raise RuntimeError(
            "could not list EOS result directory; check klist and the URL: "
            f"{detail}") from error
    return {Path(line.strip()).name for line in result.stdout.splitlines()
            if line.strip()}


def remote_file_sizes(directory: str) -> Dict[str, int]:
    """Return file sizes for one EOS directory without trusting names alone."""
    normalized = normalize_remote_directory(directory)
    endpoint, path = split_remote(normalized)

    # lxplus exposes EOS through FUSE. This is both faster and less ambiguous
    # than parsing the human-readable long-listing formats of different client
    # versions. Condor itself still uses the root:// transfer plug-in.
    local_directory = Path(path)
    if local_directory.is_dir():
        return {item.name: item.stat().st_size
                for item in local_directory.iterdir() if item.is_file()}

    names = remote_basenames(normalized)
    sizes: Dict[str, int] = {}
    for name in names:
        remote_path = path.rstrip("/") + "/" + name
        if shutil.which("xrdfs"):
            command = ["xrdfs",endpoint,"stat",remote_path]
        elif shutil.which("gfal-stat"):
            command = ["gfal-stat",normalized+name]
        else:
            raise RuntimeError(
                "neither xrdfs nor gfal-stat is available; run this on lxplus")
        try:
            result = subprocess.run(command,check=True,text=True,
                                    stdout=subprocess.PIPE,
                                    stderr=subprocess.PIPE)
        except subprocess.CalledProcessError as error:
            detail = error.stderr.strip() or error.stdout.strip()
            raise RuntimeError(
                f"could not stat EOS result {name}: {detail}") from error
        match = re.search(r"(?m)^\s*Size:\s*(\d+)\s*$",result.stdout)
        if not match:
            raise RuntimeError(
                f"could not parse file size for {name}: {result.stdout.strip()}")
        sizes[name] = int(match.group(1))
    return sizes


def ensure_remote_directory(directory: str) -> None:
    normalized = normalize_remote_directory(directory)
    endpoint, path = split_remote(normalized)
    if shutil.which("xrdfs"):
        subprocess.run(["xrdfs",endpoint,"mkdir","-p",path],check=True)
    else:
        subprocess.run(["gfal-mkdir","-p",normalized],check=True)


def remove_remote_files(urls: List[str]) -> None:
    for url in urls:
        if shutil.which("xrdfs"):
            endpoint, path = split_remote(url)
            subprocess.run(["xrdfs", endpoint, "rm", path], check=True)
        else:
            subprocess.run(["gfal-rm", url], check=True)
