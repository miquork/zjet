#!/usr/bin/env python3
"""Convert one or more {"files": [...]} JSON documents to a ROOT file list."""

import argparse
import json
from pathlib import Path
from urllib.parse import urlsplit


def logical_file_name(value: str) -> str:
    """Return the /store/... logical file name from an LFN or URL."""
    path = urlsplit(value).path if "://" in value else value
    path = "/" + path.lstrip("/")
    if not path.startswith("/store/"):
        raise ValueError(f"file does not contain a /store LFN: {value}")
    return path


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("json", nargs="+", type=Path,
                        help="input JSON file(s); duplicates are removed")
    parser.add_argument("-o", "--output", required=True, type=Path,
                        help="output text file for TChain/TFileCollection")
    parser.add_argument(
        "--redirector", default="root://cms-xrd-global.cern.ch/",
        help="XRootD redirector (default: CMS global redirector)")
    parser.add_argument("--keep-urls", action="store_true",
                        help="keep the URLs from JSON instead of making XRootD URLs")
    parser.add_argument("--max-files", type=int, default=-1,
                        help="write at most this many unique files")
    args = parser.parse_args()

    files = []
    seen = set()
    for input_path in args.json:
        with input_path.open(encoding="utf-8") as handle:
            payload = json.load(handle)
        values = payload.get("files") if isinstance(payload, dict) else None
        if not isinstance(values, list) or not all(isinstance(x, str) for x in values):
            raise ValueError(f"{input_path}: expected a string array named 'files'")
        for value in values:
            if args.max_files >= 0 and len(files) >= args.max_files:
                break
            output_value = value
            if not args.keep_urls:
                output_value = args.redirector.rstrip("/") + "/" + logical_file_name(value)
            if output_value in seen:
                continue
            seen.add(output_value)
            files.append(output_value)
        if args.max_files >= 0 and len(files) >= args.max_files:
            break

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text("".join(f"{value}\n" for value in files), encoding="utf-8")
    print(f"Wrote {len(files)} unique file(s) to {args.output}")


if __name__ == "__main__":
    main()
