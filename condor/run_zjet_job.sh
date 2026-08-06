#!/usr/bin/env bash
set -euo pipefail

sample="${1:?missing sample argument}"
input_list="${2:?missing input-list argument}"
output_file="${3:?missing output-file argument}"
golden_json="${4:-}"
lumi_pileup="${5:-}"
pileup_weights="${6:-}"

preserve_failed_status() {
  status=$?
  trap - EXIT
  if [[ ${status} -ne 0 ]]; then
    # A declared transfer output must exist even when the payload fails.
    # The empty marker is overwritten by a successful retry, and the merge
    # tool rejects it if all retries are exhausted.
    rm -f -- "${output_file}"
    : > "${output_file}"
    echo "Job failed with status ${status}; returning an empty failure marker." >&2
  fi
  exit "${status}"
}
trap preserve_failed_status EXIT

case "${sample}" in
  mc) is_mc=true ;;
  data) is_mc=false ;;
  *) echo "ERROR: sample must be mc or data, got ${sample}" >&2; exit 2 ;;
esac

echo "Job started at $(date -u '+%Y-%m-%dT%H:%M:%SZ') on $(hostname)."
echo "Working directory: $(pwd)"
echo "Sample: ${sample}; input list: ${input_list}; output: ${output_file}"

if [[ -z "${X509_USER_PROXY:-}" || ! -r "${X509_USER_PROXY}" ]]; then
  echo "ERROR: HTCondor did not provide a readable X509_USER_PROXY." >&2
  exit 10
fi

if command -v voms-proxy-info >/dev/null 2>&1; then
  proxy_seconds="$(voms-proxy-info -timeleft)"
  proxy_fqan="$(voms-proxy-info -fqan)"
  echo "Proxy lifetime at job start: ${proxy_seconds} seconds."
  case "${proxy_fqan}" in
    *"/cms/"*) echo "CMS VOMS attribute found." ;;
    *) echo "ERROR: delegated proxy has no CMS VOMS attribute." >&2; exit 11 ;;
  esac
fi

if ! command -v root >/dev/null 2>&1; then
  echo "ERROR: ROOT is not available in the inherited lxplus environment." >&2
  exit 12
fi

root-config --version
root -l -b -q mk_compile.C

root_macro="run_zjet_job.C(\"${input_list}\",${is_mc},\"${output_file}\",\"${golden_json}\",\"${lumi_pileup}\",\"${pileup_weights}\")"
root -l -b -q "${root_macro}"

if [[ ! -s "${output_file}" ]]; then
  echo "ERROR: expected output ${output_file} is missing or empty." >&2
  exit 13
fi

echo "Job finished successfully at $(date -u '+%Y-%m-%dT%H:%M:%SZ')."
