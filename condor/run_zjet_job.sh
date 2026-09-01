#!/usr/bin/env bash
set -euo pipefail

sample="${1:?missing sample argument}"
input_list="${2:?missing input-list argument}"
output_file="${3:?missing output-file argument}"
golden_json="${4:-}"
lumi_pileup="${5:-}"
pileup_weights="${6:-}"
jec_l2="${7:-}"
jec_residual="${8:-}"
jer_resolution="${9:-}"
jer_scale_factor="${10:-}"
muon_corrections="${11:-}"
jet_veto_map="${12:-}"

[[ "${golden_json}" == "-" ]] && golden_json=""
[[ "${lumi_pileup}" == "-" ]] && lumi_pileup=""
[[ "${pileup_weights}" == "-" ]] && pileup_weights=""
[[ "${jec_l2}" == "-" ]] && jec_l2=""
[[ "${jec_residual}" == "-" ]] && jec_residual=""
[[ "${jer_resolution}" == "-" ]] && jer_resolution=""
[[ "${jer_scale_factor}" == "-" ]] && jer_scale_factor=""
[[ "${muon_corrections}" == "-" ]] && muon_corrections=""
[[ "${jet_veto_map}" == "-" ]] && jet_veto_map=""

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

# getenv=True is needed for the lxplus ROOT runtime, but it can also inherit a
# ccache directory below the user's small AFS home quota.  Keep every compiler
# temporary and cache access on the execution node instead.
export CCACHE_DISABLE=1
export CCACHE_DIR="${PWD}/.ccache"
export CCACHE_TEMPDIR="${PWD}/.ccache/tmp"
export XDG_CACHE_HOME="${PWD}/.cache"
export TMPDIR="${PWD}/tmp"
mkdir -p "${CCACHE_TEMPDIR}" "${XDG_CACHE_HOME}" "${TMPDIR}"
echo "Compiler cache disabled; temporary files are local to ${PWD}."

if command -v sha256sum >/dev/null 2>&1; then
  echo "Transferred analysis source SHA256 values:"
  sha256sum zjet.C zjet.h FlavorMatrixTools.h ZJetJerResolution.h \
    ZJetMuonCorrections.h \
    data/MuonCorrections/2024_Summer24_generated.h \
    mk_compile.C run_zjet_job.C validateFlavorMatrix.C
fi

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
echo "Compiling the transferred zjet.C and JetMETObjects sources."
root -l -b -q mk_compile.C
for library in \
  CondFormats/JetMETObjects/src/Utilities_cc.so \
  CondFormats/JetMETObjects/src/JetCorrectorParameters_cc.so \
  CondFormats/JetMETObjects/src/SimpleJetCorrector_cc.so \
  CondFormats/JetMETObjects/src/FactorizedJetCorrector_cc.so \
  zjet_C.so; do
  if [[ ! -s "${library}" ]]; then
    echo "ERROR: compilation did not produce ${library}." >&2
    exit 14
  fi
done
echo "Compilation finished; starting the analysis payload."

root_macro="run_zjet_job.C(\"${input_list}\",${is_mc},\"${output_file}\",\"${golden_json}\",\"${lumi_pileup}\",\"${pileup_weights}\",\"${jec_l2}\",\"${jec_residual}\",\"${jer_resolution}\",\"${jer_scale_factor}\",\"${muon_corrections}\",\"${jet_veto_map}\")"
root -l -b -q "${root_macro}"

if [[ ! -s "${output_file}" ]]; then
  echo "ERROR: expected output ${output_file} is missing or empty." >&2
  exit 13
fi

validation_macro="validateFlavorMatrix.C(\"${output_file}\",${is_mc})"
root -l -b -q "${validation_macro}"

echo "Job finished successfully at $(date -u '+%Y-%m-%dT%H:%M:%SZ')."
