# Running the all-pairs Z+jet analysis

This note documents the local workflow and the intended small-scale lxplus
validation. It deliberately contains no grid proxy, certificate, private file
list, or other credential material.

## Local build and regression run

Compile the analysis after changing `zjet.C`, `zjet.h`, or `ZJetLumi.h`:

```bash
root -l -b -q mk_compile.C
```

Run the default local Summer24 MC and Run2024I data files:

```bash
root -l -b -q mk_zjet.C
```

The outputs are written directly to `rootfiles/zjet_MC.root` and
`rootfiles/zjet_DATA.root`. The macro no longer creates a temporary output and
renames it with an interactive shell command.

The optional arguments of `mk_zjet` are, in order:

1. MC ROOT file
2. data ROOT file
3. golden-JSON file for data
4. lumisection pileup file for the data `mu` control
5. ROOT file containing the MC pileup-weight histogram
6. maximum number of MC files (`-1` means all)
7. maximum number of data files (`-1` means all)

The pileup-weight histogram must be named `pileup_ratio` (preferred) or
`pileup`. The lumisection pileup reader accepts a brilcalc CSV containing an
`avgpu` column, or a whitespace-separated `run ls mu` file. Empty optional
arguments disable only the corresponding correction or selection and print no
credential information.

Example:

```bash
root -l -b -q 'mk_zjet.C("mc.root","data.root","Cert.json","lumi.csv","pu_weights.root",-1,-1)'
```

The nominal dimuon selection scans all opposite-sign pairs and keeps the one
closest to the Z mass. It requires one `pT > 27 GeV`, tight-ID, tight-PF-
isolation tag and one `pT > 10 GeV`, medium-ID, loose-PF-isolation probe. The
single-muon `HLT_IsoMu24` requirement constrains the tag, not the probe. The
`control/h_muon_selection` histogram evaluates alternative ID, isolation, and
threshold choices on the same HLT-and-filter input events.

The two transverse directions are independent half-weight
signal-sideband pairs. If a lepton overlaps only one sideband, that sideband
and only the corresponding half of the signal are vetoed. The
`control/h_probe_veto` and `control/h_probe_pair_state` histograms expose this
acceptance explicitly.

An input ending in `.txt` is interpreted as a newline-separated ROOT file
list. A single `.root` input continues to work as before.

## Converting the sample JSON files

The received sample descriptions contain a top-level `files` string array.
Convert the MC JSON to one text list and all non-overlapping data JSON files to
a second list:

```bash
python3 scripts/json_to_filelist.py /path/to/mc.json \
  -o textfiles/generated/summer24_mc.txt

python3 scripts/json_to_filelist.py /path/to/data_part1.json \
  /path/to/data_part2.json /path/to/data_part3.json /path/to/data_part4.json \
  -o textfiles/generated/run2024i_data.txt
```

The converter removes duplicate URLs and, by default, maps every `/store/...`
logical file name to the CMS global XRootD redirector. If a file is not
published to the federation, pass the storage site's XRootD endpoint explicitly
with `--redirector root://<site-host>:1094/`. Generated lists are ignored by
Git and must not be committed.

The global redirector can return XRootD error 3011 even when a file exists at
its original storage site. In that case, preserve the direct URLs supplied in
the JSON documents:

```bash
python3 scripts/json_to_filelist.py /path/to/mc.json --keep-urls \
  -o textfiles/generated/summer24_mc.txt

python3 scripts/json_to_filelist.py /path/to/data_part1.json \
  /path/to/data_part2.json /path/to/data_part3.json /path/to/data_part4.json \
  --keep-urls -o textfiles/generated/run2024i_data.txt
```

Test the first direct URL before starting ROOT:

```bash
gfal-ls "$(sed -n '1p' textfiles/generated/summer24_mc.txt)"
gfal-ls "$(sed -n '1p' textfiles/generated/run2024i_data.txt)"
```

The storage site used for these samples supports both XRootD and HTTPS/WebDAV.
If direct HTTPS works with `gfal-ls` but the installed ROOT cannot open it, find
the site's current XRootD prefix in
`/cvmfs/cms.cern.ch/SITECONF/<site>/storage.json` and regenerate the lists with
`--redirector`. Do not guess the endpoint or commit a site-specific file list.

## Control plots

The original control plots are produced with:

```bash
root -l -b -q drawControl.C
```

Pileup and truth-matching controls are produced with:

```bash
root -l -b -q drawPileupControl.C
```

The latter writes DB and hybrid-MPF profiles versus `PV_npvs`, `rho`, and
`mu`, together with selection efficiencies, truth-matched versus unmatched
jet controls, matching-definition comparisons, explicit signed-yield ratios,
and eta-split closure plots under `pdf/drawPileupControl/`. The primary MC
pileup-jet definition is a missing valid `Jet_genJetIdx`. Extra generator-pT
and DeltaR cuts are retained only as a diagnostic comparison because the
NanoAOD index already encodes the reco-to-gen match.

## jecsys3 compatibility file

Create the combined directory hierarchy expected by `jecsys3/L2Res.C`:

```bash
root -l -b -q writeJecsys3.C
```

The default output is `rootfiles/zjet_JMENANO_compat.root`, containing
`data/l2res`, `data/l2res1`, `mc/l2res`, and `mc/l2res1`. Point both the
`ZMM_<era>_DATA` and `ZMM_<era>_MC` entries in `jecsys3/Config.C` to this same
file when testing the new method. Keeping the old file paths in a separate
configuration makes the old/new comparison reversible.

## lxplus small-file validation

Log in with the CERN account and verify the proxy:

```bash
ssh <cern-user>@lxplus.cern.ch
voms-proxy-info -timeleft
voms-proxy-info -fqan
export X509_USER_PROXY="$(voms-proxy-info -path)"
```

If either proxy command fails or the lifetime is too short, create a fresh CMS
proxy and export its path explicitly:

```bash
voms-proxy-init --voms cms --valid 192:00
export X509_USER_PROXY="$(voms-proxy-info -path)"
voms-proxy-info -timeleft
```

The PEM pass phrase belongs only in the interactive `voms-proxy-init` step. If
ROOT itself asks for it, stop the job: the XRootD client did not find a usable
proxy and is trying to create one from the long-lived certificate. Never put
the certificate password in a command, script, file list, or batch job.

After cloning or pulling this repository on lxplus, source a ROOT environment
compatible with the JMENANO files, compile, and first run one MC and one data
file through XRootD:

```bash
root -l -b -q mk_compile.C
root -l -b -q \
  'mk_zjet.C("textfiles/generated/summer24_mc.txt","textfiles/generated/run2024i_data.txt","","","",1,1)'
```

Then repeat with a few files, for example by replacing the final `1,1` with
`3,3`. Check that:

- both remote files open before starting the event loop;
- the cutflow is non-zero through the paired probe-lepton veto;
- the `h_muon_selection` and `h_probe_veto` diagnostics are populated;
- `l2res/p2m0tc` and `l2res/p2restc` both have entries;
- the `l2res` eta underflow is empty;
- the residual profile is exactly one in every populated bin;
- the two half-weight transverse windows have the expected normalization.

Do not submit the full sample or HTCondor jobs before the one- and few-file
checks pass.

At startup, the macro reports that it is opening the remote files before
calling `GetEntries()`. It then reports the first successfully read event,
early rate estimates at 1k, 10k, and 100k events, and updated completion-time
estimates at least once per minute. Output directories are created
automatically. A zero-entry chain or a read failure is treated as an error and
does not leave an apparently valid output file behind.

## CERN HTCondor campaign

Do not leave the full sample in an interactive SSH session. The campaign tool
splits each input list into small lists, submits MC and data as independent
jobs, delegates the CMS VOMS proxy, and returns only histogram ROOT files. The
default chunk size is ten remote files per job. For the current lists this
means 77 MC jobs and 69 data jobs instead of one multi-hour process.

Run all commands below in the same bash environment in which both `root` and
`condor_submit` are available. Check the Kerberos ticket and CMS proxy first:

```bash
klist
voms-proxy-info -timeleft
voms-proxy-info -fqan
```

The FQAN must contain `/cms/`. Create a fresh long-lived proxy before an
overnight campaign if necessary:

```bash
voms-proxy-init --rfc --voms cms --valid 192:00
```

The default proxy is under `/tmp` on one lxplus host and is not visible to the
HTCondor scheduler on another host. Copy it to a private AFS directory before
preparing any campaign, then point `X509_USER_PROXY` to that shared copy:

```bash
mkdir -p "$HOME/private"
chmod 700 "$HOME/private"
fs setacl -dir "$HOME/private" -acl system:anyuser none
fs setacl -dir "$HOME/private" -acl "$(whoami)" all

afs_proxy="$HOME/private/x509up_u$(id -u)"
cp "$(voms-proxy-info -path)" "$afs_proxy"
chmod 600 "$afs_proxy"
export X509_USER_PROXY="$afs_proxy"

voms-proxy-info -file "$X509_USER_PROXY" -timeleft
voms-proxy-info -file "$X509_USER_PROXY" -fqan
```

Never copy the proxy into the repository. The campaign generator rejects
non-AFS proxies and proxy files located under the public repository. Refresh
the AFS copy whenever a new VOMS proxy is created.

First prepare a two-job worker-node smoke test. Campaign names are immutable;
use a new name if a campaign directory already exists:

```bash
python3 scripts/prepare_condor.py \
  --mc-list textfiles/generated/summer24_mc.txt \
  --data-list textfiles/generated/run2024i_data.txt \
  --campaign run2024i_worker_test \
  --files-per-job 1 \
  --max-mc-files 1 \
  --max-data-files 1

condor_submit condor/jobs/run2024i_worker_test/zjet.sub
```

Monitor the jobs and inspect their logs:

```bash
condor_q -nobatch
condor_wait condor/jobs/run2024i_worker_test/logs/condor.log
python3 scripts/status_condor.py run2024i_worker_test
find condor/jobs/run2024i_worker_test/results -name '*.root' -type f
sed -n '1,160p' condor/jobs/run2024i_worker_test/logs/mc_0000.out
sed -n '1,160p' condor/jobs/run2024i_worker_test/logs/data_0000.out
```

Both logs must end in `Job finished successfully`, and both ROOT files must be
present. The worker validates the proxy FQAN and output cutflow. Transient
worker failures are retried twice on a different worker.

After the smoke test succeeds, prepare and submit the full campaign:

```bash
python3 scripts/prepare_condor.py \
  --mc-list textfiles/generated/summer24_mc.txt \
  --data-list textfiles/generated/run2024i_data.txt \
  --campaign run2024i_full \
  --files-per-job 10

condor_submit condor/jobs/run2024i_full/zjet.sub
```

`condor_submit` returns immediately; the campaign continues after the SSH
connection closes. The generated submit file requests one CPU, 2 GB of memory,
EL9, the two-hour `longlunch` flavour, and an HTCondor-delegated X.509 proxy.
The source is compiled inside each worker scratch directory, so compiled files
are never shared between concurrent jobs.

In the morning, wait until `condor_q` shows no campaign jobs and merge only
after every expected output has returned:

```bash
python3 scripts/status_condor.py run2024i_full
python3 scripts/merge_condor.py run2024i_full --force
```

The status tool prints `READY TO MERGE` only after every output is non-empty
and its worker log ends successfully. The merge tool independently refuses to
run if any expected job output is missing or empty.
Without `--force`, it also refuses to replace existing
`rootfiles/zjet_MC.root` and `rootfiles/zjet_DATA.root`. After merging, rerun
the normal control plots and compatibility writer.

Optional data certification, data pileup, and MC pileup-weight inputs can be
transferred to every job with `--golden-json`, `--lumi-pileup`, and
`--pileup-weights`. Generated campaigns, direct storage URLs, logs, results,
and credentials are ignored by Git and must remain outside the public history.
