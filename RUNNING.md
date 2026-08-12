# Running the all-pairs Z+jet analysis

This note documents the local workflow and the intended small-scale lxplus
validation. It deliberately contains no grid proxy, certificate, private file
list, or other credential material.

## Local build and regression run

Compile the analysis and the bundled standalone JetMET correction code after
changing the analysis:

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
8. L2 JEC text file (empty disables JEC recomputation)
9. data residual JEC text file (ignored for MC)
10. MC jet-resolution text file
11. MC JER scale-factor text file
12. Summer24 muon-correction JSON
13. data-only jet-veto-map ROOT file

The pileup-weight histogram must be named `pileup_ratio` (preferred) or
`pileup`. The lumisection pileup reader accepts a brilcalc CSV containing an
`avgpu` column, or a whitespace-separated `run ls mu` file. By default jets
are recomputed from raw pT with the Summer24 V2 L2Relative correction, plus
the Run2024I nib1 V11M L2L3Residual correction for data. The first three
lepton-cleaned MC jets in NanoAOD order receive the configured JER smearing.
Nominal Summer24 muon scale corrections are applied to data and MC and the
additional resolution correction is applied to MC. The 2024 veto map is
applied only to data probe jets. `RawPuppiMET` is rebuilt with all
lepton-cleaned JEC/JER-corrected jets above 15 GeV, without a Jet-ID or veto-map
requirement in that Type-I sum. Empty optional arguments disable only the
corresponding correction or selection and print no credential information.

Example:

```bash
root -l -b -q 'mk_zjet.C("mc.root","data.root","Cert.json","lumi.csv","pu_weights.root",-1,-1)'
```

The nominal synchronized dimuon selection uses
`HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8`, requires both tight isolated
muons to match trigger objects, and keeps the opposite-sign pair closest to
90 GeV. Corrected muon kinematics are used for pair ranking, the 20/10 GeV and
eta cuts, the Z boson, and jet cleaning; trigger matching and the initial 8 GeV
preselection use the stored NanoAOD muons, matching the reference analysis.
The `control/h_muon_selection` histogram evaluates alternative selections on
the same HLT-and-filter input events.

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

An optional data overlay for the four one-dimensional controls is produced
with:

```bash
root -l -b -q drawControlData.C
```

Data is normalized once to the MC parallel-region jet-pT yield over
`30 < pT(jet) < 200 GeV`; the same scale is used for every plot. The upper pad
shows data and MC separately for the parallel, transverse and subtracted
stages, while the lower pad shows their three data/MC ratios. This preserves
the relative stage yields instead of silently renormalizing each distribution.
Input paths and the normalization interval can be changed through the macro
arguments.

`drawPileupControl.C` writes DB and hybrid-MPF profiles versus `PV_npvs`, `rho`, and
`mu`, together with selection efficiencies, truth-matched versus unmatched
jet controls, matching-definition comparisons, explicit signed-yield ratios,
and eta-split closure plots under `pdf/drawPileupControl/`. The primary MC
pileup-jet definition is a missing valid `Jet_genJetIdx`. Extra generator-pT
and DeltaR cuts are retained only as a diagnostic comparison because the
NanoAOD index already encodes the reco-to-gen match.

## jecsys3 global-fit compatibility file

Create the combined directory hierarchy expected by the standalone L2/L3
macros and the main `reprocess.C` -> `softrad3.C` -> `globalFit.C` workflow:

```bash
root -l -b -q writeJecsys3.C
```

The default output is `rootfiles/zjet_JMENANO_compat.root`. It contains the raw
`data|mc/l2res*` trees plus `data|mc/eta_00_13` profiles for MPF, direct
balance, MPF recoil components, PF composition, rho, event counts, Z mass and
MC reco/gen closure. The response and count objects are written for all three
reference-pT choices used by `reprocess.C`: `zmmjet` (Z pT), `jetpt` (probe-jet
pT), and `ptave` (average Z--jet pT). Point both entries below to the same
compatibility file:

```cpp
mfile["ZMM_2024I_nib1_DATA"] =
  "/absolute/path/to/zjet/rootfiles/zjet_JMENANO_compat.root";
mfile["ZMM_2024I_nib1_MC"] =
  "/absolute/path/to/zjet/rootfiles/zjet_JMENANO_compat.root";
```

Both local paths and `root://` EOS inputs are supported through `TFile::Open`.
The EOS `/eos/user/...` namespace is independent of the AFS `/afs/cern.ch/work`
checkout location. The checkpointed driver writes to a `.part` file and only
atomically replaces the requested local compatibility file after ROOT exits
successfully; a failure preserves the previous output and the `merged`
checkpoint.

The normal jecsys3 chain can then run for `2024I_nib1`; its gamma+jet,
multijet, inclusive-jet and W inputs remain configured as before. Keeping the
old ZMM paths as commented alternatives makes the old/new comparison
reversible. Campaign metadata is copied when present.

If the merged `zjet_DATA.root` and `zjet_MC.root` files already exist, changes
to this compatibility mapping do not require rerunning the NanoAOD event loop.
Pull the updated code and rerun only `writeJecsys3.C`; the macro projects all
three pT representations from the raw `l2res` profiles already stored in the
merged files.

Current outputs fill the reco-tag and generator-flavor subsamples expected by
`reprocess.C`. The reco categories use the Bettina/Sami DeepJet definition:
`DeepFlavB > 0.7527`; otherwise
`0.5*(DeepFlavCvB+DeepFlavCvL) > 0.3985` for c; otherwise the quark/gluon split
is `DeepFlavQG = 0.5`. MC generator categories use the matched
`GenJet_partonFlavour`; an unmatched jet is unclassified. Flavor profiles use
the same signed signal-minus-sideband weights as the inclusive response.

`writeJecsys3.C` automatically writes measured flavor inputs when both merged
files contain the `flavor` directory. For older output files, the fourth
argument can create the complete Z+flavor directory and object contract as
empty placeholders for short-term integration testing only:

```bash
root -l -b -q \
  'writeJecsys3.C("rootfiles/zjet_DATA.root","rootfiles/zjet_MC.root","rootfiles/zjet_JMENANO_compat.root",true)'
```

This mode prints a warning, labels every placeholder histogram, and writes the
top-level `zjet_flavor_status` marker. The placeholders only let downstream
code traverse the expected ROOT hierarchy; they contain no flavor measurement
and must not be used as physics input. The default fourth argument is `false`.
When measured inputs are present they take precedence even if the fallback
argument is `true`; `zjet_flavor_status` then identifies them as measured.

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
- the MC residual profile is one and the data residual profile contains the
  inverse V11M residual correction;
- `control/h_jer_smear_factor` is populated only in MC;
- `control/h_muon_scale_factor` is populated in both samples and
  `control/h_muon_resolution_factor` only in MC;
- `control/h_jet_veto_map` is populated only in data;
- the output metadata names the JER, muon, veto-map and Type-I MET settings;
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
  --max-data-files 1 \
  --golden-json data/Cert_Collisions2024_378981_386951_Golden.json

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

After the smoke test succeeds, prepare and submit the full campaign. Store the
large partial ROOT files in EOS, not below the AFS checkout. CERN's batch
XRootD transfer plugin uses the active Kerberos identity for this EOS transfer;
the X.509 proxy is still needed separately by the analysis to read CMS input
files.

```bash
python3 scripts/prepare_condor.py \
  --mc-list textfiles/generated/summer24_mc.txt \
  --data-list textfiles/generated/run2024i_data.txt \
  --campaign run2024i_full \
  --files-per-job 10 \
  --job-flavour workday \
  --golden-json data/Cert_Collisions2024_378981_386951_Golden.json \
  --log-dir \
    /afs/cern.ch/work/v/voutila/zjet-condor/run2024i_full \
  --eos-results \
    root://eosuser.cern.ch//eos/user/v/voutila/zjet/run2024i_full

condor_submit condor/jobs/run2024i_full/zjet.sub
```

`condor_submit` returns immediately; the campaign continues after the SSH
connection closes. The generated submit file requests one CPU, 2 GB of memory,
EL9, the selected job flavour, and an HTCondor-delegated X.509 proxy.
The source is compiled inside each worker scratch directory, so compiled files
are never shared between concurrent jobs.

The generated `campaign.json` records the Git commit and dirty state, hashes
and basenames of the input lists, golden JSON, L2 and residual JECs, JER files,
muon corrections, the data jet veto map, and optional pileup files. It also
records the correction order and Type-I MET definition. It intentionally does
not copy the full private file URLs into the publishable provenance.

The `--log-dir` path keeps stdout, stderr and the HTCondor event log outside
the limited home AFS volume. Use a new campaign-specific directory below the
CERN AFS work area. Do not use an `/eos/...` FUSE path for these submit-file
fields; large ROOT outputs are transferred separately to the `root://`
destination by HTCondor's XRootD plug-in. The selected log directory is stored
in `campaign.json` and used automatically by `status_condor.py`.

In the morning, wait until `condor_q` shows no campaign jobs and merge only
after every expected output has returned:

```bash
python3 scripts/status_condor.py run2024i_full
python3 scripts/merge_condor.py run2024i_full \
  --output-dir \
    root://eosuser.cern.ch//eos/user/v/voutila/zjet/run2024i_full/merged
```

The status tool prints `READY TO MERGE` only after every output is non-empty
and its worker log ends successfully. The merge tool independently refuses to
run if any expected job output is missing or empty.
Without `--force`, it also refuses to replace existing output files. With an
EOS output directory, the merge uses a local temporary directory and uploads
the two merged outputs, their JSON provenance, and a human-readable merge log
to EOS. The same JSON is embedded in both ROOT files as
`zjet_campaign_metadata`.

### Checkpointed one-command driver

`runCondorAnalysis.py` wraps the same small scripts without hiding their
commands or safety checks. It currently supports only the Run 2024I preset.
List or inspect the supported plan without changing anything:

```bash
python3 runCondorAnalysis.py --list-presets
python3 runCondorAnalysis.py --preset run2024i --plan
```

Start a new, immutable workflow name after the generated MC and data lists are
available:

```bash
python3 runCondorAnalysis.py \
  --preset run2024i \
  --campaign run2024i_v11m_20260812
```

The driver validates Kerberos and a CMS VOMS proxy with at least 24 hours
remaining, copies a valid `/tmp` proxy to protected AFS storage after approval,
optionally pulls with `git pull --ff-only`, reports AFS/work space, compiles the
analysis, tests the first MC and data URLs, and prepares a one-file-per-sample
worker smoke test. The preflight also verifies that the committed nominal muon
lookup header was generated from the committed JSON. It then prepares the full
campaign only after the smoke test has completed and passed
`status_condor.py`.

Submission, merging, and replacement of the compatibility ROOT file each
require a separate `y` answer. Answering `n` stops at a recorded checkpoint;
submitted jobs continue independently. Resume later with the command printed
by the driver, for example:

```bash
python3 runCondorAnalysis.py --resume run2024i_v11m_20260811
```

Workflow state is stored below the git-ignored `condor/jobs/_workflows`
directory. `Ctrl-C` is also safe between polling calls and prints the same
resume command. The displayed CPU and wall-time estimates are deliberately
rough planning numbers based on the first Run 2024I tests, not scheduler
guarantees. Keep eras as separate workflow names; their EOS partial and merged
outputs remain separate and can be combined explicitly only after each era has
been validated.

### Recovering an AFS log quota hold

A job held with code `12/122` and a reason naming its `.out` or `.err` file has
finished the payload but could not return a standard-stream log because the
AFS quota was exhausted. Do not remove the jobs. Either free enough space and
release the cluster, or redirect the held jobs to a campaign-specific AFS work
directory before releasing them:

```bash
cluster_id=1234567
recovery_logs=/afs/cern.ch/work/u/username/zjet-condor/recovery_1234567
mkdir -p "$recovery_logs"
chmod 700 "$recovery_logs"

while read -r proc_id old_out old_err; do
  condor_qedit "${cluster_id}.${proc_id}" Out \
    "\"${recovery_logs}/${old_out##*/}\""
  condor_qedit "${cluster_id}.${proc_id}" Err \
    "\"${recovery_logs}/${old_err##*/}\""
done < <(condor_q -constraint \
  "ClusterId == ${cluster_id} && JobStatus == 5" -af ProcId Out Err)

condor_qedit -constraint \
  "ClusterId == ${cluster_id} && JobStatus == 5" UserLog \
  "\"${recovery_logs}/condor.log\""
condor_release "$cluster_id"

python3 scripts/status_condor.py run2024i_full \
  --additional-log-dir "$recovery_logs"
```

Changing `Out`, `Err` and `UserLog` preserves the EOS `output_destination` for
the ROOT files. Check EOS and `condor_q -hold` again after the release. A
released transfer hold may retry the transfer or restart the payload, so do not
delete any existing EOS partial output until the cluster has reached a stable
completed state.

To merge into the checkout instead, omit `--output-dir`. This consumes the
size of the final two files in AFS, but not an additional same-sized temporary
file. After merging, rerun the normal control plots and compatibility writer.

Inspect the space occupied by an old campaign before deleting anything:

```bash
python3 scripts/cleanup_condor.py run2024i_full
```

After a successful merge, archive its small provenance files and remove the
AFS logs/chunks. Add the second flag only when the EOS partial files are no
longer needed:

```bash
python3 scripts/cleanup_condor.py run2024i_full --delete
# Or, when the EOS partial outputs are no longer needed:
python3 scripts/cleanup_condor.py run2024i_full \
  --delete --delete-remote-results
```

The cleanup tool accepts only a directory below `condor/jobs` and deletes only
the exact EOS objects listed in that campaign's manifest. It normally requires
merge provenance and keeps the small manifest directory so that the settings
remain inspectable and remote partials can still be removed later. An
intentionally abandoned smoke test can be cleaned without merging only by
adding the explicit `--discard-unmerged` flag.

Before remote deletion, cleanup verifies that a CMS VOMS proxy has at least ten
minutes remaining. It treats already absent manifest files as successfully
cleaned, so the same command can safely resume after a partially completed
deletion. Other per-file failures are collected and reported after independent
files have been attempted.

Optional data certification, data pileup, and MC pileup-weight inputs can be
transferred to every job with `--golden-json`, `--lumi-pileup`, and
`--pileup-weights`. Generated campaigns, direct storage URLs, logs, results,
and credentials are ignored by Git and must remain outside the public history.
