# Legacy-only flavor baseline with DY and TT

Use the `run2024i_legacy_tt` preset for this production. It runs the synchronized
leading-jet analysis and **does not run the new parallel/sideband selection**.
Legacy unit MC weights, its alpha cut and its disabled JetID are retained for
baseline replication. Full-generator-weight and signed-weight controls are also
stored, but do not silently replace the reference convention.

## Tag definition

The primary result is under `TaggingControls/legacy/bT_cvlT_pnet045`:

1. B Tight on `Jet_btagUParTAK4B` (not `1-CvB` as a proxy).
2. If B fails, CvL Tight on `Jet_btagUParTAK4CvL`.
3. If both fail, `Jet_btagPNetQvG >= 0.45` selects q; below selects g.

The CvL-only choice matches the available reference and is **not** labelled as
the full official two-dimensional C WP. B-Medium and the eight B-M/T x
two-dimensional C-M/T x PNet-0.30/0.45 variants remain available for comparison.
The official values are read from the Git-ignored `data/Tagging/official_2024.txt`.
No missing charm or QvG SF is invented or automatically applied.

## Launch

Keep using the work-area checkout and its **working DY/data URL lists**. Some
old locally generated lists use the global redirector, which cannot locate these
private skims; do not overwrite the working HTTPS lists with those.

The privately transferred `legacy_run_inputs.tgz` contains only:

- `data/Tagging/official_2024.txt`
- `textfiles/generated/summer24_tt.txt`

It contains no credentials and is not tracked in Git. Copy it to your CERN work
area, extract it in the work checkout, pull the code, then run:

```bash
python3 runCondorAnalysis.py --preset run2024i_legacy_tt \
  --campaign run2024i_legacy_flavor_20260910_v1
```

Use a fresh campaign name for another attempt. The workflow checks credentials,
stages a protected AFS proxy if necessary, rebuilds `mk_compile.C`, runs three
one-file smoke jobs, validates them, and asks before submitting the full sample.
Logs use AFS work storage; partial/merged results use EOS. The final download step
copies DATA, DY (`zjet_MC.root`) and TT (`zjet_TT.root`) locally.

The workflow deliberately skips compatibility writing in this preset: the old
writer's flavor inputs are all-pairs, not these new Legacy flavor inputs. No
mislabelled all-pairs flavor file or unnormalized DY+TT combination is produced.

## Retag and reweight without another NanoAOD production

`LegacyFlavor/events` has one row for every accepted barrel Legacy probe **before
tagging**, including invalid/undefined scores. It stores all three QvG outputs,
UParT B/CvB/CvL, raw parton/hadron labels, corrected jet pT/eta, response components,
inverse residual, muon fraction, mass, MET projections and full/signed weights.
Scores of additional lepton-cleaned HF-candidate jets are retained as vectors.
Run/lumi/event IDs permit event-level resampling within each dataset; they do not
uniquely identify events across different MC datasets.

For the nominal B-Tight / CvL-Tight / PNet-0.45 reference:

```bash
root -l -b -q \
  'replayLegacyFlavor.C("rootfiles/zjet_MC.root","rootfiles/legacy_MC_nominal.root")'
```

Change to B-Medium / CvL-Tight / PNet-0.30 with no event rerun:

```bash
root -l -b -q \
  'replayLegacyFlavor.C("rootfiles/zjet_MC.root","rootfiles/legacy_MC_medium_q030.root","M","T",0.30)'
```

The sixth argument enables the two-dimensional C cut. Arguments 7--8 choose
`unit`, `gen` or `signed` MC weights and an optional explicit SF-weight table.
Arguments 9--11 override B, CvL and CvB thresholds; negative values select the
configured WP. Argument 12 selects raw `pnet`, `deep` or `upart` QvG. Verify the
score's quark/gluon orientation before interpreting a different branch.

The optional SF table uses whitespace-delimited columns:

```text
# absParton hadron recoTag ptMin ptMax absEtaMin absEtaMax factor
```

Flavor/tag `-1` is a wildcard; raw parton 21 denotes gluons. Bounds are half-open.
Rows cannot overlap. Unspecified cells have weight 1. Factors apply only to MC.
This permits independent flavor/tag efficiency variations, including failed-b
migration to charm versus q/g. It does **not** infer these factors from data or
convert a binary official B SF into all conditional transition SFs. Such a
conversion must supply an explicit efficiency/complement model.

Replays store counts and linear response profiles (plus component products),
not averages of per-event nonlinear HDM. Compute HDM from the final component
means. Use a fresh output name; existing files are not overwritten. Kinematic
selection cannot be loosened beyond the accepted Legacy sample after production.

DY and TT remain separate. Only the full-generator-weight variant can use the
stored `Runs` sumw for cross-section normalization. Raw pre-skim counts and sign
weight sums are not interchangeable. Never raw-hadd the two processes.

## Checks performed

- Full local one-file runs for data, DY and TT with production corrections.
- Empty new-method output populations in Legacy-only mode.
- Offline nominal retag counts equal production counts in every bin.
- Replayed component means agree within float storage precision, below 3e-7 in
  the tested samples (well below one per mille).
- Offline three-sample Condor manifest preparation, without submission.
- Additive pre-skim/Runs counters and metadata deduplication across merging.

The private payloads, input lists, replay trees, numerical results and PDFs are
not part of the public Git commit.
