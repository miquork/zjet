# Compact local results and full EOS replay data

The full Legacy flavor merge contains `LegacyFlavor/events`: one row per
accepted barrel leading-jet event, with raw tagger scores, linear response
components, weights and additional-jet scores. This supports arbitrary offline
working-point/SF scans. Its size grows with the event count, unlike a fixed-bin
histogram output. Keep the full merge once per campaign in EOS, not in every
Dropbox-backed analysis iteration.

## Finish or repeat the local download

From the original work-area checkout on lxplus, in a persistent tmux session:

```bash
git pull --ff-only
python3 runCondorAnalysis.py \
  --resume run2024i_legacy_flavor_20260910_v1 \
  --compact-download
```

The option is saved in the workflow state. Subsequent resumes of that workflow
keep compact mode. It also works if the previous download was skipped or the
workflow already completed, without resetting the merge checkpoint or rerunning
the analysis. The workflow still asks before replacing local outputs.

For DATA, DY (`MC`) and TT separately, `writeCompactOutput.C` reads the original
EOS file in READ mode and copies all non-tree objects. It reads the replay tree
header, but does not read/download its event baskets. No remerge or new Condor
jobs are needed when the workflow has reached `merged`.

The resulting local files retain the usual names:

- `rootfiles/zjet_DATA.root`
- `rootfiles/zjet_MC.root`
- `rootfiles/zjet_TT.root`

Their size is printed after export. New tagging controls still occupy space;
histogram-only does not imply exactly the size of older, less detailed outputs.
Normal plotting can use these compact files. The full EOS files are unchanged.

## Validation and publication

Each output is created under a unique, same-directory temporary name. Normal
FlavorMatrix and ResponseAudit validators run, followed by
`validateTaggingControls.C(path,isMC,expectedFiles,true)`. Only then does an atomic
local rename replace the previous local file. Failures leave that previous file
intact and remove only this invocation's temporary output. Do not run competing
full and compact downloads into the same destination.

Compact outputs retain campaign provenance, counters, histogram/profile/graph
names and binning. Additional metadata explicitly identifies the compact view:

- `zjet_compact_definition`
- `zjet_replay_entries`
- `zjet_replay_compressed_bytes`
- `zjet_compact_source_uuid`

`LegacyFlavor/definition` remains, but `LegacyFlavor/events` does not. The
exporter refuses pre-existing destinations, duplicate key cycles, already
compact sources, and unexpected additional TTrees rather than silently losing
data. The compact validator checks the stored replay entry count against the
nominal Legacy count. Default worker/full-merge validators still reject compact
files: the relaxed tree requirement is explicit and restricted to this view.

## Offline retagging

`replayLegacyFlavor.C` requires the original full file. Pass its EOS URL directly
on lxplus, or place a working copy outside Dropbox. The workflow records those
URLs under `compact_sources` in its private/local state JSON. No automatic remote
redirect is made when a compact file is passed by mistake: the replay macro
reports which kind of input is required.

Only the resulting small WP/SF histogram variants need to be copied to the
Mac. Do not copy the full merge into Dropbox merely to produce each variant.
The full sample remains one immutable source shared across those iterations.

Tests compare every histogram/profile cell, error, weight sum and statistical
sum against the original file, including all pre-existing metadata strings.
Interrupted-export tests verify that original local outputs and workflow
checkpoints are preserved. This change does not alter physics selection,
normalization, response calculations or the full EOS archive.
