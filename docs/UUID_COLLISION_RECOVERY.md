# Recovering a skim UUID collision

ROOT file UUIDs are not guaranteed to be unique in these skims. Two TT input
files in the September 2026 production shared a UUID but had 296726 and 296441
distinct event IDs respectively, with no intersection. Rejecting every repeated
UUID stopped one chunk even though its two event sets were disjoint.

## Counter policy

`ZJetInputCounters.h` still rejects an input path seen twice. For distinct paths
with the same UUID it now opens separate read-only handles and reads all
`(run, luminosityBlock, event)` values. Every prior file with that UUID is checked.
A disjoint pair is accepted and counted in
`configInfo/UUIDCollisionChecks` (additive through hadd). Shared IDs, repeated IDs
inside a checked file, missing branches, and read failures stop the job. Shared
MC IDs warrant inspection; they do not by themselves prove identical physics
content. There is no automatic event deletion or file skipping.

Normal files do not incur this extra scan. The active analysis tree's branch
addresses and loaded event are not touched. Physics selection, weights, random
seeds and calibration are unchanged. The existing `SkimCounter`,
`GeneratorCounters`, and `counter_definition` schema/content remain compatible
with successful jobs made before the repair. The rejected-collision diagnostic
bin remains zero for verified disjoint collisions.

This is a within-job collision guard, not an exhaustive dataset overlap audit.
Overlaps across job boundaries or files with different UUIDs are not checked.

## Prepare one repaired chunk

From the original lxplus work-area checkout, after `git pull --ff-only`:

```bash
campaign=run2024i_legacy_flavor_20260910_v1_full
python3 scripts/prepare_uuid_retry.py "$campaign" tt_0084 --cluster 9584441
retry="condor/jobs/$campaign/retries/tt_0084_uuidfix"
condor_submit -maxjobs 1 "$retry/zjet.sub"
condor_wait "$retry/condor.log"
python3 scripts/status_condor.py "$campaign"
```

Preparation does not submit. It checks that the original cluster has no queued
jobs, the selected output is missing/empty, the original protected AFS CMS proxy
has at least two hours left, and the original input chunks still match their
recorded selected-list hash. Among recorded worker sources only
`ZJetInputCounters.h` may have changed; hashed calibration inputs must match.
It copies the failed worker logs into the retry directory without removing the
originals. Submission retains the original worker stdout/stderr and ROOT result
destinations, but uses a separate Condor event log. The worker recompiles the
transferred source. No successful ROOT result is removed by preparation.

If proxy validation fails, renew the proxy **at the protected AFS path recorded
in `zjet.sub`**; changing only `X509_USER_PROXY` does not change that submit file.
If another source/configuration changed, investigate instead of bypassing the
check. Reusing a retry label is refused; for another intentional attempt pass a
new `--label` and use its printed commands.

The retry submit file has one explicit queue row, not a second appended queue.
`-maxjobs 1` also limits accidental expansion; see the
[HTCondor submit manual](https://htcondor.readthedocs.io/en/25.0/man-pages/condor_submit.html).
`condor_wait` ending does not certify physics success: require the status check
above to report all outputs ready before resuming:

```bash
python3 runCondorAnalysis.py --resume run2024i_legacy_flavor_20260910_v1
```

The original workflow checkpoint need not be edited. Once the replacement is
ready, its check passes and normal merging/download resumes. Repair preparation
hashes are retained in `retries/*/repair.json` and included in the merged
provenance under `prepared_job_repairs`. These document preparation, not proof
of execution; the worker stdout contains the actual transferred source hashes.
