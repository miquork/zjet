# Controlled tagger and unclustered-response study

## Scope and decision order

The `ResponseAudit` addition changes no production selection, calibration,
sideband estimator, legacy output, or default Ru. It records parallel probes
with `abs(eta)<1.3` alongside the existing analysis. Use the same full Run2024I
and Summer24 MC preset first; do not change JER or muon corrections concurrently.

1. Align discriminator directions and check MC-as-data inversion closure.
2. Compare taggers with identical selected probes, weights, heavy-flavor veto,
   and (for the continuous light-score controls) valid-score support.
3. Compare unmodified MC purities and the minimum-KL/IPF inferred data matrix.
   Inspect conditioning and unknown/invalid categories before interpreting JES.
4. Separate mean unclustered closure, event regression, and jet-threshold migration.
5. Only then propagate a validated model to HDM and examine fn+fu as an independent
   physics test. Do not choose Ru by forcing a desired color-factor ratio.

## Important discriminator convention

`Jet_btagDeepFlavQG` is **GvsQ**, whereas PNet and UParT QvG are **QvsG**.
This is documented in the local NanoAOD branch descriptions and the
[CMSSW NanoAOD producer](https://github.com/cms-sw/cmssw/blob/master/PhysicsTools/NanoAOD/python/jetsAK4_Puppi_cff.py).
All new aligned score controls use `1-DeepFlavQG` for DeepJet. Do not transform
invalid negative values into valid scores. Existing historical outputs are kept
unchanged for reproducibility: their q/g reco labels are inverted. A consistent
row permutation alone cannot change a matrix-inverted result.

## Payload layout

`ResponseAudit/{new,legacy,common_new,common_legacy}/{policy}` contains `counts`
(TH3D) and TProfile3D objects `m0,m2,mn,mu,db,inverse_residual,reco_gen`:

- x: the existing fine pT,Z bins; y: reco tag ID; z: Jet parton flavor ID.
- `new`: all accepted parallel barrel probes, with their usual pair weights.
- `legacy`: accepted leading barrel probe, with the historical legacy weight.
- `common_new` and `common_legacy`: exactly the same accepted jet index, both
  weighted with the new parallel-pair weight; only the response definition differs.
- `inverse_residual` is the stored inverse residual JEC, not the complete JES.
- `reco_gen` is corrected/smeared reco pT divided by the matched GenJet pT.

Policies: historical and direction-corrected DeepJet; UParT HF at 0.5 with PNet
QvG thresholds 0.2/0.3/0.4, aligned DeepJet 0.3/0.5/0.7, UParT 0.5/0.8/0.95.
These are analysis scans, **not certified tagging-SF working points**.

`{new,common_new}/{deep,pnet,upart}_light_scores` uses the common UParT HF veto
and requires all three QvG scores to be valid. Its y axis has 100 score bins.
Choose equal MC gluon efficiencies on this common sample after production.
The older `FlavorMatrix/taggerAudit` has no common HF veto and a coarser grid.

`truth_labels` compares Jet and matched GenJet parton-flavor assignments.
The older flavor analysis and the compact matrix do not use the same branch.
Truth 0 is retained separately in this audit, not silently added to gluons.
`selection` uses z codes 0=no accepted barrel legacy, 1=same jet, 2=different jet.
Numeric axes deliberately have no labels that could be reordered by ROOT hadd.

## Ru sufficient statistics

For each pT bin and truth/reco category, `Ru/<mode>_<label>_blocks` stores
20 deterministic event-hash blocks, each with weighted sums of
`1,x,y,x*x,x*y,y*y,n,n*n,x*n,y*n`. Every probe from the same event shares a block.
They add under hadd; use delete-group jackknife for correlated-ratio uncertainty.
These estimates need enough populated events/blocks and are not reliable in sparse
tails. Signed generator weights are retained.

- `native`: x=gen mu, y=reco mu, n=gen mn, matched probe and matched generator Z.
- `native_lowrho/highrho`: same, divided at rho=25 GeV (not a calibrated UE cut).
- `native_complete`: native, restricted to complete unique GenJet matching of
  all reco HT jets and inclusion of the probe in that list.
- `common`: same subset, but generator HT is built from that reco-selected
  list, replacing each reco jet with its matched GenJet. Compensating in U keeps
  gen mn+mu invariant. The matching requirement can bias this subset; compare
  to `native_complete`, **not** directly to inclusive `native`.
- `closure`: x=1-reco m2-reco mn, y=reco mu, n=reco mn (also available in data).
- `closure_matched`: closure on the native truth-matched population.
- `closure_truth`: x=gen mu, y=the reco closure proxy, n=gen mn.

The accompanying conditional profiles retain x under/overflow; x in [-1,1]
is only the plotted profile axis, not a selection on the moment sums.

For x=gen mu and y=reco mu distinguish:

```
zero-intercept slope R0 = <xy>/<x^2>
mean-closing ratio Rmean = <y>/<x>
affine slope b = Cov(x,y)/Var(x), intercept a = <y>-b<x>
mean bias after slope correction = <y>/R0 - <x>
```

R0 minimizes the forward prediction loss, not HDM mean bias. Mean ratios can
be unstable when a signed denominator is near zero. The new cross moments also
allow `y = a + b_u*x + b_n*n`; correlated n/u partition fluctuations need not
be represented by a single universal Ru.

The data closure proxy reuses the same measured recoil: `x=1-m0+mu`.
Correlated noise, recoil resolution, assumed R2/Rn and invisible generator
momentum can spoil its interpretation. Data proxy agreement is not an independent
Ru measurement. First validate/calibrate the proxy against truth in MC.

## Running the read-only study

From the repository containing the updated data/MC ROOT files:

```bash
root -l -b -q 'studyTaggerRu.C()'
latexmk -pdf -interaction=nonstopmode -halt-on-error -outdir=output/pdf taggerRuStudy.tex
```

This works on the current full production without the new controls and writes
central-value diagnostics under `output/taggerRu`. Full error correlations are
not present in the older files, so those graphs intentionally have no error bars.
The new production additionally generates `block_regressions.tsv`.
The Beamer source reads generated figures/tables: numerical results remain in
ignored local output directories, not in the public repository.

With the new audit production:

```bash
root -l -b -q 'validateResponseAudit.C("rootfiles/zjet_DATA.root")'
root -l -b -q 'validateResponseAudit.C("rootfiles/zjet_MC.root")'
root -l -b -q 'compareTaggerResponse.C()'
root -l -b -q 'compareTaggerResponse.C("rootfiles/zjet_MC.root","rootfiles/zjet_MC.root","output/taggerRu/mc_closure.tsv",true)'
```

The crossed table uses the same algorithm for every policy/cohort, no unity prior,
and Ru scans 0.5/0.92/1.1. Only full-rank solutions receive numerical scales.
Errors are **conditional** on templates/purities and omit component covariance.
The four-flavor fit excludes truth 0 (reported separately) and undefined reco
tags. It is not the central published flavor analysis.

For HDM, form N=m0-mn-mu and D=1-mn/Rn-mu/Ru. Mix N and D before taking their
ratio. The fit uses `sum(P_f N_f k_f)/sum(P_f D_f)`, with fixed MC D_f and flavor
scales k_f; this is a specific response model. Ordinary purity-weighted averaging
of per-cell HDM is different. Compare models before interpreting a residual.

## Condor preparation

The new header and validator are included in the submit payload and source hashes.
`mk_compile.C` already forces a fresh worker build; no manual compile is needed
on lxplus. Worker validation checks identical policy populations, common-pair
counts, finite moments, and common-partition N+U closure. Merge repeats the audit
validation when present and permits historical productions without it.

Use a fresh campaign, not resume of a campaign with the old source snapshot.
From the existing work-area checkout on lxplus, in bash:

```bash
git pull --ff-only
campaign="run2024i_tagger_ru_$(date +%Y%m%d_%H%M%S)"
python3 runCondorAnalysis.py --preset run2024i --campaign "$campaign" --plan
python3 runCondorAnalysis.py --preset run2024i --campaign "$campaign"
```

The second invocation is interactive and requires explicit submission approval.
Check the shown /work log quota and EOS destination; keep the smoke checkpoint.
The workflow checks the proxy, rebuilds locally and on workers, runs one MC and
one data file, waits, merges, and offers copying merged files into local rootfiles.
No jobs are submitted by the read-only study or by `--plan`.

## Tests

```bash
root -l -b -q 'tests/testResponseAudit.C()'
python3 -m py_compile scripts/prepare_condor.py scripts/merge_condor.py runCondorAnalysis.py
bash -n condor/run_zjet_job.sh
```

The local integration test must use separate output files, never overwrite the
full-production `rootfiles/zjet_DATA.root` and `rootfiles/zjet_MC.root`.
Check MC-as-data closure and equality of historical/corrected DeepJet fitted
scales (a pure row permutation). ROOT macros are run in separate processes to
avoid interpreter redefinition problems on older ROOT versions.
