# Legacy Z+jet synchronization audit

This audit compares the local `legacy/` control implementation with
`slehti/ZbAnalysis` master commit `46dbf3401250a79154129268efc9903df7db2c3d`.
The CERN GitLab HEAD was rechecked on 2026-08-12 and still pointed to this
commit.

## Already synchronized

- Trigger: `HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8`.
- Both selected muons are matched to a muon trigger object within DeltaR 0.3.
- Tight muon ID, `Muon_pfRelIso04_all < 0.15`, pT thresholds 20/10 GeV,
  and `|eta| < 2.3`.
- The opposite-sign pair closest to 90 GeV is selected, followed by
  `pT(Z) > 12 GeV` and `|m(mumu)-90 GeV| < 20 GeV`.
- The leading lepton-cleaned jet has `pT >= 12 GeV` and `|eta| <= 5`.
- The back-to-back residual requirement is below 0.44 radians.
- The second jet is taken from the pT-ordered, lepton-cleaned collection after
  the 10 GeV jet threshold. Alpha is `pT(jet2)/pT(Z)`, is set to zero when the
  second jet is below 15 GeV, and the `a100` profiles require `alpha < 1`.
- The central response profiles use the strict interval `0 < |eta| < 1.3`.
- MC events with more than 100 true pileup interactions are excluded before
  JER smearing. This also keeps the stochastic-smearing random sequence aligned
  with the reference and gives the legacy and all-pairs methods the same event
  sample.
- Summer24 V2 MC L2Relative and 2024I V11M data L2L3Residual corrections are
  recomputed from raw jet pT.
- Pileup reweighting is disabled in both nominal workflows. It is not the
  source of the current legacy difference.

## Corrections made in this audit

- The legacy MPF1 and MPFn four-vectors are now explicitly flattened to the
  transverse plane before projection. Without this step their longitudinal
  components produced unphysical values as large as O(100), although the
  components cancelled in total MPF.
- Jet ID is disabled only in the `legacy/` synchronization control, matching
  the current `ZbAnalysis` reference as confirmed by its maintainer. The normal
  all-pairs analysis continues to require reconstructed Run-3 Tight Jet ID for
  probe jets. Type-I MET and non-leading recoil sums do not use Jet ID.
- The 2024 forward spike veto used by `ZbAnalysis` is reproduced.
- The selected-muon multiplicity is restricted to two or three.
- The closest-mass pair is chosen from all tight, isolated, trigger-matched
  muons above 8 GeV before applying the selected pair's eta and 20/10 GeV pT
  cuts. Jets are cleaned against all two or three selected muons, not only the
  pair assigned to the Z boson.
- The `ZbAnalysis` MET-filter list and its data-only application are used.
- Legacy MC profiles use unit event weights, matching `ZbAnalysis`. The new
  all-pairs method retains signed generator weights.
- Summer24 nominal muon scale corrections are applied before pair ranking and
  kinematic cuts; the deterministic additional resolution is applied to MC.
- The first three lepton-cleaned MC jets in NanoAOD order are smeared with the
  same JER resolution, scale-factor parametrization, and random-engine seed as
  `ZbAnalysis`.
- `jetvetoReReco2024_V9M.root:jetvetomap` is applied only to data analysis
  jets. It is deliberately absent from MC and from the Type-I MET sum.
- Type-I PUPPI MET is rebuilt from `RawPuppiMET` with the same lepton-cleaned,
  JEC/JER-corrected jets above 15 GeV used in the recoil sum.
- The missing `alpha < 1` requirement was restored. Earlier local files used
  the `a100` name without applying the cut, which selected severely unbalanced
  multijet events in high jet-pT bins. In those events a large positive MPF1
  contribution was partly cancelled by a negative MPFn contribution, hiding
  the problem in total MPF and HDM.
- The reference analysis's orientation-dependent forward spike veto is
  reproduced literally for synchronization. It should be revisited together
  with Jet ID after the legacy result has been matched.
- `legacy/control` stores alpha, second-jet pT, alpha versus leading-jet pT,
  DB/MPF1 profiles before the alpha cut, and rho versus Z pT before the alpha
  cut for direct validation.
- Generator balance is projected using only transverse x and y components.
  Earlier files accidentally included the generator jet and Z longitudinal
  components in the dot product; those files require a new event-processing
  pass before their generator-response profile can be used.

These changes require a new event-processing pass before they can appear in a
compatibility file.

## Previous residual correction audit

Both implementations recompute the Run 2024I V11M data JEC from raw jet pT.
`ZbAnalysis::applyJEC` stores the final residual factor as the ratio of the
last two cumulative subcorrections and passes its inverse to
`L2ResHistograms::fill`. The local analysis stores the same inverse quantity
directly as `jetInverseResidual`. The basic correction convention is therefore
the same.

The previous compatibility writer nevertheless introduced a method mismatch:
legacy response profiles were selected from `legacy/profiles1d`, but
`data/l2res/p2res` was always copied from the signed all-pairs analysis. The
legacy response was consequently multiplied later by a previous-JEC profile
formed from a different jet population and with transverse sideband
subtraction. New event output now contains `legacy/l2res/p2res`, `p2respf` and
`p2restc` (and the corresponding full-JEC profiles and counts).
`writeJecsys3.C(..., useLegacyMethod=true)` overlays these method-specific
objects in the compatibility file. Old event files remain readable but carry
an explicit warning metadata object when this directory is absent.

For exact synchronization, the legacy `p2jes*` controls use the original
inverse NanoAOD JEC, `1 - Jet_rawFactor`, saved before the local raw-pT JEC
recomputation. The normal all-pairs controls continue to describe the locally
recomputed correction. This follows the reference implementation and prevents
the legacy diagnostic itself from changing when the input JEC is recomputed.

The all-pairs `p2res*` profiles intentionally use the same signed
parallel-minus-transverse weights as the response. This is the correct
first-order average correction for the subtracted population. The new truth
controls also store the event-correlated products DB/residual and
MPF1/residual, which permit testing the approximation made by multiplying a
non-linear HDM result by an independently averaged correction.

Two limitations remain in the downstream JEC chain:

- The historical `ZbAnalysis` unsuffixed `p2res` is filled versus Z pT even
  though its title says average-projection pT. The local legacy profile
  reproduces that behavior for synchronization. The all-pairs unsuffixed
  profile is genuinely binned in average pT.
- `reprocess.C` currently reduces `p2res` to one central `presz` and uses the
  resulting `herr_l2l3res` for all Z+jet reference-pT variants. It therefore
  cannot simultaneously represent the Z-pT, jet-pT and average-pT conditional
  previous correction. In addition, multiplying the final non-linear HDM
  response by a mean inverse residual neglects correlations with MPF1, MPFn
  and MPFu.

The compatibility output now retains `residual_zmmjet_a100`,
`residual_jetpt_a100`, and `residual_ptave_a100` from the native profiles.
Only the first is consumed by the current `reprocess.C`; the other two make a
future axis-consistent iterative treatment possible without rerunning the
event analysis.

Iteration is now method-consistent for the legacy Z-pT synchronization path,
provided the correction derived from one iteration is used as the input JEC
of the next event pass. Per-axis convergence at the per-mille level still
requires carrying `p2restc`, `p2respf` and `p2res` separately through
`reprocess.C` and applying the previous correction at component level before
the HDM equation.

## Intentional analysis difference

The normal all-pairs analysis reconstructs and requires Run-3 Tight Jet ID for
probe jets because JMENANOv15 does not store `Jet_jetId`. The synchronized
legacy control deliberately omits this requirement. This difference should be
kept explicit when legacy and nominal outputs are compared.

## Statistics normalization

The historical `jecdata2024I_nib1.root` MC statistics integral is about 258.
The pre-alpha legacy file contained about 20.5 million raw selected events;
the synchronized alpha requirement reduces that number. This normalization
difference does not affect profile means, but it makes raw statistics overlays
and the copied `ratio/counts_*` objects misleading. `compareJECdata.C`
therefore uses unit-area count shapes, derives a normalized data/MC shape
ratio, and reports raw integrals separately.
