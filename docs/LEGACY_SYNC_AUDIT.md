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
- The leading lepton-cleaned jet has pT above 12 GeV and `|eta| < 5`.
- The back-to-back residual requirement is below 0.44 radians.
- Summer24 V2 MC L2Relative and 2024I V11M data L2L3Residual corrections are
  recomputed from raw jet pT.
- Pileup reweighting is disabled in both nominal workflows. It is not the
  source of the current legacy difference.

## Corrections made in this audit

- The legacy MPF1 and MPFn four-vectors are now explicitly flattened to the
  transverse plane before projection. Without this step their longitudinal
  components produced unphysical values as large as O(100), although the
  components cancelled in total MPF.
- The legacy leading jet is now chosen before applying tight Jet ID, as in
  `ZbAnalysis`; an event is rejected if that leading jet fails the ID. Jets
  entering the type-1 MET and non-leading recoil sums are no longer required
  to pass the leading-jet ID.
- The 2024 forward spike veto used by `ZbAnalysis` is reproduced.
- The selected-muon multiplicity is restricted to two or three.
- The closest-mass pair is chosen from all tight, isolated, trigger-matched
  muons above 8 GeV before applying the selected pair's eta and 20/10 GeV pT
  cuts. Jets are cleaned against all two or three selected muons, not only the
  pair assigned to the Z boson.
- The `ZbAnalysis` MET-filter list and its data-only application are used.
- Legacy MC profiles use unit event weights, matching `ZbAnalysis`. The new
  all-pairs method retains signed generator weights.

These changes require a new event-processing pass before they can appear in a
compatibility file.

## Remaining differences requiring implementation

1. **JER smearing (high priority).** `ZbAnalysis` enables smearing by default
   for MC. For 2024 it uses
   `JR_Winter22Run3_V1_MC_PtResolution_AK4PFPuppi.txt` and
   `Prompt24_2024_nib_JRV11M_MC_SF_AK4PFPuppi.txt`. Only the first three
   lepton-cleaned jets in NanoAOD order are smeared before pT reordering. The
   local analysis currently records JER smearing as disabled. This affects
   reference-pT migrations and is the most direct missing correction for the
   MC response.

2. **Muon scale and resolution (high priority for exact event sync).**
   `ZbAnalysis` applies `2024_Summer24.json` to data and MC muons before the
   final pT, eta, mass, and pair-ranking decisions. The local analysis uses
   stored muon kinematics. The response effect should be smaller than the JER
   effect, but events near thresholds and pT-bin edges will differ.

3. **Jet veto map (data only).** `ZbAnalysis` applies
   `jetvetoReReco2024_V9M.root:jetvetomap` to the leading data jet and does not
   apply a 2024 map to MC. The local analysis has no jet-veto-map interface.
   The direct effect on the central `|eta| < 1.305` input may be modest, but it
   is required for an event-by-event synchronization.

4. **Jet ID source.** `ZbAnalysis` reads `Jet_jetId >= 4`; JMENANOv15 inputs
   used locally do not expose that branch, so the local code reconstructs the
   Run-3 tight PF Jet ID from jet fractions and multiplicities. This should be
   validated event by event against the stored bit in a file that contains
   both representations.

5. **Type-1 MET after JER.** Once JER is added, the smeared jet momenta must be
   propagated through the raw-PUPPI-MET correction exactly as in
   `ZbAnalysis::recalculateMET`. Adding JER only to the leading-jet response
   would leave MPF inconsistent.

## Statistics normalization

The historical `jecdata2024I_nib1.root` MC statistics integral is about 258,
whereas the new legacy file contains about 20.5 million raw selected events.
This normalization difference does not affect profile means, but it makes raw
statistics overlays and the copied `ratio/counts_*` objects misleading.
`compareJECdata.C` therefore uses unit-area count shapes, derives a normalized
data/MC shape ratio, and reports raw integrals separately.
