# Compact hybrid-tagger flavor matrix

`zjet.C` writes the new all-pairs flavor study under `FlavorMatrix/`.  The
objects are filled only for probe jets with `|eta| < 1.3`.  They complement the
legacy `flavor/` directory; no existing jecsys3 input is replaced.

## Axes and tag definition

The first axis uses the requested reference-pT bin edges

```text
10, 15, 20, 25, 30, 35, 40, 50, 60, 70, 85, 100, 125, 155,
180, 210, 250, 300, 350, 400, 500, 600, 800, 1000, 1200,
1500, 1800, 2100, 2400, 2700, 3000, 3500, 4000 GeV.
```

The reconstructed category uses UParTAK4 CvB and CvL together with the more
stable ParticleNet QvG output. It is encoded on the second axis:

| ID | Category | Selection |
|---:|---|---|
| 0 | undefined | any of UParT CvB, UParT CvL or PNet QvG is negative or non-finite |
| 1 | uds | `UParT CvB >= 0.5`, `UParT CvL < 0.5`, `PNet QvG >= 0.3` |
| 2 | reserved d | empty in the initial combined-uds definition |
| 3 | reserved s | empty in the initial combined-uds definition |
| 4 | c | `UParT CvB >= 0.5`, `UParT CvL >= 0.5` |
| 5 | b | `UParT CvB < 0.5` |
| 6 | g | `UParT CvB >= 0.5`, `UParT CvL < 0.5`, `PNet QvG < 0.3` |

ParticleNet QvG is quark (`udsbc`) versus gluon, whereas the UParTAK4 node is
documented as `uds` versus gluon. The b and c regions are applied before QvG,
so this difference affects only b/c leakage into the residual light-flavor
split. The parallel control `h3_upartqvg_pnetqvg_trueflavor` retains both
outputs for direct validation.

The third axis is the absolute NanoAOD `Jet_partonFlavour`.  The d and u
partons are combined in ID 1, while s, c, b and g use IDs 3, 4, 5 and 6.
Unclassified MC jets use ID 0.  Data have no generator flavor and therefore
always use ID 0.

## Stored objects

- `h3counts_flavormatrix` is the signed parallel-minus-transverse estimator in
  Z-pT bins.  Each of the two transverse windows carries weight `-0.5`.
- `h3counts_parallel_flavormatrix` and
  `h3counts_transverse_flavormatrix` retain the corresponding un-subtracted
  region populations.  Negative generator weights can still occur in MC.
  These are the stable initial inputs for the tagging-efficiency inference.
- `p3{m0,m2,mn,mu,mnu,hdm}{ab,ad,tc,pf}_flavormatrix` stores all six HDM
  components in four reference-pT conventions: bisector, arithmetic average,
  Z pT and probe-jet pT.  The response projection is unchanged; only the
  x-axis binning variable changes.
- `p3{m0,m2,mn,mu,mnu,hdm}{ab,ad,tc,pf}_parallel_flavormatrix` stores the
  corresponding pure-parallel barrel response. The local flavor-response fit
  uses this payload while the signed profiles remain available for pileup and
  forward-region studies.
- `FlavorMatrix/controls/h3_cvb_cvl_trueflavor`, `h3_cvb_qvg_trueflavor` and
  `h3_cvl_qvg_trueflavor` store pairwise discriminator densities for the
  un-subtracted parallel sample with probe pT above 30 GeV. The QvG coordinate
  in these objects is ParticleNet QvG.
- `FlavorMatrix/controls/h3_cvb_cvl_qvg_true<ID>` stores the corresponding
  three-dimensional discriminator cube separately for every truth ID.
- `h3counts_heavytopology` and the three
  `FlavorMatrix/controls/h3_*_heavytopology` objects split the same population
  by matched-GenJet heavy-hadron topology. The codes are no-heavy=0,
  single-c=1, double-c=2, reserved-other=3, single-b=4, double-b=5 and no
  GenJet match=6. Bottom has precedence over charm because charm hadrons from
  bottom decays are common.
- `FlavorMatrix/controls/h3_genjet_nc_nb_trueflavor` stores the underlying
  `GenJet_nCHadrons` versus `GenJet_nBHadrons` multiplicities. Double-heavy
  jets are gluon-splitting-enriched controls, not exclusive production-origin
  labels. A direct-versus-splitting measurement needs an additional GenPart or
  LHE ancestry refinement and must retain a category for truncated ancestry.

The `hdm` profile uses the diagnostic constants `Rn = 1.00` and `Ru = 0.92`.
It is never averaged event by event.  Each worker derives it from that file's
component means for local inspection, and `embedCampaignMetadata.C` replaces
it after `hadd` with `HDM(<m0>,<mn>,<mu>)` from the merged components.  This is
required because the HDM equation is nonlinear.  The raw `m0`, `mn` and `mu`
components remain authoritative when alternative recoil responses or a full
covariance treatment are needed.

## Data tagging inference

Data reco-category fractions do not uniquely determine every truth-to-reco
tagging efficiency.  `analyzeFlavorMatrix.C` therefore reports a constrained,
model-dependent estimate:

- the empty reconstructed-undefined tag is omitted;
- generator d/u and s are combined into `uds`;
- the generator-undefined population is combined with gluons because its
  discriminator shapes are gluon-like.

This leaves four truth groups for four measured reconstructed tags and removes
the prior-only null modes of the earlier six-column smoke fit. The inference
then proceeds as follows:

1. build the MC joint distribution `W_MC(t,f) = pi_f A_MC(t|f)`;
2. keep the MC truth fractions `pi_f` fixed;
3. constrain the reconstructed marginal to the observed data fractions;
4. find the non-negative joint distribution nearest to `W_MC` in KL distance
   with iterative proportional fitting;
5. define the inferred data efficiency `A_data(t|f) = W_data(t,f)/pi_f` and
   the transition scale factor `SF(t,f) = A_data/A_MC`.

The resulting data matrix is an inference under the MC-odds and truth-fraction
assumptions, not a direct measurement.  Per-flavor response residuals are then
obtained from the tagged data responses and corrected MC cell responses with a
regularized SVD solve.  The output records the condition number so an
underconstrained or sparse smoke test is visible rather than silently treated
as a precision result.

The IPF fit and the current response fit both use the un-subtracted parallel
barrel population. IPF requires non-negative cell probabilities, and the
barrel pileup-jet contamination is small enough for this to be the more stable
one-file starting point. The signed parallel-minus-transverse profiles remain
available for pileup controls and forward-region studies. A precision
iteration should still compare this choice with a simultaneous non-negative
signal/background likelihood. The ROOT metadata and JSON summary record the
choice explicitly.

## Local validation

After producing one data and one MC output, run

```bash
root -l -b -q 'analyzeFlavorMatrix.C("rootfiles/zjet_DATA.root","rootfiles/zjet_MC.root")'
root -l -b -q 'drawFlavorMatrix.C("rootfiles/zjet_MC.root","rootfiles/zjet_DATA.root","output/flavorMatrix/flavorMatrixAnalysis.root")'
```

The quantitative ROOT output and machine-readable tables are written below
`output/flavorMatrix/`.  The plotting macro writes fixed-axis PDF and PNG
figures in the same directory.  A one-file run is a schema and numerical
stability test only; final efficiencies and response residuals require the
full sample and uncertainty propagation.

The vector validation deck is maintained as `flavorMatrixAnalysis.tex`. After
regenerating the ROOT plots, rebuild it with

```bash
latexmk -pdf -interaction=nonstopmode -halt-on-error \
  -outdir=output/pdf flavorMatrixAnalysis.tex
```

The deck deliberately rejects pT slices whose four-flavor response system is
rank deficient or has condition number at least 100. The TSV tables retain
all slices for diagnosis even when they are omitted from the summary graph.

## Metadata key cycles

Flavor axes intentionally have no `SetBinLabel` values in worker files.  ROOT
otherwise performs label-aware histogram merging and can append another seven
bins for every input file.  Plotting code adds human-readable labels only to
detached display copies.

Worker output writes each top-level `TObjString` only once.  After `hadd`,
`embedCampaignMetadata.C` requires every cycle of a given metadata key to have
an identical value, aborts on conflicts, removes all duplicate key cycles and
writes one copy before adding the campaign record.  This keeps the TBrowser
top level readable without hiding inconsistent worker configurations.
