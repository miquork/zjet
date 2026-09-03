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

The analysis also writes a diagnostic HDM variant using the merge-safe
zero-intercept MC regression

```text
Ru_slope = <mu_reco * fu_gen> / <fu_gen^2>.
```

It is evaluated independently in every reference-pT interval and applied to
both data and MC component means. It does not replace the nominal constant:
the fitted flavor shifts are written to `response_ru_slope_impact.tsv` and to
separate ROOT graphs.

## Near-cone and wide-angle radiation controls

New productions retain a jet-axis decomposition intended to test whether the
apparent pT dependence of the gluon-to-quark recoil ratio contains a growing
wide-angle ISR contribution. For each accepted probe jet, every other stored
jet with positive pT is lepton cleaned and assigned by its distance from the
probe axis:

- `near`: `0.4 <= DeltaR(other jet, probe) < 1.0`, an out-of-cone-FSR-enriched
  annulus;
- `wide`: `DeltaR(other jet, probe) >= 1.0`, a wide-angle/ISR-enriched control;
- `hard`: `pT > 15 GeV`, matching the reconstructed Type-I/HT boundary;
- `soft`: `0 < pT <= 15 GeV`, a resolved-jet proxy for the part of the
  unclustered recoil below that boundary.

The selected probe is excluded. Each other-jet momentum is projected onto the
same axis used by the parallel or transverse response estimator and normalized
to pT,Z. The following profiles are stored for every reference-pT convention:

```text
p3rad{near,wide}{hard,soft}{raw,ue}<variant>_flavormatrix
p3rad{near,wide}{hard,soft}<variant>_flavormatrix
p3rad{near,wide}{raw,ue}<variant>_flavormatrix
p3rad{near,wide}<variant>_flavormatrix
p3rad{near,wide}{hard,soft}count<variant>_flavormatrix
```

`raw` is the measured projected jet-axis recoil, `ue` is the corresponding
projected `rho * Jet_area` estimate, and the name without a suffix is
`raw - ue`. The combined near/wide quantities are the sums of their hard and
soft parts. Keeping all three forms is deliberate: PUPPI already suppresses
pileup and a global FastJet rho need not be an unbiased estimate of the UE
inside a selected jet population. The subtraction must therefore be tested
against data and particle-level closure rather than assumed.

For MC, matching particle-level profiles are written as

```text
p3genrad{near,wide}{hard,soft}<variant>_flavormatrix
p3genrad{near,wide}<variant>_flavormatrix
```

using the generator Z projection, excluding the selected GenJet and generator
jets overlapping the selected generator muons. The current custom JMENANO
stores GenJets down to 3 GeV. These objects expose the 3--15 GeV resolved
component but cannot recover truly diffuse particles below the GenJet
threshold.

The near and wide regions are operational enrichments, not truth-level ISR and
FSR labels. Initial--final color interference, Born-channel composition,
acceptance, and the finite jet threshold all remain. A useful closure sequence
is therefore raw reconstructed activity, rho-area-subtracted activity, and the
direct GenJet analogue, followed by the near/wide pT dependence and their
gluon-to-quark ratios.

For a conditional data constraint, new productions store the event-level
closure proxy

```text
fu_closure = 1 - m2/R2 - mn/Rn,
mu * fu_closure,
fu_closure^2,
```

with `R2 = Rn = 1`. The corresponding zero-intercept slope can be formed in
data, but is not identifiable without these response assumptions. The same
proxy and the true generator regression are saved in MC so resolution
dilution and flavor-dependent bias can be measured explicitly.

`FlavorMatrix/taggerAudit` stores counts and response-component profiles as a
continuous function of the DeepJet, ParticleNet, and UParT QvG scores. This
allows tagger response sculpting to be compared without committing to a
single working point. The current NanoAOD reader has no legacy QGL branch, so
that comparison requires a compatible input production or an optional branch
extension. Heavy-topology and jet-muon-fraction response profiles provide the
first separation of the low-pT b-jet feature into double-heavy production,
semileptonic, and detector/selection components.

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
