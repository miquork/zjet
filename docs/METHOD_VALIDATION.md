# Z+jet method validation

`compareMethods.C` is the quantitative validation entry point for the
synchronized leading-jet method and the pileup-subtracted all-pairs method. It
uses the three copied `jecdata` files by default:

```text
rootfiles/jecdata2024I_nib1.root
rootfiles/jecdata2024I_nix_legacy.root
rootfiles/jecdata2024I_nix_newmethod.root
```

Run the ROOT comparison and build the fixed-layout Beamer deck from the
repository root:

```bash
root -l -b -q 'compareMethods.C()'
latexmk -pdf -interaction=nonstopmode -halt-on-error \
  -output-directory=output/pdf compareMethods.tex
```

The macro writes its plots and `method_metrics.tsv` to
`output/compareMethods/`. The metrics file contains the direct HDM difference
and an exact, order-independent Shapley decomposition into data MPF, MPFn,
MPFu, and the corresponding three MC inputs. The six contributions close to
the direct difference bin by bin; this is an algebraic attribution and not a
claim that the inputs are statistically independent.

## HDM definition and closure

The active `softrad3.C` Z+jet equation is

```text
R_HDM = (R_MPF - r_n - r_u) / (1 - r_n/R_n - r_u/R_u),
R_n = 1.00, R_u = 0.92.
```

`writeJecsys3.C` evaluates this equation separately for data and MC and stores
their ratio. The generated synchronized-legacy and new-method compatibility
histograms were checked against their copied 2024I `jecdata` HDM histograms in
all three reference-pT binnings. The largest absolute numerical difference was
below `2.3e-10`.

## Current quantitative result

For average-pT binning, the residual legacy-minus-baseline HDM differences
below 50 GeV are between `-0.533` and `+0.918` per mille. In the first bin
(17.5 GeV), the exact decomposition is dominated by the MC total-MPF input:
`+0.934` per mille, while the five remaining signed contributions sum to about
`-0.016` per mille. This localizes the remaining synchronization discrepancy
to the selected MC MPF population rather than to the HDM soft-recoil algebra.

The new-minus-legacy HDM shift is much larger and pT dependent. In average-pT
bins it is `+11.3`, `-3.64`, `+5.38`, `+9.29`, `+10.54`, `+11.16`, and `+8.61`
per mille from 17.5 through 47.5 GeV. Between 30 and 50 GeV the dominant term
is the data total-MPF change, with a smaller, partly compensating MC total-MPF
change and soft-recoil terms at the few-per-mille level. This demonstrates
that the observed method difference is not created by the HDM combination,
but it does not by itself prove that the all-pairs estimator is unbiased.

## Generator-recoil controls for the next event pass

The analysis now books, for all, truth-matched, and pileup jets:

- parallel, transverse, and background-subtracted DB, MPF, MPF1, MPFn, MPFu,
  and MPFnu versus Z pT for `|eta(jet)| < 1.305`;
- the same controls versus absolute jet eta for `15 < pT,Z < 30 GeV`;
- legacy rho profiles at alpha maxima 0.10, 0.15, 0.20, and 0.30.

`compareMethods.C` derives HDM bin by bin from the MPF, MPFn, and MPFu truth
controls. It also derives a diagnostic linear alpha-maximum-to-zero intercept
from the four nested alpha selections. Because those selections are strongly
correlated, the displayed intercept error is the fit spread and not an
independent-sample statistical uncertainty.

The `truth_hdm/` directory now stores the following independently for the
parallel, transverse and signed parallel-minus-transverse samples, in Z-pT,
jet-pT and average-pT bins:

- selected, reco-gen matched and complete generator-recoil counts;
- match fraction, match DeltaR, generator-jet pT, reco/gen pT, gen/bin pT and
  reco/bin pT;
- generator-Z over reconstructed-Z pT;
- reconstructed MPF1, MPFn and MPFu for all selected pairs and separately
  for the same complete truth-matched population used in the ratios;
- generator MPF1, MPFn and MPFu projected both on the reconstructed-Z axis
  and on the independently reconstructed generator-Z axis;
- reconstructed-times-generator and generator-squared component moments;
- the inverse previous residual correction, plus DB and MPF1 with the selected
  jet's previous residual removed.

The generator Z is built from status-one generator muons matched to the two
selected reconstructed muons within DeltaR 0.1. Generator HT uses generator
jets above 15 GeV cleaned from those muons, and the invisible component comes
from GenMET. The generator decomposition otherwise mirrors the reconstructed
`met1`, `metn` and `metu` definitions.

Derived `TGraphErrors` named `response_r1_reco_axis`,
`response_rn_reco_axis` and `response_ru_reco_axis` contain ratios of profile
means over an identical truth-matched event population in numerator and
denominator. This is deliberate: an event-wise `reco/gen` recoil-component
ratio is singular and biased when the signed generator projection is close to zero.
The corresponding `closure_*_gen_axis` graphs expose generator-Z resolution
and final-state-radiation effects.

The additional `slope_r1_reco_axis`, `slope_rn_reco_axis`, and
`slope_ru_reco_axis` graphs evaluate
`<reco component * gen component>/<gen component squared>`. They provide a
zero-intercept response estimate that remains meaningful when the signed mean
of a soft component is close to zero. Both the ratio-of-means and moment-slope
inputs are retained so their model dependence can be tested explicitly. The
simple graph errors do not include the numerator-denominator covariance; use
the stored profiles or resampling for precision uncertainties.

The synchronized leading-jet controls are stored in
`legacy/truth_hdm/parallel/`. A new event pass is required before these inputs
can be used in the comparison deck.
