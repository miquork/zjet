# JEC data comparison

`compareJECdata.C` compares the historical baseline, the synchronized legacy
method, and the pileup-subtracted method in a reproducible fixed-layout
validation deck. It covers the three Z+jet reference-pT choices (`jetz`,
`zjav`, and `zjet`) for data, simulation, and their derived result.

Generate the control plots, comparison plots, and Beamer frame list from the
repository root:

```bash
root -l -b -q 'drawPileupControl.C()'
root -l -b -q 'drawControl.C()'
root -l -b -q 'drawControlData.C()'
root -l -b -q 'drawZjetSummary.C()'
root -l -b -q 'compareJECdata.C()'
```

The comparison macro appends whichever generated control PDFs are present.
Running the four control macros first therefore produces the complete deck;
missing optional control PDFs do not prevent a jecdata-only comparison.

The defaults read:

```text
../jecsys3/rootfiles/jecdata2024I_nib1.root
../jecsys3/rootfiles/jecdata2024I_nix_legacy.root
../jecsys3/rootfiles/jecdata2024I_nix_newmethod.root
```

Explicit inputs and labels can be supplied as the first seven macro arguments.
The optional eighth and ninth arguments set `Rn` and `Ru`; their defaults are
`1.00` and `0.92`, matching the active `softrad3.C` settings. Because jecdata
does not store the MPFn--MPFu covariance, the displayed uncertainty for the
derived

```text
MPFnu = MPFn/Rn + MPFu/Ru
```

is propagated without that covariance.
Build the deck with:

```bash
latexmk -pdf -interaction=nonstopmode \
  -output-directory=output/pdf compareJECdata.tex
```

The lower pads contain a difference only when the baseline and candidate graph
points have exactly the same pT center. This is intentional: interpolation
would hide an accidental binning change. `summary.tsv` records the point and
matching counts for every panel.

The historical MC statistics histogram uses a normalization that is not
comparable to the raw event-count convention of the new files. Statistics
panels therefore compare unit-area shapes. Their third column is derived as
the ratio of separately normalized data and MC distributions; it does not use
the copied `ratio/counts_*` object from `reprocess.C`. Raw integrals are retained
in `normalization.tsv` and on the deck's normalization-diagnostics page.

For DB and total MPF, the `ratio` directory contains data/MC. For MPF1, MPFn,
MPFu, PF fractions, and rho, `reprocess.C` stores data minus MC. The plot labels
and fixed y-axis ranges follow those definitions. The final section overlays
jet-pT, average-pT (HDM), and Z-pT binnings in the same panels for each input.
`rjet` and `gjet` are MC-only closure objects by construction in `reprocess.C`,
so the deck no longer draws empty data and data/MC panels for them. The data,
MC, and data/MC direct-balance response is the `ptchs` observable.

The compatibility writer now has two method controls after the historical
placeholder argument:

```cpp
writeJecsys3(dataFile, mcFile, outputFile,
             addFlavorPlaceholders, useLegacyMethod,
             preferOneDimensional);
```

The default uses the pileup-subtracted all-pairs method and native 1D pT bins.
Set `useLegacyMethod=true` to select the synchronized leading-jet control stored
under `legacy/`. Set `preferOneDimensional=false` only for reproducing the old
central-eta projection of the L2Res 2D histograms.
