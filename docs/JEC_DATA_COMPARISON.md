# JEC data comparison

`compareJECdata.C` compares two `jecdata*.root` outputs in a reproducible,
fixed-layout validation deck. It covers the three Z+jet reference-pT choices
(`jetz`, `zjav`, and `zjet`) for data, simulation, and their data/MC result.

Generate the plots and the Beamer frame list from the repository root:

```bash
root -l -b -q 'compareJECdata.C()'
```

The defaults compare `2024I_nib1` with `2024I_nix`. Explicit inputs and labels
can be supplied as the five macro arguments. Build the deck with:

```bash
latexmk -pdf -interaction=nonstopmode \
  -output-directory=output/pdf compareJECdata.tex
```

The lower pads contain a difference only when the old and new graph points
have exactly the same pT center. This is intentional: interpolation would hide
an accidental binning change. `summary.tsv` records the number of old, new,
matched, and unmatched bins for every panel.

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
