# Running the all-pairs Z+jet analysis

This note documents the local workflow and the intended small-scale lxplus
validation. It deliberately contains no grid proxy, certificate, private file
list, or other credential material.

## Local build and regression run

Compile the analysis after changing `zjet.C`, `zjet.h`, or `ZJetLumi.h`:

```bash
root -l -b -q mk_compile.C
```

Run the default local Summer24 MC and Run2024I data files:

```bash
root -l -b -q mk_zjet.C
```

The outputs are written directly to `rootfiles/zjet_MC.root` and
`rootfiles/zjet_DATA.root`. The macro no longer creates a temporary output and
renames it with an interactive shell command.

The optional arguments of `mk_zjet` are, in order:

1. MC ROOT file
2. data ROOT file
3. golden-JSON file for data
4. lumisection pileup file for the data `mu` control
5. ROOT file containing the MC pileup-weight histogram

The pileup-weight histogram must be named `pileup_ratio` (preferred) or
`pileup`. The lumisection pileup reader accepts a brilcalc CSV containing an
`avgpu` column, or a whitespace-separated `run ls mu` file. Empty optional
arguments disable only the corresponding correction or selection and print no
credential information.

Example:

```bash
root -l -b -q 'mk_zjet.C("mc.root","data.root","Cert.json","lumi.csv","pu_weights.root")'
```

## Control plots

The original control plots are produced with:

```bash
root -l -b -q drawControl.C
```

Pileup and truth-matching controls are produced with:

```bash
root -l -b -q drawPileupControl.C
```

The latter writes DB and hybrid-MPF profiles versus `PV_npvs`, `rho`, and
`mu`, together with truth-matched versus unmatched jet controls, under
`pdf/drawPileupControl/`.

## jecsys3 compatibility file

Create the combined directory hierarchy expected by `jecsys3/L2Res.C`:

```bash
root -l -b -q writeJecsys3.C
```

The default output is `rootfiles/zjet_JMENANO_compat.root`, containing
`data/l2res`, `data/l2res1`, `mc/l2res`, and `mc/l2res1`. Point both the
`ZMM_<era>_DATA` and `ZMM_<era>_MC` entries in `jecsys3/Config.C` to this same
file when testing the new method. Keeping the old file paths in a separate
configuration makes the old/new comparison reversible.

## lxplus small-file validation

Log in with the CERN account and verify the proxy:

```bash
ssh <cern-user>@lxplus.cern.ch
voms-proxy-info -timeleft
```

After cloning or pulling this repository on lxplus, source a ROOT environment
compatible with the JMENANO files, compile, and first run one MC and one data
file through XRootD. Check that:

- both remote files open before starting the event loop;
- the cutflow is non-zero through the probe-lepton veto;
- `l2res/p2m0tc` and `l2res/p2restc` both have entries;
- the `l2res` eta underflow is empty;
- the residual profile is exactly one in every populated bin;
- the two half-weight transverse windows have the expected normalization.

The exact JSON-to-XRootD file-list conversion will be added after the sample
JSON schema has been inspected. Do not submit the full sample or HTCondor jobs
before the one- and few-file checks pass.
