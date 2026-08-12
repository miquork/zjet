# zjet
Z+jet response measurement at low pT that is robust against pileup 

Run as follows:

```bash
root -l -b -q mk_compile.C
root -l -b -q mk_zjet.C
root -l -b -q drawControl.C
root -l -b -q drawControlData.C
root -l -b -q writeJecsys3.C
```

The default Run2024I/Summer24 setup applies the nominal Summer24 muon scale and
MC resolution, recomputes JEC from raw jet pT, applies MC JER smearing and the
data-only 2024 jet veto map, and propagates the corrected jets to Type-I PUPPI
MET. It writes both the raw response profiles and the central-eta input
contract for the `jecsys3` reprocess/soft-radiation/global-fit chain. See
`RUNNING.md` for local, XRootD and HTCondor workflows.
