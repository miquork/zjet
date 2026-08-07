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

The default Run2024I/Summer24 setup recomputes JEC from raw jet pT and writes
both the raw response profiles and the central-eta input contract for the
`jecsys3` reprocess/soft-radiation/global-fit chain. See `RUNNING.md` for local,
XRootD and HTCondor workflows.
