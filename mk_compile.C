// Compile and load the exact sources transferred to a worker.  The explicit
// library checks turn ACLiC diagnostics into a non-zero batch-job status.
{
  const std::pair<const char*,const char*> builds[] = {
    {"CondFormats/JetMETObjects/src/Utilities.cc",
     "CondFormats/JetMETObjects/src/Utilities_cc.so"},
    {"CondFormats/JetMETObjects/src/JetCorrectorParameters.cc",
     "CondFormats/JetMETObjects/src/JetCorrectorParameters_cc.so"},
    {"CondFormats/JetMETObjects/src/SimpleJetCorrector.cc",
     "CondFormats/JetMETObjects/src/SimpleJetCorrector_cc.so"},
    {"CondFormats/JetMETObjects/src/FactorizedJetCorrector.cc",
     "CondFormats/JetMETObjects/src/FactorizedJetCorrector_cc.so"},
    {"zjet.C", "zjet_C.so"},
  };
  for (const auto &build : builds) {
    const char *source = build.first;
    const char *library = build.second;
    cout << "Compiling " << source << endl << flush;
    // Remove the previous product so a failed rebuild cannot look successful.
    gSystem->Unlink(library);
    gROOT->ProcessLine(Form(".L %s++g",source));
    FileStat_t status;
    if (gSystem->GetPathInfo(library,status)!=0 || status.fSize<=0) {
      cerr << "ERROR: compilation failed for " << source
           << "; missing output " << library << endl;
      gSystem->Exit(20);
      return;
    }
  }
}
