#include "zjet.h"

#include <TChain.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TSystem.h>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <string>

R__LOAD_LIBRARY(CondFormats/JetMETObjects/src/Utilities_cc.so)
R__LOAD_LIBRARY(CondFormats/JetMETObjects/src/JetCorrectorParameters_cc.so)
R__LOAD_LIBRARY(CondFormats/JetMETObjects/src/SimpleJetCorrector_cc.so)
R__LOAD_LIBRARY(CondFormats/JetMETObjects/src/FactorizedJetCorrector_cc.so)
R__LOAD_LIBRARY(zjet_C.so)

namespace {

int addInputList(TChain *chain, const char *inputList) {
  std::ifstream stream(inputList);
  if (!stream.is_open()) {
    std::cerr << "Could not open input list " << inputList << std::endl;
    return 0;
  }

  int added = 0;
  std::string line;
  while (std::getline(stream,line)) {
    const std::string::size_type first = line.find_first_not_of(" \t\r");
    if (first==std::string::npos || line[first]=='#') continue;
    const std::string::size_type last = line.find_last_not_of(" \t\r");
    added += chain->AddFile(line.substr(first,last-first+1).c_str());
  }
  return added;
}

} // namespace

void run_zjet_job(const char *inputList, bool isMC, const char *outputFile,
                  const char *goldenJson="", const char *lumiPileup="",
                  const char *pileupWeights="", const char *jecL2="",
                  const char *jecResidual="") {
  TChain *chain = new TChain("Events","Events");
  const int added = addInputList(chain,inputList);
  if (added<=0) {
    std::cerr << "ERROR: no input files were added from " << inputList
              << std::endl;
    gSystem->Exit(2);
    return;
  }

  std::cout << "Added " << added << " file(s) for the "
            << (isMC ? "MC" : "data") << " batch job." << std::endl;
  zjet analysis(chain,isMC,outputFile,goldenJson,lumiPileup,pileupWeights,
                jecL2,jecResidual);
  analysis.Loop();

  TFile check(outputFile,"READ");
  TH2 *inclusiveCounts =
    dynamic_cast<TH2*>(check.Get("l2res/h2ptetatc"));
  TH1 *flavorCounts =
    dynamic_cast<TH1*>(check.Get("flavor/counts_gii"));
  TH1 *nativeCounts =
    dynamic_cast<TH1*>(check.Get("profiles1d/zmmjet/statistics_rmpf"));
  TH1 *nativeMpfNu =
    dynamic_cast<TH1*>(check.Get("profiles1d/zmmjet/rmpfjetnu"));
  TH1 *legacyCounts = dynamic_cast<TH1*>(
    check.Get("legacy/profiles1d/zmmjet/statistics_rmpf"));
  TH1 *legacyMpfNu =
    dynamic_cast<TH1*>(check.Get("legacy/profiles1d/zmmjet/rmpfjetnu"));
  bool flavorCountsClose = false;
  if (inclusiveCounts && flavorCounts) {
    const int firstEtaBin =
      inclusiveCounts->GetXaxis()->FindFixBin(0.+1.e-6);
    const int lastEtaBin =
      inclusiveCounts->GetXaxis()->FindFixBin(1.305-1.e-6);
    const double inclusiveIntegral = inclusiveCounts->Integral(
      firstEtaBin,lastEtaBin,1,inclusiveCounts->GetNbinsY());
    const double flavorIntegral = flavorCounts->Integral();
    flavorCountsClose =
      (std::fabs(inclusiveIntegral-flavorIntegral) <=
       1.e-9*std::max(1.,std::fabs(inclusiveIntegral)));
  }
  if (check.IsZombie() || !check.Get("control/h_cutflow") ||
      !flavorCounts || !check.Get("zjet_flavor_definition") ||
      !nativeCounts || !nativeMpfNu || !legacyCounts || !legacyMpfNu ||
      !check.Get("legacy/control/h_cutflow") ||
      !check.Get("zjet_synchronized_selection") ||
      !flavorCountsClose) {
    std::cerr << "ERROR: output validation failed for " << outputFile
              << std::endl;
    gSystem->Exit(3);
    return;
  }
  std::cout << "Validated output " << outputFile
            << ": flavor profiles, native 1D all-pairs profiles, and "
            << "synchronized legacy profiles are present." << std::endl;
}
