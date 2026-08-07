#include "zjet.h"

#include <TChain.h>
#include <TFile.h>
#include <TSystem.h>

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
  if (check.IsZombie() || !check.Get("control/h_cutflow")) {
    std::cerr << "ERROR: output validation failed for " << outputFile
              << std::endl;
    gSystem->Exit(3);
    return;
  }
  std::cout << "Validated output " << outputFile << std::endl;
}
