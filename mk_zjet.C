#include "zjet.h"

#include <fstream>
#include <string>

R__LOAD_LIBRARY(zjet_C.so)

namespace {

int addInput(TChain *chain, const char *input, int maximumFiles) {
  const std::string name(input ? input : "");
  if (name.size()<4 || name.substr(name.size()-4)!=".txt")
    return chain->AddFile(input);

  std::ifstream stream(input);
  if (!stream.is_open()) {
    std::cerr << "Could not open input list " << input << std::endl;
    return 0;
  }
  int added = 0;
  std::string line;
  while (std::getline(stream,line)) {
    const std::string::size_type first = line.find_first_not_of(" \t\r");
    if (first==std::string::npos || line[first]=='#') continue;
    const std::string::size_type last = line.find_last_not_of(" \t\r");
    if (maximumFiles>=0 && added>=maximumFiles) break;
    added += chain->AddFile(line.substr(first,last-first+1).c_str());
  }
  std::cout << "Added " << added << " file(s) from " << input << std::endl;
  return added;
}

} // namespace

void mk_zjet(
  const char *mcFile="../data/zjet/events_1_DYto2L_4Jets_M_50_Summer24.root",
  const char *dataFile="../data/zjet/events_1_Muon0_Run2024I_JMENANOv15_v2_v1.root",
  const char *goldenJson="",
  const char *lumiPileup="",
  const char *pileupWeights="",
  int maximumMcFiles=-1,
  int maximumDataFiles=-1) {
  
  bool isMC = false;
  TChain *cm = new TChain("Events","Events");
  //cm->AddFile("../data/zjet/DYto2Mu-4Jets_Bin-MLL-50_TuneCP5_13p6TeV_madgraphMLM-pythia8_Summer24.root"); //isMC =  true;
  addInput(cm,mcFile,maximumMcFiles);

  TChain *cd = new TChain("Events","Events");
  //cd->AddFile("../data/zjet/Muon0_Run2025G_PromptReco_v1.root"); //isMC = false;
  addInput(cd,dataFile,maximumDataFiles);

  zjet zm(cm, isMC = true, "rootfiles/zjet_MC.root", "", "",
          pileupWeights);
  zm.Loop();
  
  zjet zd(cd, isMC = false, "rootfiles/zjet_DATA.root", goldenJson,
          lumiPileup, "");
  zd.Loop();

}
