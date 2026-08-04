#include "zjet.h"

R__LOAD_LIBRARY(zjet_C.so)

void mk_zjet() {
  
  bool isMC = false;
  TChain *cm = new TChain("Events","Events");
  //cm->AddFile("../data/zjet/DYto2Mu-4Jets_Bin-MLL-50_TuneCP5_13p6TeV_madgraphMLM-pythia8_Summer24.root"); //isMC =  true;
  cm->AddFile("../data/zjet/events_1_DYto2L_4Jets_M_50_Summer24.root");

  TChain *cd = new TChain("Events","Events");
  //cd->AddFile("../data/zjet/Muon0_Run2025G_PromptReco_v1.root"); //isMC = false;
  cd->AddFile("../data/zjet/events_1_Muon0_Run2024I_JMENANOv15_v2_v1.root");

  //TFileCollection fc;
  //fc.AddFromFile("textfiles/2025G_Muon0.txt"); isMC = false;
  //fc.AddFromFile("textfiles/2025G_Muon1.txt"); isMC = false;
  //fc.AddFromFile("textfiles/Summer24MC.txt"); isMC = true;
  //c->AddFileInfoList(fc.GetList());
  
  zjet zm(cm, isMC = true);
  zm.Loop();
  gROOT->ProcessLine(".! mv -i rootfiles/zjet.root rootfiles/zjet_MC.root");
  
  zjet zd(cd, isMC = false);
  zd.Loop();
  gROOT->ProcessLine(".! mv -i rootfiles/zjet.root rootfiles/zjet_DATA.root");

}
