#ifndef ZJET_LEGACY_REPLAY_H
#define ZJET_LEGACY_REPLAY_H
#include <TTree.h>
#include <TDirectory.h>
#include <TObjString.h>
#include <vector>
// One accepted barrel leading jet per event. No tag selection is applied here.
// Keep raw scores and linear response components for arbitrary offline WPs/SFs.
struct ZJetLegacyRecord {
  UInt_t run=0,lumi=0;ULong64_t event=0;Bool_t isMC=false;
  Int_t parton=0,hadron=0,jetIndex=-1;
  Float_t ptz=0,jetpt=0,eta=0,phi=0,mass=0,rho=0,npv=0,
    b=0,cvb=0,cvl=0,pnet=0,deep=0,upart=0,m0=0,m2=0,mn=0,mu=0,db=0,
    inverseResidual=0,muEF=0,recoGen=0,metParallel=0,metPerp=0;
  Double_t genWeight=1,signedWeight=1;
  std::vector<float> otherB,otherCvB,otherCvL,otherQ;
  void bind(TTree *t,bool writing) {
#define BR(name,type) if(writing)t->Branch(#name,&name,#name "/" type);else t->SetBranchAddress(#name,&name);
    BR(run,"i") BR(lumi,"i") BR(event,"l") BR(isMC,"O")
    BR(parton,"I") BR(hadron,"I") BR(jetIndex,"I")
    BR(ptz,"F") BR(jetpt,"F") BR(eta,"F") BR(phi,"F") BR(mass,"F") BR(rho,"F") BR(npv,"F")
    BR(b,"F") BR(cvb,"F") BR(cvl,"F") BR(pnet,"F") BR(deep,"F") BR(upart,"F")
    BR(m0,"F") BR(m2,"F") BR(mn,"F") BR(mu,"F") BR(db,"F") BR(inverseResidual,"F")
    BR(muEF,"F") BR(recoGen,"F") BR(metParallel,"F") BR(metPerp,"F")
    BR(genWeight,"D") BR(signedWeight,"D")
#undef BR
    // Reading variable-length vectors uses TTreeReader in postprocessing.
    if(writing) {
      t->Branch("otherB",&otherB);t->Branch("otherCvB",&otherCvB);
      t->Branch("otherCvL",&otherCvL);t->Branch("otherQ",&otherQ);
    }
  }
};
struct ZJetLegacyReplay {
  ZJetLegacyRecord row;TTree *tree;
  explicit ZJetLegacyReplay(TDirectory *out) {
    TDirectory::TContext context;out->mkdir("LegacyFlavor")->cd();
    TObjString("v1: one lepton-cleaned accepted legacy barrel leading jet/event before tagging; "
      "B=Jet_btagUParTAK4B, CvB/CvL=UParT, PNet/DeepFlav/UParT raw QvG branches; "
      "nominal reference B Tight, CvL Tight only, PNet>=0.45; no SF applied. "
      "Legacy unit weights reproduced by weight=1; full genWeight and signedWeight kept separately. "
      "other* scores: additional lepton-cleaned jets pT>30, abs(eta)<2.5, no JetID. "
      "Truth labels are raw signed parton flavor and hadron flavor. "
      "Only accepted Legacy kinematics are retained; event/jet selection cannot be loosened offline.").Write("definition");
    tree=new TTree("events","Retaggable Legacy flavor probes");tree->SetAutoFlush(-16000000);
    tree->SetAutoSave(0);row.bind(tree,true);
  }
  ZJetLegacyReplay(const ZJetLegacyReplay&)=delete;
};
#endif
