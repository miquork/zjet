// Retag accepted Legacy probes and change MC weights without rerunning NanoAOD.
// SF table columns: absParton hadron tag ptMin ptMax absEtaMin absEtaMax factor.
// flavor/tag=-1 are wildcards; gluon parton=21. Unspecified cells have factor 1.
// These are final per-jet weights, NOT an automatic pass/fail conversion of a B SF.
#include "ZJetLegacyReplay.h"
#include "ZJetTaggingControls.h"
#include <TFile.h>
#include <TH3D.h>
#include <TProfile3D.h>
#include <fstream>
#include <sstream>
#include <memory>
#include <iostream>
#include <stdexcept>

void replayLegacyFlavor(const char *input,const char *output,const char *bwp="T",
    const char *cwp="T",double qcut=.45,bool officialTwoDimensionalC=false,
    const char *weights="unit",const char *sfTable="",double bcut=-1,double cvlcut=-1,
    double cvbcut=-1,const char *qbranch="pnet") {
  ZJetTaggingControls::Config config;if(!config.enabled)throw std::runtime_error("Missing official WP configuration");
  if((std::string(bwp)!="M"&&std::string(bwp)!="T")||
     (std::string(cwp)!="M"&&std::string(cwp)!="T"))throw std::runtime_error("WP must be M or T");
  const std::string mode=weights,qtype=qbranch;
  if(mode!="unit"&&mode!="gen"&&mode!="signed")throw std::runtime_error("Weight mode: unit, gen or signed");
  if(qtype!="pnet"&&qtype!="deep"&&qtype!="upart")throw std::runtime_error("QvG branch: pnet, deep or upart");
  if(bcut<0)bcut=config.values.at(std::string("b")+bwp);
  if(cvlcut<0)cvlcut=config.values.at(std::string("c")+cwp+"CvL");
  if(cvbcut<0)cvbcut=config.values.at(std::string("c")+cwp+"CvB");
  for(double cut:{bcut,cvlcut,cvbcut,qcut})if(!ZJetResponseAudit::valid(cut))throw std::runtime_error("Invalid tag cut");
  struct SF {int parton,hadron,tag;double pl,ph,el,eh,w;};std::vector<SF> factors;
  std::string tableText;
  if(sfTable&&*sfTable) {
    std::ifstream f(sfTable);if(!f)throw std::runtime_error("Cannot open SF table");
    std::string line;while(std::getline(f,line)) {
      tableText+=line+"\n";if(line.empty()||line[0]=='#')continue;
      SF s;std::istringstream in(line);
      if(!(in>>s.parton>>s.hadron>>s.tag>>s.pl>>s.ph>>s.el>>s.eh>>s.w)||
         !std::isfinite(s.w)||s.w<0||s.pl>=s.ph||s.el<0||s.el>=s.eh)
        throw std::runtime_error("Invalid SF-table row");
      factors.push_back(s);
    }
  }
  std::unique_ptr<TFile> in(TFile::Open(input));if(!in||in->IsZombie())throw std::runtime_error("Cannot open input");
  auto *tree=dynamic_cast<TTree*>(in->Get("LegacyFlavor/events"));
  if(!tree)throw std::runtime_error(in->Get("zjet_compact_definition")
    ? "This is a compact histogram file. For arbitrary retagging, pass the original full EOS merge URL."
    : "Input lacks LegacyFlavor/events");
  ZJetLegacyRecord row;row.bind(tree,false);
  TFile out(output,"CREATE");if(out.IsZombie())throw std::runtime_error("Output exists or cannot be created; choose a fresh filename");
  out.mkdir("LegacyFlavor")->cd();
  std::ostringstream metadata;metadata<<"Legacy replay v1; weight="<<mode<<"; B>="<<bcut<<"; CvL>="<<cvlcut
    <<"; twoDimensionalC="<<officialTwoDimensionalC<<"; CvB>="<<cvbcut<<"; raw "<<qtype<<" QvG>="<<qcut
    <<"; unspecified SF cells=1; no cross-section normalization; no SF is inferred from data.\n"<<tableText;
  TObjString(metadata.str().c_str()).Write("definition");
  const double edges[]={10,15,20,25,30,35,40,50,60,70,85,100,125,155,180,210,250,300,350,400,500,600,800,1000,1200,1500,1800,2100,2400,2700,3000,3500,4000};
  double ids[]={-.5,.5,1.5,2.5,3.5,4.5,5.5,6.5};const int nx=sizeof(edges)/sizeof(double)-1;
  TH3D counts("counts",";Z pT;reco tag;truth parton",nx,edges,7,ids,7,ids);counts.Sumw2();
  std::map<std::string,TProfile3D*> profiles;
  for(const std::string name:{"m0","m2","mn","mu","db","inverseResidual","muEF","recoGen",
      "m0m0","m0mn","m0mu","mnmn","mnmu","mumu"})
    profiles[name]=new TProfile3D(name.c_str(),";Z pT;reco tag;truth parton",nx,edges,7,ids,7,ids);
  Long64_t reweighted=0;
  for(Long64_t i=0;i<tree->GetEntries();++i) {
    if(tree->GetEntry(i)<=0)throw std::runtime_error("Failed replay entry read");
    const double q=qtype=="pnet"?row.pnet:qtype=="deep"?row.deep:row.upart;
    int tag=0;
    if(ZJetResponseAudit::valid(row.b)) {
      if(row.b>=bcut)tag=5;
      else if(ZJetResponseAudit::valid(row.cvl)&&(!officialTwoDimensionalC||ZJetResponseAudit::valid(row.cvb))) {
        if(row.cvl>=cvlcut&&(!officialTwoDimensionalC||row.cvb>=cvbcut))tag=4;
        else if(ZJetResponseAudit::valid(q))tag=q>=qcut?1:6;
      }
    }
    int truth=std::abs(row.parton);if(truth==21)truth=6;
    else if(truth>5)truth=0;if(truth==2)truth=1;
    double weight=!row.isMC||mode=="unit"?1.:mode=="gen"?row.genWeight:row.signedWeight;
    int matches=0;
    if(row.isMC)for(const auto &s:factors)if((s.parton<0||s.parton==std::abs(row.parton))&&(s.hadron<0||s.hadron==row.hadron)&&
        (s.tag<0||s.tag==tag)&&row.jetpt>=s.pl&&row.jetpt<s.ph&&std::abs(row.eta)>=s.el&&std::abs(row.eta)<s.eh) {
      weight*=s.w;++matches;
    }
    if(matches>1)throw std::runtime_error("Overlapping SF-table rows");
    if(matches)++reweighted;
    if(!std::isfinite(weight))throw std::runtime_error("Non-finite replay weight");
    counts.Fill(row.ptz,tag,truth,weight);
    const std::map<std::string,double> values={{"m0",row.m0},{"m2",row.m2},{"mn",row.mn},{"mu",row.mu},
      {"db",row.db},{"inverseResidual",row.inverseResidual},{"muEF",row.muEF},{"recoGen",row.recoGen},
      {"m0m0",row.m0*row.m0},{"m0mn",row.m0*row.mn},{"m0mu",row.m0*row.mu},
      {"mnmn",row.mn*row.mn},{"mnmu",row.mn*row.mu},{"mumu",row.mu*row.mu}};
    for(const auto &v:values)if(std::isfinite(v.second))profiles[v.first]->Fill(row.ptz,tag,truth,v.second,weight);
  }
  tree->ResetBranchAddresses();
  // Carry independent generator normalization counters; response components remain mergeable.
  out.mkdir("configInfo")->cd();
  for(const char *name:{"SkimCounter","GeneratorCounters"})if(auto *h=in->Get((std::string("configInfo/")+name).c_str()))h->Write(name);
  out.Write();out.Close();
  std::cout<<"Replayed "<<tree->GetEntries()<<" accepted Legacy jets; "<<reweighted<<" MC rows received explicit SF-table weights.\n";
}
