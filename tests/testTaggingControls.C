#include "../ZJetTaggingControls.h"
#include "../ZJetInputCounters.h"
#include "../drawTaggingControls.C"
#include <TMemFile.h>
#include <cassert>
#include <iostream>
void testTaggingControls() {
  using namespace ZJetTaggingControls;
  Config c("/nonexistent/zjet-synthetic-wp-test");c.enabled=true;
  c.values={{"bM",.1},{"bT",.6},{"cMCvB",.4},{"cMCvL",.3},{"cTCvB",.5},{"cTCvL",.7}};
  assert(c.tags(1,-1,-1,-1).at("bM_cM_pnet030")==5);
  assert(c.tags(-1,1,1,1).at("bM_cM_pnet030")==0);
  assert(c.tags(0,1,1,1).at("bM_cM_pnet030")==4);
  assert(c.tags(0,0,0,.3).at("bM_cM_pnet030")==1);
  assert(c.tags(0,0,0,.3).at("bM_cM_pnet045")==6);
  const double score=.5*(c.values.at("bM")+c.values.at("bT"));
  assert(c.tags(score,1,1,1).at("bM_cM_pnet030")==5);
  assert(c.tags(score,1,1,1).at("bT_cM_pnet030")==4);
  TMemFile out("tagtest.root","RECREATE");double x[]={30,600};auto h=book(&out,1,x,c);
  std::map<std::string,int> other=c.tags(0,0,0,0);other["hybrid_pnet030"]=0;for(auto &v:other)v.second=0;
  ZJetResponseAudit::Scores s{0,0,0,0,1,1,1,1};
  fill(h,"new",50,55,0,5,5,s,score,{1,1,0,0,1,1},.1,91,0,10,other,.5);
  assert(h.counts.at("new/bM_cM_pnet030")->Integral()==.5);
  assert(h.counts.at("new/bT_cM_pnet030")->GetBinContent(1,5,6)==.5);
  // Profile projection must preserve weighted first/second moments and effective entries.
  double px[]={30,50,100},py[]={0,1,2},pz[]={-.5,.5};
  TProfile3D joint("projection_test","",2,px,2,py,1,pz);
  joint.Fill(40,.5,0,2,.5);joint.Fill(60,.5,0,4,1.5);
  std::unique_ptr<TProfile> projection(ZJetControlPlots::project(&joint,"projected",1,2,1));
  TProfile direct("direct","",2,py);direct.Fill(.5,2,.5);direct.Fill(.5,4,1.5);
  assert(std::abs(projection->GetBinContent(1)-direct.GetBinContent(1))<1e-12);
  assert(std::abs(projection->GetBinError(1)-direct.GetBinError(1))<1e-12);
  ZJetInputCounters counters(&out);
  TMemFile input("inputtest.root","RECREATE");input.mkdir("configInfo")->cd();
  TH1D skim("SkimCounter","",1,0,1);skim.GetXaxis()->SetBinLabel(1,"Skim: All events");skim.SetBinContent(1,1234);skim.Write();
  counters.add(&input,false);assert(counters.skim->GetBinContent(1)==1234);
  TH1D merged(*counters.skim);merged.SetDirectory(nullptr);
  assert(merged.Add(counters.skim));assert(merged.GetBinContent(1)==2468);
  TList list;list.Add(counters.skim);
  merged.Merge(&list);assert(merged.GetBinContent(1)==3702);
  bool caught=false;try {counters.add(&input,false);}catch(const std::runtime_error&){caught=true;}assert(caught);
  std::cout<<"Tagging and input-counter tests passed\n";
}
