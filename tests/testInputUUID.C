#undef NDEBUG
#include "../ZJetInputCounters.h"
#include <TMemFile.h>
#include <TSystem.h>
#include <cassert>
#include <functional>

namespace {
// Deliberately reproduce the skim-production collision, without private inputs.
struct UUIDFile : TFile {
  UUIDFile(const char *path):TFile(path,"RECREATE") {
    fUUID=TUUID("11111111-2222-1333-8444-555555555555");
  }
};
void input(const std::string &path,std::vector<ULong64_t> ids,bool missing=false) {
  UUIDFile f(path.c_str());
  TTree t("Events","");UInt_t run=1,lumi=2;ULong64_t event=0;float payload=42;
  t.Branch("run",&run);t.Branch("luminosityBlock",&lumi);
  if(!missing)t.Branch("event",&event);
  t.Branch("payload",&payload);
  for(auto id:ids){event=id;t.Fill();}t.Write();
  f.mkdir("configInfo")->cd();TH1D h("SkimCounter","",1,0,1);
  h.GetXaxis()->SetBinLabel(1,"Skim: All events");h.SetBinContent(1,100);h.Write();
}
void mustFail(const std::function<void()> &f) {
  bool failed=false;try{f();}catch(const std::runtime_error &e){failed=true;std::cout<<"Expected: "<<e.what()<<"\n";}
  assert(failed);
}
}
void testInputUUID() {
  const std::string dir="tmp/uuid_test_"+std::to_string(gSystem->GetPid());
  assert(gSystem->mkdir(dir.c_str(),true)==0);
  input(dir+"/a.root",{1,2,3});input(dir+"/b.root",{4,5});
  input(dir+"/c.root",{6});input(dir+"/overlap.root",{3,7});
  input(dir+"/missing.root",{8},true);input(dir+"/internal.root",{9,9});
  TMemFile out("uuid_output.root","RECREATE");ZJetInputCounters counts(&out);
  TFile a((dir+"/a.root").c_str()),b((dir+"/b.root").c_str()),c((dir+"/c.root").c_str());
  assert(std::string(a.GetUUID().AsString())==b.GetUUID().AsString());
  // Audit must not alter the event already loaded by the caller.
  auto *tree=static_cast<TTree*>(b.Get("Events"));ULong64_t active=0;
  tree->SetBranchAddress("event",&active);tree->GetEntry(1);assert(active==5);
  counts.add(&a,false);counts.add(&b,false);assert(active==5);
  tree->GetEntry(0);assert(active==4);tree->ResetBranchAddresses();
  counts.add(&c,false);
  assert(counts.skim->GetBinContent(1)==300&&counts.skim->GetBinContent(2)==3);
  assert(counts.runs->GetBinContent(5)==0&&counts.uuidChecks->GetBinContent(1)==3);
  mustFail([&]{counts.add(&a,false);});
  for(const char *name:{"overlap","missing","internal"}) {
    TFile f((dir+"/"+name+".root").c_str());mustFail([&]{counts.add(&f,false);});
  }
  assert(counts.skim->GetBinContent(2)==3);
  // Plain copies must still fail even with different paths.
  assert(gSystem->CopyFile((dir+"/a.root").c_str(),(dir+"/copy.root").c_str())==0);
  TFile copy((dir+"/copy.root").c_str());mustFail([&]{counts.add(&copy,false);});
  TH1D sum(*counts.skim);sum.SetDirectory(nullptr);sum.Add(counts.skim);
  assert(sum.GetBinContent(1)==600&&sum.GetBinContent(2)==6);
  std::cout<<"UUID tests passed: disjoint collisions, overlap/copies, missing IDs, active-tree preservation, additive counters\n";
}
