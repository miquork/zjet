#ifndef ZJET_INPUT_COUNTERS_H
#define ZJET_INPUT_COUNTERS_H
#include <TDirectory.h>
#include <TFile.h>
#include <TH1D.h>
#include <TTree.h>
#include <TObjString.h>
#include <set>
#include <stdexcept>
// Numeric bins with fixed labels: additive through hadd, no duplicate strings.
struct ZJetInputCounters {
  TH1D *skim, *runs;
  std::set<std::string> files;
  explicit ZJetInputCounters(TDirectory *out) {
    TDirectory::TContext context;out->mkdir("configInfo")->cd();
    skim=new TH1D("SkimCounter","Input counters before any analysis cut",3,.5,3.5);
    skim->GetXaxis()->SetBinLabel(1,"Skim: All events");
    skim->GetXaxis()->SetBinLabel(2,"Input files");
    skim->GetXaxis()->SetBinLabel(3,"Missing SkimCounter");
    runs=new TH1D("GeneratorCounters","Do not interchange counts and weight sums",5,.5,5.5);
    const char *labels[]={"genEventCount","genEventSumw","genEventSumw2","Missing Runs","Duplicate file UUID"};
    for(int i=0;i<5;++i)runs->GetXaxis()->SetBinLabel(i+1,labels[i]);
    TObjString("SkimCounter bin 1 sums each input's configInfo/SkimCounter bin labelled Skim: All events once. "
      "Runs totals kept separately. Neither is a pre-skim sign-weight sum. Missing counters explicitly flagged; "
      "no fallback to selected Events entries. Check dataset disjointness across jobs.").Write("counter_definition");
  }
  void add(TFile *file,bool mc) {
    if(!file)throw std::runtime_error("No current input file for counters");
    const std::string id=file->GetUUID().AsString();
    if(!files.insert(id).second) {runs->AddBinContent(5,1);throw std::runtime_error("Duplicate input file UUID");}
    skim->AddBinContent(2,1);
    auto *h=dynamic_cast<TH1*>(file->Get("configInfo/SkimCounter"));int bin=0;
    if(h)for(int i=1;i<=h->GetNbinsX();++i) {
      std::string label=h->GetXaxis()->GetBinLabel(i);
      const auto first=label.find_first_not_of(" \t\r\n"),last=label.find_last_not_of(" \t\r\n");
      if(first!=std::string::npos)label=label.substr(first,last-first+1);
      if(label=="Skim: All events")bin=i;
    }
    if(bin)skim->AddBinContent(1,h->GetBinContent(bin));else skim->AddBinContent(3,1);
    if(!mc)return;
    auto *tree=dynamic_cast<TTree*>(file->Get("Runs"));
    if(!tree||!tree->GetBranch("genEventCount")||!tree->GetBranch("genEventSumw")||!tree->GetBranch("genEventSumw2")) {runs->AddBinContent(4,1);return;}
    Long64_t count=0;double sumw=0,sumw2=0;
    tree->SetBranchAddress("genEventCount",&count);tree->SetBranchAddress("genEventSumw",&sumw);
    tree->SetBranchAddress("genEventSumw2",&sumw2);
    for(Long64_t i=0;i<tree->GetEntries();++i) {
      if(tree->GetEntry(i)<=0)throw std::runtime_error("Cannot read Runs input counters");
      runs->AddBinContent(1,count);runs->AddBinContent(2,sumw);runs->AddBinContent(3,sumw2);
    }
    tree->ResetBranchAddresses();
  }
};
#endif
