#ifndef ZJET_INPUT_COUNTERS_H
#define ZJET_INPUT_COUNTERS_H
#include <TDirectory.h>
#include <TFile.h>
#include <TH1D.h>
#include <TTree.h>
#include <TObjString.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <map>
#include <memory>
#include <set>
#include <tuple>
#include <vector>
#include <iostream>
#include <stdexcept>
// Numeric bins with fixed labels: additive through hadd, no duplicate strings.
struct ZJetInputCounters {
  TH1D *skim, *runs, *uuidChecks;
  std::set<std::string> paths;
  std::map<std::string,std::vector<std::string>> files;
  using EventID=std::tuple<UInt_t,UInt_t,ULong64_t>;
  // Open a separate handle: never change the active TChain's branch addresses,
  // branch status or current event while auditing a UUID collision.
  static std::set<EventID> eventIDs(const std::string &path) {
    TDirectory::TContext context;
    std::unique_ptr<TFile> input(TFile::Open(path.c_str(),"READ"));
    if(!input||input->IsZombie())throw std::runtime_error("UUID audit cannot open "+path);
    auto *events=dynamic_cast<TTree*>(input->Get("Events"));
    if(!events)throw std::runtime_error("UUID audit missing Events: "+path);
    for(const char *name:{"run","luminosityBlock","event"})
      if(!events->GetBranch(name))throw std::runtime_error("UUID audit missing event-ID branch: "+path);
    TTreeReader reader(events);
    TTreeReaderValue<UInt_t> run(reader,"run"),lumi(reader,"luminosityBlock");
    TTreeReaderValue<ULong64_t> event(reader,"event");
    std::set<EventID> ids;
    Long64_t n=0;
    while(reader.Next()) {
      if(run.GetSetupStatus()<0||lumi.GetSetupStatus()<0||event.GetSetupStatus()<0)
        throw std::runtime_error("UUID audit invalid event-ID types: "+path);
      if(!ids.emplace(*run,*lumi,*event).second)
        throw std::runtime_error("UUID audit repeated event ID within "+path);
      ++n;
    }
    if(n!=events->GetEntries()||reader.GetEntryStatus()!=TTreeReader::kEntryBeyondEnd)
      throw std::runtime_error("UUID audit incomplete event-ID read: "+path);
    return ids;
  }
  explicit ZJetInputCounters(TDirectory *out) {
    TDirectory::TContext context;out->mkdir("configInfo")->cd();
    skim=new TH1D("SkimCounter","Input counters before any analysis cut",3,.5,3.5);
    skim->GetXaxis()->SetBinLabel(1,"Skim: All events");
    skim->GetXaxis()->SetBinLabel(2,"Input files");
    skim->GetXaxis()->SetBinLabel(3,"Missing SkimCounter");
    runs=new TH1D("GeneratorCounters","Do not interchange counts and weight sums",5,.5,5.5);
    const char *labels[]={"genEventCount","genEventSumw","genEventSumw2","Missing Runs","Duplicate file UUID"};
    for(int i=0;i<5;++i)runs->GetXaxis()->SetBinLabel(i+1,labels[i]);
    // Additive diagnostic; absent in older outputs. Keep the existing counter
    // schema/definition unchanged so repaired jobs can merge with those outputs.
    uuidChecks=new TH1D("UUIDCollisionChecks","Full event-ID checks for repeated UUIDs",1,.5,1.5);
    uuidChecks->GetXaxis()->SetBinLabel(1,"Disjoint file pairs verified");
    TObjString("SkimCounter bin 1 sums each input's configInfo/SkimCounter bin labelled Skim: All events once. "
      "Runs totals kept separately. Neither is a pre-skim sign-weight sum. Missing counters explicitly flagged; "
      "no fallback to selected Events entries. Check dataset disjointness across jobs.").Write("counter_definition");
  }
  void add(TFile *file,bool mc) {
    if(!file)throw std::runtime_error("No current input file for counters");
    const std::string id=file->GetUUID().AsString();
    const std::string path=file->GetName();
    if(paths.count(path)) {
      runs->AddBinContent(5,1);
      throw std::runtime_error("Repeated input file path: "+path);
    }
    auto &previous=files[id];
    if(!previous.empty()) {
      std::cout<<"Repeated input UUID "<<id<<"; checking full event IDs for "<<path<<std::endl;
      try {
        const auto currentIDs=eventIDs(path);
        for(const auto &other:previous) {
          const auto priorIDs=eventIDs(other);
          for(const auto &key:currentIDs)if(priorIDs.count(key))
            throw std::runtime_error("Repeated UUID with shared event ID between "+other+" and "+path+
              "; inspect possible overlapping input (MC IDs alone do not prove identical content)");
          uuidChecks->AddBinContent(1,1);
          std::cout<<"UUID collision verified disjoint: "<<other<<" and "<<path
                   <<" ("<<priorIDs.size()<<", "<<currentIDs.size()<<" unique event IDs)"<<std::endl;
        }
      } catch(...) {runs->AddBinContent(5,1);throw;}
    }
    previous.push_back(path);paths.insert(path);
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
