#include <TFile.h>
#include <TTree.h>
#include <TObjString.h>
#include <TH1.h>
#include <fstream>
#include <memory>
#include <stdexcept>
#include <iostream>

// Read-only inspection of a downloaded merge. Write only its provenance JSON
// to a separate scratch file for structural comparison by Python.
void inspectMergedFile(const char *path,const char *metadataOutput,bool mc,
                       int expectedFiles) {
  TFile::SetOpenTimeout(60000);
  std::unique_ptr<TFile> f(TFile::Open(path,"READ"));
  if(!f||f->IsZombie()||f->TestBit(TFile::kRecovered)||f->GetEND()>f->GetSize())
    throw std::runtime_error("Merged file is unreadable, truncated or ROOT-recovered");
  auto *metadata=dynamic_cast<TObjString*>(f->Get("zjet_campaign_metadata"));
  auto *counts=dynamic_cast<TH1*>(f->Get("configInfo/SkimCounter"));
  if(!metadata||!counts||counts->GetBinContent(2)!=expectedFiles)
    throw std::runtime_error("Missing merged provenance or unexpected input-file count");
  auto *mode=dynamic_cast<TObjString*>(f->Get("zjet_analysis_mode"));
  auto *tree=dynamic_cast<TTree*>(f->Get("LegacyFlavor/events"));
  if(mode&&mode->GetString()=="legacy"&&!tree)
    throw std::runtime_error("Missing Legacy replay tree");
  if(tree) {
    // Validate baskets, not just the tree header or the first event. GetEntry
    // reads all branches, including the variable-length offline-retag vectors.
    Bool_t isMC=false;
    if(!tree->GetBranch("isMC")||tree->SetBranchAddress("isMC",&isMC)<0)
      throw std::runtime_error("Replay sample-type branch missing or invalid");
    tree->SetBranchStatus("*",1);
    for(Long64_t i=0;i<tree->GetEntries();++i) {
      if(tree->GetEntry(i)<=0||isMC!=mc)
        throw std::runtime_error("Replay event read failed or sample type disagrees");
      if(i&&i%1000000==0)std::cout<<"Checked "<<i<<" replay entries"<<std::endl;
    }
    tree->ResetBranchAddresses();
  }
  std::ofstream out(metadataOutput);
  out<<metadata->GetString().Data()<<"\n";
  out.close();
  if(!out)throw std::runtime_error("Cannot export merge provenance to scratch");
  std::cout<<"Merged-file integrity check passed: "<<path<<std::endl;
}
