#include <TFile.h>
#include <TKey.h>
#include <TClass.h>
#include <TTree.h>
#include <TH1.h>
#include <TObjString.h>
#include <TParameter.h>
#include <TSystem.h>
#include <functional>
#include <memory>
#include <set>
#include <stdexcept>
#include <iostream>

// Derive a histogram-only view without reading/copying replay event baskets.
// The source (possibly on EOS) is always READ-only. The caller publishes the
// newly-created output only after running the physics validators.
void writeCompactOutput(const char *source,const char *destination) {
  TDirectory::TContext context;
  if(!gSystem->AccessPathName(destination))
    throw std::runtime_error("Compact output already exists; choose a new temporary path");
  std::unique_ptr<TFile> input(TFile::Open(source,"READ"));
  if(!input||input->IsZombie()||input->TestBit(TFile::kRecovered)||input->GetEND()>input->GetSize())
    throw std::runtime_error("Cannot compact an unreadable, truncated or recovered file");
  if(input->GetKey("zjet_compact_definition"))
    throw std::runtime_error("Input is already compact; use the original merged source");
  TFile output(destination,"CREATE");
  if(output.IsZombie())throw std::runtime_error("Cannot create compact output");
  Long64_t replayEntries=-1,replayBytes=0,copied=0;
  std::function<void(TDirectory*,TDirectory*,std::string)> copy;
  copy=[&](TDirectory *from,TDirectory *to,std::string prefix) {
    std::set<std::string> names;
    TIter next(from->GetListOfKeys());
    while(auto *key=dynamic_cast<TKey*>(next())) {
      const std::string name=key->GetName(),path=prefix+name;
      if(!names.insert(name).second)
        throw std::runtime_error("Duplicate key cycles in source: "+path+"; finalize the merge first");
      auto *type=TClass::GetClass(key->GetClassName());
      if(!type)throw std::runtime_error("Unknown stored class at "+path);
      if(type->InheritsFrom(TDirectory::Class())) {
        auto *child=from->GetDirectory(name.c_str());
        auto *target=to->mkdir(name.c_str(),key->GetTitle());
        if(!child||!target)throw std::runtime_error("Cannot copy directory "+path);
        copy(child,target,path+"/");
      } else if(type->InheritsFrom(TTree::Class())) {
        if(path!="LegacyFlavor/events")
          throw std::runtime_error("Unexpected TTree "+path+"; refusing to silently drop or copy it");
        // Read only the tree header for the audit summary, never GetEntry().
        std::unique_ptr<TTree> tree(dynamic_cast<TTree*>(key->ReadObj()));
        if(!tree)throw std::runtime_error("Cannot read replay header");
        tree->SetDirectory(nullptr);
        replayEntries=tree->GetEntries();replayBytes=tree->GetZipBytes();
      } else {
        std::unique_ptr<TObject> object(key->ReadObj());
        if(!object)throw std::runtime_error("Cannot read "+path);
        if(auto *h=dynamic_cast<TH1*>(object.get()))h->SetDirectory(nullptr);
        to->cd();
        if(object->Write(name.c_str(),TObject::kOverwrite)<=0)
          throw std::runtime_error("Cannot write "+path);
        ++copied;
      }
    }
  };
  copy(input.get(),&output,"");
  if(replayEntries<0)throw std::runtime_error("Full source lacks LegacyFlavor/events");
  output.cd();
  TObjString("v1: histogram-only derivative; all non-tree objects preserved; LegacyFlavor/events omitted. "
             "Arbitrary offline retagging requires the original full EOS merge, not this file.")
    .Write("zjet_compact_definition",TObject::kOverwrite);
  TParameter<Long64_t>("zjet_replay_entries",replayEntries).Write();
  TParameter<Long64_t>("zjet_replay_compressed_bytes",replayBytes).Write();
  TObjString(input->GetUUID().AsString()).Write("zjet_compact_source_uuid",TObject::kOverwrite);
  output.Close();
  if(output.TestBit(TFile::kWriteError))throw std::runtime_error("Compact ROOT output write failed");
  TFile check(destination,"READ");
  if(check.IsZombie()||check.TestBit(TFile::kRecovered)||check.GetEND()>check.GetSize()||
      check.Get("LegacyFlavor/events")||!check.Get("zjet_compact_definition"))
    throw std::runtime_error("Compact output integrity check failed");
  std::cout<<"Compact output: "<<check.GetSize()<<" bytes; copied "<<copied
           <<" non-tree objects; replay retained only in source: "<<replayEntries
           <<" entries, "<<replayBytes<<" compressed bytes"<<std::endl;
}
