#ifndef ZJET_VALIDATE_TAGGING_CONTROLS_C
#define ZJET_VALIDATE_TAGGING_CONTROLS_C
#include <TFile.h>
#include <TH1.h>
#include <TKey.h>
#include <TObjString.h>
#include <TTree.h>
#include <stdexcept>
#include <iostream>
#include <memory>
void validateTaggingControls(const char *path,bool mc=true,int expectedFiles=-1) {
  std::unique_ptr<TFile> f(TFile::Open(path));
  if(!f||f->IsZombie())throw std::runtime_error("Cannot open tagging-control output");
  auto *n=dynamic_cast<TH1*>(f->Get("configInfo/SkimCounter"));
  auto *r=dynamic_cast<TH1*>(f->Get("configInfo/GeneratorCounters"));
  if(!n||!r||n->GetBinContent(2)<=0||n->GetBinContent(3)!=0||r->GetBinContent(5)!=0||
      (expectedFiles>=0&&n->GetBinContent(2)!=expectedFiles)||
      (mc&&(r->GetBinContent(1)<=0||r->GetBinContent(4)!=0)))
    throw std::runtime_error("Missing, duplicated or incomplete pre-skim input counters");
  for(const char *co:{"legacy","new","legacy_genweight","new_genweight"}) {
    const std::string base=std::string("TaggingControls/")+co+"/hybrid_pnet030/";
    for(const char *name:{"counts","h_mass","h_met_perp","h_second_hf","p_fnu_vs_muef","calibration_tag5"})
      if(!f->Get((base+name).c_str()))throw std::runtime_error("Missing TT/HF control");
  }
  auto *definition=dynamic_cast<TObjString*>(f->Get("TaggingControls/definition"));
  auto *mode=dynamic_cast<TObjString*>(f->Get("zjet_analysis_mode"));
  if(mode&&mode->GetString()=="legacy") {
    auto *tree=dynamic_cast<TTree*>(f->Get("LegacyFlavor/events"));
    auto *counts=dynamic_cast<TH1*>(f->Get("TaggingControls/legacy/bT_cvlT_pnet045/counts"));
    if(!tree||!counts||tree->GetEntries()!=counts->GetEntries())
      throw std::runtime_error("Legacy replay tree and selected-probe count disagree");
    if(tree->GetEntries()>0&&(!tree->GetBranch("b")||!tree->GetBranch("genWeight")))
      throw std::runtime_error("Incomplete Legacy replay schema");
    auto *newCounts=dynamic_cast<TH1*>(f->Get("TaggingControls/new/hybrid_pnet030/counts"));
    if(!newCounts||newCounts->GetEntries()!=0)throw std::runtime_error("Legacy-only run filled new-method controls");
  }
  if(!definition)throw std::runtime_error("Missing tagging definition");
  if(definition->GetString().Contains("btag_sha256")) {
    for(const char *b:{"M","T"})for(const char *c:{"M","T"})for(const char *q:{"030","045"})
      if(!f->Get((std::string("TaggingControls/new/b")+b+"_c"+c+"_pnet"+q+"/counts").c_str()))
        throw std::runtime_error("Missing official WP policy");
  }
  std::cout<<"Validated TT/HF controls and "<<n->GetBinContent(2)<<" input-file counters: "<<path<<"\n";
}
#endif
