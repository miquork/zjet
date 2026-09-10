#include <TFile.h>
#include <TTree.h>
#include <TH3D.h>
#include <TProfile3D.h>
#include <stdexcept>
#include <iostream>
void testLegacyReplay(const char *source,const char *replay) {
  TFile a(source),b(replay);if(a.IsZombie()||b.IsZombie())throw std::runtime_error("Missing test input");
  auto *x=dynamic_cast<TH3D*>(a.Get("TaggingControls/legacy/bT_cvlT_pnet045/counts"));
  auto *y=dynamic_cast<TH3D*>(b.Get("LegacyFlavor/counts"));
  if(!x||!y||x->GetNcells()!=y->GetNcells())throw std::runtime_error("Replay counts schema mismatch");
  double maxDelta=0;
  for(int i=0;i<x->GetNcells();++i)if(x->GetBinContent(i)!=y->GetBinContent(i))throw std::runtime_error("Replay changed nominal tag counts");
  for(const char *obs:{"m0","m2","mn","mu","muEF"}) {
    auto *p=dynamic_cast<TProfile3D*>(a.Get((std::string("TaggingControls/legacy/bT_cvlT_pnet045/")+obs).c_str()));
    auto *q=dynamic_cast<TProfile3D*>(b.Get((std::string("LegacyFlavor/")+obs).c_str()));
    if(!p||!q)throw std::runtime_error("Missing response profile");
    for(int i=0;i<p->GetNcells();++i) {
      if(p->GetBinEntries(i)!=q->GetBinEntries(i))throw std::runtime_error("Replay changed response entries");
      const double delta=std::abs(p->GetBinContent(i)-q->GetBinContent(i));maxDelta=std::max(delta,maxDelta);
      if(delta>1e-6*std::max(1.,std::abs(p->GetBinContent(i))))throw std::runtime_error("Replay response mismatch");
    }
  }
  std::cout<<"Replay nominal counts exact; maximum component difference "<<maxDelta<<" (float storage)\n";
}
