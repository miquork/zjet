#ifndef ZJET_VALIDATE_RESPONSE_AUDIT_C
#define ZJET_VALIDATE_RESPONSE_AUDIT_C
#include "ZJetResponseAudit.h"
#include <TFile.h>
#include <TKey.h>
#include <stdexcept>
#include <iostream>

void validateResponseAudit(const char *filename,bool required=true) {
  TFile f(filename);
  if(f.IsZombie()) throw std::runtime_error("Cannot validate ResponseAudit input");
  if(!f.Get("ResponseAudit/definition")) {
    if(required) throw std::runtime_error("Missing ResponseAudit definition");
    return;
  }
  int definitions=0;
  TIter next(f.GetDirectory("ResponseAudit")->GetListOfKeys());
  while(auto *key=dynamic_cast<TKey*>(next()))
    if(std::string(key->GetName())=="definition")++definitions;
  if(definitions!=1)throw std::runtime_error("Duplicate ResponseAudit metadata cycles");
  auto get=[&](const std::string &path) {
    auto *h=dynamic_cast<TH3D*>(f.Get(path.c_str()));
    if(!h)throw std::runtime_error("Missing audit histogram "+path);
    return h;
  };
  auto close=[](double a,double b) {
    return std::isfinite(a)&&std::isfinite(b)&&std::abs(a-b)<1.e-8*std::max(1.,std::max(std::abs(a),std::abs(b)));
  };
  for(const std::string sample:{"new","legacy","common_new","common_legacy"}) {
    auto *reference=get("ResponseAudit/"+sample+"/hybrid_pnet030/counts");
    for(const auto &policy:ZJetResponseAudit::categories({})) {
      auto *h=get("ResponseAudit/"+sample+"/"+policy.first+"/counts");
      for(int x=0;x<=h->GetNbinsX()+1;++x)for(int z=1;z<=7;++z) {
        double a=0.,b=0.;for(int y=1;y<=7;++y) {
          a+=h->GetBinContent(x,y,z);b+=reference->GetBinContent(x,y,z);
        }
        if(!close(a,b))throw std::runtime_error("Tag policy changed the probe population");
      }
    }
  }
  for(const auto &policy:ZJetResponseAudit::categories({})) {
    auto *a=get("ResponseAudit/common_new/"+policy.first+"/counts");
    auto *b=get("ResponseAudit/common_legacy/"+policy.first+"/counts");
    for(int i=0;i<a->GetNcells();++i)
      if(!close(a->GetBinContent(i),b->GetBinContent(i)))
        throw std::runtime_error("Common-probe counts differ");
  }
  for(const std::string mode:{"native","native_complete","native_lowrho",
                            "native_highrho","common","closure",
                            "closure_matched","closure_truth"})
    for(const std::string label:{"truth","reco"}) {
      auto *h=get("ResponseAudit/Ru/"+mode+"_"+label+"_blocks");
      if(h->GetNbinsZ()!=200)throw std::runtime_error("Wrong moment-block schema");
      for(int i=0;i<h->GetNcells();++i)if(!std::isfinite(h->GetBinContent(i)))
        throw std::runtime_error("Non-finite regression moment");
    }
  // Common/native-complete have identical samples and identical N+U.
  auto *a=get("ResponseAudit/Ru/native_complete_truth_blocks");
  auto *b=get("ResponseAudit/Ru/common_truth_blocks");
  for(int x=0;x<=a->GetNbinsX()+1;++x)for(int y=1;y<=7;++y)for(int k=0;k<20;++k) {
    const int base=10*k+1;
    if(!close(a->GetBinContent(x,y,base),b->GetBinContent(x,y,base))||
       !close(a->GetBinContent(x,y,base+1)+a->GetBinContent(x,y,base+6),
              b->GetBinContent(x,y,base+1)+b->GetBinContent(x,y,base+6)))
      throw std::runtime_error("Common-partition population or N+U changed");
  }
  std::cout<<"ResponseAudit validated: "<<filename<<std::endl;
}
#endif
