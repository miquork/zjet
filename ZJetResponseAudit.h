#ifndef ZJET_RESPONSE_AUDIT_H
#define ZJET_RESPONSE_AUDIT_H

// Additive controls only: never change the production selection/calibration.
#include <TDirectory.h>
#include <TH3D.h>
#include <TProfile3D.h>
#include <TObjString.h>
#include <array>
#include <cmath>
#include <cstdint>
#include <map>
#include <string>

namespace ZJetResponseAudit {
constexpr int blocks = 20;
constexpr int nMoments = 10;
inline unsigned block(std::uint64_t event, unsigned run, unsigned lumi) {
  // Same event (including different accepted probes) always in the same block.
  std::uint64_t v = event ^ (std::uint64_t(run)<<32) ^ lumi;
  v += 0x9e3779b97f4a7c15ULL;
  v = (v^(v>>30))*0xbf58476d1ce4e5b9ULL;
  v = (v^(v>>27))*0x94d049bb133111ebULL;
  return unsigned((v^(v>>31))%blocks);
}
inline bool valid(double v) { return std::isfinite(v) && v>=0. && v<=1.; }
struct Scores { double b, cvb, cvl, qdeep, ucvb, ucvl, qpnet, qupart; };
inline int hybrid(const Scores &s, double q, double cut) {
  if (!valid(s.ucvb) || !valid(s.ucvl) || !valid(q)) return 0;
  if (s.ucvb<.5) return 5;
  if (s.ucvl>=.5) return 4;
  return q>=cut ? 1 : 6;
}
inline std::map<std::string,int> categories(const Scores &s) {
  int deep = 0;
  // Match the historical priority and validity behavior exactly.
  if (s.b>.7527) deep=5;
  else if((s.cvb+s.cvl)*.5>.3985) deep=4;
  else if(s.qdeep>=.5) deep=1;
  else if(s.qdeep>=0.) deep=6;
  const int corrected=deep==1 ? 6 : (deep==6 ? 1 : deep);
  const double q=valid(s.qdeep) ? 1.-s.qdeep : -1.;
  return {{"deepjet_asstored",deep},{"deepjet",corrected},
          {"hybrid_pnet020",hybrid(s,s.qpnet,.2)},
          {"hybrid_pnet030",hybrid(s,s.qpnet,.3)},
          {"hybrid_pnet040",hybrid(s,s.qpnet,.4)},
          {"hybrid_deep030",hybrid(s,q,.3)},
          {"hybrid_deep050",hybrid(s,q,.5)},
          {"hybrid_deep070",hybrid(s,q,.7)},
          {"hybrid_upart050",hybrid(s,s.qupart,.5)},
          {"hybrid_upart080",hybrid(s,s.qupart,.8)},
          {"hybrid_upart095",hybrid(s,s.qupart,.95)}};
}
struct Response { double m0=0., m2=0., mn=0., mu=0., db=0., inverseResidual=0.; };
struct Histograms {
  std::map<std::string,TH3D*> counts;
  std::map<std::string,TProfile3D*> profiles;
  std::map<std::string,TH3D*> moments;
  TH3D *truthLabels=nullptr;
  TH3D *selection=nullptr;
};
inline Histograms book(TDirectory *parent, int nx, const double *xbins) {
  TDirectory::TContext context;
  Histograms h;
  auto *top = parent->mkdir("ResponseAudit");
  top->cd();
  TObjString("v1; parallel barrel |eta|<1.3; all policies use identical probes/weights; "
             "common uses the same accepted legacy jet; common_legacy is reweighted "
             "to the new pair weight. Truth is Jet_partonFlavour; gen labels separate. "
             "Ru moments: 1,x,y,x*x,x*y,y*y,n,n*n,x*n,y*n; "
             "20 event-hash blocks; x=genmu/native,common-partition or reco closure; "
             "y=reco mu,n=genmn (closure: reco mn); signed weights retained. "
             "DeepJet QvG=1-btagDeepFlavQG; deepjet_asstored retains historical "
             "inverted labels. rho split at 25 GeV. "
             "No SF or central calibration is modified.").Write("definition");
  double ids[8]={-.5,.5,1.5,2.5,3.5,4.5,5.5,6.5};
  double score[101]; for (int i=0;i<=100;++i) score[i]=i*.01;
  h.truthLabels = new TH3D("truth_labels",";p_{T,Z};Jet parton flavor;GenJet parton flavor",
                          nx,xbins,7,ids,7,ids);
  h.truthLabels->Sumw2();
  h.selection = new TH3D("selection",";p_{T,Z};Jet parton flavor;Legacy overlap",
                         nx,xbins,7,ids,7,ids);
  // Numeric codes, not labels, avoid ROOT label-axis remapping during hadd:
  // 0=no accepted barrel legacy probe, 1=same jet, 2=different jet.
  h.selection->Sumw2();
  const Scores dummy{};
  for (const std::string sample : {"new","legacy","common_new","common_legacy"}) {
    auto *dir=top->mkdir(sample.c_str());
    for (const auto &policy:categories(dummy)) {
      dir->mkdir(policy.first.c_str())->cd();
      const std::string key=sample+"/"+policy.first;
      auto *c=new TH3D("counts",";p_{T,Z};Reco flavor;Jet parton flavor",
                       nx,xbins,7,ids,7,ids); c->Sumw2(); h.counts[key]=c;
      for (const std::string obs : {"m0","m2","mn","mu","db","inverse_residual","reco_gen"})
        h.profiles[key+"/"+obs]=new TProfile3D(obs.c_str(),
          ";p_{T,Z};Reco flavor;Jet parton flavor",nx,xbins,7,ids,7,ids);
    }
    if (sample!="new" && sample!="common_new") continue;
    for (const std::string tag : {"deep","pnet","upart"}) {
      dir->mkdir((tag+"_light_scores").c_str())->cd();
      const std::string key=sample+"/"+tag+"_light_scores";
      auto *c=new TH3D("counts",";p_{T,Z};QvG (common valid support, UParT HF veto);Jet parton flavor",
                       nx,xbins,100,score,7,ids); c->Sumw2(); h.counts[key]=c;
      for (const std::string obs : {"m0","m2","mn","mu","db","inverse_residual","reco_gen"})
        h.profiles[key+"/"+obs]=new TProfile3D(obs.c_str(),
          ";p_{T,Z};QvG;Jet parton flavor",nx,xbins,100,score,7,ids);
    }
  }
  auto *ru=top->mkdir("Ru");
  double slots[blocks*nMoments+1];
  for (int i=0;i<=blocks*nMoments;++i) slots[i]=i-.5;
  double genBins[101]; for(int i=0;i<=100;++i) genBins[i]=-1.+.02*i;
  for(const std::string mode:{"native","native_complete","native_lowrho",
                              "native_highrho","common","closure",
                              "closure_matched","closure_truth"}) {
    for (const std::string label:{"truth","reco"}) {
      const std::string key=mode+"_"+label;
      ru->cd();
      h.moments[key]=new TH3D((key+"_blocks").c_str(),
        ";p_{T,Z};Flavor;10*event block + moment",nx,xbins,7,ids,
        blocks*nMoments,slots);
      h.moments[key]->Sumw2();
      h.profiles[key]=new TProfile3D((key+"_conditional").c_str(),
        ";p_{T,Z};x (under/overflow retained);Flavor",nx,xbins,100,genBins,7,ids);
    }
  }
  return h;
}
inline void fill(Histograms &h,const std::string &sample,double pt,int truth,
                 const Scores &scores,const Response &v,double recoGen,double w) {
  const std::map<std::string,double> obs={{"m0",v.m0},{"m2",v.m2},
    {"mn",v.mn},{"mu",v.mu},{"db",v.db},{"inverse_residual",v.inverseResidual},{"reco_gen",recoGen}};
  auto cell=[&](const std::string &key,double y) {
    h.counts.at(key)->Fill(pt,y,truth,w);
    for(const auto &o:obs) if(std::isfinite(o.second))
      h.profiles.at(key+"/"+o.first)->Fill(pt,y,truth,o.second,w);
  };
  for(const auto &p:categories(scores)) cell(sample+"/"+p.first,p.second);
  if (sample!="new" && sample!="common_new") return;
  // Identical heavy veto and validity support: compare QvG alone. Exact
  // equal-efficiency thresholds can be chosen later on this 0.01 score grid.
  if(!valid(scores.ucvb)||!valid(scores.ucvl)||scores.ucvb<.5||scores.ucvl>=.5||
     !valid(scores.qdeep)||!valid(scores.qpnet)||!valid(scores.qupart)) return;
  for(const auto &s:std::map<std::string,double>{{"deep",1.-scores.qdeep},
      {"pnet",scores.qpnet},{"upart",scores.qupart}})
    cell(sample+"/"+s.first+"_light_scores",std::min(.999999,s.second));
}
inline void fillRu(Histograms &h,const std::string &mode,double pt,int truth,int reco,
                   double x,double y,double n,double w,unsigned eventBlock) {
  if(!std::isfinite(x)||!std::isfinite(y)||!std::isfinite(n)) return;
  const std::array<double,nMoments> m={{1.,x,y,x*x,x*y,y*y,n,n*n,x*n,y*n}};
  for (const auto &label:std::map<std::string,int>{{"truth",truth},{"reco",reco}}) {
    const std::string key=mode+"_"+label.first;
    for(int i=0;i<nMoments;++i)
      h.moments.at(key)->Fill(pt,label.second,eventBlock*nMoments+i,w*m[i]);
    h.profiles.at(key)->Fill(pt,x,label.second,y,w);
  }
}
} // namespace
#endif
