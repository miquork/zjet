#ifndef ZJET_TAGGING_CONTROLS_H
#define ZJET_TAGGING_CONTROLS_H
#include "ZJetResponseAudit.h"
#include <TH2D.h>
#include <TProfile2D.h>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <vector>

namespace ZJetTaggingControls {
struct Config {
  bool enabled=false;
  std::map<std::string,double> values;
  std::string source;
  explicit Config(const char *path="data/Tagging/official_2024.txt") {
    std::ifstream in(path); if(!in) return;
    std::string line;
    while(std::getline(in,line)) {
      source+=line+"\n";
      if(line.empty()||line[0]=='#')continue;
      std::istringstream s(line);std::string key;double v;
      if(!(s>>key>>v)||!ZJetResponseAudit::valid(v)||values.count(key))
        throw std::runtime_error("Invalid official tagging WP configuration");
      values[key]=v;
    }
    for(const auto *key:{"bM","bT","cMCvB","cMCvL","cTCvB","cTCvL"})
      if(!values.count(key))throw std::runtime_error("Missing official WP key");
    enabled=true;
  }
  std::map<std::string,int> tags(double b,double cvb,double cvl,double q) const {
    std::map<std::string,int> out;
    if(!enabled)return out;
    // Synchronized reference: CvL-only charm after the official B veto.
    // This is not claimed to be the two-dimensional official C working point.
    for(const std::string bw:{"M","T"}) {
      int tag=0;
      if(ZJetResponseAudit::valid(b)) {
        if(b>=values.at("b"+bw))tag=5;
        else if(ZJetResponseAudit::valid(cvl)) {
          if(cvl>=values.at("cTCvL"))tag=4;
          else if(ZJetResponseAudit::valid(q))tag=q>=.45?1:6;
        }
      }
      out["b"+bw+"_cvlT_pnet045"]=tag;
    }
    for(const std::string bw:{"M","T"})for(const std::string cw:{"M","T"})
      for(const double cut:{.30,.45}) {
        int tag=0;
        // Do not discard a valid B decision because an unrelated score is invalid.
        if(ZJetResponseAudit::valid(b)) {
          if(b>=values.at("b"+bw))tag=5;
          else if(ZJetResponseAudit::valid(cvb)&&ZJetResponseAudit::valid(cvl)) {
            if(cvb>=values.at("c"+cw+"CvB")&&cvl>=values.at("c"+cw+"CvL"))tag=4;
            else if(ZJetResponseAudit::valid(q))tag=q>=cut?1:6;
          }
        }
        out["b"+bw+"_c"+cw+(cut<.4?"_pnet030":"_pnet045")]=tag;
      }
    return out;
  }
};
struct Histograms {
  Config config;
  std::map<std::string,TH3D*> counts, distributions;
  std::map<std::string,TProfile3D*> response;
  std::map<std::string,TH3D*> calibration;
};
inline const std::vector<std::string>& cohorts() {
  static const std::vector<std::string> v={"legacy","legacy_signed","legacy_tightphi",
    "legacy_ptratio","legacy_jetid","legacy_pairveto","new","common_unit","common_event","common_pair",
    "legacy_genweight","new_genweight"};
  return v;
}
inline Histograms book(TDirectory *parent,int nx,const double *xbins,Config config=Config()) {
  TDirectory::TContext context;Histograms h;h.config=config;
  auto *top=parent->mkdir("TaggingControls");top->cd();
  TObjString((std::string("v1; additive controls; no SF applied. Baseline hybrid plus official BvsAll/CvB/CvL WPs. ")+
    "X=Z pT except calibration X=corrected jet pT, Y=eta, Z=hadron flavor. "+
    "Mass window inherited: 70<mumu<110 GeV. Components are signed projections. "
    "legacy_genweight=full genWeight without PU; new_genweight=new signed pair weight times abs(genWeight). "
    "Use Runs genEventSumw only for the full-generator-weight cohorts, not nominal unit/sign cohorts. "+h.config.source).c_str()).Write("definition");
  std::map<std::string,int> policies=h.config.tags(0,0,0,0);policies["hybrid_pnet030"]=0;
  double ids[8]={-.5,.5,1.5,2.5,3.5,4.5,5.5,6.5};
  for(const auto &co:cohorts()) {
    auto *d=top->mkdir(co.c_str());
    for(const auto &pol:policies) {
      d->mkdir(pol.first.c_str())->cd();const auto key=co+"/"+pol.first;
      h.counts[key]=new TH3D("counts",";Z pT;reco tag;parton flavor",nx,xbins,7,ids,7,ids);h.counts[key]->Sumw2();
      for(const std::string obs:{"m0","m2","mn","mu","fnu","muEF","reco_gen","genmn","genmu"})
        h.response[key+"/"+obs]=new TProfile3D(obs.c_str(),";Z pT;reco tag;parton flavor",nx,xbins,7,ids,7,ids);
      if(co!="new"&&co!="legacy"&&co!="new_genweight"&&co!="legacy_genweight")continue;
      // Control axes: X=Z pT, Y=observable, Z=reco category.
      for(const auto &spec:std::map<std::string,std::array<double,3>>{
        {"mass",{80,70,110}},{"met_perp",{60,0,180}},{"met_parallel",{100,-250,250}},
        {"second_hf",{2,-.5,1.5}},{"n_hf",{6,-.5,5.5}},
        {"muef",{50,0,.5}},{"fnu",{100,-1,1}},{"mn",{100,-1,1}},{"mu",{100,-1,1}}}) {
        std::vector<double> y(int(spec.second[0])+1);
        for(size_t i=0;i<y.size();++i)y[i]=spec.second[1]+(spec.second[2]-spec.second[1])*i/(y.size()-1);
        h.distributions[key+"/"+spec.first]=new TH3D(("h_"+spec.first).c_str(),";Z pT;observable;reco tag",
          nx,xbins,int(spec.second[0]),y.data(),7,ids);
        h.distributions[key+"/"+spec.first]->Sumw2();
      }
      // Joint controls allow conditional means, not merely separate marginals.
      for(const std::string control:{"mass","met_perp","second_hf","muef"})
        for(const std::string obs:{"mn","mu","fnu"}) {
          auto *axes=h.distributions.at(key+"/"+control);
          h.response[key+"/"+obs+"_vs_"+control]=new TProfile3D(
            ("p_"+obs+"_vs_"+control).c_str(),";Z pT;control;reco tag",nx,xbins,
            axes->GetNbinsY(),axes->GetYaxis()->GetXbins()->GetArray(),7,ids);
        }
      for(int r:{0,1,4,5,6}) {
        const auto name="calibration_tag"+std::to_string(r);
        double eta[14];for(int i=0;i<14;++i)eta[i]=-1.3+.2*i;
        h.calibration[key+"/"+name]=new TH3D(name.c_str(),";corrected jet pT;eta;hadron flavor",nx,xbins,13,eta,7,ids);
        h.calibration[key+"/"+name]->Sumw2();
      }
    }
  }
  return h;
}
inline void fill(Histograms &h,const std::string &co,double pt,double jetpt,double eta,int truth,int hadron,
    const ZJetResponseAudit::Scores &s,double b,const ZJetResponseAudit::Response &v,double muEF,
    double mass,double metx,double mety,const std::map<std::string,int> &otherHF,double weight,
    double recoGen=std::nan(""),double genmn=std::nan(""),double genmu=std::nan("")) {
  auto tags=h.config.tags(b,s.ucvb,s.ucvl,s.qpnet);tags["hybrid_pnet030"]=ZJetResponseAudit::hybrid(s,s.qpnet,.3);
  // Derive HDM only from merged component means, never average event ratios.
  const double fnu=v.mn+v.mu/.92;
  const std::map<std::string,double> values={{"m0",v.m0},{"m2",v.m2},{"mn",v.mn},{"mu",v.mu},
    {"fnu",fnu},
    {"muEF",muEF},{"reco_gen",recoGen},{"genmn",genmn},{"genmu",genmu}};
  for(const auto &tag:tags) {
    const auto key=co+"/"+tag.first;const int r=tag.second;
    h.counts.at(key)->Fill(pt,r,truth,weight);
    for(const auto &v:values)if(std::isfinite(v.second))h.response.at(key+"/"+v.first)->Fill(pt,r,truth,v.second,weight);
    if(co!="new"&&co!="legacy"&&co!="new_genweight"&&co!="legacy_genweight")continue;
    const int n=otherHF.at(tag.first);
    const std::map<std::string,double> controls={{"mass",mass},{"met_perp",std::abs(mety)},
      {"met_parallel",metx},{"second_hf",n>0?1.:0.},{"n_hf",double(n)},
      {"muef",muEF},{"mn",v.mn},{"mu",v.mu},{"fnu",fnu}};
    for(const auto &v:controls)if(std::isfinite(v.second))h.distributions.at(key+"/"+v.first)->Fill(pt,v.second,r,weight);
    for(const std::string ctl:{"mass","met_perp","second_hf","muef"})for(const std::string obs:{"mn","mu","fnu"})
      h.response.at(key+"/"+obs+"_vs_"+ctl)->Fill(pt,controls.at(ctl),r,values.at(obs),weight);
    h.calibration.at(key+"/calibration_tag"+std::to_string(r))->Fill(jetpt,eta,hadron,weight);
  }
}
}
#endif
