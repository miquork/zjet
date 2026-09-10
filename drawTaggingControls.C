// Standalone TT/HF controls. No absolute normalization or DY+TT mixing is implicit.
#include "tdrstyle_mod22.C"
#include <TFile.h>
#include <TH3D.h>
#include <TProfile3D.h>
#include <TProfile.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TSystem.h>
#include <memory>
#include <stdexcept>
#include <vector>
#include <map>
#include <iostream>

namespace ZJetControlPlots {
// Merge stored profile moments exactly across pT bins, including sum(w^2).
TProfile *project(const TProfile3D *h,const char *name,int first,int last,int tag) {
  auto *axis=h->GetYaxis();
  auto *p=new TProfile(name,"",axis->GetNbins(),axis->GetXbins()->GetArray());p->SetDirectory(nullptr);p->Sumw2();
  const auto *ww=h->GetBinSumw2();
  for(int y=0;y<=axis->GetNbins()+1;++y) {
    double w=0,wy=0,wyy=0,w2=0;
    for(int x=first;x<=last;++x) {
      const int bin=h->GetBin(x,y,tag);const double v=h->GetBinEntries(bin);
      w+=v;wy+=h->GetArray()[bin];wyy+=h->GetSumw2()->At(bin);
      w2+=ww->GetSize()?ww->At(bin):v;
    }
    p->SetBinEntries(y,w);p->SetBinContent(y,wy);p->GetSumw2()->SetAt(wyy,y);p->GetBinSumw2()->SetAt(w2,y);
  }
  return p;
}
}
void drawTaggingControls(const char *data="rootfiles/zjet_DATA.root",
    const char *dy="rootfiles/zjet_MC.root",const char *tt="",
    const char *policy="bM_cM_pnet045",const char *cohort="new_genweight",
    int tag=5,double ptmin=30,double ptmax=100,const char *directory="output/taggingControls") {
  if(tag!=1&&tag!=4&&tag!=5&&tag!=6)throw std::runtime_error("Use tag 1,4,5,6");
  std::vector<std::unique_ptr<TFile>> files;
  for(const auto *path:{data,dy,tt})if(path&&*path) {
    files.emplace_back(TFile::Open(path));if(!files.back()||files.back()->IsZombie())throw std::runtime_error("Cannot open TT control input");
  }
  const std::string base=std::string("TaggingControls/")+cohort+"/"+policy+"/";
  gSystem->mkdir(directory,true);setTDRStyle();extraText="Work in progress";
  lumi_136TeV="2024I + Summer24";
  const char *labels[]={"Data","DY","t#bar{t} (separate)"};int colors[]={1,601,633};
  const std::map<std::string,std::string> axes={{"mass","m_{#mu#mu} (GeV)"},
    {"met_perp","|MET perpendicular to Z| (GeV)"},{"second_hf","Additional HF tag (0 / 1)"},{"muef","Jet PF muon fraction"}};
  for(const auto &axis:axes)for(const std::string moment:{"","mn","mu","fnu"}) {
    const std::string object=moment.empty()?"h_"+axis.first:"p_"+moment+"_vs_"+axis.first;
    std::vector<std::unique_ptr<TH1>> histograms;
    double maximum=0;
    for(size_t s=0;s<files.size();++s) {
      auto *h=dynamic_cast<TH3*>(files[s]->Get((base+object).c_str()));
      if(!h)throw std::runtime_error("Missing new production control: "+base+object);
      int first=h->GetXaxis()->FindFixBin(ptmin+1e-6),last=h->GetXaxis()->FindFixBin(ptmax-1e-6);
      if(std::abs(h->GetXaxis()->GetBinLowEdge(first)-ptmin)>1e-5||
         std::abs(h->GetXaxis()->GetBinUpEdge(last)-ptmax)>1e-5)
        throw std::runtime_error("Choose native pT boundaries; partial bins are not split");
      const auto name=object+std::to_string(s);TH1 *p=nullptr;
      if(moment.empty()) {
        p=h->ProjectionY(name.c_str(),first,last,tag+1,tag+1,"e");p->SetDirectory(nullptr);
        const double n=p->Integral(0,p->GetNbinsX()+1);
        if(n<=0)throw std::runtime_error("Nonpositive shape normalization");
        p->Scale(1/n);
        std::cout<<object<<" "<<labels[s]<<" fraction outside visible axis="
                 <<p->GetBinContent(0)+p->GetBinContent(p->GetNbinsX()+1)<<"\n";
      } else {
        auto *prof=dynamic_cast<TProfile3D*>(h);if(!prof)throw std::runtime_error("Expected profile");
        p=ZJetControlPlots::project(prof,name.c_str(),first,last,tag+1);
      }
      maximum=std::max(maximum,p->GetMaximum());histograms.emplace_back(p);
    }
    const std::string name=std::string(cohort)+"_"+policy+"_tag"+std::to_string(tag)+"_"+object+
        "_pt"+std::to_string(int(ptmin))+"to"+std::to_string(int(ptmax));
    auto *x=histograms.front()->GetXaxis();TH1D frame((name+"frame").c_str(),"",1,x->GetXmin(),x->GetXmax());frame.SetDirectory(nullptr);
    frame.GetXaxis()->SetTitle(axis.second.c_str());frame.GetYaxis()->SetTitle(moment.empty()?"Fraction / bin (each sample unit area)":
       moment=="fnu"?"<m_{n}+m_{u}/0.92>":("<m_{"+moment.substr(1)+"}>").c_str());
    frame.SetMinimum(moment.empty()?0:-.2);frame.SetMaximum(moment.empty()?1.55*maximum:std::max(.4,1.25*maximum));
    auto *c=tdrCanvas(name.c_str(),&frame,8,11,kSquare);c->SetLeftMargin(.18);
    frame.GetYaxis()->SetTitleOffset(1.9);frame.GetYaxis()->SetTitleSize(.04);
    TLegend leg(.62,.73,.93,.88);leg.SetFillStyle(0);leg.SetBorderSize(0);leg.SetTextSize(.032);
    for(size_t s=0;s<histograms.size();++s) {
      auto *h=histograms[s].get();h->SetLineColor(colors[s]);h->SetMarkerColor(colors[s]);h->SetMarkerStyle(s?24:20);
      h->Draw("E1 SAME");leg.AddEntry(h,labels[s],"lp");
    }
    leg.Draw();TLatex text;text.SetNDC();text.SetTextSize(.025);
    text.DrawLatex(.22,.65,Form("%g < p_{T,Z} < %g GeV, |#eta^{jet}| < 1.3",ptmin,ptmax));
    text.DrawLatex(.22,.61,(std::string(cohort)+", "+policy+", reco tag "+std::to_string(tag)).c_str());
    c->SaveAs((std::string(directory)+"/"+name+".pdf").c_str());delete c;
  }
}
