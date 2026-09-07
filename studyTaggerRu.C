// Read-only physics audit. Never alters either input or the central calibration.
#include "analyzeFlavorMatrix.C"
#include "tdrstyle_mod22.C"
#include <TLegend.h>
#include <TLine.h>
#include <array>

namespace TaggerRuStudy {
using namespace FlavorMatrixAnalysis;
const std::vector<std::pair<std::string,std::vector<int>>> flavors={
  {"uds",{1,3}},{"c",{4}},{"b",{5}},{"g",{6}},{"unknown",{0}}};
const double edges[]={30,50,70,100,150,210,300,400,600};
const int colors[]={kBlue+1,kGreen+2,kRed+1,kMagenta+1,kGray+2};
struct Estimates {
  double slope=NAN,mean=NAN,affine=NAN,intercept=NAN,bias=NAN,mixing=NAN;
};
Estimates estimate(double x,double y,double xx,double xy,double n=0.,
                   double nn=0.,double xn=0.,double yn=0.) {
  Estimates r;
  if(xx>1.e-15) r.slope=xy/xx;
  if(std::abs(x)>1.e-12) r.mean=y/x;
  const double v=xx-x*x;
  if(v>1.e-15) {
    r.affine=(xy-x*y)/v; r.intercept=y-r.affine*x;
  }
  if(std::isfinite(r.slope)&&std::abs(r.slope)>1.e-12) r.bias=y/r.slope-x;
  const double vn=nn-n*n,c=xn-x*n,det=v*vn-c*c;
  if(det>1.e-15) r.mixing=((xy-x*y)*vn-(yn-y*n)*c)/det;
  return r;
}
void point(TGraph &g,double x,double y) {
  if(std::isfinite(y)) g.SetPoint(g.GetN(),x,y);
}
void plot(const std::string &out,const std::string &name,const char *ytitle,
          double low,double high,const std::vector<TGraph*> &graphs,
          const std::vector<std::string> &labels) {
  TH1D frame((name+"_frame").c_str(),"",1,30,600);
  frame.SetDirectory(nullptr); frame.SetMinimum(low); frame.SetMaximum(high);
  frame.GetXaxis()->SetTitle("p_{T,Z} (GeV)");
  frame.GetYaxis()->SetTitle(ytitle);
  auto *c=tdrCanvas(name.c_str(),&frame,8,11,kSquare); c->SetLogx();
  frame.GetXaxis()->SetMoreLogLabels(); frame.GetXaxis()->SetNoExponent();
  TLegend legend(.49,.65,.91,.88); legend.SetBorderSize(0); legend.SetFillStyle(0);
  legend.SetTextSize(.032);
  for(size_t i=0;i<graphs.size();++i) {
    graphs[i]->SetLineColor(colors[i%5]); graphs[i]->SetMarkerColor(colors[i%5]);
    graphs[i]->SetMarkerStyle(20+i%5); graphs[i]->SetLineWidth(2);
    graphs[i]->Draw("LP SAME"); legend.AddEntry(graphs[i],labels[i].c_str(),"lp");
  }
  legend.Draw(); drawCustomLogXLabels(&frame,{30,100,300,600});
  c->SaveAs((out+"/"+name+".pdf").c_str()); delete c;
}
// Aggregate ordinary score histograms; never interpret data bins as truth.
double counts(TH3D *h,double lo,double hi,const std::vector<int> &ids,double cut,bool reverse=false) {
  double s=0.;
  for(int x=1;x<=h->GetNbinsX();++x) {
    const double center=h->GetXaxis()->GetBinCenter(x);
    if(center<lo||center>=hi) continue;
    for(int y=1;y<=h->GetNbinsY();++y) {
      const double score=h->GetYaxis()->GetBinCenter(y);
      if((reverse ? 1.-score : score)>=cut) continue;
      for(int id:ids) s+=h->GetBinContent(x,y,h->GetZaxis()->FindBin(id));
    }
  } return s;
}
double scoreMean(TProfile3D *h,double lo,double hi,const std::vector<int> &ids,double cut,bool reverse=false) {
  double w=0.,s=0.;
  for(int x=1;x<=h->GetNbinsX();++x) {
    const double center=h->GetXaxis()->GetBinCenter(x);
    if(center<lo||center>=hi) continue;
    for(int y=1;y<=h->GetNbinsY();++y) {
      const double score=h->GetYaxis()->GetBinCenter(y);
      if((reverse ? 1.-score : score)>=cut) continue;
      for(int id:ids) {
        const int b=h->GetBin(x,y,h->GetZaxis()->FindBin(id));
        const double v=h->GetBinEntries(b); w+=v; s+=v*h->GetBinContent(b);
      }
    }
  } return std::abs(w)>1.e-12?s/w:NAN;
}
// Block delete-group estimates preserve event correlations across accepted
// probes and across all regression moments; merge by adding the raw TH3D.
void blocks(TFile &file,const std::string &sample,std::ofstream &out) {
  for(const std::string mode:{"native","native_complete","native_lowrho",
                              "native_highrho","common","closure",
                              "closure_matched","closure_truth"}) {
    for(const std::string label:{"truth","reco"}) {
      auto *h=dynamic_cast<TH3D*>(file.Get(("ResponseAudit/Ru/"+mode+"_"+label+"_blocks").c_str()));
      if(!h) continue;
      for(const auto &fl:flavors) for(int p=0;p<8;++p) {
        double sums[20][10]={{0.}};
        for(int ix=1;ix<=h->GetNbinsX();++ix) {
          const double pt=h->GetXaxis()->GetBinCenter(ix);
          if(pt<edges[p]||pt>=edges[p+1]) continue;
          for(int id:fl.second) for(int b=0;b<20;++b) for(int m=0;m<10;++m)
            sums[b][m]+=h->GetBinContent(ix,h->GetYaxis()->FindBin(id),b*10+m+1);
        }
        double all[10]={0.}; for(auto &row:sums) for(int m=0;m<10;++m) all[m]+=row[m];
        if(all[0]<=0.) continue;
        auto calc=[](const double *a) {
          if(a[0]<=0.) return Estimates();
          return estimate(a[1]/a[0],a[2]/a[0],a[3]/a[0],a[4]/a[0],
                          a[6]/a[0],a[7]/a[0],a[8]/a[0],a[9]/a[0]);
        };
        auto r=calc(all); double avg=0.,sq=0.; int good=0;
        for(auto &row:sums) {
          double a[10]; for(int m=0;m<10;++m) a[m]=all[m]-row[m];
          const double b=calc(a).slope;
          if(std::isfinite(b)) {avg+=b;sq+=b*b;++good;}
        }
        const double err=good==20?std::sqrt(std::max(0.,19./20.*(sq-avg*avg/20.))):NAN;
        out<<sample<<'\t'<<mode<<'\t'<<label<<'\t'<<fl.first<<'\t'<<edges[p]<<'\t'
           <<edges[p+1]<<'\t'<<all[0]<<'\t'<<r.slope<<'\t'<<err<<'\t'<<r.mean<<'\t'
           <<r.affine<<'\t'<<r.intercept<<'\t'<<r.bias<<'\t'<<r.mixing<<'\n';
      }
    }
  }
}
} // namespace

void studyTaggerRu(const char *dataName="rootfiles/zjet_DATA.root",
                   const char *mcName="rootfiles/zjet_MC.root",
                   const char *directory="output/taggerRu") {
  using namespace TaggerRuStudy;
  std::unique_ptr<TFile> data(TFile::Open(dataName)),mc(TFile::Open(mcName));
  if(!data||!mc||data->IsZombie()||mc->IsZombie()) throw std::runtime_error("Missing inputs");
  gSystem->mkdir(directory,true); const std::string out=directory;
  setTDRStyle(); extraText="Work in progress"; extraText2="";
  lumi_136TeV="Run 2024I + Summer24 DY";
  TFile results((out+"/study.root").c_str(),"RECREATE");
  TObjString((std::string("data=")+dataName+"; mc="+mcName+
    "; parallel, barrel; truth 0 kept separate; old-profile ratios have no covariance errors").c_str()).Write("inputs");
  std::ofstream table(out+"/ru_estimators.tsv");
  std::ofstream summary(out+"/summary_table.tex");
  summary<<"\\begin{tabular}{lrrrr}\\toprule Flavor & $R_0$ & $R_{\\rm mean}$ & $b$ & $a$ \\\\\n\\midrule\n";
  table<<"flavor\tlow\thigh\tsumw\tslope0\tmean_ratio\taffine_slope\tintercept\tmean_fu_bias\n";
  std::vector<TGraph*> slope,means,bias;
  std::vector<std::string> labels;
  for(const auto &fl:flavors) {
    auto *s=new TGraph(),*m=new TGraph(),*b=new TGraph();
    s->SetName(("slope_"+fl.first).c_str());m->SetName(("mean_"+fl.first).c_str());
    b->SetName(("bias_"+fl.first).c_str());
    for(int p=0;p<8;++p) {
      auto get=[&](const char *name){return readComponentAcrossReco(mc.get(),name,"tc",fl.second,edges[p],edges[p+1]);};
      auto x=get("genmu"),y=get("recoumatched"),xx=get("genmu2"),xy=get("recogenmu");
      if(!x.valid||!y.valid||!xx.valid||!xy.valid) continue;
      auto r=estimate(x.mean,y.mean,xx.mean,xy.mean);
      table<<fl.first<<'\t'<<edges[p]<<'\t'<<edges[p+1]<<'\t'<<x.sumWeights<<'\t'
           <<r.slope<<'\t'<<r.mean<<'\t'<<r.affine<<'\t'<<r.intercept<<'\t'<<r.bias<<'\n';
      if(p==1 && fl.first!="unknown") summary<<std::fixed<<std::setprecision(3)<<fl.first
        <<" & "<<r.slope<<" & "<<r.mean<<" & "<<r.affine<<" & "<<r.intercept<<" \\\\\n";
      const double pt=std::sqrt(edges[p]*edges[p+1]);
      point(*s,pt,r.slope);point(*m,pt,r.mean);point(*b,pt,100*r.bias);
    }
    s->Write();m->Write();b->Write();
    if(fl.first!="unknown") {slope.push_back(s);means.push_back(m);bias.push_back(b);labels.push_back(fl.first);}
  }
  summary<<"\\bottomrule\\end{tabular}\n";summary.close();
  plot(out,"ru_slope","R_{u}: <yx>/<x^{2}> (MC)",.2,.8,slope,labels);
  plot(out,"ru_mean","R_{u}: <y>/<x> (MC)",0,5.,means,labels);
  plot(out,"ru_bias","100 [<y>/R_{u,slope} - <x>]",-5,20,bias,labels);
  std::ofstream closure(out+"/data_closure_proxy.tsv");
  closure<<"sample\tlow\thigh\tslope0\taffine_slope\tintercept\n";
  std::vector<TGraph*> proxies;
  for(auto *file:{data.get(),mc.get()}) {
    auto *g=new TGraph();
    for(int p=0;p<8;++p) {
      auto get=[&](const char *name){return readComponentAcrossReco(file,name,"tc",
        file==data.get()?std::vector<int>{0}:std::vector<int>{0,1,3,4,5,6},edges[p],edges[p+1]);};
      auto x=get("fuclosure"),y=get("mu"),xx=get("fuclosure2"),xy=get("mufuclosure");
      if(!x.valid||!y.valid||!xx.valid||!xy.valid)continue;
      auto r=estimate(x.mean,y.mean,xx.mean,xy.mean);
      closure<<(file==data.get()?"data":"mc")<<'\t'<<edges[p]<<'\t'<<edges[p+1]<<'\t'
             <<r.slope<<'\t'<<r.affine<<'\t'<<r.intercept<<'\n';
      point(*g,std::sqrt(edges[p]*edges[p+1]),r.slope);
    }
    proxies.push_back(g);
  }
  plot(out,"closure_proxy","Closure-proxy slope (not R_{u})",0,.4,
       proxies,{"Data proxy (not R_{u})","MC proxy (not R_{u})"});
  std::ofstream scores(out+"/score_audit.tsv");
  scores<<"tagger\tlow\thigh\tcut\tg_eff\tq_eff\tg_mpf1_shift_percent\tdata_mpf1\tmc_mpf1\n";
  std::vector<TGraph*> sculpt; std::vector<std::string> tagLabels;
  auto *reference=dynamic_cast<TH3D*>(mc->Get("FlavorMatrix/taggerAudit/h3counts_pnetqvg"));
  if(!reference) throw std::runtime_error("Missing current score audit");
  for(const std::string tag:{"deepjet","pnet","upart"}) {
    const bool reverse=tag=="deepjet";
    auto *c=dynamic_cast<TH3D*>(mc->Get(("FlavorMatrix/taggerAudit/h3counts_"+tag+"qvg").c_str()));
    auto *m=dynamic_cast<TProfile3D*>(mc->Get(("FlavorMatrix/taggerAudit/p3m2_"+tag+"qvg").c_str()));
    auto *d=dynamic_cast<TProfile3D*>(data->Get(("FlavorMatrix/taggerAudit/p3m2_"+tag+"qvg").c_str()));
    if(!c||!m||!d) throw std::runtime_error("Missing tagger scores");
    auto *g=new TGraph();g->SetName(("sculpt_"+tag).c_str());
    for(int p=0;p<8;++p) {
      const double lo=edges[p],hi=edges[p+1];
      const double ref=counts(reference,lo,hi,{6},.3)/counts(reference,lo,hi,{6},1.);
      double best=.025,dist=1.e99;
      for(int k=1;k<=40;++k) {
        const double e=counts(c,lo,hi,{6},k*.025,reverse)/counts(c,lo,hi,{6},1.,reverse);
        if(std::abs(e-ref)<dist) {dist=std::abs(e-ref);best=k*.025;}
      }
      const double base=scoreMean(m,lo,hi,{6},1.);
      const double shift=100*(scoreMean(m,lo,hi,{6},best,reverse)-base);
      point(*g,std::sqrt(lo*hi),shift);
      scores<<tag<<'\t'<<lo<<'\t'<<hi<<'\t'<<best<<'\t'
        <<counts(c,lo,hi,{6},best,reverse)/counts(c,lo,hi,{6},1.)<<'\t'
        <<counts(c,lo,hi,{1,3},best,reverse)/counts(c,lo,hi,{1,3},1.)<<'\t'<<shift<<'\t'
        <<scoreMean(d,lo,hi,{0},best,reverse)<<'\t'<<scoreMean(m,lo,hi,{0,1,3,4,5,6},best,reverse)<<'\n';
    }
    g->Write();sculpt.push_back(g);tagLabels.push_back(tag);
  }
  plot(out,"score_sculpting","100 [<m_{2}>_{g,tag} - <m_{2}>_{g,all}]",-2,3,sculpt,tagLabels);
  std::ofstream regress(out+"/block_regressions.tsv");
  regress<<"sample\tmode\tlabel\tflavor\tlow\thigh\tsumw\tslope0\tblock_error\tmean_ratio\taffine_slope\tintercept\tmean_fu_bias\tpartial_u_slope\n";
  blocks(*mc,"mc",regress);blocks(*data,"data",regress);
  std::cout<<"Wrote "<<out<<"; old-production graphs are central-value diagnostics, not precision fits.\n";
}
