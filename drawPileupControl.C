#include "TCanvas.h"
#include "TFile.h"
#include "TLegend.h"
#include "TLine.h"
#include "TProfile.h"
#include "TStyle.h"
#include "TSystem.h"

#include <cassert>
#include <string>

namespace {

TProfile *getProfile(TFile *file, const std::string &name, bool required=true) {
  TProfile *profile = dynamic_cast<TProfile*>(file->Get(name.c_str()));
  if (required) assert(profile);
  return profile;
}

void style(TProfile *profile, Color_t color, Style_t marker,
           Style_t lineStyle=1) {
  profile->SetLineColor(color);
  profile->SetMarkerColor(color);
  profile->SetMarkerStyle(marker);
  profile->SetMarkerSize(0.8);
  profile->SetLineStyle(lineStyle);
}

void drawResponse(TFile *mc, TFile *data, const char *response,
                  const char *observable, const char *xTitle,
                  const char *outputDirectory) {
  TProfile *mcParallel = getProfile(mc,Form("control/p_%s_vs_%s_parallel",
                                            response,observable));
  TProfile *mcTransverse = getProfile(mc,Form("control/p_%s_vs_%s_transverse",
                                              response,observable));
  TProfile *mcSubtracted = getProfile(mc,Form("control/p_%s_vs_%s_subtracted",
                                              response,observable));
  TProfile *dataSubtracted = getProfile(data,Form("control/p_%s_vs_%s_subtracted",
                                                  response,observable));

  style(mcParallel,kGray+2,kOpenCircle,2);
  style(mcTransverse,kRed+1,kOpenSquare,2);
  style(mcSubtracted,kBlue+1,kFullCircle);
  style(dataSubtracted,kBlack,kFullSquare);

  TCanvas canvas(Form("c_%s_%s",response,observable),"",700,650);
  canvas.SetBottomMargin(0.13);
  mcParallel->SetTitle("");
  mcParallel->GetXaxis()->SetTitle(xTitle);
  mcParallel->GetYaxis()->SetTitle(std::string(response)=="db" ?
                                   "Direct balance" : "Hybrid MPF");
  mcParallel->GetXaxis()->CenterTitle();
  mcParallel->GetXaxis()->SetTitleOffset(1.15);
  mcParallel->GetYaxis()->CenterTitle();
  mcParallel->GetYaxis()->SetRangeUser(0.5,1.5);
  mcParallel->Draw("E1");
  mcTransverse->Draw("E1 SAME");
  mcSubtracted->Draw("E1 SAME");
  if (dataSubtracted->GetEntries()>0) dataSubtracted->Draw("E1 SAME");

  TLine unity(mcParallel->GetXaxis()->GetXmin(),1.,
              mcParallel->GetXaxis()->GetXmax(),1.);
  unity.SetLineStyle(3);
  unity.Draw();

  TLegend legend(0.53,0.70,0.88,0.88);
  legend.SetBorderSize(0);
  legend.AddEntry(mcParallel,"MC parallel","lp");
  legend.AddEntry(mcTransverse,"MC transverse average","lp");
  legend.AddEntry(mcSubtracted,"MC subtracted","lp");
  if (dataSubtracted->GetEntries()>0)
    legend.AddEntry(dataSubtracted,"Data subtracted","lp");
  legend.Draw();

  canvas.SaveAs(Form("%s/%s_vs_%s.pdf",outputDirectory,response,observable));
}

void drawPileupFraction(TFile *mc, const char *observable, const char *xTitle,
                        const char *outputDirectory) {
  TProfile *parallel = getProfile(mc,Form("control/p_pujet_fraction_vs_%s_parallel",
                                          observable));
  TProfile *transverse = getProfile(mc,Form("control/p_pujet_fraction_vs_%s_transverse",
                                            observable));
  TProfile *subtracted = getProfile(mc,Form("control/p_pujet_fraction_vs_%s_subtracted",
                                            observable));
  style(parallel,kBlue+1,kFullCircle);
  style(transverse,kRed+1,kFullSquare);
  style(subtracted,kBlack,kOpenCircle);

  TCanvas canvas(Form("c_pujet_%s",observable),"",700,650);
  canvas.SetBottomMargin(0.13);
  parallel->SetTitle("");
  parallel->GetXaxis()->SetTitle(xTitle);
  parallel->GetYaxis()->SetTitle("Unmatched (pileup-jet) fraction");
  parallel->GetXaxis()->CenterTitle();
  parallel->GetXaxis()->SetTitleOffset(1.15);
  parallel->GetYaxis()->CenterTitle();
  parallel->GetYaxis()->SetRangeUser(0.,1.05);
  parallel->Draw("E1");
  transverse->Draw("E1 SAME");
  subtracted->Draw("E1 SAME");

  TLegend legend(0.58,0.72,0.88,0.88);
  legend.SetBorderSize(0);
  legend.AddEntry(parallel,"Parallel","lp");
  legend.AddEntry(transverse,"Transverse average","lp");
  legend.AddEntry(subtracted,"Subtracted","lp");
  legend.Draw();
  canvas.SaveAs(Form("%s/pujet_fraction_vs_%s.pdf",outputDirectory,observable));
}

void drawTruthResponse(TFile *mc, const char *response, const char *region,
                       const char *outputDirectory) {
  TProfile *matched = getProfile(mc,Form("control/p_%s_vs_ptz_matched_%s",
                                         response,region));
  TProfile *pileup = getProfile(mc,Form("control/p_%s_vs_ptz_pileup_%s",
                                        response,region));
  style(matched,kBlue+1,kFullCircle);
  style(pileup,kRed+1,kFullSquare);

  TCanvas canvas(Form("c_%s_truth_%s",response,region),"",700,650);
  canvas.SetBottomMargin(0.13);
  matched->SetTitle("");
  matched->GetXaxis()->SetTitle("p_{T,Z} (GeV)");
  matched->GetYaxis()->SetTitle(std::string(response)=="db" ?
                                "Direct balance" : "Hybrid MPF");
  matched->GetXaxis()->CenterTitle();
  matched->GetXaxis()->SetTitleOffset(1.15);
  matched->GetYaxis()->CenterTitle();
  matched->GetYaxis()->SetRangeUser(0.,3.);
  matched->Draw("E1");
  pileup->Draw("E1 SAME");

  TLegend legend(0.58,0.75,0.88,0.88);
  legend.SetBorderSize(0);
  legend.AddEntry(matched,"Truth matched","lp");
  legend.AddEntry(pileup,"Pileup / unmatched","lp");
  legend.Draw();
  canvas.SaveAs(Form("%s/%s_truth_%s.pdf",outputDirectory,response,region));
}

} // namespace

void drawPileupControl(
  const char *mcFile="rootfiles/zjet_MC.root",
  const char *dataFile="rootfiles/zjet_DATA.root",
  const char *outputDirectory="pdf/drawPileupControl") {
  gStyle->SetOptStat(0);
  gSystem->mkdir(outputDirectory,true);

  TFile mc(mcFile,"READ");
  TFile data(dataFile,"READ");
  assert(!mc.IsZombie());
  assert(!data.IsZombie());

  const char *observables[] = {"npvs","rho","mu"};
  const char *titles[] = {"N_{PV}","#rho (GeV)","#mu"};
  for (int io = 0; io != 3; ++io) {
    drawResponse(&mc,&data,"db",observables[io],titles[io],outputDirectory);
    drawResponse(&mc,&data,"mpf",observables[io],titles[io],outputDirectory);
  }

  const char *truthObservables[] = {"ptz","npvs","rho","mu"};
  const char *truthTitles[] = {"p_{T,Z} (GeV)","N_{PV}","#rho (GeV)","#mu"};
  for (int io = 0; io != 4; ++io)
    drawPileupFraction(&mc,truthObservables[io],truthTitles[io],outputDirectory);

  for (const char *region : {"parallel","transverse"})
    for (const char *response : {"db","mpf"})
      drawTruthResponse(&mc,response,region,outputDirectory);
}
