#include "TCanvas.h"
#include "TFile.h"
#include "TGraphErrors.h"
#include "TH1D.h"
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

TH1D *getHistogram(TFile *file, const std::string &name) {
  TH1D *histogram = dynamic_cast<TH1D*>(file->Get(name.c_str()));
  assert(histogram);
  return histogram;
}

void style(TProfile *profile, Color_t color, Style_t marker,
           Style_t lineStyle=1) {
  profile->SetLineColor(color);
  profile->SetMarkerColor(color);
  profile->SetMarkerStyle(marker);
  profile->SetMarkerSize(0.8);
  profile->SetLineStyle(lineStyle);
}

void fillGraph(TGraphErrors &graph, TProfile *profile, Color_t color,
               Style_t marker) {
  graph.SetLineColor(color);
  graph.SetMarkerColor(color);
  graph.SetMarkerStyle(marker);
  graph.SetMarkerSize(0.8);
  int point = 0;
  for (int bin = 1; bin <= profile->GetNbinsX(); ++bin) {
    // A non-positive signed denominator does not define a useful fraction.
    if (profile->GetBinEntries(bin)<=0.) continue;
    graph.SetPoint(point,profile->GetBinCenter(bin),profile->GetBinContent(bin));
    graph.SetPointError(point,0.,profile->GetBinError(bin));
    ++point;
  }
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

void drawMatchDefinition(TFile *mc, const char *outputDirectory) {
  TProfile *analysis = getProfile(
    mc,"control/p_pujet_fraction_vs_ptz_subtracted");
  TProfile *extraCuts = getProfile(
    mc,"control/p_extra_match_cuts_unmatched_fraction_vs_ptz_subtracted");
  style(analysis,kBlack,kFullCircle);
  style(extraCuts,kBlue+1,kOpenSquare);

  TCanvas canvas("c_match_definition","",700,650);
  canvas.SetBottomMargin(0.13);
  analysis->SetTitle("");
  analysis->GetXaxis()->SetTitle("p_{T,Z} (GeV)");
  analysis->GetYaxis()->SetTitle("Subtracted unmatched fraction");
  analysis->GetYaxis()->SetRangeUser(-0.1,1.05);
  analysis->Draw("E1");
  extraCuts->Draw("E1 SAME");
  TLine zero(0.,0.,200.,0.);
  zero.SetLineStyle(3);
  zero.Draw();
  TLegend legend(0.49,0.74,0.88,0.88);
  legend.SetBorderSize(0);
  legend.AddEntry(analysis,"No valid genJetIdx","lp");
  legend.AddEntry(extraCuts,"Also require gen p_{T}>8, #DeltaR<0.4","lp");
  legend.Draw();
  canvas.SaveAs(Form("%s/match_definition_vs_ptz.pdf",outputDirectory));
}

void drawEtaDependence(TFile *mc, const char *outputDirectory) {
  TProfile *central = getProfile(
    mc,"control/p_pujet_fraction_vs_ptz_subtracted_central");
  TProfile *endcap = getProfile(
    mc,"control/p_pujet_fraction_vs_ptz_subtracted_endcap");
  TProfile *forward = getProfile(
    mc,"control/p_pujet_fraction_vs_ptz_subtracted_forward");
  TGraphErrors centralGraph, endcapGraph, forwardGraph;
  fillGraph(centralGraph,central,kBlack,kFullCircle);
  fillGraph(endcapGraph,endcap,kBlue+1,kOpenSquare);
  fillGraph(forwardGraph,forward,kRed+1,kOpenTriangleUp);

  TCanvas canvas("c_eta_dependence","",700,650);
  canvas.SetBottomMargin(0.13);
  TH1D frame("h_eta_frame",";p_{T,Z} (GeV);Subtracted unmatched fraction",
             100,0.,200.);
  frame.GetYaxis()->SetRangeUser(-0.25,0.5);
  frame.Draw("AXIS");
  centralGraph.Draw("P E SAME");
  endcapGraph.Draw("P E SAME");
  forwardGraph.Draw("P E SAME");
  TLine zero(0.,0.,200.,0.);
  zero.SetLineStyle(3);
  zero.Draw();
  TLegend legend(0.58,0.70,0.88,0.88);
  legend.SetBorderSize(0);
  legend.AddEntry(&centralGraph,"|#eta| < 1.3","lp");
  legend.AddEntry(&endcapGraph,"1.3 #leq |#eta| < 2.5","lp");
  legend.AddEntry(&forwardGraph,"|#eta| #geq 2.5","lp");
  legend.Draw();
  canvas.SaveAs(Form("%s/pujet_fraction_eta_vs_ptz.pdf",outputDirectory));
}

void drawCountCheck(TFile *mc, const char *outputDirectory) {
  TH1D *unmatched = getHistogram(
    mc,"control/h_unmatched_jets_vs_ptz_subtracted");
  TH1D *all = getHistogram(mc,"control/h_all_jets_vs_ptz_subtracted");
  TH1D *ratio = dynamic_cast<TH1D*>(unmatched->Clone("h_count_ratio"));
  ratio->SetDirectory(0);
  ratio->Divide(all);
  TProfile *profile = getProfile(
    mc,"control/p_pujet_fraction_vs_ptz_subtracted");
  ratio->SetLineColor(kBlue+1);
  ratio->SetMarkerColor(kBlue+1);
  ratio->SetMarkerStyle(kOpenSquare);
  style(profile,kBlack,kFullCircle);

  TCanvas canvas("c_count_check","",700,650);
  canvas.SetBottomMargin(0.13);
  ratio->SetTitle("");
  ratio->GetXaxis()->SetTitle("p_{T,Z} (GeV)");
  ratio->GetYaxis()->SetTitle("Subtracted unmatched fraction");
  ratio->GetYaxis()->SetRangeUser(-0.2,1.05);
  ratio->Draw("E1");
  profile->Draw("E1 SAME");
  TLine zero(0.,0.,200.,0.);
  zero.SetLineStyle(3);
  zero.Draw();
  TLegend legend(0.50,0.74,0.88,0.88);
  legend.SetBorderSize(0);
  legend.AddEntry(profile,"TProfile with signed weights","lp");
  legend.AddEntry(ratio,"Explicit unmatched/all yields","lp");
  legend.Draw();
  canvas.SaveAs(Form("%s/pujet_fraction_count_check_vs_ptz.pdf",
                     outputDirectory));
  delete ratio;
}

void drawSelectionEfficiency(TFile *mc, TFile *data, const char *histogram,
                             const char *outputDirectory) {
  TH1D *mcSource = getHistogram(mc,Form("control/%s",histogram));
  TH1D *dataSource = getHistogram(data,Form("control/%s",histogram));
  TH1D *mcEfficiency = dynamic_cast<TH1D*>(mcSource->Clone(
    Form("%s_mc_efficiency",histogram)));
  TH1D *dataEfficiency = dynamic_cast<TH1D*>(dataSource->Clone(
    Form("%s_data_efficiency",histogram)));
  mcEfficiency->SetDirectory(0);
  dataEfficiency->SetDirectory(0);
  if (mcEfficiency->GetBinContent(1)!=0.)
    mcEfficiency->Scale(1./mcEfficiency->GetBinContent(1));
  if (dataEfficiency->GetBinContent(1)!=0.)
    dataEfficiency->Scale(1./dataEfficiency->GetBinContent(1));
  TGraphErrors mcGraph, dataGraph;
  for (int bin = 1; bin <= mcEfficiency->GetNbinsX(); ++bin) {
    const int point = bin-1;
    mcGraph.SetPoint(point,mcEfficiency->GetBinCenter(bin),
                     mcEfficiency->GetBinContent(bin));
    mcGraph.SetPointError(point,0.,mcEfficiency->GetBinError(bin));
    dataGraph.SetPoint(point,dataEfficiency->GetBinCenter(bin),
                       dataEfficiency->GetBinContent(bin));
    dataGraph.SetPointError(point,0.,dataEfficiency->GetBinError(bin));
  }
  mcGraph.SetLineColor(kBlue+1);
  mcGraph.SetMarkerColor(kBlue+1);
  mcGraph.SetMarkerStyle(kOpenSquare);
  dataGraph.SetLineColor(kBlack);
  dataGraph.SetMarkerColor(kBlack);
  dataGraph.SetMarkerStyle(kFullCircle);

  TCanvas canvas(Form("c_%s",histogram),"",850,650);
  canvas.SetBottomMargin(0.32);
  canvas.SetLeftMargin(0.12);
  TH1D frame(Form("%s_frame",histogram),"",mcEfficiency->GetNbinsX(),
             mcEfficiency->GetXaxis()->GetXmin(),
             mcEfficiency->GetXaxis()->GetXmax());
  for (int bin = 1; bin <= frame.GetNbinsX(); ++bin)
    frame.GetXaxis()->SetBinLabel(bin,mcEfficiency->GetXaxis()->GetBinLabel(bin));
  frame.GetYaxis()->SetTitle("Fraction of first bin");
  frame.GetYaxis()->SetRangeUser(0.,1.1);
  frame.GetXaxis()->LabelsOption("v");
  frame.Draw("AXIS");
  mcGraph.Draw("P E SAME");
  dataGraph.Draw("P E SAME");
  TLegend legend(0.71,0.76,0.88,0.88);
  legend.SetBorderSize(0);
  legend.AddEntry(&mcGraph,"MC","lp");
  legend.AddEntry(&dataGraph,"Data","lp");
  legend.Draw();
  canvas.SaveAs(Form("%s/%s.pdf",outputDirectory,histogram));
  delete mcEfficiency;
  delete dataEfficiency;
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

  drawMatchDefinition(&mc,outputDirectory);
  drawEtaDependence(&mc,outputDirectory);
  drawCountCheck(&mc,outputDirectory);
  drawSelectionEfficiency(&mc,&data,"h_cutflow",outputDirectory);
  drawSelectionEfficiency(&mc,&data,"h_muon_selection",outputDirectory);
  drawSelectionEfficiency(&mc,&data,"h_probe_veto",outputDirectory);
}
