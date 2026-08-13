// Compact, publication-style summaries of the all-pairs Z+jet construction.
#include "TCanvas.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLine.h"
#include "TStyle.h"
#include "TSystem.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <string>

namespace {

TH2D *getMap(TFile &file, const char *name, const char *suffix) {
  TH2D *source = dynamic_cast<TH2D*>(file.Get(Form("control/%s",name)));
  assert(source);
  TH2D *copy = dynamic_cast<TH2D*>(source->Clone(Form("%s_%s",name,suffix)));
  assert(copy);
  copy->SetDirectory(nullptr);
  return copy;
}

TH1D *slice(TH2D *source, const char *name, double minimum, double maximum) {
  const int first = source->GetXaxis()->FindFixBin(minimum+1.e-6);
  const int last = source->GetXaxis()->FindFixBin(maximum-1.e-6);
  TH1D *projection = source->ProjectionY(name,first,last,"e");
  projection->SetDirectory(nullptr);
  return projection;
}

void annotate(const char *sample, const char *region) {
  TLatex text;
  text.SetNDC();
  text.SetTextFont(42);
  text.SetTextSize(0.040);
  text.DrawLatex(0.13,0.93,"CMS #it{Work in progress}");
  text.SetTextSize(0.034);
  text.DrawLatex(0.13,0.87,sample);
  text.SetTextSize(0.040);
  text.DrawLatex(0.58,0.87,region);
}

void drawRegions(TFile &data, const char *outputDirectory) {
  TH2D *parallel = getMap(data,"h2_dbp","data_parallel");
  TH2D *transverse = getMap(data,"h2_dbt","data_transverse");
  TH2D *subtracted = getMap(data,"h2_db","data_subtracted");
  TCanvas canvas("c_regions_db","",1800,600);
  canvas.Divide(3,1,0.002,0.002);
  TH2D *histograms[] = {parallel,transverse,subtracted};
  const char *labels[] = {"Parallel","Average of #pm#pi/2","Parallel - transverse"};
  for (int index = 0; index != 3; ++index) {
    canvas.cd(index+1);
    gPad->SetRightMargin(0.13);
    gPad->SetLeftMargin(0.12);
    gPad->SetBottomMargin(0.12);
    TH2D *histogram = histograms[index];
    histogram->SetTitle("");
    histogram->GetXaxis()->SetTitle("p_{T,Z} (GeV)");
    histogram->GetYaxis()->SetTitle("p_{T,jet} (GeV)");
    histogram->GetZaxis()->SetTitle("Pairs / bin");
    histogram->GetXaxis()->SetRangeUser(0.,100.);
    histogram->GetYaxis()->SetRangeUser(0.,100.);
    if (index==2) {
      const double extent = std::max(std::abs(histogram->GetMinimum()),
                                     std::abs(histogram->GetMaximum()));
      histogram->SetMinimum(-extent);
      histogram->SetMaximum(extent);
    }
    histogram->Draw("COLZ");
    annotate("Run2024I data",labels[index]);
  }
  canvas.SaveAs(Form("%s/regions_db.pdf",outputDirectory));
  delete parallel;
  delete transverse;
  delete subtracted;
}

void drawSampleSlice(TFile &file, const char *sample, int pad,
                     TCanvas &canvas, double ptMin, double ptMax) {
  TH2D *parallelMap = getMap(file,"h2_dbp",Form("%s_parallel",sample));
  TH2D *transverseMap = getMap(file,"h2_dbt",Form("%s_transverse",sample));
  TH2D *subtractedMap = getMap(file,"h2_db",Form("%s_subtracted",sample));
  TH1D *parallel = slice(parallelMap,Form("parallel_%s",sample),ptMin,ptMax);
  TH1D *transverse = slice(transverseMap,Form("transverse_%s",sample),ptMin,ptMax);
  TH1D *subtracted = slice(subtractedMap,Form("subtracted_%s",sample),ptMin,ptMax);
  const double normalization = parallel->Integral();
  if (normalization!=0.) {
    parallel->Scale(1./normalization);
    transverse->Scale(1./normalization);
    subtracted->Scale(1./normalization);
  }
  parallel->SetLineColor(kGray+2);
  parallel->SetLineWidth(2);
  transverse->SetLineColor(kRed+1);
  transverse->SetLineWidth(2);
  subtracted->SetLineColor(kBlue+1);
  subtracted->SetLineWidth(3);
  const double maximum = std::max({parallel->GetMaximum(),transverse->GetMaximum(),
                                   subtracted->GetMaximum()});
  const double minimum = std::min({0.,parallel->GetMinimum(),
                                   transverse->GetMinimum(),subtracted->GetMinimum()});

  canvas.cd(pad);
  gPad->SetLeftMargin(0.13);
  gPad->SetBottomMargin(0.13);
  parallel->SetTitle("");
  parallel->GetXaxis()->SetTitle("p_{T,jet} (GeV)");
  parallel->GetYaxis()->SetTitle("Pairs / parallel yield");
  parallel->GetXaxis()->SetRangeUser(0.,100.);
  parallel->GetYaxis()->SetRangeUser(minimum-0.08*maximum,1.25*maximum);
  parallel->Draw("HIST");
  transverse->Draw("HIST SAME");
  subtracted->Draw("HIST SAME");
  TLine zero(0.,0.,100.,0.);
  zero.SetLineStyle(3);
  zero.Draw();
  TLegend legend(0.56,0.68,0.88,0.86);
  legend.SetBorderSize(0);
  legend.SetFillStyle(0);
  legend.AddEntry(parallel,"Parallel","l");
  legend.AddEntry(transverse,"Average of #pm#pi/2","l");
  legend.AddEntry(subtracted,"Parallel - transverse","l");
  legend.Draw();
  annotate(sample,Form("%.0f < p_{T,Z} < %.0f GeV",ptMin,ptMax));

  delete parallelMap;
  delete transverseMap;
  delete subtractedMap;
  // The canvas owns references to the projected histograms until SaveAs in
  // drawSlices. Keep them alive for that draw; ROOT cleans them up at exit.
}

void drawSlices(TFile &data, TFile &mc, const char *outputDirectory) {
  TCanvas canvas("c_balance_slices","",1400,650);
  canvas.Divide(2,1,0.002,0.002);
  drawSampleSlice(data,"Run2024I data",1,canvas,28.,32.);
  drawSampleSlice(mc,"Summer24 DY MC",2,canvas,28.,32.);
  canvas.SaveAs(Form("%s/balance_slices.pdf",outputDirectory));
}

} // namespace

void drawZjetSummary(
  const char *dataFile="rootfiles/zjet_DATA.root",
  const char *mcFile="rootfiles/zjet_MC.root",
  const char *outputDirectory="pdf/drawZjetSummary") {
  gStyle->SetOptStat(0);
  gStyle->SetPalette(kBird);
  gSystem->mkdir(outputDirectory,true);
  TFile data(dataFile,"READ");
  TFile mc(mcFile,"READ");
  assert(!data.IsZombie() && !mc.IsZombie());
  drawRegions(data,outputDirectory);
  drawSlices(data,mc,outputDirectory);
}
