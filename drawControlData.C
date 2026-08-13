// Compare every pileup-subtraction stage between data and simulation.
#include "TFile.h"
#include "TH1D.h"
#include "TLatex.h"
#include "TLine.h"
#include "TSystem.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <vector>

#include "tdrstyle_mod22.C"

namespace {

TH1D *getHistogram(TFile *file, const char *name, const char *suffix) {
  TH1D *source = dynamic_cast<TH1D*>(file->Get(Form("control/%s",name)));
  assert(source);
  TH1D *copy = dynamic_cast<TH1D*>(source->Clone(Form("%s_%s",name,suffix)));
  assert(copy);
  copy->SetDirectory(nullptr);
  return copy;
}

TH1D *makeRatio(const TH1D *data, const TH1D *mc, const char *name) {
  assert(data->GetNbinsX()==mc->GetNbinsX());
  TH1D *ratio = dynamic_cast<TH1D*>(data->Clone(name));
  assert(ratio);
  ratio->Reset();
  ratio->SetDirectory(nullptr);
  for (int bin = 1; bin != ratio->GetNbinsX()+1; ++bin) {
    const double d = data->GetBinContent(bin);
    const double m = mc->GetBinContent(bin);
    const double ed = data->GetBinError(bin);
    const double em = mc->GetBinError(bin);
    if (m==0. || !std::isfinite(d) || !std::isfinite(m)) continue;
    ratio->SetBinContent(bin,d/m);
    ratio->SetBinError(bin,std::sqrt(ed*ed/(m*m) +
                                    d*d*em*em/(m*m*m*m)));
  }
  return ratio;
}

void drawComparison(const char *canvasName, const char *parallelName,
                    const char *transverseName, const char *subtractedName,
                    const char *xTitle, double xMin, double xMax,
                    bool logarithmic, double dataScale,
                    double ratioMin, double ratioMax,
                    double normalizationPtMin, double normalizationPtMax,
                    const char *mcFile, const char *dataFile) {
  TFile mc(mcFile,"READ");
  TFile data(dataFile,"READ");
  assert(!mc.IsZombie() && !data.IsZombie());

  TH1D *mcParallel = getHistogram(&mc,parallelName,"mc_parallel");
  TH1D *mcTransverse = getHistogram(&mc,transverseName,"mc_transverse");
  TH1D *mcSubtracted = getHistogram(&mc,subtractedName,"mc_subtracted");
  TH1D *dataParallel = getHistogram(&data,parallelName,"data_parallel");
  TH1D *dataTransverse = getHistogram(&data,transverseName,"data_transverse");
  TH1D *dataSubtracted = getHistogram(&data,subtractedName,"data_subtracted");
  for (TH1D *histogram : {dataParallel,dataTransverse,dataSubtracted})
    histogram->Scale(dataScale);

  double maximum = 0.;
  for (TH1D *histogram : {mcParallel,mcTransverse,mcSubtracted,
                          dataParallel,dataTransverse,dataSubtracted})
    maximum = std::max(maximum,histogram->GetMaximum());
  const double yMin = logarithmic ? std::max(0.1,1.e-5*maximum) : 0.;
  const double yMax = (logarithmic ? 20. : 1.45)*maximum;
  TH1D *upperFrame = tdrHist(Form("upper_%s",canvasName),"Events / bin",
                             yMin,yMax,xTitle,xMin,xMax);
  TH1D *lowerFrame = tdrHist(Form("lower_%s",canvasName),"Data / MC",
                             ratioMin,ratioMax,xTitle,xMin,xMax);
  lumi_136TeV = "Run2024I data and Summer24 DY";
  TCanvas *canvas = tdrDiCanvas(canvasName,upperFrame,lowerFrame,8,11);

  canvas->cd(1);
  if (logarithmic) gPad->SetLogy();
  tdrDraw(mcParallel,"HIST",kNone,kGray+2,kSolid,kGray+2,kNone,0,1.,2.);
  tdrDraw(mcTransverse,"HIST",kNone,kRed-7,kSolid,kRed-7,kNone,0,1.,2.);
  tdrDraw(mcSubtracted,"HIST",kNone,kGreen+2,kSolid,kGreen+2,kNone,0,1.,2.);
  tdrDraw(dataParallel,"Pz",kOpenCircle,kGray+2,kSolid,kGray+2,kNone,0,0.7);
  tdrDraw(dataTransverse,"Pz",kOpenTriangleUp,kRed+1,kSolid,kRed+1,kNone,0,0.7);
  tdrDraw(dataSubtracted,"Pz",kOpenSquare,kGreen+3,kSolid,kGreen+3,kNone,0,0.8);

  TLegend *legend = tdrLeg(0.46,0.87-6*0.045,0.88,0.87);
  legend->SetHeader("|#Delta#phi|<#pi/16, 0.5<p_{T,jet}/p_{T,Z}<2");
  legend->AddEntry(mcParallel,"MC parallel","L");
  legend->AddEntry(dataParallel,"Data parallel","P");
  legend->AddEntry(mcTransverse,"MC transverse","L");
  legend->AddEntry(dataTransverse,"Data transverse","P");
  legend->AddEntry(mcSubtracted,"MC parallel - transverse","L");
  legend->AddEntry(dataSubtracted,"Data parallel - transverse","P");

  TLatex text;
  text.SetNDC();
  text.SetTextSize(0.032);
  text.DrawLatex(0.46,0.52,Form("Common data scale: %.4g",dataScale));
  text.DrawLatex(0.46,0.48,Form("Normalized for %.0f < p_{T} < %.0f GeV",
                                normalizationPtMin,normalizationPtMax));
  gPad->RedrawAxis();

  TH1D *ratioParallel = makeRatio(dataParallel,mcParallel,
                                  Form("ratio_parallel_%s",canvasName));
  TH1D *ratioTransverse = makeRatio(dataTransverse,mcTransverse,
                                    Form("ratio_transverse_%s",canvasName));
  TH1D *ratioSubtracted = makeRatio(dataSubtracted,mcSubtracted,
                                    Form("ratio_subtracted_%s",canvasName));
  canvas->cd(2);
  TLine unity(xMin,1.,xMax,1.);
  unity.SetLineStyle(kDashed);
  unity.SetLineColor(kGray+2);
  unity.Draw("SAME");
  tdrDraw(ratioParallel,"Pz",kOpenCircle,kGray+2,kSolid,kGray+2,kNone,0,0.55);
  tdrDraw(ratioTransverse,"Pz",kOpenTriangleUp,kRed+1,kSolid,kRed+1,kNone,0,0.55);
  tdrDraw(ratioSubtracted,"Pz",kOpenSquare,kGreen+3,kSolid,kGreen+3,kNone,0,0.65);
  gPad->RedrawAxis();
  canvas->SaveAs(Form("pdf/drawControl/%s.pdf",canvasName));
}

} // namespace

// Data is scaled once using the parallel pT yield in the requested interval.
// The same factor is used at all three stages and for every observable.
void drawControlData(const char *mcFile="rootfiles/zjet_MC.root",
                     const char *dataFile="rootfiles/zjet_DATA.root",
                     double normalizationPtMin=30.,
                     double normalizationPtMax=200.) {
  gSystem->mkdir("pdf/drawControl",true);
  setTDRStyle();

  TFile mc(mcFile,"READ");
  TFile data(dataFile,"READ");
  assert(!mc.IsZombie() && !data.IsZombie());
  TH1D *mcPt = getHistogram(&mc,"h_parpt","mc_norm");
  TH1D *dataPt = getHistogram(&data,"h_parpt","data_norm");
  const int first = mcPt->GetXaxis()->FindFixBin(normalizationPtMin+1.e-6);
  const int last = mcPt->GetXaxis()->FindFixBin(normalizationPtMax-1.e-6);
  const double mcYield = mcPt->Integral(first,last);
  const double dataYield = dataPt->Integral(first,last);
  assert(mcYield>0 && dataYield>0);
  const double dataScale = mcYield/dataYield;

  drawComparison("drawControlData_c1_1_pt","h_parpt","h_tranpt","h_mixpt",
                 "p_{T,jet} (GeV)",0,200,true,dataScale,0.5,1.5,
                 normalizationPtMin,normalizationPtMax,mcFile,dataFile);
  drawComparison("drawControlData_c1_2_eta","h_pareta","h_traneta","h_mixeta",
                 "#eta_{jet}",-5,5,false,dataScale,0.5,1.5,
                 normalizationPtMin,normalizationPtMax,mcFile,dataFile);
  drawComparison("drawControlData_c1_3_db","h_dbp","h_dbt","h_db",
                 "DB",0,2,false,dataScale,0.75,1.25,
                 normalizationPtMin,normalizationPtMax,mcFile,dataFile);
  drawComparison("drawControlData_c1_4_mpf","h_mpfp","h_mpft","h_mpf",
                 "MPF",-3,4,false,dataScale,0.75,1.25,
                 normalizationPtMin,normalizationPtMax,mcFile,dataFile);
}
