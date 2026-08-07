// Optional shape comparison of data with the one-dimensional MC controls.
#include "TFile.h"
#include "TH1D.h"
#include "TLatex.h"
#include "TSystem.h"

#include <algorithm>
#include <cassert>

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

void drawComparison(const char *canvasName, const char *parallelName,
                    const char *transverseName, const char *subtractedName,
                    const char *xTitle, double xMin, double xMax,
                    bool logarithmic, double dataScale,
                    const char *mcFile, const char *dataFile) {
  TFile mc(mcFile,"READ");
  TFile data(dataFile,"READ");
  assert(!mc.IsZombie() && !data.IsZombie());

  TH1D *parallel = getHistogram(&mc,parallelName,"mc");
  TH1D *transverse = getHistogram(&mc,transverseName,"mc");
  TH1D *subtracted = getHistogram(&mc,subtractedName,"mc");
  TH1D *dataSubtracted = getHistogram(&data,subtractedName,"data");
  dataSubtracted->Scale(dataScale);

  const double maximum = std::max(parallel->GetMaximum(),
                                  dataSubtracted->GetMaximum());
  const double yMin = logarithmic ? std::max(0.1,1.e-5*maximum) : 0.;
  const double yMax = (logarithmic ? 20. : 1.35)*maximum;
  TH1D *frame = tdrHist(Form("frame_%s",canvasName),"Events / bin",
                        yMin,yMax,xTitle,xMin,xMax);
  lumi_136TeV = "Run2024I data and Summer24 DY";
  TCanvas *canvas = tdrCanvas(canvasName,frame,8,11,kSquare);
  if (logarithmic) gPad->SetLogy();

  tdrDraw(parallel,"HIST",kNone,kGray+1,kSolid,-1,1001,kGray);
  tdrDraw(transverse,"HIST",kNone,kRed-8,kSolid,-1,1001,kRed-9);
  tdrDraw(subtracted,"Pz",kFullCircle,kGreen+2,kSolid,-1,kNone,0,0.5);
  tdrDraw(dataSubtracted,"Pz",kOpenSquare,kBlack,kSolid,-1,kNone,0,0.8);

  TLegend *legend = tdrLeg(0.48,0.87-5*0.05,0.86,0.87);
  legend->SetHeader("|#Delta#phi|<#pi/16, 0.5<p_{T,jet}/p_{T,Z}<2");
  legend->AddEntry(parallel,"MC parallel");
  legend->AddEntry(transverse,"MC transverse");
  legend->AddEntry(subtracted,"MC parallel - transverse");
  legend->AddEntry(dataSubtracted,"Data parallel - transverse");

  TLatex text;
  text.SetNDC();
  text.SetTextSize(0.035);
  text.DrawLatex(0.48,0.56,Form("Data scale factor: %.4g",dataScale));
  gPad->RedrawAxis();
  canvas->SaveAs(Form("pdf/drawControl/%s.pdf",canvasName));
}

} // namespace

// Data is scaled once using the parallel pT yield in the requested interval.
// The same factor is then used for every observable so that the relative
// parallel, transverse, and subtracted yields are preserved.
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
                 "p_{T,Z} (GeV)",0,200,true,dataScale,mcFile,dataFile);
  drawComparison("drawControlData_c1_2_eta","h_pareta","h_traneta","h_mixeta",
                 "#eta_{jet}",-5,5,false,dataScale,mcFile,dataFile);
  drawComparison("drawControlData_c1_3_db","h_dbp","h_dbt","h_db",
                 "DB",0,2,false,dataScale,mcFile,dataFile);
  drawComparison("drawControlData_c1_4_mpf","h_mpfp","h_mpft","h_mpf",
                 "MPF",-3,4,false,dataScale,mcFile,dataFile);
}
