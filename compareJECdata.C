#include "TCanvas.h"
#include "TFile.h"
#include "TGraphErrors.h"
#include "TH1.h"
#include "TH1D.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TPad.h"
#include "TStyle.h"
#include "TSystem.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace {

struct ComparisonPoint {
  double x = 0.;
  double ex = 0.;
  double y = 0.;
  double ey = 0.;
};

std::vector<ComparisonPoint> readPoints(TFile &file, const std::string &path) {
  TObject *object = file.Get(path.c_str());
  std::vector<ComparisonPoint> points;
  if (TGraphErrors *graph = dynamic_cast<TGraphErrors*>(object)) {
    for (int i=0; i<graph->GetN(); ++i) {
      double x = 0.;
      double y = 0.;
      graph->GetPoint(i,x,y);
      points.push_back({x,graph->GetErrorX(i),y,graph->GetErrorY(i)});
    }
  }
  else if (TH1 *histogram = dynamic_cast<TH1*>(object)) {
    for (int i=1; i<=histogram->GetNbinsX(); ++i) {
      const double y = histogram->GetBinContent(i);
      const double ey = histogram->GetBinError(i);
      if (y==0. && ey==0.) continue;
      points.push_back({histogram->GetBinCenter(i),
                        0.5*histogram->GetBinWidth(i),y,ey});
    }
  }
  return points;
}

TGraphErrors *makeGraph(const std::vector<ComparisonPoint> &points, int color,
                        int marker) {
  TGraphErrors *graph = new TGraphErrors(points.size());
  for (size_t i=0; i<points.size(); ++i) {
    graph->SetPoint(i,points[i].x,points[i].y);
    graph->SetPointError(i,points[i].ex,points[i].ey);
  }
  graph->SetMarkerStyle(marker);
  graph->SetMarkerSize(0.85);
  graph->SetMarkerColor(color);
  graph->SetLineColor(color);
  return graph;
}

std::pair<double,double> yRange(const std::string &observable,
                                const std::string &sample) {
  if (observable=="counts") return {0.5,1.e9};
  if (observable=="rho") return {0.,60.};
  if (observable=="chf" || observable=="nef" || observable=="nhf" ||
      observable=="cef" || observable=="muf") return {0.,1.};
  if (observable=="mpfn") return {-0.6,0.3};
  if (observable=="mpfu" || observable=="mpfnu") return {-0.2,0.8};
  if (sample=="ratio") return {0.75,1.25};
  return {0.5,1.5};
}

double differenceRange(const std::string &observable) {
  if (observable=="counts") return 1.;
  if (observable=="rho") return 20.;
  return 0.25;
}

std::string objectPath(const std::string &sample,
                       const std::string &observable,
                       const std::string &channel) {
  const std::string base = sample+"/eta00-13/";
  if (observable=="counts")
    return base+"counts_"+channel+"_a100";
  return base+"orig/"+observable+"_"+channel+"_a100";
}

std::string observableLabel(const std::string &observable) {
  if (observable=="counts") return "Statistics";
  if (observable=="ptchs") return "Direct balance";
  if (observable=="mpfchs1") return "MPF";
  if (observable=="mpf1") return "MPF1";
  if (observable=="mpfn") return "MPFn";
  if (observable=="mpfu") return "MPFu";
  if (observable=="mpfnu") return "MPFnu";
  if (observable=="rjet") return "DB response";
  if (observable=="gjet") return "Generator response";
  return observable;
}

std::string channelLabel(const std::string &channel) {
  if (channel=="jetz") return "jet pT binning";
  if (channel=="zjav") return "average pT (HDM) binning";
  return "Z pT binning";
}

std::vector<ComparisonPoint> matchedDifferences(
                                      const std::vector<ComparisonPoint> &oldPoints,
                                      const std::vector<ComparisonPoint> &newPoints,
                                      bool relative, int &matched) {
  std::vector<ComparisonPoint> differences;
  matched = 0;
  for (const ComparisonPoint &newPoint : newPoints) {
    auto oldPoint = std::find_if(
      oldPoints.begin(),oldPoints.end(),[&](const ComparisonPoint &candidate) {
        return fabs(candidate.x-newPoint.x)<1.e-5*std::max(1.,newPoint.x);
      });
    if (oldPoint==oldPoints.end()) continue;
    ++matched;
    double value = newPoint.y-oldPoint->y;
    double error = std::hypot(newPoint.ey,oldPoint->ey);
    if (relative) {
      if (oldPoint->y==0.) continue;
      value /= oldPoint->y;
      error /= fabs(oldPoint->y);
    }
    differences.push_back({newPoint.x,newPoint.ex,value,error});
  }
  return differences;
}

void drawMissing(const std::string &message) {
  TLatex latex;
  latex.SetNDC();
  latex.SetTextAlign(22);
  latex.SetTextSize(0.045);
  latex.DrawLatex(0.5,0.52,message.c_str());
}

} // namespace

void compareJECdata(
  const char *oldFile="../jecsys3/rootfiles/jecdata2024I_nib1.root",
  const char *newFile="../jecsys3/rootfiles/jecdata2024I_nix.root",
  const char *outputDirectory="output/compareJECdata",
  const char *oldLabel="2024I_nib1",
  const char *newLabel="2024I_nix") {
  TFile oldInput(oldFile,"READ");
  TFile newInput(newFile,"READ");
  if (oldInput.IsZombie() || newInput.IsZombie())
    throw std::runtime_error("Could not open one or both jecdata inputs");
  gSystem->mkdir(outputDirectory,true);
  gStyle->SetOptStat(0);

  const std::vector<std::string> samples = {"data","mc","ratio"};
  const std::vector<std::string> channels = {"jetz","zjav","zjet"};
  const std::vector<std::string> observables = {
    "counts", "ptchs", "mpfchs1", "mpf1", "mpfn", "mpfu", "mpfnu",
    "rjet", "gjet", "chf", "nef", "nhf", "cef", "muf", "rho",
  };

  std::ofstream summary(std::string(outputDirectory)+"/summary.tsv");
  summary << "sample\tobservable\tchannel\told_points\tnew_points\tmatched_bins"
          << "\told_only\tnew_only\n";
  std::ofstream frames(std::string(outputDirectory)+
                       "/compareJECdata_frames.tex");
  frames << "% Generated by compareJECdata.C. Do not edit by hand.\n";

  for (const std::string &observable : observables) {
    const bool hasAllReferenceAxes =
      (observable=="counts" || observable=="ptchs" ||
       observable=="mpfchs1" || observable=="mpf1" ||
       observable=="mpfn" || observable=="mpfu" ||
       observable=="mpfnu");
    const std::vector<std::string> observableChannels =
      (hasAllReferenceAxes ? channels : std::vector<std::string>{"zjet"});
    for (const std::string &channel : observableChannels) {
      frames << "\\begin{frame}{" << observableLabel(observable) << " -- "
             << channelLabel(channel) << "}\n"
             << "\\begin{columns}[T,onlytextwidth]\n";
      for (const std::string &sample : samples) {
        const std::string path = objectPath(sample,observable,channel);
        const std::vector<ComparisonPoint> oldPoints = readPoints(oldInput,path);
        const std::vector<ComparisonPoint> newPoints = readPoints(newInput,path);
        int matched = 0;
        const std::vector<ComparisonPoint> differences = matchedDifferences(
          oldPoints,newPoints,observable=="counts",matched);
        summary << sample << '\t' << observable << '\t' << channel << '\t'
                << oldPoints.size() << '\t' << newPoints.size() << '\t'
                << matched << '\t' << int(oldPoints.size())-matched << '\t'
                << int(newPoints.size())-matched << '\n';

        const std::string stem =
          sample+"_"+observable+"_"+channel;
        TCanvas canvas(stem.c_str(),stem.c_str(),720,720);
        TPad upper("upper","upper",0.,0.30,1.,1.);
        TPad lower("lower","lower",0.,0.,1.,0.30);
        upper.SetBottomMargin(0.02);
        upper.SetLeftMargin(0.14);
        upper.SetRightMargin(0.04);
        lower.SetTopMargin(0.03);
        lower.SetBottomMargin(0.31);
        lower.SetLeftMargin(0.14);
        lower.SetRightMargin(0.04);
        upper.Draw();
        lower.Draw();

        upper.cd();
        upper.SetLogx();
        if (observable=="counts") upper.SetLogy();
        const auto range = yRange(observable,sample);
        TH1D upperFrame("upperFrame","",100,12.,1500.);
        upperFrame.SetMinimum(range.first);
        upperFrame.SetMaximum(range.second);
        upperFrame.GetYaxis()->SetTitle(observableLabel(observable).c_str());
        upperFrame.GetYaxis()->SetTitleSize(0.055);
        upperFrame.GetYaxis()->SetLabelSize(0.045);
        upperFrame.GetXaxis()->SetLabelSize(0.);
        upperFrame.Draw("AXIS");
        TGraphErrors *oldGraph = makeGraph(oldPoints,kBlack,20);
        TGraphErrors *newGraph = makeGraph(newPoints,kRed+1,24);
        if (!oldPoints.empty()) oldGraph->Draw("P SAME");
        if (!newPoints.empty()) newGraph->Draw("P SAME");
        if (oldPoints.empty() && newPoints.empty())
          drawMissing(("Missing in both files: "+path).c_str());
        TLegend legend(0.55,0.73,0.92,0.90);
        legend.SetBorderSize(0);
        legend.SetFillStyle(0);
        legend.SetTextSize(0.040);
        legend.AddEntry(oldGraph,oldLabel,"pl");
        legend.AddEntry(newGraph,newLabel,"pl");
        legend.Draw();
        TLatex label;
        label.SetNDC();
        label.SetTextSize(0.045);
        label.DrawLatex(0.16,0.93,(sample+"; "+channelLabel(channel)).c_str());

        lower.cd();
        lower.SetLogx();
        const double deltaRange = differenceRange(observable);
        TH1D lowerFrame("lowerFrame","",100,12.,1500.);
        lowerFrame.SetMinimum(-deltaRange);
        lowerFrame.SetMaximum(deltaRange);
        lowerFrame.GetXaxis()->SetTitle("p_{T} (GeV)");
        lowerFrame.GetYaxis()->SetTitle(
          observable=="counts" ? "(new-old)/old" : "new-old");
        lowerFrame.GetXaxis()->SetTitleSize(0.12);
        lowerFrame.GetXaxis()->SetLabelSize(0.10);
        lowerFrame.GetYaxis()->SetTitleSize(0.10);
        lowerFrame.GetYaxis()->SetTitleOffset(0.64);
        lowerFrame.GetYaxis()->SetLabelSize(0.085);
        lowerFrame.GetYaxis()->SetNdivisions(505);
        lowerFrame.Draw("AXIS");
        TGraphErrors *differenceGraph = makeGraph(differences,kBlue+1,20);
        if (!differences.empty()) differenceGraph->Draw("P SAME");
        TLatex matchLabel;
        matchLabel.SetNDC();
        matchLabel.SetTextSize(0.085);
        matchLabel.DrawLatex(0.17,0.84,
          Form("exactly matched bins: %d",matched));

        const std::string plotPath =
          std::string(outputDirectory)+"/"+stem+".pdf";
        canvas.SaveAs(plotPath.c_str());
        delete oldGraph;
        delete newGraph;
        delete differenceGraph;

        frames << "\\begin{column}{0.327\\textwidth}\n"
               << "\\centering\\textbf{" << sample << "}\\par\n"
               << "\\includegraphics[width=\\linewidth]{"
               << outputDirectory << "/" << stem << ".pdf}\n"
               << "\\end{column}\n";
      }
      frames << "\\end{columns}\n"
             << "\\vspace{0.3ex}{\\tiny Top: old and new inputs. Bottom: "
             << (observable=="counts" ? "relative" : "absolute")
             << " bin-by-bin difference at exactly common bin centers.}\n"
             << "\\end{frame}\n";
    }
  }
  frames.close();
  summary.close();
  std::cout << "Wrote comparison plots, Beamer frames and summary to "
            << outputDirectory << std::endl;
}
