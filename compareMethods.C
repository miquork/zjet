// Quantitative validation of the synchronized legacy and all-pairs methods.
#include "TCanvas.h"
#include "TDirectory.h"
#include "TFile.h"
#include "TGraphErrors.h"
#include "TH1.h"
#include "TH1D.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLine.h"
#include "TProfile.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TSystem.h"
#include "tdrstyle_mod22.C"

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <numeric>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

struct Point {
  double x = 0.;
  double ex = 0.;
  double y = 0.;
  double ey = 0.;
};

struct MethodInput {
  TFile *file = nullptr;
  std::string label;
  int color = kBlack;
  int marker = 20;
};

const double kRn = 1.00;
const double kRu = 0.92;
const std::vector<std::string> kChannels = {"jetz","zjav","zjet"};

void applyMethodTDRStyle() {
  setTDRStyle();
  writeExtraText = true;
  extraText = "Work in progress";
  extraText2 = "";
  // CMS_lumi appends the collision energy for iPeriod=8.
  lumi_136TeV = "Run2024I";
}

bool sameCenter(double first, double second) {
  return std::fabs(first-second)<1.e-6*std::max({1.,std::fabs(first),
                                                 std::fabs(second)});
}

std::vector<Point> readPoints(TFile &file, const std::string &path,
                              bool required=true) {
  TObject *object = file.Get(path.c_str());
  if (!object) {
    if (required)
      throw std::runtime_error("Missing "+path+" in "+file.GetName());
    return {};
  }
  std::vector<Point> points;
  if (TGraphErrors *graph = dynamic_cast<TGraphErrors*>(object)) {
    for (int index = 0; index != graph->GetN(); ++index) {
      double x = 0.;
      double y = 0.;
      graph->GetPoint(index,x,y);
      if (!std::isfinite(x) || !std::isfinite(y)) continue;
      points.push_back({x,graph->GetErrorX(index),y,graph->GetErrorY(index)});
    }
    return points;
  }
  if (TH1 *histogram = dynamic_cast<TH1*>(object)) {
    for (int bin = 1; bin <= histogram->GetNbinsX(); ++bin) {
      const double value = histogram->GetBinContent(bin);
      if (value==0. || !std::isfinite(value)) continue;
      points.push_back({histogram->GetBinCenter(bin),
                        0.5*histogram->GetBinWidth(bin),value,
                        histogram->GetBinError(bin)});
    }
    return points;
  }
  throw std::runtime_error(path+" has unsupported class "+object->ClassName());
}

std::string originalPath(const std::string &sample,
                         const std::string &observable,
                         const std::string &channel) {
  return sample+"/eta00-13/orig/"+observable+"_"+channel+"_a100";
}

std::string hdmPath(const std::string &sample, const std::string &channel) {
  return sample+"/eta00-13/hdm_mpfchs1_"+channel;
}

std::vector<Point> combineSoftRecoil(TFile &file, const std::string &sample,
                                     const std::string &channel) {
  const std::vector<Point> neutral = readPoints(
    file,originalPath(sample,"mpfn",channel));
  const std::vector<Point> unclustered = readPoints(
    file,originalPath(sample,"mpfu",channel));
  std::vector<Point> result;
  for (const Point &n : neutral) {
    const auto u = std::find_if(
      unclustered.begin(),unclustered.end(),[&](const Point &point) {
        return sameCenter(n.x,point.x);
      });
    if (u==unclustered.end()) continue;
    result.push_back({n.x,n.ex,n.y/kRn+u->y/kRu,
                      std::hypot(n.ey/kRn,u->ey/kRu)});
  }
  return result;
}

std::vector<Point> differences(const std::vector<Point> &first,
                               const std::vector<Point> &second) {
  std::vector<Point> result;
  for (const Point &left : first) {
    const auto right = std::find_if(
      second.begin(),second.end(),[&](const Point &point) {
        return sameCenter(left.x,point.x);
      });
    if (right==second.end()) continue;
    result.push_back({left.x,left.ex,left.y-right->y,
                      std::hypot(left.ey,right->ey)});
  }
  return result;
}

std::vector<Point> ratios(const std::vector<Point> &numerator,
                          const std::vector<Point> &denominator) {
  std::vector<Point> result;
  for (const Point &top : numerator) {
    const auto bottom = std::find_if(
      denominator.begin(),denominator.end(),[&](const Point &point) {
        return sameCenter(top.x,point.x);
      });
    if (bottom==denominator.end() || bottom->y==0.) continue;
    const double value = top.y/bottom->y;
    const double error = std::hypot(top.ey/bottom->y,
                                    top.y*bottom->ey/(bottom->y*bottom->y));
    result.push_back({top.x,top.ex,value,error});
  }
  return result;
}

double interpolate(const std::vector<Point> &points, double x);

std::vector<Point> series(TFile &file, const std::string &sample,
                          const std::string &observable,
                          const std::string &channel) {
  if (observable=="hdm") return readPoints(file,hdmPath(sample,channel));
  if (observable=="mpfnu") {
    if (sample=="data" || sample=="mc")
      return combineSoftRecoil(file,sample,channel);
    return differences(combineSoftRecoil(file,"data",channel),
                       combineSoftRecoil(file,"mc",channel));
  }
  return readPoints(file,originalPath(sample,observable,channel));
}

std::vector<Point> previousJEC(TFile &file) {
  return readPoints(file,"ratio/eta00-13/herr_l2l3res");
}

double globalFitAxisScale(const std::string &channel) {
  // Current globalFitSettings.h keeps scaleJZperEra disabled. The average-pT
  // branch currently enables its Run2024I factor internally.
  if (channel=="jetz") return 1.0000;
  if (channel=="zjav") return 1.0025;
  return 1.0000;
}

std::vector<Point> globalFitHDM(TFile &file, const std::string &channel) {
  std::vector<Point> result = series(file,"ratio","hdm",channel);
  const std::vector<Point> correction = previousJEC(file);
  const double axisScale = globalFitAxisScale(channel);
  for (Point &point : result) {
    const double factor = interpolate(correction,point.x);
    if (!std::isfinite(factor)) {
      point.y = std::numeric_limits<double>::quiet_NaN();
      continue;
    }
    point.y *= factor*axisScale;
    // This reproduces scaleGraph in globalFit.C: the central JEC is treated as
    // an input scale, so its uncertainty is not added to the graph error.
    point.ey *= std::fabs(factor*axisScale);
  }
  result.erase(std::remove_if(result.begin(),result.end(),[](const Point &point) {
    return !std::isfinite(point.y);
  }),result.end());
  return result;
}

std::unique_ptr<TGraphErrors> graph(const std::vector<Point> &points,
                                    int color, int marker,
                                    bool line=false) {
  std::unique_ptr<TGraphErrors> result(new TGraphErrors(points.size()));
  for (size_t index = 0; index != points.size(); ++index) {
    result->SetPoint(index,points[index].x,points[index].y);
    result->SetPointError(index,points[index].ex,points[index].ey);
  }
  result->SetLineColor(color);
  result->SetMarkerColor(color);
  result->SetMarkerStyle(marker);
  result->SetMarkerSize(0.8);
  result->SetLineWidth(line ? 2 : 1);
  return result;
}

void configureLogAxis(TH1 *frame) {
  frame->GetXaxis()->SetMoreLogLabels();
  frame->GetXaxis()->SetNoExponent();
}

void drawMethodLogLabels(TH1 *frame, double offset=0.018,
                         double textSize=-1.) {
  if (!frame || !gPad || !gPad->GetLogx()) return;
  gPad->Update();
  const double xmin = frame->GetXaxis()->GetXmin();
  const double xmax = frame->GetXaxis()->GetXmax();
  frame->GetXaxis()->SetLabelOffset(999.);
  gPad->Modified();
  gPad->Update();
  TLatex label;
  label.SetNDC();
  label.SetTextAlign(23);
  label.SetTextFont(frame->GetXaxis()->GetLabelFont());
  label.SetTextSize(textSize>0. ? textSize : frame->GetXaxis()->GetLabelSize());
  for (double value : {30.,100.,300.,1000.,3000.}) {
    if (value<xmin || value>xmax) continue;
    const double fraction = (std::log10(value)-std::log10(xmin))/
                            (std::log10(xmax)-std::log10(xmin));
    const double x = gPad->GetLeftMargin()+fraction*
      (1.-gPad->GetLeftMargin()-gPad->GetRightMargin());
    label.DrawLatex(x,gPad->GetBottomMargin()-offset,Form("%g",value));
  }
}

void zeroLine(double minimum, double maximum) {
  TLine line(minimum,0.,maximum,0.);
  line.SetLineStyle(kDashed);
  line.SetLineColor(kGray+2);
  line.DrawClone();
}

double symmetricRange(const std::vector<std::vector<Point> > &collections,
                      double minimum) {
  double maximum = minimum;
  for (const auto &points : collections)
    for (const Point &point : points)
      maximum = std::max(maximum,std::fabs(point.y)+point.ey);
  return 1.15*maximum;
}

double focusedSymmetricRange(
  const std::vector<std::vector<Point> > &collections, double minimum,
  double focusMaximum=500.) {
  double maximum = minimum;
  for (const auto &points : collections)
    for (const Point &point : points)
      if (point.x<=focusMaximum) maximum = std::max(maximum,std::fabs(point.y));
  return 1.20*maximum;
}

std::vector<Point> scaleDifference(const std::vector<Point> &candidate,
                                   const std::vector<Point> &reference,
                                   double scale=1000.) {
  std::vector<Point> result = differences(candidate,reference);
  for (Point &point : result) {
    point.y *= scale;
    point.ey *= scale;
  }
  return result;
}

double interpolate(const std::vector<Point> &points, double x) {
  if (points.empty()) return std::numeric_limits<double>::quiet_NaN();
  for (const Point &point : points)
    if (sameCenter(point.x,x)) return point.y;
  for (size_t index = 1; index != points.size(); ++index) {
    if (points[index-1].x<x && x<points[index].x) {
      const double fraction =
        (std::log(x)-std::log(points[index-1].x))/
        (std::log(points[index].x)-std::log(points[index-1].x));
      return points[index-1].y+fraction*(points[index].y-points[index-1].y);
    }
  }
  return std::numeric_limits<double>::quiet_NaN();
}

void writeFrame(std::ofstream &frames, const std::string &title,
                const std::vector<std::pair<std::string,std::string> > &plots,
                const std::string &note="") {
  frames << "\\begin{frame}{" << title << "}\n"
         << "\\begin{columns}[T,onlytextwidth]\n";
  const double width = 0.98/plots.size();
  for (const auto &plot : plots) {
    frames << "\\begin{column}{" << std::fixed << std::setprecision(3)
           << width << "\\textwidth}\\centering\n"
           << "\\textbf{" << plot.second << "}\\par\n"
           << "\\includegraphics[width=\\linewidth,height=0.76\\textheight,"
           << "keepaspectratio]{" << plot.first << "}\n"
           << "\\end{column}\n";
  }
  frames << "\\end{columns}\n";
  if (!note.empty()) frames << "\\vspace{0.3ex}{\\tiny " << note << "}\n";
  frames << "\\end{frame}\n";
}

std::vector<Point> transformedResult(TFile &file,
                                     const std::string &observable,
                                     const std::string &channel) {
  std::vector<Point> points = series(file,"ratio",observable,channel);
  for (Point &point : points) {
    point.y = 1000.*(point.y-((observable=="mpfchs1" || observable=="hdm")
                              ? 1. : 0.));
    point.ey *= 1000.;
  }
  return points;
}

double solveHDM(double r0, double rn, double ru) {
  return (r0-rn-ru)/(1.-rn/kRn-ru/kRu);
}

double hdmRatio(const std::array<double,6> &values) {
  const double data = solveHDM(values[0],values[1],values[2]);
  const double mc = solveHDM(values[3],values[4],values[5]);
  return data/mc;
}

bool valuesAt(TFile &file, const std::string &channel, double x,
              std::array<double,6> &values) {
  const char *samples[] = {"data","mc"};
  const char *observables[] = {"mpfchs1","mpfn","mpfu"};
  int position = 0;
  for (const char *sample : samples) {
    for (const char *observable : observables) {
      const std::vector<Point> points = series(file,sample,observable,channel);
      const auto point = std::find_if(
        points.begin(),points.end(),[&](const Point &candidate) {
          return sameCenter(candidate.x,x);
        });
      if (point==points.end()) return false;
      values[position++] = point->y;
    }
  }
  return true;
}

std::array<double,6> shapley(const std::array<double,6> &legacy,
                             const std::array<double,6> &modern) {
  constexpr int count = 6;
  const double factorial[] = {1.,1.,2.,6.,24.,120.,720.};
  std::array<double,count> contributions{};
  for (int variable = 0; variable != count; ++variable) {
    for (int mask = 0; mask != (1<<count); ++mask) {
      if (mask&(1<<variable)) continue;
      const int size = __builtin_popcount(static_cast<unsigned>(mask));
      const double weight = factorial[size]*factorial[count-size-1]/
                            factorial[count];
      std::array<double,count> before = legacy;
      for (int index = 0; index != count; ++index)
        if (mask&(1<<index)) before[index] = modern[index];
      std::array<double,count> after = before;
      after[variable] = modern[variable];
      contributions[variable] += weight*(hdmRatio(after)-hdmRatio(before));
    }
  }
  return contributions;
}

void drawHDMDifference(TFile &reference, TFile &candidate,
                       const std::string &referenceLabel,
                       const std::string &candidateLabel,
                       const std::string &stem,
                       const std::string &outputDirectory,
                       std::ofstream &metrics) {
  const std::vector<Point> difference = scaleDifference(
    series(candidate,"ratio","hdm","zjav"),
    series(reference,"ratio","hdm","zjav"));
  const double range = symmetricRange({difference},2.);
  TH1D frame(("h_"+stem).c_str(),"",100,12.,1500.);
  frame.SetDirectory(nullptr);
  frame.SetMinimum(-range);
  frame.SetMaximum(range);
  frame.GetXaxis()->SetTitle("p_{T,ave} (GeV)");
  frame.GetYaxis()->SetTitle("HDM candidate - reference (per mille)");
  configureLogAxis(&frame);
  std::unique_ptr<TCanvas> canvas(
    tdrCanvas(("c_"+stem).c_str(),&frame,8,11,kSquare));
  canvas->SetLogx();
  zeroLine(12.,1500.);
  auto result = graph(difference,kRed+1,kFullCircle);
  result->Draw("P SAME");
  TLegend legend(0.50,0.78,0.91,0.87);
  legend.SetBorderSize(0);
  legend.AddEntry(result.get(),(candidateLabel+" - "+referenceLabel).c_str(),"p");
  legend.Draw();
  drawMethodLogLabels(&frame);
  canvas->SaveAs((outputDirectory+"/"+stem+".pdf").c_str());

  metrics << "# HDM " << candidateLabel << " minus " << referenceLabel
          << ", average-pT binning\n"
          << "pt\tdifference_per_mille\terror_per_mille\n";
  for (const Point &point : difference)
    metrics << point.x << '\t' << point.y << '\t' << point.ey << '\n';
}

void drawGlobalFitDifference(TFile &reference, TFile &candidate,
                             const std::string &referenceLabel,
                             const std::string &candidateLabel,
                             const std::string &stem,
                             const std::string &outputDirectory,
                             std::ofstream &metrics) {
  const std::string channel = "zjet";
  const std::vector<Point> referenceRaw =
    series(reference,"ratio","hdm",channel);
  const std::vector<Point> candidateRaw =
    series(candidate,"ratio","hdm",channel);
  const std::vector<Point> referenceJEC = previousJEC(reference);
  const std::vector<Point> candidateJEC = previousJEC(candidate);
  const std::vector<Point> total = scaleDifference(
    globalFitHDM(candidate,channel),globalFitHDM(reference,channel));
  std::vector<Point> rawContribution;
  std::vector<Point> jecContribution;
  metrics << "\n# globalFit.C input difference: " << candidateLabel
          << " minus " << referenceLabel << ", Z-pT binning\n"
          << "pt\ttotal_per_mille\traw_hdm_contribution"
          << "\tprevious_jec_contribution\tclosure_per_mille\n";
  for (const Point &base : referenceRaw) {
    const double rawNew = interpolate(candidateRaw,base.x);
    const double jecOld = interpolate(referenceJEC,base.x);
    const double jecNew = interpolate(candidateJEC,base.x);
    const double exact = interpolate(total,base.x);
    if (!std::isfinite(rawNew) || !std::isfinite(jecOld) ||
        !std::isfinite(jecNew) || !std::isfinite(exact)) continue;
    // Exact two-variable Shapley decomposition of raw HDM times previous JEC.
    const double raw = 500.*(rawNew-base.y)*(jecOld+jecNew);
    const double jec = 500.*(jecNew-jecOld)*(base.y+rawNew);
    rawContribution.push_back({base.x,base.ex,raw,0.});
    jecContribution.push_back({base.x,base.ex,jec,0.});
    metrics << base.x << '\t' << exact << '\t' << raw << '\t' << jec
            << '\t' << exact-raw-jec << '\n';
  }
  const double range = symmetricRange({total,rawContribution,jecContribution},5.);
  TH1D frame(("h_"+stem).c_str(),"",100,12.,1500.);
  frame.SetDirectory(nullptr);
  frame.SetMinimum(-range);
  frame.SetMaximum(range);
  frame.GetXaxis()->SetTitle("p_{T,Z} (GeV)");
  frame.GetYaxis()->SetTitle("Candidate - reference (per mille)");
  configureLogAxis(&frame);
  std::unique_ptr<TCanvas> canvas(
    tdrCanvas(("c_"+stem).c_str(),&frame,8,11,kSquare));
  canvas->SetLogx();
  zeroLine(12.,1500.);
  auto direct = graph(total,kBlack,kFullSquare,true);
  auto raw = graph(rawContribution,kRed+1,kFullCircle,true);
  auto jec = graph(jecContribution,kBlue+1,kFullTriangleUp,true);
  direct->SetLineWidth(3);
  direct->Draw("LP SAME");
  raw->Draw("LP SAME");
  jec->Draw("LP SAME");
  TLegend legend(0.43,0.66,0.91,0.88);
  legend.SetBorderSize(0);
  legend.AddEntry(direct.get(),
                  (candidateLabel+" - "+referenceLabel).c_str(),"lp");
  legend.AddEntry(raw.get(),"raw HDM contribution","lp");
  legend.AddEntry(jec.get(),"previous-JEC contribution","lp");
  legend.Draw();
  drawMethodLogLabels(&frame);
  canvas->SaveAs((outputDirectory+"/"+stem+".pdf").c_str());
}

void drawHDMDecomposition(TFile &reference, TFile &candidate,
                          const std::string &referenceLabel,
                          const std::string &candidateLabel,
                          const std::string &stem,
                          const std::string &outputDirectory,
                          std::ofstream &metrics) {
  const std::vector<Point> legacyHDM = series(reference,"ratio","hdm","zjav");
  const std::vector<Point> modernHDM = series(candidate,"ratio","hdm","zjav");
  std::array<std::vector<Point>,6> componentPoints;
  std::vector<Point> totalPoints;
  std::vector<Point> closurePoints;
  metrics << "\n# Exact Shapley decomposition: " << candidateLabel
          << " minus " << referenceLabel << "\n"
          << "pt\ttotal_per_mille\tdata_mpf\tdata_mpfn\tdata_mpfu"
          << "\tmc_mpf\tmc_mpfn\tmc_mpfu\tclosure_per_mille\n";
  for (const Point &referencePoint : legacyHDM) {
    const auto candidatePoint = std::find_if(
      modernHDM.begin(),modernHDM.end(),[&](const Point &point) {
        return sameCenter(referencePoint.x,point.x);
      });
    if (candidatePoint==modernHDM.end()) continue;
    std::array<double,6> oldValues{};
    std::array<double,6> newValues{};
    if (!valuesAt(reference,"zjav",referencePoint.x,oldValues) ||
        !valuesAt(candidate,"zjav",referencePoint.x,newValues)) continue;
    const std::array<double,6> contribution = shapley(oldValues,newValues);
    const double total = candidatePoint->y-referencePoint.y;
    const double sum = std::accumulate(contribution.begin(),contribution.end(),0.);
    totalPoints.push_back({referencePoint.x,referencePoint.ex,1000.*total,0.});
    closurePoints.push_back({referencePoint.x,referencePoint.ex,
                             1000.*(total-sum),0.});
    metrics << referencePoint.x << '\t' << 1000.*total;
    for (int index = 0; index != 6; ++index) {
      componentPoints[index].push_back(
        {referencePoint.x,referencePoint.ex,1000.*contribution[index],0.});
      metrics << '\t' << 1000.*contribution[index];
    }
    metrics << '\t' << 1000.*(total-sum) << '\n';
  }
  std::vector<std::vector<Point> > ranges(componentPoints.begin(),
                                           componentPoints.end());
  ranges.push_back(totalPoints);
  const double range = symmetricRange(ranges,1.);
  TH1D frame(("h_"+stem).c_str(),"",100,12.,1500.);
  frame.SetDirectory(nullptr);
  frame.SetMinimum(-range);
  frame.SetMaximum(range);
  frame.GetXaxis()->SetTitle("p_{T,ave} (GeV)");
  frame.GetYaxis()->SetTitle("HDM contribution (per mille)");
  configureLogAxis(&frame);
  std::unique_ptr<TCanvas> canvas(
    tdrCanvas(("c_"+stem).c_str(),&frame,8,11,kSquare));
  canvas->SetLogx();
  zeroLine(12.,1500.);
  const int colors[] = {kBlue+1,kAzure+7,kCyan+2,kRed+1,kOrange+7,kMagenta+1};
  const int styles[] = {1,2,3,1,2,3};
  const char *labels[] = {"data MPF","data MPFn","data MPFu",
                          "MC MPF","MC MPFn","MC MPFu"};
  std::vector<std::unique_ptr<TGraphErrors> > graphs;
  TLegend legend(0.53,0.60,0.92,0.88);
  legend.SetBorderSize(0);
  legend.SetNColumns(2);
  for (int index = 0; index != 6; ++index) {
    graphs.push_back(graph(componentPoints[index],colors[index],kFullCircle,true));
    graphs.back()->SetLineStyle(styles[index]);
    graphs.back()->SetMarkerSize(0.45);
    graphs.back()->Draw("LP SAME");
    legend.AddEntry(graphs.back().get(),labels[index],"lp");
  }
  auto total = graph(totalPoints,kBlack,kFullSquare,true);
  total->SetLineWidth(3);
  total->Draw("LP SAME");
  legend.AddEntry(total.get(),"exact total","lp");
  legend.Draw();
  drawMethodLogLabels(&frame);
  canvas->SaveAs((outputDirectory+"/"+stem+".pdf").c_str());
}

void drawWReference(TFile &baseline, TFile &legacy, TFile &modern,
                    const std::string &outputDirectory) {
  // W->qq' is an independent mass-peak response, not a Z+jet HDM channel.
  const std::vector<Point> w = readPoints(
    baseline,originalPath("ratio","mpfchs1","wqq"));
  const std::vector<Point> zBaseline = series(baseline,"ratio","hdm","zjav");
  const std::vector<Point> zLegacy = series(legacy,"ratio","hdm","zjav");
  const std::vector<Point> zModern = series(modern,"ratio","hdm","zjav");
  std::vector<std::vector<Point> > residuals;
  for (const auto &input : {zBaseline,zLegacy,zModern}) {
    std::vector<Point> residual;
    for (const Point &point : input) {
      const double reference = interpolate(w,point.x);
      if (!std::isfinite(reference)) continue;
      residual.push_back({point.x,point.ex,1000.*(point.y-reference),
                          1000.*point.ey});
    }
    residuals.push_back(residual);
  }
  const double lowerRange = symmetricRange(residuals,5.);
  TH1D upperFrame("h_w_upper","",100,12.,1500.);
  TH1D lowerFrame("h_w_lower","",100,12.,1500.);
  upperFrame.SetDirectory(nullptr);
  lowerFrame.SetDirectory(nullptr);
  upperFrame.SetMinimum(0.93); upperFrame.SetMaximum(1.07);
  upperFrame.GetYaxis()->SetTitle("HDM data / MC");
  upperFrame.GetXaxis()->SetLabelSize(0.);
  configureLogAxis(&upperFrame);
  lowerFrame.SetMinimum(-lowerRange); lowerFrame.SetMaximum(lowerRange);
  lowerFrame.GetXaxis()->SetTitle("reference p_{T} (GeV)");
  lowerFrame.GetYaxis()->SetTitle("Z - W (10^{-3})");
  configureLogAxis(&lowerFrame);
  std::unique_ptr<TCanvas> canvas(
    tdrDiCanvas("c_w_reference",&upperFrame,&lowerFrame,8,11));
  canvas->cd(1); gPad->SetLogx();
  TLine unity(12.,1.,1500.,1.); unity.SetLineStyle(kDashed); unity.DrawClone();
  auto gw = graph(w,kBlack,kFullSquare,true);
  auto gb = graph(zBaseline,kGray+2,kOpenCircle,true);
  auto gl = graph(zLegacy,kBlue+1,kOpenSquare,true);
  auto gn = graph(zModern,kRed+1,kFullCircle,true);
  for (auto *g : {gw.get(),gb.get(),gl.get(),gn.get()}) g->Draw("LP SAME");
  TLegend legend(0.50,0.66,0.92,0.89);
  legend.SetBorderSize(0);
  legend.AddEntry(gw.get(),"W#rightarrowqq' reference","lp");
  legend.AddEntry(gb.get(),"Z baseline","lp");
  legend.AddEntry(gl.get(),"Z synchronized legacy","lp");
  legend.AddEntry(gn.get(),"Z all-pairs","lp");
  legend.Draw();
  canvas->cd(2); gPad->SetLogx();
  zeroLine(12.,1500.);
  auto grb = graph(residuals[0],kGray+2,kOpenCircle,true);
  auto grl = graph(residuals[1],kBlue+1,kOpenSquare,true);
  auto grn = graph(residuals[2],kRed+1,kFullCircle,true);
  for (auto *g : {grb.get(),grl.get(),grn.get()}) g->Draw("LP SAME");
  drawMethodLogLabels(&lowerFrame,0.045,0.10);
  canvas->SaveAs((outputDirectory+"/hdm_wqq_reference.pdf").c_str());
}

void drawMCTruthClosure(const std::vector<MethodInput> &inputs,
                        const std::string &outputDirectory,
                        std::ofstream &metrics) {
  std::array<std::vector<Point>,3> hdm;
  std::array<std::vector<Point>,3> truth;
  std::array<std::vector<Point>,3> closure;
  metrics << "\n# MC HDM / stored generator balance, Z-pT binning\n"
          << "method\tpt\thdm\tgen_balance\tclosure_per_mille\tstatus\n";
  std::vector<std::vector<Point> > closureRanges;
  for (int method=0; method!=3; ++method) {
    hdm[method] = series(*inputs[method].file,"mc","hdm","zjet");
    truth[method] = series(*inputs[method].file,"mc","gjet","zjet");
    for (const Point &response : hdm[method]) {
      const double gen = interpolate(truth[method],response.x);
      const bool physical = std::isfinite(gen) && gen>0.5 && gen<1.5;
      metrics << inputs[method].label << '\t' << response.x << '\t'
              << response.y << '\t' << gen << '\t';
      if (physical) {
        const double value = 1000.*(response.y/gen-1.);
        const double error = 1000.*response.ey/std::fabs(gen);
        closure[method].push_back({response.x,response.ex,value,error});
        metrics << value << "\tused\n";
      } else {
        metrics << "nan\trejected_unphysical_gen_balance\n";
      }
    }
    closureRanges.push_back(closure[method]);
  }
  const double closureRange = 25.*std::ceil(
    focusedSymmetricRange(closureRanges,25.)/25.);
  TH1D upperFrame("h_mc_truth_closure_upper","",100,12.,1500.);
  TH1D lowerFrame("h_mc_truth_closure_lower","",100,12.,1500.);
  upperFrame.SetDirectory(nullptr);
  lowerFrame.SetDirectory(nullptr);
  upperFrame.SetMinimum(-0.10);
  upperFrame.SetMaximum(1.55);
  upperFrame.GetYaxis()->SetTitle("MC response");
  upperFrame.GetXaxis()->SetLabelSize(0.);
  configureLogAxis(&upperFrame);
  lowerFrame.SetMinimum(-closureRange);
  lowerFrame.SetMaximum(closureRange);
  lowerFrame.GetXaxis()->SetTitle("p_{T,Z} (GeV)");
  lowerFrame.GetYaxis()->SetTitle("Closure (10^{-3})");
  configureLogAxis(&lowerFrame);
  std::unique_ptr<TCanvas> canvas(
    tdrDiCanvas("c_mc_truth_closure",&upperFrame,&lowerFrame,8,11));
  std::array<std::unique_ptr<TGraphErrors>,3> hdmGraphs;
  std::array<std::unique_ptr<TGraphErrors>,3> truthGraphs;
  std::array<std::unique_ptr<TGraphErrors>,3> closureGraphs;
  canvas->cd(1);
  gPad->SetLogx();
  TLegend legend(0.43,0.55,0.91,0.86);
  legend.SetBorderSize(0);
  legend.SetNColumns(2);
  for (int method=0; method!=3; ++method) {
    hdmGraphs[method] = graph(hdm[method],inputs[method].color,
                              inputs[method].marker,true);
    truthGraphs[method] = graph(truth[method],inputs[method].color,
                                inputs[method].marker,true);
    truthGraphs[method]->SetLineStyle(kDashed);
    truthGraphs[method]->SetMarkerStyle(inputs[method].marker+4);
    hdmGraphs[method]->Draw("LP SAME");
    truthGraphs[method]->Draw("LP SAME");
    legend.AddEntry(hdmGraphs[method].get(),
                    (inputs[method].label+" HDM").c_str(),"lp");
    legend.AddEntry(truthGraphs[method].get(),
                    (inputs[method].label+" gen").c_str(),"lp");
  }
  legend.Draw();
  canvas->cd(2);
  gPad->SetLogx();
  zeroLine(12.,1500.);
  for (int method=0; method!=3; ++method) {
    closureGraphs[method] = graph(closure[method],inputs[method].color,
                                  inputs[method].marker,true);
    closureGraphs[method]->Draw("LP SAME");
  }
  TLatex warning;
  warning.SetNDC();
  warning.SetTextSize(0.070);
  warning.DrawLatex(0.26,0.82,"gen outside 0.5--1.5: no division");
  drawMethodLogLabels(&lowerFrame,0.045,0.10);
  canvas->SaveAs((outputDirectory+"/mc_hdm_truth_closure_zpt.pdf").c_str());
}

void drawComponentAxes(const std::vector<MethodInput> &inputs,
                       const std::string &observable,
                       const std::string &outputDirectory) {
  std::array<std::array<std::vector<Point>,3>,3> points;
  std::vector<std::vector<Point> > ranges;
  for (int method = 0; method != 3; ++method)
    for (int channel = 0; channel != 3; ++channel) {
      points[method][channel] = transformedResult(
        *inputs[method].file,observable,kChannels[channel]);
      ranges.push_back(points[method][channel]);
    }
  const double minimum = (observable=="hdm" ? 5. :
                          observable=="mpfchs1" ? 10. : 5.);
  // Preserve the full pT reach but choose the vertical zoom from the region
  // that has useful precision.  Sparse TeV-bin errors must not hide permille
  // structure below 500 GeV.
  const double range = focusedSymmetricRange(ranges,minimum);
  TCanvas canvas(("c_axes_"+observable).c_str(),"",1800,600);
  canvas.Divide(3,1,0.002,0.002);
  // Keep pad primitives alive until SaveAs.  TGraph::Draw and TLegend::Draw do
  // not clone their objects, so pad-local smart pointers would leave an empty
  // canvas once the channel loop advances.
  std::array<std::vector<std::unique_ptr<TGraphErrors> >,3> drawnGraphs;
  std::array<std::unique_ptr<TLegend>,3> drawnLegends;
  const char *channelLabels[] = {"jet p_{T} binning","average p_{T} binning",
                                 "Z p_{T} binning"};
  for (int channel = 0; channel != 3; ++channel) {
    canvas.cd(channel+1);
    gPad->SetLogx();
    gPad->SetLeftMargin(0.14);
    gPad->SetRightMargin(0.03);
    gPad->SetBottomMargin(0.13);
    TH1D *frame = new TH1D(Form("h_axes_%s_%d",observable.c_str(),channel),
                           "",100,12.,1500.);
    frame->SetMinimum(-range); frame->SetMaximum(range);
    frame->GetXaxis()->SetTitle("reference p_{T} (GeV)");
    frame->GetYaxis()->SetTitle(
      (observable=="hdm" || observable=="mpfchs1")
        ? "data/MC - 1 (per mille)" : "data - MC (per mille)");
    configureLogAxis(frame);
    frame->Draw("AXIS");
    zeroLine(12.,1500.);
    for (int method = 0; method != 3; ++method) {
      drawnGraphs[channel].push_back(
        graph(points[method][channel],inputs[method].color,
              inputs[method].marker,true));
      if (method==0) drawnGraphs[channel].back()->SetLineStyle(kDashed);
      drawnGraphs[channel].back()->Draw("LP SAME");
    }
    drawnLegends[channel].reset(new TLegend(0.48,0.70,0.91,0.88));
    drawnLegends[channel]->SetBorderSize(0);
    for (int method = 0; method != 3; ++method)
      drawnLegends[channel]->AddEntry(drawnGraphs[channel][method].get(),
                                      inputs[method].label.c_str(),"lp");
    drawnLegends[channel]->Draw();
    TLatex title;
    title.SetNDC(); title.SetTextSize(0.042);
    title.DrawLatex(0.16,0.93,channelLabels[channel]);
    drawMethodLogLabels(frame);
  }
  canvas.SaveAs((outputDirectory+"/axes_"+observable+".pdf").c_str());
}

void drawRhoStability(const std::vector<MethodInput> &inputs,
                      const std::string &outputDirectory) {
  std::vector<std::vector<Point> > rhoRatios;
  for (const MethodInput &input : inputs)
    rhoRatios.push_back(ratios(series(*input.file,"data","rho","zjet"),
                               series(*input.file,"mc","rho","zjet")));
  std::vector<std::vector<Point> > lower = {
    scaleDifference(rhoRatios[1],rhoRatios[0]),
    scaleDifference(rhoRatios[2],rhoRatios[0]),
  };
  const double lowerRange = symmetricRange(lower,5.);
  TH1D upperFrame("h_rho_upper","",100,12.,1500.);
  TH1D lowerFrame("h_rho_lower","",100,12.,1500.);
  upperFrame.SetDirectory(nullptr);
  lowerFrame.SetDirectory(nullptr);
  upperFrame.SetMinimum(0.90); upperFrame.SetMaximum(1.10);
  upperFrame.GetYaxis()->SetTitle("#rho(data) / #rho(MC)");
  upperFrame.GetXaxis()->SetLabelSize(0.);
  configureLogAxis(&upperFrame);
  lowerFrame.SetMinimum(-lowerRange); lowerFrame.SetMaximum(lowerRange);
  lowerFrame.GetXaxis()->SetTitle("p_{T,Z} (GeV)");
  lowerFrame.GetYaxis()->SetTitle("#Delta (10^{-3})");
  configureLogAxis(&lowerFrame);
  std::unique_ptr<TCanvas> canvas(
    tdrDiCanvas("c_rho_stability",&upperFrame,&lowerFrame,8,11));
  canvas->cd(1); gPad->SetLogx();
  TLine unity(12.,1.,1500.,1.); unity.SetLineStyle(kDashed); unity.DrawClone();
  std::vector<std::unique_ptr<TGraphErrors> > graphs;
  for (int method = 0; method != 3; ++method) {
    graphs.push_back(graph(rhoRatios[method],inputs[method].color,
                           inputs[method].marker,true));
    graphs.back()->Draw("LP SAME");
  }
  TLegend legend(0.51,0.70,0.91,0.88);
  legend.SetBorderSize(0);
  for (int method = 0; method != 3; ++method)
    legend.AddEntry(graphs[method].get(),inputs[method].label.c_str(),"lp");
  legend.Draw();
  canvas->cd(2); gPad->SetLogx(); zeroLine(12.,1500.);
  auto legacy = graph(lower[0],inputs[1].color,inputs[1].marker,true);
  auto modern = graph(lower[1],inputs[2].color,inputs[2].marker,true);
  legacy->Draw("LP SAME"); modern->Draw("LP SAME");
  drawMethodLogLabels(&lowerFrame,0.045,0.10);
  canvas->SaveAs((outputDirectory+"/rho_data_mc_stability.pdf").c_str());
}

TProfile *profile(TFile &file, const std::string &name, bool required=false) {
  TProfile *result = dynamic_cast<TProfile*>(file.Get(name.c_str()));
  if (required && !result)
    throw std::runtime_error("Missing "+name+" in "+file.GetName());
  return result;
}

std::vector<Point> profilePoints(TProfile *input) {
  std::vector<Point> points;
  if (!input) return points;
  for (int bin = 1; bin <= input->GetNbinsX(); ++bin) {
    if (input->GetBinEntries(bin)==0.) continue;
    const double value = input->GetBinContent(bin);
    const double error = input->GetBinError(bin);
    if (!std::isfinite(value) || !std::isfinite(error)) continue;
    points.push_back({input->GetBinCenter(bin),
                      0.5*input->GetBinWidth(bin),value,error});
  }
  return points;
}

std::vector<Point> truthResponse(TFile &file, const std::string &component,
                                 const std::string &axis,
                                 const std::string &category,
                                 const std::string &region,
                                 const std::string &suffix) {
  const auto name = [&](const std::string &observable) {
    return "control/p_"+observable+"_vs_"+axis+"_"+category+"_"+
      region+"_"+suffix;
  };
  if (component!="hdm") return profilePoints(profile(file,name(component)));

  const std::vector<Point> total = profilePoints(profile(file,name("mpf")));
  const std::vector<Point> neutral = profilePoints(profile(file,name("mpfn")));
  const std::vector<Point> unclustered = profilePoints(profile(file,name("mpfu")));
  if (total.empty() || neutral.empty() || unclustered.empty()) return {};
  std::vector<Point> result;
  for (const Point &r0 : total) {
    const auto rn = std::find_if(neutral.begin(),neutral.end(),[&](const Point &p) {
      return sameCenter(p.x,r0.x);
    });
    const auto ru = std::find_if(unclustered.begin(),unclustered.end(),[&](const Point &p) {
      return sameCenter(p.x,r0.x);
    });
    if (rn==neutral.end() || ru==unclustered.end()) continue;
    const double denominator = 1.-rn->y/kRn-ru->y/kRu;
    const double numerator = r0.y-rn->y-ru->y;
    if (std::fabs(denominator)<1.e-8) continue;
    const double dR0 = 1./denominator;
    const double dRn = (-denominator+numerator/kRn)/
                       (denominator*denominator);
    const double dRu = (-denominator+numerator/kRu)/
                       (denominator*denominator);
    const double error = std::sqrt(
      std::pow(dR0*r0.ey,2)+std::pow(dRn*rn->ey,2)+std::pow(dRu*ru->ey,2));
    result.push_back({r0.x,r0.ex,solveHDM(r0.y,rn->y,ru->y),error});
  }
  return result;
}

bool drawTruthComponent(TFile &mc, const std::string &component,
                        const std::string &axis,
                        const std::string &outputDirectory) {
  const bool ptAxis = axis=="ptz";
  const std::string suffix = ptAxis ? "central" : "ptz15to30";
  std::vector<std::vector<Point> > curves;
  for (const char *region : {"parallel","transverse"})
    for (const char *category : {"all","matched","pileup"})
      curves.push_back(truthResponse(mc,component,axis,category,region,suffix));
  const std::vector<Point> subtracted =
    truthResponse(mc,component,axis,"all","subtracted",suffix);
  if (subtracted.empty()) return false;
  for (const auto &curve : curves) if (curve.empty()) return false;
  const int colors[] = {kBlack,kBlue+1,kRed+1,kBlack,kBlue+1,kRed+1};
  const int markers[] = {kOpenCircle,kOpenSquare,kOpenTriangleUp,
                         kFullCircle,kFullSquare,kFullTriangleUp};
  const int styles[] = {1,1,1,2,2,2};
  double minimum = std::numeric_limits<double>::max();
  double maximum = -std::numeric_limits<double>::max();
  for (const auto &curve : curves)
    for (const Point &point : curve) {
      minimum = std::min(minimum,point.y-point.ey);
      maximum = std::max(maximum,point.y+point.ey);
    }
  for (const Point &point : subtracted) {
    minimum = std::min(minimum,point.y-point.ey);
    maximum = std::max(maximum,point.y+point.ey);
  }
  if (!(maximum>minimum)) { minimum = -0.5; maximum = 1.5; }
  const double margin = 0.12*(maximum-minimum);
  TH1D frame(Form("h_truth_%s_%s",component.c_str(),axis.c_str()),"",100,
             ptAxis ? 12. : 0.,ptAxis ? 200. : 5.);
  frame.SetDirectory(nullptr);
  frame.SetMinimum(minimum-margin); frame.SetMaximum(maximum+margin);
  frame.GetXaxis()->SetTitle(ptAxis ? "p_{T,Z} (GeV)" : "|#eta_{jet}|");
  frame.GetYaxis()->SetTitle(component.c_str());
  if (ptAxis) configureLogAxis(&frame);
  const TString savedLumi = lumi_136TeV;
  lumi_136TeV = ptAxis ? "Summer24 MC, |#eta|<1.305" :
                         "Summer24 MC, 15<p_{T,Z}<30 GeV";
  std::unique_ptr<TCanvas> canvas(tdrCanvas(
    Form("c_truth_%s_%s",component.c_str(),axis.c_str()),&frame,8,11,kSquare));
  lumi_136TeV = savedLumi;
  if (ptAxis) canvas->SetLogx();
  if (minimum<1. && maximum>1.) {
    TLine unity(frame.GetXaxis()->GetXmin(),1.,frame.GetXaxis()->GetXmax(),1.);
    unity.SetLineStyle(kDotted); unity.DrawClone();
  }
  std::vector<std::unique_ptr<TGraphErrors> > graphs;
  for (int index = 0; index != 6; ++index) {
    graphs.push_back(graph(curves[index],colors[index],markers[index],true));
    graphs.back()->SetLineStyle(styles[index]);
    graphs.back()->SetMarkerSize(0.55);
    graphs.back()->Draw("LP SAME");
  }
  auto subtractedGraph = graph(subtracted,kGreen+2,kFullDiamond,true);
  subtractedGraph->Draw("LP SAME");
  TLegend legend(0.46,0.56,0.92,0.89);
  legend.SetBorderSize(0); legend.SetNColumns(2);
  const char *labels[] = {"parallel all","parallel matched","parallel pileup",
                          "transverse all","transverse matched","transverse pileup"};
  for (int index = 0; index != 6; ++index)
    legend.AddEntry(graphs[index].get(),labels[index],"lp");
  legend.AddEntry(subtractedGraph.get(),"background-subtracted","lp");
  legend.Draw();
  if (ptAxis) drawMethodLogLabels(&frame);
  canvas->SaveAs((outputDirectory+"/truth_"+component+"_"+axis+".pdf").c_str());
  return true;
}

bool drawLegacyRhoAlpha(TFile &data, TFile &mc,
                        const std::string &outputDirectory) {
  const int cuts[] = {10,15,20,30};
  std::vector<std::vector<Point> > curves;
  for (int cut : cuts) {
    TProfile *dataProfile = profile(
      data,Form("legacy/control/p_rho_vs_zpt_alpha%03d",cut));
    TProfile *mcProfile = profile(
      mc,Form("legacy/control/p_rho_vs_zpt_alpha%03d",cut));
    if (!dataProfile || !mcProfile) return false;
    std::vector<Point> dataPoints;
    std::vector<Point> mcPoints;
    for (int bin = 1; bin <= dataProfile->GetNbinsX(); ++bin) {
      if (dataProfile->GetBinEntries(bin)>0.)
        dataPoints.push_back({dataProfile->GetBinCenter(bin),0.,
                              dataProfile->GetBinContent(bin),
                              dataProfile->GetBinError(bin)});
      if (mcProfile->GetBinEntries(bin)>0.)
        mcPoints.push_back({mcProfile->GetBinCenter(bin),0.,
                            mcProfile->GetBinContent(bin),
                            mcProfile->GetBinError(bin)});
    }
    curves.push_back(ratios(dataPoints,mcPoints));
  }
  std::vector<Point> extrapolated;
  for (const Point &anchor : curves.front()) {
    std::array<double,4> values{};
    bool complete = true;
    for (int index = 0; index != 4; ++index) {
      const auto point = std::find_if(
        curves[index].begin(),curves[index].end(),[&](const Point &candidate) {
          return sameCenter(candidate.x,anchor.x);
        });
      complete = complete && point!=curves[index].end();
      if (point!=curves[index].end()) values[index] = point->y;
    }
    if (!complete) continue;
    double sx = 0., sy = 0., sxx = 0., sxy = 0.;
    for (int index = 0; index != 4; ++index) {
      const double x = 0.01*cuts[index];
      sx += x; sy += values[index]; sxx += x*x; sxy += x*values[index];
    }
    const double denominator = 4.*sxx-sx*sx;
    if (denominator==0.) continue;
    const double intercept = (sxx*sy-sx*sxy)/denominator;
    const double slope = (4.*sxy-sx*sy)/denominator;
    double residual2 = 0.;
    for (int index = 0; index != 4; ++index) {
      const double x = 0.01*cuts[index];
      residual2 += std::pow(values[index]-intercept-slope*x,2);
    }
    // The alpha selections are nested, so this is a diagnostic fit-spread
    // uncertainty, not an independent-sample statistical uncertainty.
    const double variance = residual2/2.;
    const double interceptError = std::sqrt(variance*sxx/denominator);
    extrapolated.push_back({anchor.x,anchor.ex,intercept,interceptError});
  }
  TH1D frame("h_rho_alpha","",100,12.,1500.);
  frame.SetDirectory(nullptr);
  frame.SetMinimum(0.90); frame.SetMaximum(1.10);
  frame.GetXaxis()->SetTitle("p_{T,Z} (GeV)");
  frame.GetYaxis()->SetTitle("#rho(data) / #rho(MC)");
  configureLogAxis(&frame);
  std::unique_ptr<TCanvas> canvas(
    tdrCanvas("c_rho_alpha",&frame,8,11,kSquare));
  canvas->SetLogx();
  TLine unity(12.,1.,1500.,1.); unity.SetLineStyle(kDashed); unity.DrawClone();
  const int colors[] = {kBlue+2,kBlue,kOrange+7,kRed+1};
  std::vector<std::unique_ptr<TGraphErrors> > graphs;
  TLegend legend(0.56,0.65,0.91,0.88); legend.SetBorderSize(0);
  for (int index = 0; index != 4; ++index) {
    graphs.push_back(graph(curves[index],colors[index],20+index,true));
    graphs.back()->Draw("LP SAME");
    legend.AddEntry(graphs.back().get(),Form("#alpha < %.2f",0.01*cuts[index]),"lp");
  }
  auto zeroAlpha = graph(extrapolated,kBlack,kFullStar,true);
  zeroAlpha->SetLineWidth(3);
  zeroAlpha->Draw("LP SAME");
  legend.AddEntry(zeroAlpha.get(),"linear #alpha_{max}#rightarrow0","lp");
  legend.Draw(); drawMethodLogLabels(&frame);
  canvas->SaveAs((outputDirectory+"/legacy_rho_alpha_scan.pdf").c_str());
  return true;
}

} // namespace

void compareMethods(
  const char *baselineFile="rootfiles/jecdata2024I_nib1.root",
  const char *legacyFile="rootfiles/jecdata2024I_nix_legacy.root",
  const char *newMethodFile="rootfiles/jecdata2024I_nix_newmethod.root",
  const char *mcEventFile="rootfiles/zjet_MC.root",
  const char *dataEventFile="rootfiles/zjet_DATA.root",
  const char *outputDirectory="output/compareMethods") {
  applyMethodTDRStyle();
  gStyle->SetOptStat(0);
  gSystem->mkdir(outputDirectory,true);
  TFile baseline(baselineFile,"READ");
  TFile legacy(legacyFile,"READ");
  TFile modern(newMethodFile,"READ");
  TFile mcEvents(mcEventFile,"READ");
  TFile dataEvents(dataEventFile,"READ");
  if (baseline.IsZombie() || legacy.IsZombie() || modern.IsZombie() ||
      mcEvents.IsZombie() || dataEvents.IsZombie())
    throw std::runtime_error("Could not open one or more method inputs");
  const std::vector<MethodInput> inputs = {
    {&baseline,"baseline",kGray+2,kOpenCircle},
    {&legacy,"legacy",kBlue+1,kOpenSquare},
    {&modern,"new method",kRed+1,kFullCircle},
  };
  std::ofstream metrics(std::string(outputDirectory)+"/method_metrics.tsv");
  drawHDMDifference(baseline,legacy,"baseline","legacy",
                    "hdm_baseline_legacy_difference",outputDirectory,metrics);
  drawHDMDecomposition(baseline,legacy,"baseline","legacy",
                       "hdm_baseline_legacy_decomposition",outputDirectory,
                       metrics);
  drawGlobalFitDifference(baseline,legacy,"baseline","legacy",
                          "globalfit_baseline_legacy_difference",
                          outputDirectory,metrics);
  drawHDMDifference(legacy,modern,"legacy","new method",
                    "hdm_legacy_new_difference",outputDirectory,metrics);
  drawHDMDecomposition(legacy,modern,"legacy","new method",
                       "hdm_legacy_new_decomposition",outputDirectory,metrics);
  drawGlobalFitDifference(legacy,modern,"legacy","new method",
                          "globalfit_legacy_new_difference",
                          outputDirectory,metrics);
  drawMCTruthClosure(inputs,outputDirectory,metrics);
  drawWReference(baseline,legacy,modern,outputDirectory);
  for (const std::string &observable :
       {"mpfchs1","mpf1","mpfn","mpfu","mpfnu","hdm"})
    drawComponentAxes(inputs,observable,outputDirectory);
  drawRhoStability(inputs,outputDirectory);

  std::vector<std::string> truthPtPlots;
  std::vector<std::string> truthEtaPlots;
  for (const std::string &component :
       {"db","mpf","mpf1","mpfn","mpfu","mpfnu","hdm"}) {
    if (drawTruthComponent(mcEvents,component,"ptz",outputDirectory))
      truthPtPlots.push_back(component);
    if (drawTruthComponent(mcEvents,component,"abseta",outputDirectory))
      truthEtaPlots.push_back(component);
  }
  const bool alphaAvailable =
    drawLegacyRhoAlpha(dataEvents,mcEvents,outputDirectory);

  std::ofstream frames(std::string(outputDirectory)+"/compareMethods_frames.tex");
  frames << "% Generated by compareMethods.C. Do not edit.\n";
  writeFrame(frames,"Residual legacy synchronization difference",{
    {std::string(outputDirectory)+"/hdm_baseline_legacy_difference.pdf",
     "Direct legacy - baseline reference"},
    {std::string(outputDirectory)+"/hdm_baseline_legacy_decomposition.pdf",
     "Exact component decomposition"}},
    "This is the sub-50 GeV synchronization target. Component contributions "
    "sum exactly to the direct HDM difference bin by bin.");
  writeFrame(frames,"Legacy difference seen by globalFit.C",{
    {std::string(outputDirectory)+"/globalfit_baseline_legacy_difference.pdf",
     "Z-pT HDM times the stored previous JEC"}},
    "The black curve reproduces the observable passed to globalFit.C. The raw "
    "HDM and previous-JEC curves are an exact two-input Shapley decomposition; "
    "their sum equals the black curve bin by bin.");
  writeFrame(frames,"All-pairs method difference",{
    {std::string(outputDirectory)+"/hdm_legacy_new_difference.pdf",
     "Direct new - legacy reference"},
    {std::string(outputDirectory)+"/hdm_legacy_new_decomposition.pdf",
     "Exact component decomposition"}},
    "All differences are in per mille. The Shapley decomposition is exact, "
    "order independent, and sums to the direct HDM difference bin by bin.");
  writeFrame(frames,"All-pairs difference seen by globalFit.C",{
    {std::string(outputDirectory)+"/globalfit_legacy_new_difference.pdf",
     "Z-pT HDM times the stored previous JEC"}},
    "Legacy and new-method files share the same previous JEC in this sample, "
    "so their difference is dominated by the raw HDM term. The decomposition "
    "also guards against future input-JEC changes.");
  writeFrame(frames,"MC truth closure and bin-migration diagnostic",{
    {std::string(outputDirectory)+"/mc_hdm_truth_closure_zpt.pdf",
     "HDM MC response versus stored gen balance"}},
    "The lower pad is the requested combined resolution/flavor closure. A "
    "division is shown only where the stored generator balance is in 0.5--1.5; "
    "the legacy and new-method gjet inputs fail this guard over most of the "
    "range and therefore cannot yet define a physical correction.");
  writeFrame(frames,"Independent W-mass reference",{
    {std::string(outputDirectory)+"/hdm_wqq_reference.pdf",
     "Z HDM compared with $W\\to qq'$"}},
    "The lower pad uses log-pT interpolation of the W reference only for the "
    "visual Z-W residual; the measured points in the upper pad are untouched.");
  const std::vector<std::pair<std::string,std::string> > components = {
    {"mpfchs1","Total MPF"}, {"mpf1","MPF1"}, {"mpfn","MPFn"},
    {"mpfu","MPFu"}, {"mpfnu","MPFnu"}, {"hdm","HDM"},
  };
  for (const auto &component : components)
    writeFrame(frames,"Reference-axis comparison - "+component.second,{
      {std::string(outputDirectory)+"/axes_"+component.first+".pdf",
       "jet pT, average pT, and Z pT"}},
      "Legacy is shown with open markers and a line; the all-pairs method uses "
      "filled red markers. The vertical scale is explicitly in per mille.");
  writeFrame(frames,"Pileup-environment stability",{
    {std::string(outputDirectory)+"/rho_data_mc_stability.pdf",
     "rho(data) / rho(MC)"}},
    "The ratio is derived directly from the stored data and MC rho graphs; "
    "the copied jecdata ratio object stores a difference and is not used.");
  for (size_t index = 0; index < truthPtPlots.size(); index += 2) {
    std::vector<std::pair<std::string,std::string> > plots;
    for (size_t item = index; item < std::min(index+2,truthPtPlots.size()); ++item)
      plots.push_back({std::string(outputDirectory)+"/truth_"+
                       truthPtPlots[item]+"_ptz.pdf",truthPtPlots[item]});
    writeFrame(frames,"Truth-separated central response",plots,
      "All curves use |eta(jet)|<1.305. Parallel and transverse totals are "
      "decomposed into matched and pileup jets; the background-subtracted "
      "estimator is overlaid. HDM is derived bin by bin from MPF, MPFn, and MPFu.");
  }
  for (size_t index = 0; index < truthEtaPlots.size(); index += 2) {
    std::vector<std::pair<std::string,std::string> > plots;
    for (size_t item = index; item < std::min(index+2,truthEtaPlots.size()); ++item)
      plots.push_back({std::string(outputDirectory)+"/truth_"+
                       truthEtaPlots[item]+"_abseta.pdf",truthEtaPlots[item]});
    writeFrame(frames,"Truth-separated eta dependence",plots,
      "This stress test uses 15 < pT,Z < 30 GeV and is available after the "
      "new control profiles have been produced.");
  }
  if (alphaAvailable)
    writeFrame(frames,"Legacy alpha scan",{
      {std::string(outputDirectory)+"/legacy_rho_alpha_scan.pdf",
       "rho(data) / rho(MC) versus alpha max"}},
      "The black curve is the binwise linear alpha-max to zero intercept. Its "
      "error is the fit-spread diagnostic; nested selections are correlated.");
  else
    frames << "\\begin{frame}{Controls queued for the next event pass}\n"
           << "\\begin{itemize}\n"
           << "\\item Central truth-separated DB, MPF, MPF1, MPFn, MPFu, "
           << "MPFnu, and derived HDM for parallel, transverse, and "
           << "background-subtracted regions.\n"
           << "\\item The same decomposition versus jet eta for "
           << "$15<p_{T,Z}<30$ GeV.\n"
           << "\\item Legacy rho profiles for $\\alpha_{max}=0.10,0.15,0.20,0.30$ "
           << "and the subsequent $\\alpha_{max}\\rightarrow0$ extrapolation.\n"
           << "\\end{itemize}\n"
           << "These objects are now booked in \\texttt{zjet.C}; the current "
           << "ROOT files predate that booking.\n\\end{frame}\n";
  frames.close();
  metrics.close();
  std::cout << "Wrote method plots, exact HDM decomposition, and Beamer frames to "
            << outputDirectory << std::endl;
}
