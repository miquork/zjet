#include "TCanvas.h"
#include "TFile.h"
#include "TGraphErrors.h"
#include "TH1.h"
#include "TH1D.h"
#include "TLegend.h"
#include "TLine.h"
#include "TLatex.h"
#include "TPad.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TSystem.h"
#include "tdrstyle_mod22.C"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <sstream>
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

struct InputSpec {
  TFile *file = nullptr;
  std::string label;
  int color = 1;
  int marker = 20;
};

const std::vector<std::string> kSamples = {"data","mc","ratio"};
const std::vector<std::string> kChannels = {"jetz","zjav","zjet"};
double kSoftRecoilRn = 1.00;
double kSoftRecoilRu = 0.92;
double kGlobalFitJetPtScale = 1.0000;
double kGlobalFitAveragePtScale = 1.0025;

void applyComparisonTDRStyle() {
  setTDRStyle();
  writeExtraText = true;
  extraText = "Work in progress";
  extraText2 = "";
  // CMS_lumi appends the collision energy for iPeriod=8.
  lumi_136TeV = "Run2024I";
}

std::string objectPath(const std::string &sample,
                       const std::string &observable,
                       const std::string &channel);

std::vector<ComparisonPoint> readPoints(TFile &file, const std::string &path) {
  TObject *object = file.Get(path.c_str());
  std::vector<ComparisonPoint> points;
  if (TGraphErrors *graph = dynamic_cast<TGraphErrors*>(object)) {
    for (int i=0; i<graph->GetN(); ++i) {
      double x = 0.;
      double y = 0.;
      graph->GetPoint(i,x,y);
      if (!std::isfinite(x) || !std::isfinite(y)) continue;
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

double countIntegral(TFile &file, const std::string &sample,
                     const std::string &channel) {
  const std::string path = sample+"/eta00-13/counts_"+channel+"_a100";
  TH1 *histogram = dynamic_cast<TH1*>(file.Get(path.c_str()));
  if (!histogram) return 0.;
  return histogram->Integral(0,histogram->GetNbinsX()+1);
}

std::vector<ComparisonPoint> normalized(
  const std::vector<ComparisonPoint> &input) {
  double integral = 0.;
  for (const ComparisonPoint &point : input) integral += point.y;
  if (!std::isfinite(integral) || integral==0.) return {};
  std::vector<ComparisonPoint> result = input;
  for (ComparisonPoint &point : result) {
    point.y /= integral;
    point.ey /= fabs(integral);
  }
  return result;
}

bool sameCenter(double first, double second) {
  return fabs(first-second)<1.e-5*std::max(1.,fabs(second));
}

std::vector<ComparisonPoint> pointRatios(
  const std::vector<ComparisonPoint> &numerator,
  const std::vector<ComparisonPoint> &denominator) {
  std::vector<ComparisonPoint> ratios;
  for (const ComparisonPoint &num : numerator) {
    const auto den = std::find_if(
      denominator.begin(),denominator.end(),[&](const ComparisonPoint &point) {
        return sameCenter(point.x,num.x);
      });
    if (den==denominator.end() || den->y==0.) continue;
    const double value = num.y/den->y;
    double relativeVariance = 0.;
    if (num.y!=0.) relativeVariance += pow(num.ey/num.y,2);
    relativeVariance += pow(den->ey/den->y,2);
    ratios.push_back({num.x,num.ex,value,
                      fabs(value)*sqrt(relativeVariance)});
  }
  return ratios;
}

std::vector<ComparisonPoint> pointDifferences(
  const std::vector<ComparisonPoint> &first,
  const std::vector<ComparisonPoint> &second) {
  std::vector<ComparisonPoint> differences;
  for (const ComparisonPoint &left : first) {
    const auto right = std::find_if(
      second.begin(),second.end(),[&](const ComparisonPoint &point) {
        return sameCenter(point.x,left.x);
      });
    if (right==second.end()) continue;
    differences.push_back({left.x,left.ex,left.y-right->y,
                           std::hypot(left.ey,right->ey)});
  }
  return differences;
}

std::vector<ComparisonPoint> correctedSoftRecoil(
  TFile &file, const std::string &sample, const std::string &channel) {
  const std::vector<ComparisonPoint> neutral = readPoints(
    file,objectPath(sample,"mpfn",channel));
  const std::vector<ComparisonPoint> unclustered = readPoints(
    file,objectPath(sample,"mpfu",channel));
  std::vector<ComparisonPoint> result;
  for (const ComparisonPoint &n : neutral) {
    const auto u = std::find_if(
      unclustered.begin(),unclustered.end(),[&](const ComparisonPoint &point) {
        return sameCenter(point.x,n.x);
      });
    if (u==unclustered.end()) continue;
    result.push_back({
      n.x,n.ex,n.y/kSoftRecoilRn+u->y/kSoftRecoilRu,
      std::hypot(n.ey/kSoftRecoilRn,u->ey/kSoftRecoilRu)
    });
  }
  return result;
}

double interpolatePoints(const std::vector<ComparisonPoint> &points,
                         double x) {
  for (const ComparisonPoint &point : points)
    if (sameCenter(point.x,x)) return point.y;
  for (size_t index=1; index<points.size(); ++index) {
    if (points[index-1].x<x && x<points[index].x) {
      const double fraction =
        (std::log(x)-std::log(points[index-1].x))/
        (std::log(points[index].x)-std::log(points[index-1].x));
      return points[index-1].y+fraction*(points[index].y-points[index-1].y);
    }
  }
  return std::numeric_limits<double>::quiet_NaN();
}

double globalFitAxisScale(const std::string &channel) {
  if (channel=="jetz") return kGlobalFitJetPtScale;
  if (channel=="zjav") return kGlobalFitAveragePtScale;
  return 1.0000;
}

std::vector<ComparisonPoint> globalFitHDMPoints(
  TFile &file, const std::string &channel) {
  std::vector<ComparisonPoint> result = readPoints(
    file,"ratio/eta00-13/hdm_mpfchs1_"+channel);
  const std::vector<ComparisonPoint> correction = readPoints(
    file,"ratio/eta00-13/herr_l2l3res");
  const double axisScale = globalFitAxisScale(channel);
  for (ComparisonPoint &point : result) {
    const double factor = interpolatePoints(correction,point.x);
    if (!std::isfinite(factor)) {
      point.y = std::numeric_limits<double>::quiet_NaN();
      continue;
    }
    point.y *= factor*axisScale;
    point.ey *= std::fabs(factor*axisScale);
  }
  result.erase(std::remove_if(result.begin(),result.end(),
    [](const ComparisonPoint &point) { return !std::isfinite(point.y); }),
    result.end());
  return result;
}

std::string objectPath(const std::string &sample,
                       const std::string &observable,
                       const std::string &channel) {
  const std::string base = sample+"/eta00-13/";
  if (observable=="counts") return base+"counts_"+channel+"_a100";
  if (observable=="hdm") return base+"hdm_mpfchs1_"+channel;
  return base+"orig/"+observable+"_"+channel+"_a100";
}

std::vector<ComparisonPoint> comparisonPoints(
  TFile &file, const std::string &sample, const std::string &observable,
  const std::string &channel) {
  if (observable=="hdm_globalfit")
    return sample=="ratio" ? globalFitHDMPoints(file,channel)
                           : std::vector<ComparisonPoint>{};
  if (observable=="mpfnu") {
    if (sample=="data" || sample=="mc")
      return correctedSoftRecoil(file,sample,channel);
    return pointDifferences(correctedSoftRecoil(file,"data",channel),
                            correctedSoftRecoil(file,"mc",channel));
  }
  if (observable!="counts")
    return readPoints(file,objectPath(sample,observable,channel));
  if (sample=="data" || sample=="mc")
    return normalized(readPoints(file,objectPath(sample,observable,channel)));
  const std::vector<ComparisonPoint> data = normalized(
    readPoints(file,objectPath("data",observable,channel)));
  const std::vector<ComparisonPoint> mc = normalized(
    readPoints(file,objectPath("mc",observable,channel)));
  return pointRatios(data,mc);
}

TGraphErrors *makeGraph(const std::vector<ComparisonPoint> &points, int color,
                        int marker, double markerSize=0.72) {
  TGraphErrors *graph = new TGraphErrors(points.size());
  for (size_t i=0; i<points.size(); ++i) {
    graph->SetPoint(i,points[i].x,points[i].y);
    graph->SetPointError(i,points[i].ex,points[i].ey);
  }
  graph->SetMarkerStyle(marker);
  graph->SetMarkerSize(markerSize);
  graph->SetMarkerColor(color);
  graph->SetLineColor(color);
  graph->SetLineWidth(2);
  return graph;
}

bool isComponent(const std::string &observable) {
  return observable=="mpf1" || observable=="mpfn" ||
         observable=="mpfu" || observable=="mpfnu";
}

bool resultIsDifference(const std::string &observable) {
  return isComponent(observable) || observable=="rho" ||
         observable=="chf" || observable=="nef" ||
         observable=="nhf" || observable=="cef" || observable=="muf";
}

std::pair<double,double> yRange(const std::string &observable,
                                const std::string &sample) {
  if (observable=="counts")
    return sample=="ratio" ? std::make_pair(0.45,1.55)
                           : std::make_pair(1.e-7,1.);
  if (sample=="ratio" && resultIsDifference(observable)) {
    if (observable=="rho") return {-12.,12.};
    return {-0.15,0.15};
  }
  if (observable=="rho") return {0.,60.};
  if (observable=="chf" || observable=="nef" || observable=="nhf" ||
      observable=="cef" || observable=="muf") return {0.,1.};
  if (observable=="mpfn") return {-0.6,0.3};
  if (observable=="mpfu" || observable=="mpfnu") return {-0.2,0.8};
  if (observable=="hdm")
    return sample=="ratio" ? std::make_pair(0.94,1.06)
                           : std::make_pair(0.75,1.25);
  if (observable=="hdm_globalfit") return {0.98,1.18};
  if (observable=="rjet") return {0.75,1.45};
  if (observable=="gjet") return {-0.15,1.15};
  if (sample=="ratio") return {0.75,1.25};
  return {0.5,1.5};
}

double comparisonRange(const std::string &observable,
                       const std::string &sample) {
  if (observable=="counts") return 1.;
  if (observable=="rho") return 20.;
  if (observable=="hdm") return 0.03;
  if (observable=="hdm_globalfit") return 0.04;
  if (observable=="gjet") return 1.0;
  if (observable=="rjet") return 0.5;
  if (sample=="ratio" && resultIsDifference(observable)) return 0.12;
  if (observable=="chf" || observable=="nef" || observable=="nhf" ||
      observable=="cef" || observable=="muf") return 0.20;
  return 0.25;
}

std::string observableLabel(const std::string &observable) {
  if (observable=="counts") return "Statistics shape";
  if (observable=="ptchs") return "Direct balance";
  if (observable=="mpfchs1") return "MPF";
  if (observable=="mpf1") return "MPF1";
  if (observable=="mpfn") return "MPFn";
  if (observable=="mpfu") return "MPFu";
  if (observable=="mpfnu") return "Response-corrected MPFnu";
  if (observable=="hdm") return "HDM-corrected MPF";
  if (observable=="hdm_globalfit") return "globalFit input: HDM x previous JEC";
  if (observable=="rjet") return "Reconstructed balance closure";
  if (observable=="gjet") return "Generator balance closure";
  return observable;
}

std::string channelLabel(const std::string &channel) {
  if (channel=="jetz") return "jet pT binning";
  if (channel=="zjav") return "average pT (HDM) binning";
  return "Z pT binning";
}

std::string sampleLabel(const std::string &sample,
                        const std::string &observable) {
  if (sample=="data") return "data";
  if (sample=="mc") return "MC";
  if (observable=="counts") return "normalized data/MC";
  return resultIsDifference(observable) ? "data-MC" : "data/MC";
}

std::string texSampleLabel(const std::string &sample,
                           const std::string &observable) {
  if (sample=="data") return "data";
  if (sample=="mc") return "MC";
  if (observable=="counts") return "normalized data/MC";
  return resultIsDifference(observable) ? "data--MC" : "data/MC";
}

std::string texEscape(std::string value) {
  size_t position = 0;
  while ((position=value.find('_',position))!=std::string::npos) {
    value.replace(position,1,"\\_");
    position += 2;
  }
  return value;
}

std::string scientific(double value) {
  std::ostringstream stream;
  stream << std::scientific << std::setprecision(3) << value;
  return stream.str();
}

std::vector<ComparisonPoint> matchedComparisons(
  const std::vector<ComparisonPoint> &reference,
  const std::vector<ComparisonPoint> &candidate, bool relative, int &matched) {
  std::vector<ComparisonPoint> differences;
  matched = 0;
  for (const ComparisonPoint &point : candidate) {
    const auto base = std::find_if(
      reference.begin(),reference.end(),[&](const ComparisonPoint &other) {
        return sameCenter(other.x,point.x);
      });
    if (base==reference.end()) continue;
    ++matched;
    double value = point.y-base->y;
    double error = std::hypot(point.ey,base->ey);
    if (relative) {
      if (base->y==0.) continue;
      value /= base->y;
      error /= fabs(base->y);
    }
    differences.push_back({point.x,point.ex,value,error});
  }
  return differences;
}

void drawMissing(const std::string &message, double size=0.040) {
  TLatex latex;
  latex.SetNDC();
  latex.SetTextAlign(22);
  latex.SetTextSize(size);
  latex.DrawLatex(0.5,0.52,message.c_str());
}

void drawZeroLine(double xmin=12., double xmax=1500.) {
  TLine line(xmin,0.,xmax,0.);
  line.SetLineStyle(kDashed);
  line.SetLineColor(kGray+2);
  line.DrawClone();
}

void configureLogAxis(TAxis *axis) {
  if (!axis) return;
  axis->SetMoreLogLabels();
  axis->SetNoExponent();
}

void drawComparisonLogXLabels(
  TH1 *frame, const std::vector<double> &values={30.,100.,300.,1000.},
  double offset=0.014, double textSize=-1.) {
  if (!frame || !gPad || !gPad->GetLogx()) return;
  gPad->Update();
  const double xmin = frame->GetXaxis()->GetXmin();
  const double xmax = frame->GetXaxis()->GetXmax();
  if (xmin<=0. || xmax<=xmin) return;
  frame->GetXaxis()->SetLabelOffset(999.);
  gPad->Modified();
  gPad->Update();
  const double left = gPad->GetLeftMargin();
  const double right = gPad->GetRightMargin();
  const double bottom = gPad->GetBottomMargin();
  TLatex label;
  label.SetNDC();
  label.SetTextAlign(23);
  label.SetTextFont(frame->GetXaxis()->GetLabelFont());
  label.SetTextSize(textSize>0. ? textSize : frame->GetXaxis()->GetLabelSize());
  for (double value : values) {
    if (value<xmin || value>xmax) continue;
    const double fraction = (std::log10(value)-std::log10(xmin))/
                            (std::log10(xmax)-std::log10(xmin));
    const double x = left+fraction*(1.-left-right);
    label.DrawLatex(x,bottom-offset,Form("%g",value));
  }
}

std::string standardPanel(
  const std::vector<InputSpec> &inputs, const std::string &sample,
  const std::string &observable, const std::string &channel,
  const std::string &outputDirectory, std::ofstream &summary) {
  std::vector<std::vector<ComparisonPoint> > points;
  for (const InputSpec &input : inputs)
    points.push_back(comparisonPoints(*input.file,sample,observable,channel));

  int legacyMatched = 0;
  int newMatched = 0;
  const bool relative = observable=="counts";
  const std::vector<ComparisonPoint> legacyDifference = matchedComparisons(
    points[0],points[1],relative,legacyMatched);
  const std::vector<ComparisonPoint> newDifference = matchedComparisons(
    points[0],points[2],relative,newMatched);
  summary << sample << '\t' << observable << '\t' << channel;
  for (const auto &collection : points) summary << '\t' << collection.size();
  summary << '\t' << legacyMatched << '\t' << newMatched << '\n';

  const std::string stem = "compare_"+sample+"_"+observable+"_"+channel;
  const auto range = yRange(observable,sample);
  TH1D upperFrame((stem+"_upper_frame").c_str(),"",100,12.,1500.);
  upperFrame.SetDirectory(nullptr);
  upperFrame.SetMinimum(range.first);
  upperFrame.SetMaximum(range.second);
  upperFrame.GetYaxis()->SetTitle(
    (observableLabel(observable)+" ("+sampleLabel(sample,observable)+")").c_str());
  upperFrame.GetYaxis()->SetTitleSize(0.050);
  upperFrame.GetYaxis()->SetTitleOffset(1.25);
  upperFrame.GetYaxis()->SetLabelSize(0.043);
  upperFrame.GetXaxis()->SetLabelSize(0.);
  configureLogAxis(upperFrame.GetXaxis());
  const double deltaRange = comparisonRange(observable,sample);
  TH1D lowerFrame((stem+"_lower_frame").c_str(),"",100,12.,1500.);
  lowerFrame.SetDirectory(nullptr);
  lowerFrame.SetMinimum(-deltaRange);
  lowerFrame.SetMaximum(deltaRange);
  lowerFrame.GetXaxis()->SetTitle("p_{T} (GeV)");
  lowerFrame.GetYaxis()->SetTitle(
    observable=="counts" ? "cand./base - 1" : "cand. - base");
  configureLogAxis(lowerFrame.GetXaxis());
  std::unique_ptr<TCanvas> canvas(
    tdrDiCanvas(stem.c_str(),&upperFrame,&lowerFrame,8,11));

  canvas->cd(1);
  gPad->SetLogx();
  if (observable=="counts" && sample!="ratio") gPad->SetLogy();

  std::vector<std::unique_ptr<TGraphErrors> > graphs;
  bool anyPoints = false;
  for (size_t index=0; index<inputs.size(); ++index) {
    graphs.emplace_back(makeGraph(points[index],inputs[index].color,
                                  inputs[index].marker));
    if (!points[index].empty()) {
      graphs.back()->Draw("P SAME");
      anyPoints = true;
    }
  }
  if (!anyPoints)
    drawMissing("Object unavailable in all three files");
  TLegend legend(0.53,0.69,0.92,0.90);
  legend.SetBorderSize(0);
  legend.SetFillStyle(0);
  legend.SetTextSize(0.036);
  for (size_t index=0; index<inputs.size(); ++index)
    legend.AddEntry(graphs[index].get(),inputs[index].label.c_str(),"pl");
  legend.Draw();
  TLatex label;
  label.SetNDC();
  label.SetTextSize(0.040);
  label.DrawLatex(0.61,0.64,channelLabel(channel).c_str());

  canvas->cd(2);
  gPad->SetLogx();
  drawZeroLine();
  std::unique_ptr<TGraphErrors> legacyGraph(
    makeGraph(legacyDifference,inputs[1].color,inputs[1].marker,0.68));
  std::unique_ptr<TGraphErrors> newGraph(
    makeGraph(newDifference,inputs[2].color,inputs[2].marker,0.68));
  if (!legacyDifference.empty()) legacyGraph->Draw("P SAME");
  if (!newDifference.empty()) newGraph->Draw("P SAME");
  TLatex matchLabel;
  matchLabel.SetNDC();
  matchLabel.SetTextSize(0.072);
  matchLabel.DrawLatex(0.17,0.84,
    Form("matched bins: legacy %d, new %d",legacyMatched,newMatched));
  drawComparisonLogXLabels(&lowerFrame,{30.,100.,300.,1000.},0.045,0.10);

  const std::string plotPath = outputDirectory+"/"+stem+".pdf";
  canvas->SaveAs(plotPath.c_str());
  return stem;
}

std::string axisOverlayPanel(
  const InputSpec &input, const std::string &sample,
  const std::string &observable, const std::string &outputDirectory) {
  const std::vector<int> colors = {kBlue+1,kGreen+2,kRed+1};
  const std::vector<int> markers = {20,21,24};
  std::vector<std::vector<ComparisonPoint> > points;
  for (const std::string &channel : kChannels)
    points.push_back(comparisonPoints(*input.file,sample,observable,channel));

  std::string labelStem = input.label;
  std::replace(labelStem.begin(),labelStem.end(),' ','_');
  const std::string stem =
    "axes_"+labelStem+"_"+sample+"_"+observable;
  const auto range = yRange(observable,sample);
  TH1D frame((stem+"_frame").c_str(),"",100,12.,1500.);
  frame.SetDirectory(nullptr);
  frame.SetMinimum(range.first);
  frame.SetMaximum(range.second);
  frame.GetXaxis()->SetTitle("reference p_{T} (GeV)");
  frame.GetYaxis()->SetTitle(
    (observableLabel(observable)+" ("+sampleLabel(sample,observable)+")").c_str());
  frame.GetXaxis()->SetTitleSize(0.050);
  frame.GetXaxis()->SetLabelSize(0.042);
  configureLogAxis(frame.GetXaxis());
  frame.GetYaxis()->SetTitleSize(0.050);
  frame.GetYaxis()->SetTitleOffset(1.25);
  frame.GetYaxis()->SetLabelSize(0.042);
  std::unique_ptr<TCanvas> canvas(
    tdrCanvas(stem.c_str(),&frame,8,11,kSquare));
  canvas->SetLogx();
  if (observable=="counts" && sample!="ratio") canvas->SetLogy();

  std::vector<std::unique_ptr<TGraphErrors> > graphs;
  bool anyPoints = false;
  for (size_t index=0; index<kChannels.size(); ++index) {
    graphs.emplace_back(makeGraph(points[index],colors[index],markers[index]));
    if (!points[index].empty()) {
      graphs.back()->Draw("P SAME");
      anyPoints = true;
    }
  }
  if (!anyPoints) drawMissing("All reference-axis objects unavailable");
  TLegend legend(0.53,0.70,0.92,0.90);
  legend.SetBorderSize(0);
  legend.SetFillStyle(0);
  legend.SetTextSize(0.037);
  for (size_t index=0; index<kChannels.size(); ++index)
    legend.AddEntry(graphs[index].get(),channelLabel(kChannels[index]).c_str(),"pl");
  legend.Draw();
  TLatex label;
  label.SetNDC();
  label.SetTextSize(0.043);
  label.DrawLatex(0.56,0.64,input.label.c_str());
  drawComparisonLogXLabels(&frame,{30.,100.,300.,1000.},0.022,0.042);
  const std::string plotPath = outputDirectory+"/"+stem+".pdf";
  canvas->SaveAs(plotPath.c_str());
  return stem;
}

void writeNormalizationFrame(std::ofstream &frames,
                             const std::vector<InputSpec> &inputs,
                             std::ofstream &normalization) {
  frames << "\\begin{frame}{Statistics normalization diagnostics}\n"
         << "\\centering\\scriptsize\n"
         << "\\begin{tabular}{llrrr}\n"
         << "Reference axis & input & data yield & MC yield & MC/data \\\\\n"
         << "\\hline\n";
  normalization << "channel\tinput\tdata_integral\tmc_integral\tmc_over_data\n";
  for (const std::string &channel : kChannels) {
    for (const InputSpec &input : inputs) {
      const double data = countIntegral(*input.file,"data",channel);
      const double mc = countIntegral(*input.file,"mc",channel);
      const double ratio = data!=0. ? mc/data : 0.;
      frames << texEscape(channelLabel(channel)) << " & "
             << texEscape(input.label) << " & " << scientific(data) << " & "
             << scientific(mc) << " & " << scientific(ratio) << " \\\\\n";
      normalization << channel << '\t' << input.label << '\t' << data << '\t'
                    << mc << '\t' << ratio << '\n';
    }
  }
  frames << "\\end{tabular}\n"
         << "\\vspace{1ex}\\parbox{0.94\\linewidth}{\\tiny Raw yields are "
         << "reported only as a diagnostic. Statistics plots use a separate "
         << "unit-area normalization for each data or MC distribution; their "
         << "third column is the derived normalized data/MC shape ratio.}\n"
         << "\\end{frame}\n";
}

double maximumAbsoluteValue(TFile &file, const std::string &observable) {
  double maximum = 0.;
  for (const std::string &sample : std::vector<std::string>{"data","mc"}) {
    for (const std::string &channel : kChannels) {
      const std::vector<ComparisonPoint> points = readPoints(
        file,objectPath(sample,observable,channel));
      for (const ComparisonPoint &point : points)
        maximum = std::max(maximum,fabs(point.y));
    }
  }
  return maximum;
}

void writeInputSanityFrame(std::ofstream &frames,
                           const std::vector<InputSpec> &inputs) {
  frames << "\\begin{frame}{Input sanity checks}\n"
         << "\\centering\\small\n"
         << "\\begin{tabular}{lrrl}\n"
         << "Input & max $|\\mathrm{MPF1}|$ & max $|\\mathrm{MPFn}|$ & status \\\\\n"
         << "\\hline\n";
  bool failed = false;
  for (const InputSpec &input : inputs) {
    const double mpf1 = maximumAbsoluteValue(*input.file,"mpf1");
    const double mpfn = maximumAbsoluteValue(*input.file,"mpfn");
    const bool pass = mpf1<5. && mpfn<5.;
    failed |= !pass;
    frames << texEscape(input.label) << " & " << scientific(mpf1) << " & "
           << scientific(mpfn) << " & " << (pass ? "PASS" : "\\textbf{FAIL}")
           << " \\\\\n";
  }
  frames << "\\end{tabular}\n"
         << "\\vspace{1.5ex}\\parbox{0.92\\linewidth}{\\small ";
  if (failed)
    frames << "At least one input contains unphysical longitudinal MPF "
           << "components. Those component panels may be empty because their "
           << "values lie outside the fixed validation range. The total MPF "
           << "can nevertheless remain finite when MPF1 and MPFn cancel. "
           << "Regenerate that input after the transverse-projection fix.";
  else
    frames << "All inputs pass the broad MPF-component magnitude check. This "
           << "is a corruption guard, not a physics-quality criterion.";
  frames << "}\n\\end{frame}\n";
}

void writeDerivedObservableFrame(std::ofstream &frames) {
  frames << "\\begin{frame}{Derived observables and MC-only closure}\n"
         << "\\small\n"
         << "\\begin{itemize}\n"
         << "\\item Response-corrected soft recoil is derived bin by bin as "
         << "$\\mathrm{MPF}_{nu}=\\mathrm{MPF}_{n}/R_n+"
         << "\\mathrm{MPF}_{u}/R_u$, with $R_n=" << kSoftRecoilRn
         << "$ and $R_u=" << kSoftRecoilRu
         << "$ from the current \\texttt{softrad3.C} defaults. "
         << "Input errors are propagated without an unavailable covariance.\n"
         << "\\item \\texttt{rjet} and \\texttt{gjet} are written by "
         << "\\texttt{reprocess.C} only for MC. They are therefore shown as "
         << "MC closure comparisons, without meaningless empty data or "
         << "data/MC panels.\n"
         << "\\item Direct balance in data, MC, and data/MC remains the "
         << "\\textbf{Direct balance} (\\texttt{ptchs}) page.\n"
         << "\\item Generator-balance curves from event outputs made before "
         << "the transverse-projection correction are invalid. In the current "
         << "legacy and new-method files the stored generator balance still "
         << "falls far below the physical response, so it is shown as a "
         << "diagnostic but is not used to correct HDM.\n"
         << "\\end{itemize}\n"
         << "\\end{frame}\n";
}

bool available(const std::string &path) {
  return !gSystem->AccessPathName(path.c_str());
}

void writeImageFrame(
  std::ofstream &frames, const std::string &title,
  const std::vector<std::pair<std::string,std::string> > &images) {
  std::vector<std::pair<std::string,std::string> > present;
  for (const auto &image : images)
    if (available(image.first)) present.push_back(image);
  if (present.empty()) return;
  frames << "\\begin{frame}{" << title << "}\n"
         << "\\begin{columns}[T,onlytextwidth]\n";
  const double width = 0.98/present.size();
  for (const auto &image : present) {
    frames << "\\begin{column}{" << std::fixed << std::setprecision(3)
           << width << "\\textwidth}\n"
           << "\\centering\\textbf{" << image.second << "}\\par\n"
           << "\\includegraphics[width=\\linewidth,height=0.77\\textheight,"
           << "keepaspectratio]{" << image.first << "}\n"
           << "\\end{column}\n";
  }
  frames << "\\end{columns}\n\\end{frame}\n";
}

void writeControlAppendix(std::ofstream &frames) {
  frames << "\\section{Analysis controls}\n";
  writeImageFrame(frames,"Pileup stability of the subtracted response",{
    {"pdf/drawPileupControl/db_vs_rho.pdf","Direct balance vs. rho"},
    {"pdf/drawPileupControl/mpf_vs_rho.pdf","Hybrid MPF vs. rho"},
  });
  writeImageFrame(frames,"Pileup-jet subtraction controls",{
    {"pdf/drawPileupControl/pujet_fraction_vs_ptz.pdf",
     "Unmatched fraction vs. Z pT"},
    {"pdf/drawPileupControl/pujet_fraction_eta_vs_ptz.pdf",
     "Eta dependence after subtraction"},
  });
  writeImageFrame(frames,"Truth-separated response controls",{
    {"pdf/drawPileupControl/db_truth_parallel.pdf",
     "Parallel direct balance"},
    {"pdf/drawPileupControl/match_definition_vs_ptz.pdf",
     "Reco-to-gen matching definition"},
  });
  writeImageFrame(frames,"Data and MC control distributions",{
    {"pdf/drawControl/drawControlData_c1_1_pt.pdf","Probe-jet pT"},
    {"pdf/drawControl/drawControlData_c1_3_db.pdf","Direct balance"},
  });
  writeImageFrame(frames,"Response maps in the signal and sidebands",{
    {"pdf/drawControl/drawControl_c2_2_db.pdf","Direct balance map"},
    {"pdf/drawControl/drawControl_c2_3_mpf.pdf","MPF map"},
  });
  writeImageFrame(frames,"All-pairs method geometry",{
    {"pdf/drawZjetSummary/regions_db.pdf","Parallel, transverse, subtracted"},
  });
  writeImageFrame(frames,"All-pairs balance slices",{
    {"pdf/drawZjetSummary/balance_slices.pdf",
     "Data and MC at approximately 30 GeV"},
  });
}

} // namespace

void compareJECdata(
  const char *baselineFile=
    "rootfiles/jecdata2024I_nib1.root",
  const char *legacyFile=
    "rootfiles/jecdata2024I_nix_legacy.root",
  const char *newMethodFile=
    "rootfiles/jecdata2024I_nix_newmethod.root",
  const char *outputDirectory="output/compareJECdata",
  const char *baselineLabel="2024I_nib1 baseline",
  const char *legacyLabel="2024I_nix legacy",
  const char *newMethodLabel="2024I_nix new method",
  double softRecoilRn=1.00,
  double softRecoilRu=0.92,
  double globalFitJetPtScale=1.0000,
  double globalFitAveragePtScale=1.0025) {
  if (softRecoilRn<=0. || softRecoilRu<=0.)
    throw std::runtime_error("Soft-recoil responses Rn and Ru must be positive");
  if (globalFitJetPtScale<=0. || globalFitAveragePtScale<=0.)
    throw std::runtime_error("globalFit axis scales must be positive");
  kSoftRecoilRn = softRecoilRn;
  kSoftRecoilRu = softRecoilRu;
  kGlobalFitJetPtScale = globalFitJetPtScale;
  kGlobalFitAveragePtScale = globalFitAveragePtScale;
  applyComparisonTDRStyle();
  std::unique_ptr<TFile> baseline(TFile::Open(baselineFile,"READ"));
  std::unique_ptr<TFile> legacy(TFile::Open(legacyFile,"READ"));
  std::unique_ptr<TFile> newMethod(TFile::Open(newMethodFile,"READ"));
  if (!baseline || baseline->IsZombie() || !legacy || legacy->IsZombie() ||
      !newMethod || newMethod->IsZombie())
    throw std::runtime_error("Could not open one or more jecdata inputs");
  gSystem->mkdir(outputDirectory,true);
  gStyle->SetOptStat(0);
  gStyle->SetErrorX(0.5);

  const std::vector<InputSpec> inputs = {
    {baseline.get(),baselineLabel,kBlack,20},
    {legacy.get(),legacyLabel,kBlue+1,21},
    {newMethod.get(),newMethodLabel,kRed+1,24},
  };
  const std::vector<std::string> observables = {
    "counts", "ptchs", "mpfchs1", "mpf1", "mpfn", "mpfu", "mpfnu",
    "hdm", "hdm_globalfit", "rjet", "gjet", "chf", "nef", "nhf", "cef",
    "muf", "rho",
  };

  std::ofstream summary(std::string(outputDirectory)+"/summary.tsv");
  summary << "sample\tobservable\tchannel\tbaseline_points\tlegacy_points"
          << "\tnew_points\tlegacy_matched\tnew_matched\n";
  std::ofstream normalization(
    std::string(outputDirectory)+"/normalization.tsv");
  std::ofstream frames(
    std::string(outputDirectory)+"/compareJECdata_frames.tex");
  frames << "% Generated by compareJECdata.C. Do not edit by hand.\n";
  writeNormalizationFrame(frames,inputs,normalization);
  writeInputSanityFrame(frames,inputs);
  writeDerivedObservableFrame(frames);

  for (const std::string &observable : observables) {
    const bool hasAllReferenceAxes =
      observable=="counts" || observable=="ptchs" ||
      observable=="mpfchs1" || observable=="mpf1" ||
      observable=="mpfn" || observable=="mpfu" || observable=="mpfnu" ||
      observable=="hdm" || observable=="hdm_globalfit";
    const std::vector<std::string> channels =
      hasAllReferenceAxes ? kChannels : std::vector<std::string>{"zjet"};
    for (const std::string &channel : channels) {
      if (observable=="hdm_globalfit") {
        const std::string stem = standardPanel(
          inputs,"ratio",observable,channel,outputDirectory,summary);
        frames << "\\begin{frame}{" << observableLabel(observable) << " -- "
               << channelLabel(channel) << "}\n"
               << "\\centering\\includegraphics[width=0.68\\linewidth,"
               << "height=0.73\\textheight,keepaspectratio]{"
               << outputDirectory << "/" << stem << ".pdf}\n"
               << "\\vspace{0.5ex}\\par{\\tiny This reproduces the quantity "
               << "passed to globalFit.C: the raw HDM data/MC ratio is "
               << "multiplied by the input file's own "
               << "\\texttt{herr\\_l2l3res}. The configured jet-pT and "
               << "average-pT axis factors (" << kGlobalFitJetPtScale
               << " and " << kGlobalFitAveragePtScale
               << ") are also included.}\n\\end{frame}\n";
        continue;
      }
      if (observable=="rjet" || observable=="gjet") {
        const std::string stem = standardPanel(
          inputs,"mc",observable,channel,outputDirectory,summary);
        frames << "\\begin{frame}{" << observableLabel(observable) << " -- "
               << channelLabel(channel) << "}\n"
               << "\\centering\\textbf{MC-only closure}\\par\n"
               << "\\includegraphics[width=0.66\\linewidth,height=0.70\\textheight,"
               << "keepaspectratio]{" << outputDirectory << "/" << stem
               << ".pdf}\n"
               << "\\vspace{0.5ex}\\par{\\tiny ";
        if (observable=="rjet")
          frames << "The corresponding direct-balance response in data, MC, "
                 << "and data/MC is shown on the Direct balance "
                 << "(\\texttt{ptchs}) pages.";
        else
          frames << "Generator balance is an MC closure diagnostic. The "
                 << "legacy and new-method curves fail the physical closure "
                 << "guard in the current files and must not be used as a "
                 << "resolution or flavor correction.";
        frames << "}\n\\end{frame}\n";
        continue;
      }
      frames << "\\begin{frame}{" << observableLabel(observable) << " -- "
             << channelLabel(channel) << "}\n"
             << "\\begin{columns}[T,onlytextwidth]\n";
      for (const std::string &sample : kSamples) {
        const std::string stem = standardPanel(
          inputs,sample,observable,channel,outputDirectory,summary);
        frames << "\\begin{column}{0.327\\textwidth}\n"
               << "\\centering\\textbf{" << texSampleLabel(sample,observable)
               << "}\\par\n"
               << "\\includegraphics[width=\\linewidth]{"
               << outputDirectory << "/" << stem << ".pdf}\n"
               << "\\end{column}\n";
      }
      frames << "\\end{columns}\n"
             << "\\vspace{0.3ex}{\\tiny Top: baseline, synchronized legacy, "
             << "and new method on fixed axes. Bottom: legacy--baseline and "
             << "new--baseline at exactly common bin centers. Statistics "
             << "differences are relative; response differences are absolute.}\n"
             << "\\end{frame}\n";
    }
  }

  const std::vector<std::string> axisObservables = {
    "counts","ptchs","mpfchs1","mpf1","mpfn","mpfu","mpfnu","hdm"
  };
  for (const std::string &observable : axisObservables) {
    for (const InputSpec &input : inputs) {
      frames << "\\begin{frame}{Reference-axis overlay -- "
             << observableLabel(observable) << " -- "
             << texEscape(input.label) << "}\n"
             << "\\begin{columns}[T,onlytextwidth]\n";
      for (const std::string &sample : kSamples) {
        const std::string stem = axisOverlayPanel(
          input,sample,observable,outputDirectory);
        frames << "\\begin{column}{0.327\\textwidth}\n"
               << "\\centering\\textbf{" << texSampleLabel(sample,observable)
               << "}\\par\n"
               << "\\includegraphics[width=\\linewidth]{"
               << outputDirectory << "/" << stem << ".pdf}\n"
               << "\\end{column}\n";
      }
      frames << "\\end{columns}\n"
             << "\\vspace{0.3ex}{\\tiny Jet-pT, average-pT (HDM), and Z-pT "
             << "binnings are overlaid with identical plot geometry. This "
             << "exposes low-pT pileup and high-pT resolution trade-offs "
             << "without moving axes between methods.}\n"
             << "\\end{frame}\n";
    }
  }

  for (const InputSpec &input : inputs) {
    const std::string stem = axisOverlayPanel(
      input,"ratio","hdm_globalfit",outputDirectory);
    frames << "\\begin{frame}{globalFit input across reference axes -- "
           << texEscape(input.label) << "}\n"
           << "\\centering\\includegraphics[width=0.68\\linewidth,"
           << "height=0.73\\textheight,keepaspectratio]{"
           << outputDirectory << "/" << stem << ".pdf}\n"
           << "\\vspace{0.5ex}\\par{\\tiny These curves include both the "
           << "stored previous JEC and the per-axis Run2024I factors used by "
           << "globalFit.C.}\n\\end{frame}\n";
  }

  writeControlAppendix(frames);
  frames.close();
  normalization.close();
  summary.close();
  std::cout << "Wrote three-input comparison plots, Beamer frames and TSV "
            << "diagnostics to " << outputDirectory
            << "; derived MPFnu used Rn=" << kSoftRecoilRn
            << " and Ru=" << kSoftRecoilRu << std::endl;
}
