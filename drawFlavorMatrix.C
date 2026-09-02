// Publication-style plots for the compact hybrid-tagger FlavorMatrix output.
//
// The macro deliberately draws score densities as binned heat maps.  It does
// not turn individual jets into TGraph points, which keeps vector output small
// and makes the figures deterministic in ROOT batch jobs.
#include "TBox.h"
#include "TCanvas.h"
#include "TColor.h"
#include "TFile.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLine.h"
#include "TList.h"
#include "TObjArray.h"
#include "TPolyLine.h"
#include "TPaveText.h"
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
#include <map>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace {

struct FlavorStyle {
  std::vector<int> ids;
  const char *name;
  const char *fileName;
  int color;
};

// ID 0 is the truth-unavailable bookkeeping bin (and the only data truth
// bin), not a physical flavor.  Keep it in the efficiency matrices, but omit
// it from unit-normalized MC shape/cube plots: signed sideband subtraction can
// make its tiny denominator nearly cancel and obscure every physical flavor.
const std::vector<FlavorStyle> kTruthFlavors = {
  {{1,3},"uds","uds",kBlue+1},
  {{4},"c","c",kGreen+2},
  {{5},"b","b",kRed+1},
  {{6},"g","g",kMagenta+1},
};

void configureFlavorStyle(const char *lumi) {
  setTDRStyle();
  writeExtraText = true;
  extraText = "Work in progress";
  extraText2 = "";
  lumi_136TeV = lumi;
  gStyle->SetOptStat(0);
  gStyle->SetNumberContours(50);
}

void labelFlavorAxis(TAxis *axis, bool reconstructed) {
  if (!axis) return;
  const char *reco[] = {
    "undefined", "uds", "unused", "unused", "c", "b", "g",
  };
  const char *truth[] = {
    "undefined", "d+u", "unused", "s", "c", "b", "g",
  };
  for (int bin=1; bin<=std::min(7,axis->GetNbins()); ++bin)
    axis->SetBinLabel(bin,reconstructed ? reco[bin-1] : truth[bin-1]);
}

template <typename T>
T *requireObject(TDirectory *directory, const char *path) {
  if (!directory)
    throw std::runtime_error("Cannot read " + std::string(path) +
                             " from a null ROOT directory");
  T *object = dynamic_cast<T*>(directory->Get(path));
  if (!object)
    throw std::runtime_error("Missing or wrongly typed ROOT object " +
                             std::string(path));
  return object;
}

template <typename T>
T *optionalObject(TDirectory *directory, const char *path) {
  return directory ? dynamic_cast<T*>(directory->Get(path)) : nullptr;
}

void saveBoth(TCanvas *canvas, const std::string &outputDirectory,
              const std::string &stem) {
  if (!canvas) return;
  canvas->Modified();
  canvas->Update();
  canvas->SaveAs((outputDirectory+"/"+stem+".pdf").c_str());
  canvas->SaveAs((outputDirectory+"/"+stem+".png").c_str());
}

double signedIntegral(const TH1 *histogram) {
  if (!histogram) return 0.;
  return histogram->Integral(1,histogram->GetNbinsX());
}

double signedIntegral(const TH2 *histogram) {
  if (!histogram) return 0.;
  return histogram->Integral(1,histogram->GetNbinsX(),
                             1,histogram->GetNbinsY());
}

std::unique_ptr<TH1D> scoreProjection(
  TH3D *cvbCvl, TH3D *cvbQvg, const std::string &score,
  const FlavorStyle &flavor) {
  const std::string name = "h_"+score+"_true"+flavor.fileName;
  std::unique_ptr<TH1D> combined;
  for (int id : flavor.ids) {
    const int zbin = id+1;
    TH1D *piece = nullptr;
    const std::string pieceName = name+"_id"+std::to_string(id);
    if (score=="CvB") {
      piece = cvbCvl->ProjectionX(
        pieceName.c_str(),1,cvbCvl->GetNbinsY(),zbin,zbin,"e");
    }
    else if (score=="CvL") {
      piece = cvbCvl->ProjectionY(
        pieceName.c_str(),1,cvbCvl->GetNbinsX(),zbin,zbin,"e");
    }
    else if (score=="QvG") {
      piece = cvbQvg->ProjectionY(
        pieceName.c_str(),1,cvbQvg->GetNbinsX(),zbin,zbin,"e");
    }
    else {
      throw std::runtime_error("Unknown hybrid-tagger score " + score);
    }
    if (!piece)
      throw std::runtime_error("Failed to project hybrid-tagger score " +
                               score);
    piece->SetDirectory(nullptr);
    if (!combined) {
      combined.reset(dynamic_cast<TH1D*>(piece->Clone(name.c_str())));
      combined->SetDirectory(nullptr);
    }
    else {
      combined->Add(piece);
    }
    delete piece;
  }
  if (!combined)
    throw std::runtime_error("No true-flavor bins selected for " + score);
  const double integral = signedIntegral(combined.get());
  if (std::fabs(integral)>1.e-12) combined->Scale(1./integral,"width");
  combined->SetLineColor(flavor.color);
  combined->SetMarkerColor(flavor.color);
  combined->SetLineWidth(3);
  return combined;
}

void drawScoreDistributions(TFile &mc, const std::string &outputDirectory) {
  TH3D *cvbCvl = requireObject<TH3D>(
    &mc,"FlavorMatrix/controls/h3_cvb_cvl_trueflavor");
  TH3D *cvbQvg = requireObject<TH3D>(
    &mc,"FlavorMatrix/controls/h3_cvb_qvg_trueflavor");
  for (const std::string score : {"CvB","CvL","QvG"}) {
    std::vector<std::unique_ptr<TH1D> > histograms;
    double maximum = 0.;
    double minimum = 0.;
    for (const FlavorStyle &flavor : kTruthFlavors) {
      histograms.push_back(scoreProjection(cvbCvl,cvbQvg,score,flavor));
      maximum = std::max(maximum,histograms.back()->GetMaximum());
      minimum = std::min(minimum,histograms.back()->GetMinimum());
    }

    configureFlavorStyle("Summer24 DY simulation");
    TH1D frame(("frame_"+score).c_str(),"",100,0.,1.);
    frame.SetDirectory(nullptr);
    frame.SetMinimum(minimum<0. ? 1.20*minimum : 0.);
    frame.SetMaximum(std::max(1.e-6,1.62*maximum));
    frame.GetXaxis()->SetTitle(
      (score=="QvG" ? "PNet QvG" : "UParT "+score).c_str());
    frame.GetYaxis()->SetTitle("1/N dN/d score");
    std::unique_ptr<TCanvas> canvas(tdrCanvas(
      ("c_"+score).c_str(),&frame,8,11,kSquare));
    for (const auto &histogram : histograms)
      histogram->Draw("HIST SAME");

    TLegend legend(0.54,0.76,0.90,0.89);
    legend.SetBorderSize(0);
    legend.SetFillStyle(0);
    legend.SetNColumns(2);
    legend.SetTextSize(0.038);
    for (size_t index=0; index<kTruthFlavors.size(); ++index)
      legend.AddEntry(histograms[index].get(),
                      (std::string("true ")+kTruthFlavors[index].name).c_str(),
                      "l");
    legend.Draw();
    const double cut = score=="QvG" ? 0.30 : 0.50;
    TLine cutLine(cut,frame.GetMinimum(),cut,maximum*1.05);
    cutLine.SetLineColor(kBlack);
    cutLine.SetLineStyle(kDashed);
    cutLine.SetLineWidth(2);
    cutLine.DrawClone();
    gPad->RedrawAxis();
    saveBoth(canvas.get(),outputDirectory,"tagger_"+score+"_trueflavor");
  }
}

struct PairSpec {
  const char *object;
  const char *stem;
  const char *xTitle;
  const char *yTitle;
  int cutKind;
};

std::unique_ptr<TH2D> pairProjection(
  TH3D *source, const PairSpec &spec, const FlavorStyle &flavor) {
  if (!source) return nullptr;
  // All stored pair controls use x=first discriminator and y=second
  // discriminator.  The requested figures put the second discriminator on
  // the horizontal axis and the first on the vertical axis.  Build this
  // transposition explicitly: ROOT's Project3D option ordering is easy to
  // misread and caused the previous titles and contents to be crossed.
  std::unique_ptr<TH2D> projection(new TH2D(
    (std::string("h2_")+spec.stem+"_true"+flavor.fileName).c_str(),"",
    source->GetNbinsY(),source->GetYaxis()->GetXmin(),
    source->GetYaxis()->GetXmax(),
    source->GetNbinsX(),source->GetXaxis()->GetXmin(),
    source->GetXaxis()->GetXmax()));
  projection->SetDirectory(nullptr);
  projection->Sumw2();
  for (int id : flavor.ids) {
    const int zbin = source->GetZaxis()->FindFixBin(id);
    for (int sourceX=1; sourceX<=source->GetNbinsX(); ++sourceX) {
      for (int sourceY=1; sourceY<=source->GetNbinsY(); ++sourceY) {
        const int sourceBin = source->GetBin(sourceX,sourceY,zbin);
        const int targetBin = projection->GetBin(sourceY,sourceX);
        projection->SetBinContent(
          targetBin,projection->GetBinContent(targetBin)+
                    source->GetBinContent(sourceBin));
        const double variance =
          std::pow(projection->GetBinError(targetBin),2)+
          std::pow(source->GetBinError(sourceBin),2);
        projection->SetBinError(targetBin,std::sqrt(variance));
      }
    }
  }
  projection->SetTitle("");
  projection->GetXaxis()->SetTitle(spec.xTitle);
  projection->GetYaxis()->SetTitle(spec.yTitle);
  projection->GetZaxis()->SetTitle("normalized density");
  for (int xbin=1; xbin<=projection->GetNbinsX(); ++xbin)
    for (int ybin=1; ybin<=projection->GetNbinsY(); ++ybin)
      if (!std::isfinite(projection->GetBinContent(xbin,ybin)) ||
          projection->GetBinContent(xbin,ybin)<0.)
        projection->SetBinContent(xbin,ybin,0.);
  projection->Rebin2D(2,2);
  projection->Smooth(1,"k5b");
  const double integral = signedIntegral(projection.get());
  if (std::fabs(integral)>1.e-12) projection->Scale(1./integral);
  return projection;
}

std::unique_ptr<TCanvas> heatmapCanvas(
  const std::string &name, TH2D *histogram, const char *lumi) {
  configureFlavorStyle(lumi);
  TH1D frame((name+"_axis_frame").c_str(),"",100,0.,1.);
  frame.SetDirectory(nullptr);
  frame.SetMinimum(0.);
  frame.SetMaximum(1.);
  std::unique_ptr<TCanvas> canvas(
    tdrCanvas(name.c_str(),&frame,8,11,kSquare));
  // Stock tdrCanvas reserves no room for a color palette.  Clear it and use
  // the same TDR geometry with a wider right margin before drawing CMS_lumi.
  canvas->Clear();
  canvas->SetLeftMargin(0.15);
  canvas->SetRightMargin(0.25);
  // Reserve a real header band for CMS, status and luminosity. Dense
  // heatmaps otherwise put the standard labels directly on populated cells.
  canvas->SetTopMargin(0.18);
  canvas->SetBottomMargin(0.13);
  histogram->UseCurrentStyle();
  histogram->GetXaxis()->SetTitleOffset(1.0);
  histogram->GetYaxis()->SetTitleOffset(1.25);
  histogram->GetZaxis()->SetTitleOffset(1.22);
  return canvas;
}

void drawHeatmapHeader(TCanvas *canvas, const char *qualifier="") {
  if (!canvas) return;
  canvas->cd();

  // The heat map replaces the temporary tdrCanvas frame and its header. Draw
  // one deterministic header in the explicitly reserved band instead of
  // calling CMS_lumi a second time (which can duplicate its first label).
  auto drawLabel = [](double x, double y, const TString &text, int font,
                      double size, int align) {
    TLatex *label = new TLatex(x,y,text);
    label->SetNDC();
    label->SetTextColor(kBlack);
    label->SetTextFont(font);
    label->SetTextSize(size);
    label->SetTextAlign(align);
    label->Draw();
  };
  const double labelLeft = std::min<double>(canvas->GetLeftMargin(),0.15);
  drawLabel(labelLeft,0.940,"CMS",61,0.050,11);
  const TString lumiText = lumi_136TeV + " (13.6 TeV)";
  drawLabel(1.-canvas->GetRightMargin(),0.940,lumiText,42,0.035,31);
  drawLabel(labelLeft,0.865,"Work in progress",52,0.032,11);
  if (qualifier && qualifier[0])
    drawLabel(1.-canvas->GetRightMargin(),0.865,qualifier,42,0.032,31);
}

void drawPairCuts(const PairSpec &spec) {
  TLine line;
  line.SetLineColor(kBlack);
  line.SetLineWidth(3);
  if (spec.cutKind==0) {
    line.DrawLine(0.,0.5,1.,0.5);
    line.DrawLine(0.5,0.5,0.5,1.);
    TLatex labels;
    labels.SetTextAlign(22);
    labels.SetTextSize(0.040);
    labels.DrawLatex(0.76,0.73,"c");
    labels.DrawLatex(0.27,0.73,"uds / g");
    labels.DrawLatex(0.50,0.23,"b");
  }
  else if (spec.cutKind==1) {
    line.DrawLine(0.,0.5,1.,0.5);
    line.SetLineStyle(kDashed);
    line.DrawLine(0.3,0.5,0.3,1.);
  }
  else {
    line.SetLineStyle(kDashed);
    line.DrawLine(0.,0.5,1.,0.5);
    line.DrawLine(0.3,0.,0.3,0.5);
  }
}

double highestDensityThreshold(const TH2D *histogram, double fraction) {
  if (!histogram || fraction<=0. || fraction>=1.)
    throw std::runtime_error("Invalid highest-density contour request");
  std::vector<double> cells;
  double total = 0.;
  for (int xbin=1; xbin<=histogram->GetNbinsX(); ++xbin) {
    for (int ybin=1; ybin<=histogram->GetNbinsY(); ++ybin) {
      const double value = histogram->GetBinContent(xbin,ybin);
      if (!std::isfinite(value) || value<=0.) continue;
      cells.push_back(value);
      total += value;
    }
  }
  if (cells.empty() || total<=0.) return std::numeric_limits<double>::infinity();
  std::sort(cells.begin(),cells.end(),std::greater<double>());
  double cumulative = 0.;
  for (double value : cells) {
    cumulative += value;
    if (cumulative>=fraction*total) return value;
  }
  return cells.back();
}

std::vector<std::unique_ptr<TGraph> > contourGraphs(
  const TH2D *source, double level, const std::string &name) {
  std::vector<std::unique_ptr<TGraph> > result;
  if (!source || !std::isfinite(level)) return result;
  std::unique_ptr<TH2D> work(dynamic_cast<TH2D*>(source->Clone(
    ("contour_work_"+name).c_str())));
  if (!work) return result;
  work->SetDirectory(nullptr);
  work->SetContour(1,&level);
  TCanvas scratch(("contour_canvas_"+name).c_str(),"",500,500);
  work->Draw("CONT LIST");
  scratch.Update();
  TObjArray *contours = dynamic_cast<TObjArray*>(
    gROOT->GetListOfSpecials()->FindObject("contours"));
  TList *levelGraphs = contours && contours->GetSize()>0
    ? dynamic_cast<TList*>(contours->At(0)) : nullptr;
  if (!levelGraphs) return result;
  TIter next(levelGraphs);
  TObject *object = nullptr;
  int index = 0;
  while ((object=next())) {
    TGraph *graph = dynamic_cast<TGraph*>(object);
    if (!graph || graph->GetN()<3) continue;
    std::unique_ptr<TGraph> clone(dynamic_cast<TGraph*>(graph->Clone(
      (name+"_"+std::to_string(index++)).c_str())));
    if (clone) result.push_back(std::move(clone));
  }
  return result;
}

struct FlavorContours {
  const FlavorStyle *style = nullptr;
  std::vector<std::unique_ptr<TGraph> > outer;
  std::vector<std::unique_ptr<TGraph> > core;
};

bool isClosedContour(const TGraph *graph) {
  if (!graph || graph->GetN()<3) return false;
  double firstX = 0.;
  double firstY = 0.;
  double lastX = 0.;
  double lastY = 0.;
  graph->GetPoint(0,firstX,firstY);
  graph->GetPoint(graph->GetN()-1,lastX,lastY);
  return std::hypot(firstX-lastX,firstY-lastY)<0.06;
}

void drawPairwiseDensities(TFile &mc, const std::string &outputDirectory) {
  const std::vector<PairSpec> specs = {
    {"FlavorMatrix/controls/h3_cvb_cvl_trueflavor","CvL_CvB",
     "UParT CvL","UParT CvB",0},
    {"FlavorMatrix/controls/h3_cvb_qvg_trueflavor","QvG_CvB",
     "PNet QvG","UParT CvB",1},
    {"FlavorMatrix/controls/h3_cvl_qvg_trueflavor","QvG_CvL",
     "PNet QvG","UParT CvL",2},
  };
  for (const PairSpec &spec : specs) {
    TH3D *source = requireObject<TH3D>(&mc,spec.object);
    std::vector<std::unique_ptr<TH2D> > projections;
    for (const FlavorStyle &flavor : kTruthFlavors) {
      projections.push_back(pairProjection(source,spec,flavor));
    }
    std::vector<FlavorContours> contours;
    for (size_t index=0; index<kTruthFlavors.size(); ++index) {
      FlavorContours entry;
      entry.style = &kTruthFlavors[index];
      const double outerLevel = highestDensityThreshold(
        projections[index].get(),0.90);
      const double coreLevel = highestDensityThreshold(
        projections[index].get(),0.50);
      entry.outer = contourGraphs(
        projections[index].get(),outerLevel,
        std::string(spec.stem)+"_"+entry.style->fileName+"_outer");
      entry.core = contourGraphs(
        projections[index].get(),coreLevel,
        std::string(spec.stem)+"_"+entry.style->fileName+"_core");
      contours.push_back(std::move(entry));
    }

    configureFlavorStyle("Summer24 DY simulation");
    TH1D frame((std::string("frame_")+spec.stem+"_all").c_str(),"",
               100,0.,1.);
    frame.SetDirectory(nullptr);
    frame.SetMinimum(0.);
    frame.SetMaximum(1.25);
    frame.GetXaxis()->SetTitle(spec.xTitle);
    frame.GetYaxis()->SetTitle(spec.yTitle);
    std::unique_ptr<TCanvas> canvas(tdrCanvas(
      (std::string("c_")+spec.stem+"_all").c_str(),&frame,8,11,kSquare));

    // Draw translucent 90% and 50% highest-density regions, followed by
    // their boundaries.  Negative signed sideband fluctuations are ignored
    // only when finding display contours; the stored histograms are intact.
    for (FlavorContours &entry : contours) {
      for (auto &graph : entry.outer) {
        if (!isClosedContour(graph.get())) continue;
        graph->SetFillColorAlpha(entry.style->color,0.08);
        graph->SetLineColor(entry.style->color);
        graph->Draw("F SAME");
      }
      for (auto &graph : entry.core) {
        if (!isClosedContour(graph.get())) continue;
        graph->SetFillColorAlpha(entry.style->color,0.16);
        graph->SetLineColor(entry.style->color);
        graph->Draw("F SAME");
      }
    }
    for (FlavorContours &entry : contours) {
      for (auto &graph : entry.outer) {
        graph->SetFillStyle(0);
        graph->SetLineColor(entry.style->color);
        graph->SetLineStyle(kDashed);
        graph->SetLineWidth(2);
        graph->Draw("L SAME");
      }
      for (auto &graph : entry.core) {
        graph->SetFillStyle(0);
        graph->SetLineColor(entry.style->color);
        graph->SetLineStyle(kSolid);
        graph->SetLineWidth(3);
        graph->Draw("L SAME");
      }
    }
    drawPairCuts(spec);

    TLegend flavorLegend(0.51,0.79,0.91,0.90);
    flavorLegend.SetBorderSize(0);
    flavorLegend.SetFillStyle(0);
    flavorLegend.SetNColumns(2);
    flavorLegend.SetTextSize(0.032);
    for (size_t index=0; index<contours.size(); ++index) {
      TObject *object = !contours[index].core.empty()
        ? static_cast<TObject*>(contours[index].core.front().get())
        : static_cast<TObject*>(projections[index].get());
      flavorLegend.AddEntry(object,
        (std::string("true ")+kTruthFlavors[index].name).c_str(),"l");
    }
    flavorLegend.Draw();
    TLine coreExample;
    coreExample.SetLineColor(kBlack);
    coreExample.SetLineWidth(3);
    TLine outerExample;
    outerExample.SetLineColor(kBlack);
    outerExample.SetLineStyle(kDashed);
    outerExample.SetLineWidth(2);
    TLegend densityLegend(0.51,0.70,0.91,0.79);
    densityLegend.SetBorderSize(0);
    densityLegend.SetFillStyle(0);
    densityLegend.SetTextSize(0.030);
    densityLegend.AddEntry(&coreExample,"50% highest density","l");
    densityLegend.AddEntry(&outerExample,"90% highest density","l");
    densityLegend.Draw();
    gPad->RedrawAxis();
    saveBoth(canvas.get(),outputDirectory,
      std::string("tagger2d_")+spec.stem+"_allflavors");
  }
}

struct CubeFaces {
  int bins = 0;
  std::vector<double> front; // CvL (horizontal) x CvB (vertical)
  std::vector<double> top;   // CvL (horizontal) x QvG (depth)
  std::vector<double> right; // CvB (vertical) x QvG (depth)
};

double &cell(std::vector<double> &values, int bins, int horizontal,
             int vertical) {
  return values[vertical*bins+horizontal];
}

double cell(const std::vector<double> &values, int bins, int horizontal,
            int vertical) {
  return values[vertical*bins+horizontal];
}

CubeFaces cubeFaces(TH3D *cube) {
  if (!cube) throw std::runtime_error("Cannot project a null UParT cube");
  const int bins = cube->GetNbinsX();
  if (cube->GetNbinsY()!=bins || cube->GetNbinsZ()!=bins)
    throw std::runtime_error("UParT cube must have equal score binning");
  CubeFaces faces;
  faces.bins = bins;
  faces.front.assign(bins*bins,0.);
  faces.top.assign(bins*bins,0.);
  faces.right.assign(bins*bins,0.);
  double total = 0.;
  for (int ib=1; ib<=bins; ++ib) {
    for (int il=1; il<=bins; ++il) {
      for (int iq=1; iq<=bins; ++iq) {
        const double value = cube->GetBinContent(ib,il,iq);
        total += value;
        cell(faces.front,bins,il-1,ib-1) += value;
        cell(faces.top,bins,il-1,iq-1) += value;
        cell(faces.right,bins,ib-1,iq-1) += value;
      }
    }
  }
  if (std::fabs(total)>1.e-12) {
    for (double &value : faces.front) value /= total;
    for (double &value : faces.top) value /= total;
    for (double &value : faces.right) value /= total;
  }
  return faces;
}

struct Point2D {
  double x;
  double y;
};

Point2D add(Point2D first, Point2D second) {
  return {first.x+second.x,first.y+second.y};
}

Point2D scale(Point2D value, double factor) {
  return {factor*value.x,factor*value.y};
}

void addFlavorCell(
  const Point2D &a, const Point2D &b, const Point2D &c, const Point2D &d,
  int color, double alpha,
  std::vector<std::unique_ptr<TPolyLine> > &polygons) {
  double x[] = {a.x,b.x,c.x,d.x,a.x};
  double y[] = {a.y,b.y,c.y,d.y,a.y};
  std::unique_ptr<TPolyLine> polygon(new TPolyLine(5,x,y));
  const int transparentColor = TColor::GetColorTransparent(
    color,std::max(0.,std::min(1.,alpha)));
  polygon->SetFillColor(transparentColor);
  polygon->SetFillStyle(1001);
  polygon->SetLineColor(transparentColor);
  polygon->SetLineWidth(0);
  polygon->Draw("f");
  polygons.push_back(std::move(polygon));
}

void drawSegment(Point2D start, Point2D end, int width=2,
                 int style=kSolid, int color=kBlack) {
  double x[] = {start.x,end.x};
  double y[] = {start.y,end.y};
  TPolyLine line(2,x,y);
  line.SetLineWidth(width);
  line.SetLineStyle(style);
  line.SetLineColor(color);
  line.DrawClone();
}

CubeFaces rebinCubeFaces(const CubeFaces &source, int factor) {
  if (factor<=0 || source.bins%factor!=0)
    throw std::runtime_error("Invalid hybrid-score cube display rebinning");
  CubeFaces result;
  result.bins = source.bins/factor;
  result.front.assign(result.bins*result.bins,0.);
  result.top.assign(result.bins*result.bins,0.);
  result.right.assign(result.bins*result.bins,0.);
  auto mergeFace = [factor,&source,&result](
      const std::vector<double> &input, std::vector<double> &output) {
    for (int y=0; y<result.bins; ++y) {
      for (int x=0; x<result.bins; ++x) {
        double sum = 0.;
        for (int dy=0; dy<factor; ++dy)
          for (int dx=0; dx<factor; ++dx)
            sum += cell(input,source.bins,factor*x+dx,factor*y+dy);
        cell(output,result.bins,x,y) = sum;
      }
    }
  };
  mergeFace(source.front,result.front);
  mergeFace(source.top,result.top);
  mergeFace(source.right,result.right);
  return result;
}

std::unique_ptr<TH3D> combinedFlavorCube(TFile &mc,
                                         const FlavorStyle &flavor) {
  std::unique_ptr<TH3D> result;
  for (int id : flavor.ids) {
    TH3D *source = requireObject<TH3D>(&mc,
      (std::string("FlavorMatrix/controls/h3_cvb_cvl_qvg_true")+
       std::to_string(id)).c_str());
    if (!result) {
      result.reset(dynamic_cast<TH3D*>(source->Clone(
        (std::string("cube_true")+flavor.fileName+"_combined").c_str())));
      if (result) result->SetDirectory(nullptr);
    }
    else result->Add(source);
  }
  if (!result)
    throw std::runtime_error("Failed to combine a true-flavor score cube");
  return result;
}

void drawFlavorCube(const std::vector<CubeFaces> &faces, double maximum,
                    const std::string &outputDirectory) {
  configureFlavorStyle("Summer24 DY simulation");
  TH1D frame("cube_frame_allflavors","",100,0.,1.12);
  frame.SetDirectory(nullptr);
  frame.SetMinimum(0.);
  frame.SetMaximum(1.02);
  frame.GetXaxis()->SetLabelSize(0.);
  frame.GetYaxis()->SetLabelSize(0.);
  frame.GetXaxis()->SetTickLength(0.);
  frame.GetYaxis()->SetTickLength(0.);
  std::unique_ptr<TCanvas> canvas(tdrCanvas(
    "c_cube_allflavors",&frame,8,11,kSquare));

  if (faces.size()!=kTruthFlavors.size())
    throw std::runtime_error("Flavor cube/style size mismatch");
  const Point2D origin = {0.12,0.16};
  const Point2D xAxis = {0.47,0.};
  const Point2D yAxis = {0.,0.50};
  const Point2D qAxis = {0.12,0.10};
  const int bins = faces.empty() ? 0 : faces.front().bins;
  std::vector<std::unique_ptr<TPolyLine> > polygons;
  polygons.reserve(3*bins*bins*faces.size());

  auto alphaFor = [maximum](double value) {
    if (!std::isfinite(value) || value<=0. || maximum<=0.) return 0.;
    const double fraction = std::max(0.,std::min(1.,value/maximum));
    if (fraction<0.003) return 0.;
    return std::min(0.82,0.04+0.78*std::sqrt(fraction));
  };
  auto drawCellForFlavors = [&](const Point2D &a, const Point2D &b,
                                const Point2D &c, const Point2D &d,
                                const std::vector<double> &values) {
    std::vector<size_t> order(values.size());
    for (size_t index=0; index<order.size(); ++index) order[index] = index;
    std::sort(order.begin(),order.end(),[&](size_t first, size_t second) {
      return values[first]<values[second];
    });
    for (size_t index : order) {
      const double alpha = alphaFor(values[index]);
      if (alpha<=0.) continue;
      addFlavorCell(a,b,c,d,kTruthFlavors[index].color,alpha,polygons);
    }
  };

  for (int iy=0; iy<bins; ++iy) {
    for (int ix=0; ix<bins; ++ix) {
      const double x0 = static_cast<double>(ix)/bins;
      const double x1 = static_cast<double>(ix+1)/bins;
      const double y0 = static_cast<double>(iy)/bins;
      const double y1 = static_cast<double>(iy+1)/bins;
      std::vector<double> values;
      for (const CubeFaces &flavorFaces : faces)
        values.push_back(cell(flavorFaces.front,bins,ix,iy));
      drawCellForFlavors(
        add(origin,add(scale(xAxis,x0),scale(yAxis,y0))),
        add(origin,add(scale(xAxis,x1),scale(yAxis,y0))),
        add(origin,add(scale(xAxis,x1),scale(yAxis,y1))),
        add(origin,add(scale(xAxis,x0),scale(yAxis,y1))),values);
    }
  }
  for (int iq=0; iq<bins; ++iq) {
    for (int ix=0; ix<bins; ++ix) {
      const double x0 = static_cast<double>(ix)/bins;
      const double x1 = static_cast<double>(ix+1)/bins;
      const double q0 = static_cast<double>(iq)/bins;
      const double q1 = static_cast<double>(iq+1)/bins;
      const Point2D base = add(origin,yAxis);
      std::vector<double> values;
      for (const CubeFaces &flavorFaces : faces)
        values.push_back(cell(flavorFaces.top,bins,ix,iq));
      drawCellForFlavors(
        add(base,add(scale(xAxis,x0),scale(qAxis,q0))),
        add(base,add(scale(xAxis,x1),scale(qAxis,q0))),
        add(base,add(scale(xAxis,x1),scale(qAxis,q1))),
        add(base,add(scale(xAxis,x0),scale(qAxis,q1))),values);
    }
  }
  for (int iq=0; iq<bins; ++iq) {
    for (int iy=0; iy<bins; ++iy) {
      const double y0 = static_cast<double>(iy)/bins;
      const double y1 = static_cast<double>(iy+1)/bins;
      const double q0 = static_cast<double>(iq)/bins;
      const double q1 = static_cast<double>(iq+1)/bins;
      const Point2D base = add(origin,xAxis);
      std::vector<double> values;
      for (const CubeFaces &flavorFaces : faces)
        values.push_back(cell(flavorFaces.right,bins,iy,iq));
      drawCellForFlavors(
        add(base,add(scale(yAxis,y0),scale(qAxis,q0))),
        add(base,add(scale(yAxis,y1),scale(qAxis,q0))),
        add(base,add(scale(yAxis,y1),scale(qAxis,q1))),
        add(base,add(scale(yAxis,y0),scale(qAxis,q1))),values);
    }
  }

  // Unique visible cube edges.  In particular, do not draw the old Y--XYQ
  // diagonal across the top face.
  const Point2D x = add(origin,xAxis);
  const Point2D y = add(origin,yAxis);
  const Point2D xy = add(x,yAxis);
  const Point2D xq = add(x,qAxis);
  const Point2D yq = add(y,qAxis);
  const Point2D xyq = add(xy,qAxis);
  drawSegment(origin,x,3);
  drawSegment(origin,y,3);
  drawSegment(x,xy,3);
  drawSegment(y,xy,3);
  drawSegment(y,yq,3);
  drawSegment(yq,xyq,3);
  drawSegment(xy,xyq,3);
  drawSegment(x,xq,3);
  drawSegment(xq,xyq,3);

  drawSegment(add(origin,scale(yAxis,0.5)),
              add(x,scale(yAxis,0.5)),3);
  drawSegment(add(x,scale(yAxis,0.5)),
              add(add(x,scale(yAxis,0.5)),qAxis),3);
  drawSegment(add(add(origin,scale(xAxis,0.5)),scale(yAxis,0.5)),
              add(add(origin,scale(xAxis,0.5)),yAxis),3);
  drawSegment(add(y,scale(xAxis,0.5)),
              add(add(y,scale(xAxis,0.5)),qAxis),3);
  drawSegment(add(y,scale(qAxis,0.3)),
              add(add(y,scale(xAxis,0.5)),scale(qAxis,0.3)),3);

  TLatex text;
  text.SetTextAlign(22);
  text.SetTextSize(0.034);
  text.DrawLatex(0.36,0.11,"UParT CvL");
  text.SetTextAngle(90.);
  text.DrawLatex(0.070,0.41,"UParT CvB");
  text.SetTextAngle(34.);
  text.DrawLatex(0.70,0.82,"PNet QvG");
  text.SetTextAngle(0.);
  text.SetTextSize(0.040);
  text.DrawLatex(0.36,0.29,"b");
  text.DrawLatex(0.48,0.54,"c");
  text.DrawLatex(0.25,0.54,"uds / g");
  text.SetTextSize(0.032);
  text.DrawLatex(0.27,0.70,"g");
  text.DrawLatex(0.33,0.75,"uds");

  std::vector<std::unique_ptr<TBox> > legendBoxes;
  TLegend legend(0.75,0.18,0.94,0.39);
  legend.SetBorderSize(0);
  legend.SetFillStyle(0);
  legend.SetTextSize(0.030);
  for (const FlavorStyle &flavor : kTruthFlavors) {
    std::unique_ptr<TBox> box(new TBox());
    box->SetFillColor(TColor::GetColorTransparent(flavor.color,0.58));
    box->SetLineColor(flavor.color);
    legend.AddEntry(box.get(),
      (std::string("true ")+flavor.name).c_str(),"f");
    legendBoxes.push_back(std::move(box));
  }
  legend.Draw();
  TLatex annotation;
  annotation.SetNDC();
  annotation.SetTextSize(0.027);
  annotation.DrawLatex(0.75,0.46,"darker = denser");
  annotation.DrawLatex(0.75,0.42,"unit-area shapes");

  gPad->RedrawAxis();
  saveBoth(canvas.get(),outputDirectory,"tagger_cube_allflavors");
}

void drawFlavorCubes(TFile &mc, const std::string &outputDirectory) {
  std::vector<CubeFaces> allFaces;
  double maximum = 0.;
  for (const FlavorStyle &flavor : kTruthFlavors) {
    std::unique_ptr<TH3D> cube = combinedFlavorCube(mc,flavor);
    CubeFaces faces = rebinCubeFaces(cubeFaces(cube.get()),2);
    for (const std::vector<double> *face : {
           &faces.front,&faces.top,&faces.right})
      for (double value : *face) maximum = std::max(maximum,value);
    allFaces.push_back(std::move(faces));
  }
  drawFlavorCube(allFaces,std::max(maximum,1.e-12),outputDirectory);
}

double taggedCubeCount(const TH3D *cube, int tag, double qvgThreshold) {
  if (!cube) return 0.;
  double result = 0.;
  for (int cvbBin=1; cvbBin<=cube->GetNbinsX(); ++cvbBin) {
    const double cvb = cube->GetXaxis()->GetBinCenter(cvbBin);
    for (int cvlBin=1; cvlBin<=cube->GetNbinsY(); ++cvlBin) {
      const double cvl = cube->GetYaxis()->GetBinCenter(cvlBin);
      for (int qvgBin=1; qvgBin<=cube->GetNbinsZ(); ++qvgBin) {
        const double qvg = cube->GetZaxis()->GetBinCenter(qvgBin);
        int reconstructed = 0;
        if (cvb<0.5) reconstructed = 5;
        else if (cvl>=0.5) reconstructed = 4;
        else reconstructed = (qvg>=qvgThreshold ? 1 : 6);
        if (reconstructed==tag)
          result += cube->GetBinContent(cvbBin,cvlBin,qvgBin);
      }
    }
  }
  return result;
}

void drawQvgWorkingPointScan(TFile &mc,
                             const std::string &outputDirectory) {
  const std::vector<FlavorStyle> scanFlavors = {
    {{1,3},"uds","uds_scan",kBlue+1},
    {{4},"c","c_scan",kGreen+2},
    {{5},"b","b_scan",kRed+1},
    {{0,6},"g","g_scan",kMagenta+1},
  };
  std::vector<std::unique_ptr<TH3D> > cubes;
  std::vector<double> totals;
  for (const FlavorStyle &flavor : scanFlavors) {
    cubes.push_back(combinedFlavorCube(mc,flavor));
    totals.push_back(cubes.back()->Integral());
  }
  TGraph gPurity;
  TGraph gEfficiency;
  TGraph udsPurity;
  TGraph udsEfficiency;
  for (int step=0; step<=50; ++step) {
    const double threshold = 0.02*step;
    double gTaggedTotal = 0.;
    double udsTaggedTotal = 0.;
    for (size_t flavor=0; flavor<cubes.size(); ++flavor) {
      gTaggedTotal += taggedCubeCount(cubes[flavor].get(),6,threshold);
      udsTaggedTotal += taggedCubeCount(cubes[flavor].get(),1,threshold);
    }
    const double gTagged = taggedCubeCount(cubes[3].get(),6,threshold);
    const double udsTagged = taggedCubeCount(cubes[0].get(),1,threshold);
    gPurity.SetPoint(step,threshold,
                     gTaggedTotal>0. ? gTagged/gTaggedTotal : 0.);
    gEfficiency.SetPoint(step,threshold,
                         totals[3]>0. ? gTagged/totals[3] : 0.);
    udsPurity.SetPoint(step,threshold,
                       udsTaggedTotal>0. ? udsTagged/udsTaggedTotal : 0.);
    udsEfficiency.SetPoint(step,threshold,
                           totals[0]>0. ? udsTagged/totals[0] : 0.);
  }
  std::ofstream table(
    (outputDirectory+"/qvg_working_point_scan.tsv").c_str());
  table << std::setprecision(10)
        << "threshold\tg_purity\tg_efficiency\tuds_purity"
           "\tuds_efficiency\n";
  for (int point=0; point<gPurity.GetN(); ++point) {
    double threshold = 0.;
    double gPurityValue = 0.;
    double ignored = 0.;
    double gEfficiencyValue = 0.;
    double udsPurityValue = 0.;
    double udsEfficiencyValue = 0.;
    gPurity.GetPoint(point,threshold,gPurityValue);
    gEfficiency.GetPoint(point,ignored,gEfficiencyValue);
    udsPurity.GetPoint(point,ignored,udsPurityValue);
    udsEfficiency.GetPoint(point,ignored,udsEfficiencyValue);
    table << threshold << '\t' << gPurityValue << '\t'
          << gEfficiencyValue << '\t' << udsPurityValue << '\t'
          << udsEfficiencyValue << '\n';
  }
  configureFlavorStyle("Summer24 DY simulation");
  TH1D frame("frame_qvg_working_point_scan","",100,0.,1.);
  frame.SetDirectory(nullptr);
  frame.SetMinimum(0.);
  frame.SetMaximum(1.25);
  frame.GetXaxis()->SetTitle("PNet QvG threshold");
  frame.GetYaxis()->SetTitle("efficiency or purity");
  std::unique_ptr<TCanvas> canvas(tdrCanvas(
    "c_qvg_working_point_scan",&frame,8,11,kSquare));
  for (auto graphAndStyle : {
         std::make_pair(&gPurity,std::make_pair(kMagenta+1,kSolid)),
         std::make_pair(&gEfficiency,std::make_pair(kMagenta+1,kDashed)),
         std::make_pair(&udsPurity,std::make_pair(kBlue+1,kSolid)),
         std::make_pair(&udsEfficiency,std::make_pair(kBlue+1,kDashed))}) {
    graphAndStyle.first->SetLineColor(graphAndStyle.second.first);
    graphAndStyle.first->SetLineStyle(graphAndStyle.second.second);
    graphAndStyle.first->SetLineWidth(3);
    graphAndStyle.first->Draw("L SAME");
  }
  TLine chosen(0.3,0.,0.3,1.0);
  chosen.SetLineColor(kBlack);
  chosen.SetLineWidth(3);
  chosen.DrawClone();
  TLine old(0.5,0.,0.5,1.0);
  old.SetLineColor(kGray+2);
  old.SetLineStyle(kDotted);
  old.SetLineWidth(2);
  old.DrawClone();
  TLegend legend(0.48,0.74,0.89,0.89);
  legend.SetBorderSize(0);
  legend.SetFillStyle(0);
  legend.SetNColumns(2);
  legend.AddEntry(&gPurity,"g purity","l");
  legend.AddEntry(&gEfficiency,"g efficiency","l");
  legend.AddEntry(&udsPurity,"uds purity","l");
  legend.AddEntry(&udsEfficiency,"uds efficiency","l");
  legend.AddEntry(&chosen,"chosen 0.3","l");
  legend.AddEntry(&old,"previous 0.5","l");
  legend.Draw();
  gPad->RedrawAxis();
  saveBoth(canvas.get(),outputDirectory,"qvg_working_point_scan");
}

std::unique_ptr<TH2D> compactFlavorMatrix(const TH2D *source,
                                          const std::string &name) {
  if (!source) return nullptr;
  const std::vector<int> recoIds = {1,4,5,6};
  const std::vector<int> truthIds = {1,4,5,6};
  const std::vector<std::string> recoLabels = {"uds","c","b","g"};
  const std::vector<std::string> truthLabels = {"uds","c","b","g"};
  std::unique_ptr<TH2D> result(new TH2D(
    name.c_str(),"",recoIds.size(),-0.5,recoIds.size()-0.5,
    truthIds.size(),-0.5,truthIds.size()-0.5));
  result->SetDirectory(nullptr);
  for (size_t reco=0; reco<recoIds.size(); ++reco) {
    result->GetXaxis()->SetBinLabel(reco+1,recoLabels[reco].c_str());
    for (size_t truth=0; truth<truthIds.size(); ++truth) {
      if (reco==0)
        result->GetYaxis()->SetBinLabel(truth+1,truthLabels[truth].c_str());
      const int sourceX = source->GetXaxis()->FindFixBin(recoIds[reco]);
      const int sourceY = source->GetYaxis()->FindFixBin(truthIds[truth]);
      result->SetBinContent(reco+1,truth+1,
                            source->GetBinContent(sourceX,sourceY));
      result->SetBinError(reco+1,truth+1,
                          source->GetBinError(sourceX,sourceY));
    }
  }
  return result;
}

void setWhiteSequentialPalette() {
  const int colors = 255;
  double stops[] = {0.00,0.03,0.28,0.62,1.00};
  double red[] =   {1.00,0.98,0.63,0.18,0.02};
  double green[] = {1.00,0.99,0.82,0.43,0.09};
  double blue[] =  {1.00,1.00,0.90,0.66,0.20};
  const int first = TColor::CreateGradientColorTable(
    5,stops,red,green,blue,colors);
  std::vector<int> palette(colors);
  for (int index=0; index<colors; ++index) palette[index] = first+index;
  gStyle->SetPalette(colors,palette.data());
}

void setUnityDivergingPalette() {
  const int colors = 255;
  double stops[] = {0.00,0.48,0.50,0.52,1.00};
  double red[] =   {0.05,0.86,1.00,0.98,0.55};
  double green[] = {0.18,0.94,1.00,0.88,0.03};
  double blue[] =  {0.55,0.99,1.00,0.82,0.04};
  const int first = TColor::CreateGradientColorTable(
    5,stops,red,green,blue,colors);
  std::vector<int> palette(colors);
  for (int index=0; index<colors; ++index) palette[index] = first+index;
  gStyle->SetPalette(colors,palette.data());
}

void drawMatrixText(const TH2D *histogram, double minimum, double maximum,
                    bool diverging) {
  if (!histogram || maximum<=minimum) return;
  TLatex text;
  text.SetTextAlign(22);
  text.SetTextFont(42);
  text.SetTextSize(0.031);
  for (int xbin=1; xbin<=histogram->GetNbinsX(); ++xbin) {
    for (int ybin=1; ybin<=histogram->GetNbinsY(); ++ybin) {
      const double value = histogram->GetBinContent(xbin,ybin);
      if (!std::isfinite(value) || (diverging && value==0.)) continue;
      const double darkness = diverging
        ? std::min(1.,std::fabs(value-1.)/
                        std::max(1.e-12,0.5*(maximum-minimum)))
        : std::max(0.,std::min(1.,(value-minimum)/(maximum-minimum)));
      text.SetTextColor(darkness>0.58 ? kWhite : kBlack);
      text.DrawLatex(histogram->GetXaxis()->GetBinCenter(xbin),
                     histogram->GetYaxis()->GetBinCenter(ybin),
                     Form("%.3f",value));
    }
  }
}

void drawAnalysisMatrix(
  TFile &analysis, const char *path, const std::string &stem,
  const char *zTitle, double minimum, double maximum,
  const std::string &outputDirectory, bool diverging=false) {
  TH2D *source = optionalObject<TH2D>(&analysis,path);
  if (!source) {
    std::cout << "Optional analysis object " << path
              << " is unavailable; skipping " << stem << std::endl;
    return;
  }
  std::unique_ptr<TH2D> histogram = compactFlavorMatrix(
    source,"plot_"+stem);
  histogram->SetTitle("");
  histogram->GetXaxis()->SetTitle("Reco hybrid flavor");
  histogram->GetYaxis()->SetTitle("Generator parton flavor");
  histogram->GetZaxis()->SetTitle(zTitle);
  histogram->SetMinimum(minimum);
  histogram->SetMaximum(maximum);
  configureFlavorStyle("2024I + Summer24");
  if (diverging) setUnityDivergingPalette();
  else setWhiteSequentialPalette();
  std::unique_ptr<TCanvas> canvas = heatmapCanvas(
    "c_"+stem,histogram.get(),"2024I + Summer24");
  if (diverging) setUnityDivergingPalette();
  else setWhiteSequentialPalette();
  histogram->Draw("COLZ");
  drawMatrixText(histogram.get(),minimum,maximum,diverging);
  drawHeatmapHeader(canvas.get());
  gPad->RedrawAxis();
  saveBoth(canvas.get(),outputDirectory,stem);
}

void drawDirectAnalysisMatrix(
  TFile &analysis, const char *path, const std::string &stem,
  double minimum, double maximum, const std::string &outputDirectory,
  bool percentage=false) {
  TH2D *source = optionalObject<TH2D>(&analysis,path);
  if (!source) return;
  std::unique_ptr<TH2D> histogram(dynamic_cast<TH2D*>(source->Clone(
    ("plot_"+stem).c_str())));
  if (!histogram) return;
  histogram->SetDirectory(nullptr);
  histogram->SetTitle("");
  if (percentage) histogram->Scale(100.);
  if (!(maximum>minimum))
    maximum = std::max(minimum+1.e-6,1.08*histogram->GetMaximum());
  histogram->SetMinimum(minimum);
  histogram->SetMaximum(maximum);
  configureFlavorStyle("2024I + Summer24");
  setWhiteSequentialPalette();
  std::unique_ptr<TCanvas> canvas = heatmapCanvas(
    "c_"+stem,histogram.get(),"2024I + Summer24");
  canvas->SetLeftMargin(0.24);
  // The absolute scale is stated in the slide/plot caption.  With the large
  // TDR fonts a rotated palette title obscures both the palette ticks and the
  // rightmost matrix cells, so keep the palette itself deliberately compact.
  histogram->GetYaxis()->SetTitle("");
  histogram->GetZaxis()->SetTitle("");
  setWhiteSequentialPalette();
  histogram->Draw("COLZ");
  drawMatrixText(histogram.get(),minimum,maximum,false);
  drawHeatmapHeader(canvas.get());
  gPad->RedrawAxis();
  saveBoth(canvas.get(),outputDirectory,stem);
}

std::unique_ptr<TGraphErrors> graphFromHistogram(TH1 *histogram) {
  if (!histogram) return nullptr;
  std::unique_ptr<TGraphErrors> graph(new TGraphErrors());
  graph->SetName("g_response_residual_plot");
  const std::vector<int> flavorIds = {1,4,5,6};
  for (size_t index=0; index<flavorIds.size(); ++index) {
    const int bin = histogram->GetXaxis()->FindFixBin(flavorIds[index]);
    const double value = histogram->GetBinContent(bin);
    const double error = histogram->GetBinError(bin);
    if (!std::isfinite(value) || (value==0. && error==0.)) continue;
    const int point = graph->GetN();
    graph->SetPoint(point,index,value);
    graph->SetPointError(point,0.,error);
  }
  return graph;
}

struct PtGraphStyle {
  const char *name;
  const char *label;
  int color;
  int marker;
};

const PtGraphStyle kPtGraphStyles[] = {
  {"uds","uds",kBlue+1,kFullCircle},
  {"c","c",kGreen+2,kFullSquare},
  {"b","b",kRed+1,kFullTriangleUp},
  {"g","g",kMagenta+1,kFullDiamond},
};

std::unique_ptr<TGraphErrors> clonePtGraph(
  TFile &analysis, const std::string &path, const std::string &name,
  const PtGraphStyle &style, bool openMarker=false) {
  TGraphErrors *source = optionalObject<TGraphErrors>(&analysis,path.c_str());
  if (!source) return nullptr;
  std::unique_ptr<TGraphErrors> graph(dynamic_cast<TGraphErrors*>(
    source->Clone(name.c_str())));
  if (!graph) return nullptr;
  graph->SetLineColor(style.color);
  graph->SetMarkerColor(style.color);
  graph->SetMarkerStyle(openMarker ? style.marker+4 : style.marker);
  graph->SetLineWidth(2);
  return graph;
}

void graphBounds(const std::vector<std::unique_ptr<TGraphErrors> > &graphs,
                 double &xmin, double &xmax, double &ymin, double &ymax) {
  xmin = std::numeric_limits<double>::infinity();
  xmax = 0.;
  ymin = std::numeric_limits<double>::infinity();
  ymax = -std::numeric_limits<double>::infinity();
  for (const auto &graph : graphs) {
    if (!graph) continue;
    for (int point=0; point<graph->GetN(); ++point) {
      double x = 0.;
      double y = 0.;
      graph->GetPoint(point,x,y);
      if (!std::isfinite(x) || !std::isfinite(y)) continue;
      xmin = std::min(xmin,std::max(1.e-3,x-graph->GetErrorX(point)));
      xmax = std::max(xmax,x+graph->GetErrorX(point));
      ymin = std::min(ymin,y-graph->GetErrorY(point));
      ymax = std::max(ymax,y+graph->GetErrorY(point));
    }
  }
}

void drawTagResponseRatio(TFile &analysis, bool compositionCorrected,
                          const std::string &outputDirectory) {
  std::vector<std::unique_ptr<TGraphErrors> > graphs;
  for (const PtGraphStyle &style : kPtGraphStyles) {
    const std::string qualifier = compositionCorrected
      ? "composition_corrected" : "raw";
    graphs.push_back(clonePtGraph(
      analysis,"response/g_tag_response_"+qualifier+
        "_data_over_mc_vs_pt_"+style.name,
      "plot_tag_response_"+qualifier+"_"+style.name,style));
  }
  double xmin, xmax, ymin, ymax;
  graphBounds(graphs,xmin,xmax,ymin,ymax);
  if (!std::isfinite(xmin) || !(xmax>xmin)) return;
  xmin = std::max(10.,xmin);
  xmax = std::min(1500.,xmax);
  double extent = std::max(0.03,std::max(std::fabs(ymin-1.),
                                        std::fabs(ymax-1.)));
  extent = std::min(0.35,1.12*extent);
  configureFlavorStyle("Run2024I + Summer24 DY");
  const std::string qualifier = compositionCorrected
    ? "composition_corrected" : "raw";
  TH1D frame(("frame_tag_response_"+qualifier).c_str(),"",100,xmin,xmax);
  frame.SetDirectory(nullptr);
  frame.SetMinimum(1.-extent);
  frame.SetMaximum(1.+1.45*extent);
  frame.GetXaxis()->SetTitle("p_{T,Z} (GeV)");
  frame.GetYaxis()->SetTitle("tagged data / MC HDM response");
  frame.GetXaxis()->SetMoreLogLabels();
  frame.GetXaxis()->SetNoExponent();
  std::unique_ptr<TCanvas> canvas(tdrCanvas(
    ("c_tag_response_"+qualifier).c_str(),&frame,8,11,kSquare));
  canvas->SetLogx();
  TLine unity(xmin,1.,xmax,1.);
  unity.SetLineColor(kGray+2);
  unity.SetLineStyle(kDashed);
  unity.DrawClone();
  for (auto &graph : graphs) if (graph) graph->Draw("PZ SAME");
  TLegend legend(0.50,0.76,0.88,0.88);
  legend.SetBorderSize(0);
  legend.SetFillStyle(0);
  legend.SetNColumns(2);
  for (size_t index=0; index<graphs.size(); ++index)
    if (graphs[index]) legend.AddEntry(
      graphs[index].get(),
      (std::string("reco ")+kPtGraphStyles[index].label).c_str(),"pl");
  legend.Draw();
  gPad->RedrawAxis();
  saveBoth(canvas.get(),outputDirectory,
    "tag_response_"+qualifier+"_data_over_mc_vs_pt");
}

void drawRecoilFraction(TFile &analysis, const std::string &component,
                        bool ratio, const std::string &outputDirectory) {
  std::vector<std::unique_ptr<TGraphErrors> > graphs;
  std::vector<std::string> labels;
  for (const PtGraphStyle &style : kPtGraphStyles) {
    if (ratio) {
      graphs.push_back(clonePtGraph(
        analysis,"response/g_"+component+
          "_fraction_ratio_data_over_mc_vs_pt_"+style.name,
        "plot_"+component+"_ratio_"+style.name,style));
      labels.push_back(style.label);
    }
    else {
      graphs.push_back(clonePtGraph(
        analysis,"response/g_"+component+"_fraction_data_vs_pt_"+
          style.name,"plot_"+component+"_data_"+style.name,style));
      labels.push_back(std::string("data ")+style.label);
      graphs.push_back(clonePtGraph(
        analysis,"response/g_"+component+"_fraction_mc_vs_pt_"+
          style.name,"plot_"+component+"_mc_"+style.name,style,true));
      if (graphs.back()) {
        graphs.back()->SetLineStyle(kDashed);
        graphs.back()->SetMarkerStyle(0);
      }
      labels.push_back(std::string("MC ")+style.label);
    }
  }
  double xmin, xmax, ymin, ymax;
  graphBounds(graphs,xmin,xmax,ymin,ymax);
  if (!std::isfinite(xmin) || !(xmax>xmin)) return;
  xmin = std::max(10.,xmin);
  xmax = std::min(1500.,xmax);
  configureFlavorStyle("Run2024I + Summer24 DY");
  TH1D frame(("frame_"+component+(ratio ? "_ratio" : "_absolute")).c_str(),
             "",100,xmin,xmax);
  frame.SetDirectory(nullptr);
  if (ratio) {
    double extent = std::max(0.12,std::max(std::fabs(ymin-1.),
                                          std::fabs(ymax-1.)));
    extent = std::min(0.65,1.10*extent);
    frame.SetMinimum(1.-extent);
    frame.SetMaximum(1.+1.35*extent);
    frame.GetYaxis()->SetTitle("inferred data / MC fraction");
  }
  else {
    const double span = std::max(0.02,ymax-ymin);
    frame.SetMinimum(std::min(0.,ymin-0.12*span));
    frame.SetMaximum(ymax+0.38*span);
    frame.GetYaxis()->SetTitle(component=="fsr"
      ? "effective extra-recoil fraction f_{n}"
      : "effective unclustered fraction f_{u}");
  }
  frame.GetXaxis()->SetTitle("p_{T,Z} (GeV)");
  frame.GetXaxis()->SetMoreLogLabels();
  frame.GetXaxis()->SetNoExponent();
  std::unique_ptr<TCanvas> canvas(tdrCanvas(
    ("c_"+component+(ratio ? "_ratio" : "_absolute")).c_str(),
    &frame,8,11,kSquare));
  canvas->SetLogx();
  if (ratio) {
    TLine unity(xmin,1.,xmax,1.);
    unity.SetLineColor(kGray+2);
    unity.SetLineStyle(kDashed);
    unity.DrawClone();
  }
  for (auto &graph : graphs)
    if (graph) graph->Draw(ratio ? "PZ SAME" : "PZL SAME");
  TLegend legend(0.49,0.72,0.89,0.88);
  legend.SetBorderSize(0);
  legend.SetFillStyle(0);
  legend.SetNColumns(2);
  for (size_t index=0; index<graphs.size(); ++index)
    if (graphs[index]) legend.AddEntry(
      graphs[index].get(),labels[index].c_str(),ratio ? "pl" : "pl");
  legend.Draw();
  gPad->RedrawAxis();
  saveBoth(canvas.get(),outputDirectory,component+"_fraction_"+
    (ratio ? "data_over_mc_vs_pt" : "data_mc_vs_pt"));
}

void drawResponseBinningComparison(TFile &analysis,
                                   const PtGraphStyle &flavor,
                                   const std::string &outputDirectory) {
  struct BinningStyle {
    const char *name;
    const char *label;
    int color;
    int marker;
  };
  const BinningStyle binnings[] = {
    {"tc","p_{T,Z}",kBlack,kFullCircle},
    {"ad","p_{T,ave}",kBlue+1,kFullSquare},
    {"ab","p_{T,ave} projection",kGreen+2,kFullTriangleUp},
    {"pf","p_{T,jet}",kRed+1,kFullDiamond},
  };
  std::vector<std::unique_ptr<TGraphErrors> > graphs;
  for (const BinningStyle &binning : binnings) {
    PtGraphStyle style = {
      flavor.name,flavor.label,binning.color,binning.marker};
    graphs.push_back(clonePtGraph(
      analysis,"response/g_response_residual_vs_pt_"+
        std::string(binning.name)+"_"+flavor.name,
      "plot_binning_"+std::string(binning.name)+"_"+flavor.name,style));
  }
  double xmin, xmax, ymin, ymax;
  graphBounds(graphs,xmin,xmax,ymin,ymax);
  if (!std::isfinite(xmin) || !(xmax>xmin)) return;
  xmin = std::max(10.,xmin);
  xmax = std::min(1500.,xmax);
  double extent = std::max(0.03,std::max(std::fabs(ymin-1.),
                                        std::fabs(ymax-1.)));
  extent = std::min(0.45,1.12*extent);
  configureFlavorStyle("Run2024I + Summer24 DY");
  TH1D frame(("frame_binning_"+std::string(flavor.name)).c_str(),"",
             100,xmin,xmax);
  frame.SetDirectory(nullptr);
  frame.SetMinimum(1.-extent);
  frame.SetMaximum(1.+1.35*extent);
  frame.GetXaxis()->SetTitle("reference p_{T} (GeV)");
  frame.GetYaxis()->SetTitle(Form("fitted data / MC %s response",
                                  flavor.label));
  frame.GetXaxis()->SetMoreLogLabels();
  frame.GetXaxis()->SetNoExponent();
  std::unique_ptr<TCanvas> canvas(tdrCanvas(
    ("c_binning_"+std::string(flavor.name)).c_str(),&frame,8,11,kSquare));
  canvas->SetLogx();
  TLine unity(xmin,1.,xmax,1.);
  unity.SetLineStyle(kDashed);
  unity.SetLineColor(kGray+2);
  unity.DrawClone();
  for (auto &graph : graphs) if (graph) graph->Draw("PZ SAME");
  TLegend legend(0.49,0.72,0.89,0.88);
  legend.SetBorderSize(0);
  legend.SetFillStyle(0);
  legend.SetNColumns(2);
  for (size_t index=0; index<graphs.size(); ++index)
    if (graphs[index]) legend.AddEntry(
      graphs[index].get(),binnings[index].label,"pl");
  legend.Draw();
  gPad->RedrawAxis();
  saveBoth(canvas.get(),outputDirectory,
    "flavor_response_binning_comparison_"+std::string(flavor.name));
}

void drawResponseResidual(TFile &analysis,
                          const std::string &outputDirectory) {
  TH1 *sourceHistogram = optionalObject<TH1>(
    &analysis,"response/h1_response_residual_data_over_mc");
  TGraphErrors *sourceGraph = optionalObject<TGraphErrors>(
    &analysis,"response/g_response_residual_data_over_mc");
  std::unique_ptr<TGraphErrors> graph;
  if (sourceHistogram)
    graph = graphFromHistogram(sourceHistogram);
  else if (sourceGraph)
    graph.reset(dynamic_cast<TGraphErrors*>(sourceGraph->Clone(
      "g_response_residual_data_over_mc_plot")));
  if (!graph || graph->GetN()==0) {
    std::cout << "Optional flavor-response residual is unavailable; skipping"
              << std::endl;
    return;
  }
  double extent = 0.02;
  for (int point=0; point<graph->GetN(); ++point) {
    double x = 0.;
    double y = 0.;
    graph->GetPoint(point,x,y);
    extent = std::max(extent,std::fabs(y-1.)+graph->GetErrorY(point));
  }
  extent = std::min(0.50,1.25*extent);
  configureFlavorStyle("Run2024I + Summer24 DY");
  TH1D frame("frame_response_residual","",4,-0.5,3.5);
  frame.SetDirectory(nullptr);
  frame.SetMinimum(1.-extent);
  frame.SetMaximum(1.+extent);
  frame.GetXaxis()->SetTitle("Generator parton flavor");
  frame.GetYaxis()->SetTitle("fitted data / MC flavor response");
  const char *labels[] = {"uds","c","b","g"};
  for (int bin=1; bin<=4; ++bin)
    frame.GetXaxis()->SetBinLabel(bin,labels[bin-1]);
  std::unique_ptr<TCanvas> canvas(tdrCanvas(
    "c_flavor_response_residual",&frame,8,11,kSquare));
  TLine unity(-0.5,1.,3.5,1.);
  unity.SetLineStyle(kDashed);
  unity.SetLineColor(kGray+2);
  unity.DrawClone();
  graph->SetMarkerStyle(kFullCircle);
  graph->SetMarkerColor(kBlack);
  graph->SetLineColor(kBlack);
  graph->Draw("PZ SAME");
  gPad->RedrawAxis();
  saveBoth(canvas.get(),outputDirectory,"flavor_response_residual");
}

void drawResponseResidualVsPt(TFile &analysis,
                              const std::string &outputDirectory) {
  std::vector<std::unique_ptr<TGraphErrors> > graphs;
  std::vector<const PtGraphStyle*> graphStyles;
  double extent = 0.03;
  double xmin = std::numeric_limits<double>::infinity();
  double xmax = 0.;
  for (const PtGraphStyle &style : kPtGraphStyles) {
    TGraphErrors *source = optionalObject<TGraphErrors>(
      &analysis,(std::string("response/g_response_residual_vs_pt_")+
                 style.name).c_str());
    if (!source) continue;
    std::unique_ptr<TGraphErrors> graph(dynamic_cast<TGraphErrors*>(
      source->Clone((std::string("plot_response_vs_pt_")+
                     style.name).c_str())));
    if (!graph) continue;
    graph->SetLineColor(style.color);
    graph->SetMarkerColor(style.color);
    graph->SetMarkerStyle(style.marker);
    graph->SetLineWidth(2);
    for (int point=0; point<graph->GetN(); ++point) {
      double x = 0.;
      double y = 0.;
      graph->GetPoint(point,x,y);
      if (std::isfinite(y)) {
        extent = std::max(extent,std::fabs(y-1.)+graph->GetErrorY(point));
        xmin = std::min(xmin,std::max(1.e-3,x-graph->GetErrorX(point)));
        xmax = std::max(xmax,x+graph->GetErrorX(point));
      }
    }
    graphs.push_back(std::move(graph));
    graphStyles.push_back(&style);
  }
  if (graphs.empty()) return;
  xmin = std::max(10.,xmin);
  xmax = std::min(1500.,xmax);
  extent = std::min(0.45,std::max(0.05,1.15*extent));
  configureFlavorStyle("Run2024I + Summer24 DY");
  TH1D frame("frame_response_residual_vs_pt","",100,xmin,xmax);
  frame.SetDirectory(nullptr);
  frame.SetMinimum(1.-extent);
  frame.SetMaximum(1.+extent);
  frame.GetXaxis()->SetTitle("p_{T,Z} (GeV)");
  frame.GetYaxis()->SetTitle("fitted data / MC flavor response");
  frame.GetXaxis()->SetMoreLogLabels();
  frame.GetXaxis()->SetNoExponent();
  std::unique_ptr<TCanvas> canvas(tdrCanvas(
    "c_flavor_response_residual_vs_pt",&frame,8,11,kSquare));
  canvas->SetLogx();
  TLine unity(xmin,1.,xmax,1.);
  unity.SetLineStyle(kDashed);
  unity.SetLineColor(kGray+2);
  unity.DrawClone();
  for (auto &graph : graphs) graph->Draw("PZ SAME");
  TLegend legend(0.58,0.75,0.88,0.88);
  legend.SetBorderSize(0);
  legend.SetFillStyle(0);
  legend.SetNColumns(2);
  for (size_t index=0; index<graphs.size(); ++index)
    legend.AddEntry(graphs[index].get(),graphStyles[index]->label,"pl");
  legend.Draw();
  gPad->RedrawAxis();
  saveBoth(canvas.get(),outputDirectory,
           "flavor_response_residual_vs_pt");
}

void drawOptionalAnalysis(TFile &analysis,
                          const std::string &outputDirectory) {
  drawAnalysisMatrix(analysis,"tagging/h2_efficiency_mc",
    "efficiency_mc","#epsilon_{MC}(reco | true)",0.,1.,outputDirectory);
  drawAnalysisMatrix(analysis,"tagging/h2_efficiency_data_inferred",
    "efficiency_data_inferred","inferred #epsilon_{data}(reco | true)",
    0.,1.,outputDirectory);
  drawAnalysisMatrix(analysis,"tagging/h2_transition_sf",
    "efficiency_transition_sf","#epsilon_{data} / #epsilon_{MC}",
    0.5,1.5,outputDirectory,true);
  drawAnalysisMatrix(analysis,"tagging/h2_composition_mc",
    "purity_mc","P_{MC}(true | reco)",0.,1.,outputDirectory);
  drawAnalysisMatrix(analysis,"tagging/h2_composition_data",
    "purity_data_inferred","inferred P_{data}(true | reco)",
    0.,1.,outputDirectory);
  drawAnalysisMatrix(analysis,
    "tagging/h2_composition_ratio_data_over_mc",
    "purity_ratio_data_over_mc","P_{data}/P_{MC}",
    0.5,1.5,outputDirectory,true);
  drawAnalysisMatrix(analysis,"response/h2_response_mc_by_transition",
    "response_transition_mc","MC HDM response",0.9,1.1,outputDirectory);
  drawAnalysisMatrix(
    analysis,"response/h2_response_data_estimated_by_transition",
    "response_transition_data_estimated","estimated data HDM response",
    0.9,1.1,outputDirectory);
  drawAnalysisMatrix(
    analysis,"response/h2_response_ratio_data_over_mc_by_transition",
    "response_transition_ratio_data_over_mc","estimated data / MC",
    0.97,1.03,outputDirectory,true);
  drawResponseResidual(analysis,outputDirectory);
  drawResponseResidualVsPt(analysis,outputDirectory);
  drawTagResponseRatio(analysis,false,outputDirectory);
  drawTagResponseRatio(analysis,true,outputDirectory);
  drawRecoilFraction(analysis,"fsr",false,outputDirectory);
  drawRecoilFraction(analysis,"fsr",true,outputDirectory);
  drawRecoilFraction(analysis,"ue",false,outputDirectory);
  drawRecoilFraction(analysis,"ue",true,outputDirectory);
  for (const PtGraphStyle &style : kPtGraphStyles)
    drawResponseBinningComparison(analysis,style,outputDirectory);
  drawDirectAnalysisMatrix(
    analysis,"diagnostics/h2_response_uncertainty_components",
    "response_uncertainty_components",0.,0.,outputDirectory,true);
}

} // namespace

void drawFlavorMatrix(
  const char *mcFile="rootfiles/zjet_MC.root",
  const char *dataFile="rootfiles/zjet_DATA.root",
  const char *analysisFile="",
  const char *outputDirectory="output/flavorMatrix") {
  gSystem->mkdir(outputDirectory,true);
  std::unique_ptr<TFile> mc(TFile::Open(mcFile,"READ"));
  std::unique_ptr<TFile> data(TFile::Open(dataFile,"READ"));
  if (!mc || mc->IsZombie())
    throw std::runtime_error("Failed to open FlavorMatrix MC input " +
                             std::string(mcFile));
  if (!data || data->IsZombie())
    throw std::runtime_error("Failed to open FlavorMatrix data input " +
                             std::string(dataFile));
  // Validate both inputs even though truth-separated tagger shapes come from
  // MC.  Data enter through the optional inference result.
  requireObject<TH3D>(mc.get(),"FlavorMatrix/h3counts_flavormatrix");
  requireObject<TH3D>(data.get(),"FlavorMatrix/h3counts_flavormatrix");

  drawScoreDistributions(*mc,outputDirectory);
  drawPairwiseDensities(*mc,outputDirectory);
  drawFlavorCubes(*mc,outputDirectory);
  drawQvgWorkingPointScan(*mc,outputDirectory);

  std::string resolvedAnalysis = (analysisFile ? analysisFile : "");
  if (resolvedAnalysis.empty()) {
    const std::string candidate =
      std::string(outputDirectory)+"/flavorMatrixAnalysis.root";
    if (!gSystem->AccessPathName(candidate.c_str())) resolvedAnalysis = candidate;
  }
  if (!resolvedAnalysis.empty()) {
    std::unique_ptr<TFile> analysis(TFile::Open(resolvedAnalysis.c_str(),"READ"));
    if (!analysis || analysis->IsZombie())
      throw std::runtime_error("Failed to open FlavorMatrix analysis input " +
                               resolvedAnalysis);
    drawOptionalAnalysis(*analysis,outputDirectory);
  }
  else {
    std::cout << "No FlavorMatrix analysis ROOT file was supplied; "
              << "efficiency/SF and response-residual plots were skipped."
              << std::endl;
  }
  std::cout << "Wrote FlavorMatrix PDF and PNG plots to "
            << outputDirectory << std::endl;
}
