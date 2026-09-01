// Publication-style plots for the compact UParT FlavorMatrix output.
//
// The macro deliberately draws score densities as binned heat maps.  It does
// not turn individual jets into TGraph points, which keeps vector output small
// and makes the figures deterministic in ROOT batch jobs.
#include "TBox.h"
#include "TCanvas.h"
#include "TColor.h"
#include "TFile.h"
#include "TGraphErrors.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLine.h"
#include "TPolyLine.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TSystem.h"

#include "tdrstyle_mod22.C"

#include <algorithm>
#include <cmath>
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
  int id;
  const char *name;
  const char *fileName;
  int color;
  int lineStyle;
};

// ID 0 is the truth-unavailable bookkeeping bin (and the only data truth
// bin), not a physical flavor.  Keep it in the efficiency matrices, but omit
// it from unit-normalized MC shape/cube plots: signed sideband subtraction can
// make its tiny denominator nearly cancel and obscure every physical flavor.
const std::vector<FlavorStyle> kTruthFlavors = {
  {1,"d+u","ud",kBlue+1,kSolid},
  {3,"s","s",kCyan+2,kDashed},
  {4,"c","c",kGreen+2,kSolid},
  {5,"b","b",kRed+1,kSolid},
  {6,"g","g",kOrange+7,kSolid},
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
  const int zbin = flavor.id+1;
  TH1D *projection = nullptr;
  const std::string name = "h_"+score+"_true"+flavor.fileName;
  if (score=="CvB") {
    projection = cvbCvl->ProjectionX(
      name.c_str(),1,cvbCvl->GetNbinsY(),zbin,zbin,"e");
  }
  else if (score=="CvL") {
    projection = cvbCvl->ProjectionY(
      name.c_str(),1,cvbCvl->GetNbinsX(),zbin,zbin,"e");
  }
  else if (score=="QvG") {
    projection = cvbQvg->ProjectionY(
      name.c_str(),1,cvbQvg->GetNbinsX(),zbin,zbin,"e");
  }
  else {
    throw std::runtime_error("Unknown UParT score " + score);
  }
  if (!projection)
    throw std::runtime_error("Failed to project UParT score " + score);
  projection->SetDirectory(nullptr);
  const double integral = signedIntegral(projection);
  if (std::fabs(integral)>1.e-12) projection->Scale(1./integral,"width");
  projection->SetLineColor(flavor.color);
  projection->SetMarkerColor(flavor.color);
  projection->SetLineStyle(flavor.lineStyle);
  projection->SetLineWidth(3);
  return std::unique_ptr<TH1D>(projection);
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
    frame.SetMaximum(std::max(1.e-6,1.30*maximum));
    frame.GetXaxis()->SetTitle(("UParT "+score).c_str());
    frame.GetYaxis()->SetTitle("1/N dN/d score");
    std::unique_ptr<TCanvas> canvas(tdrCanvas(
      ("c_"+score).c_str(),&frame,8,11,kSquare));
    for (const auto &histogram : histograms)
      histogram->Draw("HIST SAME");

    TLegend legend(0.58,0.42,0.90,0.63);
    legend.SetBorderSize(0);
    legend.SetFillStyle(0);
    legend.SetNColumns(2);
    legend.SetTextSize(0.038);
    for (size_t index=0; index<kTruthFlavors.size(); ++index)
      legend.AddEntry(histograms[index].get(),
                      (std::string("true ")+kTruthFlavors[index].name).c_str(),
                      "l");
    legend.Draw();
    TLatex selection;
    selection.SetNDC();
    selection.SetTextSize(0.038);
    selection.DrawLatex(0.18,0.72,"|#eta_{jet}| < 1.3, p_{T,jet} > 30 GeV");
    selection.SetTextSize(0.032);
    selection.DrawLatex(0.18,0.67,"valid UParT nodes; unit-area shapes");
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
  const int oldFirst = source->GetZaxis()->GetFirst();
  const int oldLast = source->GetZaxis()->GetLast();
  source->GetZaxis()->SetRange(flavor.id+1,flavor.id+1);
  TH2D *temporary = dynamic_cast<TH2D*>(source->Project3D("yxe"));
  source->GetZaxis()->SetRange(oldFirst,oldLast);
  if (!temporary)
    throw std::runtime_error("Failed to project " + std::string(spec.object));
  temporary->SetDirectory(nullptr);
  TH2D *projection = dynamic_cast<TH2D*>(temporary->Clone(
    (std::string("h2_")+spec.stem+"_true"+flavor.fileName).c_str()));
  delete temporary;
  if (!projection)
    throw std::runtime_error("Failed to clone pairwise score projection");
  projection->SetDirectory(nullptr);
  projection->SetTitle("");
  projection->GetXaxis()->SetTitle(spec.xTitle);
  projection->GetYaxis()->SetTitle(spec.yTitle);
  projection->GetZaxis()->SetTitle("normalized density");
  const double integral = signedIntegral(projection);
  if (std::fabs(integral)>1.e-12) projection->Scale(1./integral);
  return std::unique_ptr<TH2D>(projection);
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
  canvas->SetRightMargin(0.22);
  // Reserve a real header band for CMS, status and luminosity. Dense
  // heatmaps otherwise put the standard labels directly on populated cells.
  canvas->SetTopMargin(0.18);
  canvas->SetBottomMargin(0.13);
  histogram->UseCurrentStyle();
  histogram->GetXaxis()->SetTitleOffset(1.0);
  histogram->GetYaxis()->SetTitleOffset(1.25);
  histogram->GetZaxis()->SetTitleOffset(0.95);
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
  drawLabel(canvas->GetLeftMargin(),0.940,"CMS",61,0.050,11);
  const TString lumiText = lumi_136TeV + " (13.6 TeV)";
  drawLabel(1.-canvas->GetRightMargin(),0.940,lumiText,42,0.035,31);
  drawLabel(canvas->GetLeftMargin(),0.865,"Work in progress",52,0.032,11);
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
    line.DrawLine(0.5,0.5,0.5,1.);
  }
  else {
    line.SetLineStyle(kDashed);
    line.DrawLine(0.,0.5,1.,0.5);
    line.DrawLine(0.5,0.,0.5,0.5);
  }
}

void drawPairwiseDensities(TFile &mc, const std::string &outputDirectory) {
  const std::vector<PairSpec> specs = {
    {"FlavorMatrix/controls/h3_cvb_cvl_trueflavor","CvL_CvB",
     "UParT CvL","UParT CvB",0},
    {"FlavorMatrix/controls/h3_cvb_qvg_trueflavor","QvG_CvB",
     "UParT QvG","UParT CvB",1},
    {"FlavorMatrix/controls/h3_cvl_qvg_trueflavor","QvG_CvL",
     "UParT QvG","UParT CvL",2},
  };
  for (const PairSpec &spec : specs) {
    TH3D *source = requireObject<TH3D>(&mc,spec.object);
    std::vector<std::unique_ptr<TH2D> > projections;
    double minimum = 0.;
    double maximum = 0.;
    for (const FlavorStyle &flavor : kTruthFlavors) {
      projections.push_back(pairProjection(source,spec,flavor));
      minimum = std::min(minimum,projections.back()->GetMinimum());
      maximum = std::max(maximum,projections.back()->GetMaximum());
    }
    // A few bins can be slightly negative after sideband subtraction.  Do
    // not let a per-mille fluctuation force the full plot onto a diverging
    // scale; use one only when the negative excursion is material.
    if (minimum < -0.05*std::max(maximum,1.e-12)) {
      const double extent = std::max(std::fabs(minimum),std::fabs(maximum));
      minimum = -extent;
      maximum = extent;
      gStyle->SetPalette(kBlueRedYellow);
    }
    else {
      minimum = 0.;
      gStyle->SetPalette(kViridis);
    }
    for (size_t index=0; index<kTruthFlavors.size(); ++index) {
      TH2D *projection = projections[index].get();
      projection->SetMinimum(minimum);
      projection->SetMaximum(std::max(maximum,1.e-12));
      const std::string lumi = "Summer24 DY";
      std::unique_ptr<TCanvas> canvas = heatmapCanvas(
        std::string("c_")+spec.stem+"_"+kTruthFlavors[index].fileName,
        projection,lumi.c_str());
      // heatmapCanvas resets the style, so select the palette again here.
      gStyle->SetPalette(minimum<0. ? kBlueRedYellow : kViridis);
      projection->Draw("COLZ");
      drawPairCuts(spec);
      const std::string qualifier =
        std::string("true ")+kTruthFlavors[index].name;
      drawHeatmapHeader(canvas.get(),qualifier.c_str());
      gPad->RedrawAxis();
      saveBoth(canvas.get(),outputDirectory,
        std::string("tagger2d_")+spec.stem+"_true"+
        kTruthFlavors[index].fileName);
    }
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

int paletteColor(double value, double minimum, double maximum) {
  if (maximum<=minimum) return kWhite;
  const double fraction = std::max(0.,std::min(1.,
    (value-minimum)/(maximum-minimum)));
  const int colors = std::max(1,gStyle->GetNumberOfColors());
  const int index = std::min(colors-1,
    static_cast<int>(fraction*(colors-1)+0.5));
  return gStyle->GetColorPalette(index);
}

void addFilledCell(
  const Point2D &a, const Point2D &b, const Point2D &c, const Point2D &d,
  double value, double minimum, double maximum,
  std::vector<std::unique_ptr<TPolyLine> > &polygons) {
  double x[] = {a.x,b.x,c.x,d.x,a.x};
  double y[] = {a.y,b.y,c.y,d.y,a.y};
  std::unique_ptr<TPolyLine> polygon(new TPolyLine(5,x,y));
  const int color = paletteColor(value,minimum,maximum);
  polygon->SetFillColor(color);
  polygon->SetFillStyle(1001);
  polygon->SetLineColor(color);
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

void drawCube(
  const CubeFaces &faces, const FlavorStyle &flavor,
  double minimum, double maximum, const std::string &outputDirectory) {
  configureFlavorStyle("Summer24 DY simulation");
  gStyle->SetPalette(minimum<0. ? kBlueRedYellow : kViridis);
  TH1D frame((std::string("cube_frame_")+flavor.fileName).c_str(),"",
             100,0.,1.12);
  frame.SetDirectory(nullptr);
  frame.SetMinimum(0.);
  frame.SetMaximum(1.02);
  frame.GetXaxis()->SetLabelSize(0.);
  frame.GetYaxis()->SetLabelSize(0.);
  frame.GetXaxis()->SetTickLength(0.);
  frame.GetYaxis()->SetTickLength(0.);
  std::unique_ptr<TCanvas> canvas(tdrCanvas(
    (std::string("c_cube_")+flavor.fileName).c_str(),&frame,8,11,kSquare));

  const Point2D origin = {0.16,0.15};
  const Point2D xAxis = {0.56,0.};
  const Point2D yAxis = {0.,0.55};
  const Point2D qAxis = {0.16,0.13};
  const int bins = faces.bins;
  std::vector<std::unique_ptr<TPolyLine> > polygons;
  polygons.reserve(3*bins*bins);

  // Front face: CvL x CvB.
  for (int iy=0; iy<bins; ++iy) {
    for (int ix=0; ix<bins; ++ix) {
      const double x0 = static_cast<double>(ix)/bins;
      const double x1 = static_cast<double>(ix+1)/bins;
      const double y0 = static_cast<double>(iy)/bins;
      const double y1 = static_cast<double>(iy+1)/bins;
      addFilledCell(
        add(origin,add(scale(xAxis,x0),scale(yAxis,y0))),
        add(origin,add(scale(xAxis,x1),scale(yAxis,y0))),
        add(origin,add(scale(xAxis,x1),scale(yAxis,y1))),
        add(origin,add(scale(xAxis,x0),scale(yAxis,y1))),
        cell(faces.front,bins,ix,iy),minimum,maximum,polygons);
    }
  }
  // Top face: CvL x QvG, at CvB=1.
  for (int iq=0; iq<bins; ++iq) {
    for (int ix=0; ix<bins; ++ix) {
      const double x0 = static_cast<double>(ix)/bins;
      const double x1 = static_cast<double>(ix+1)/bins;
      const double q0 = static_cast<double>(iq)/bins;
      const double q1 = static_cast<double>(iq+1)/bins;
      const Point2D base = add(origin,yAxis);
      addFilledCell(
        add(base,add(scale(xAxis,x0),scale(qAxis,q0))),
        add(base,add(scale(xAxis,x1),scale(qAxis,q0))),
        add(base,add(scale(xAxis,x1),scale(qAxis,q1))),
        add(base,add(scale(xAxis,x0),scale(qAxis,q1))),
        cell(faces.top,bins,ix,iq),minimum,maximum,polygons);
    }
  }
  // Right face: CvB x QvG, at CvL=1.
  for (int iq=0; iq<bins; ++iq) {
    for (int iy=0; iy<bins; ++iy) {
      const double y0 = static_cast<double>(iy)/bins;
      const double y1 = static_cast<double>(iy+1)/bins;
      const double q0 = static_cast<double>(iq)/bins;
      const double q1 = static_cast<double>(iq+1)/bins;
      const Point2D base = add(origin,xAxis);
      addFilledCell(
        add(base,add(scale(yAxis,y0),scale(qAxis,q0))),
        add(base,add(scale(yAxis,y1),scale(qAxis,q0))),
        add(base,add(scale(yAxis,y1),scale(qAxis,q1))),
        add(base,add(scale(yAxis,y0),scale(qAxis,q1))),
        cell(faces.right,bins,iy,iq),minimum,maximum,polygons);
    }
  }

  // Cube outline and exact visible tag-region boundaries.
  drawSegment(origin,add(origin,xAxis),3);
  drawSegment(origin,add(origin,yAxis),3);
  drawSegment(add(origin,xAxis),add(add(origin,xAxis),yAxis),3);
  drawSegment(add(origin,yAxis),add(add(origin,yAxis),xAxis),3);
  drawSegment(add(origin,yAxis),add(add(add(origin,yAxis),xAxis),qAxis),3);
  drawSegment(add(add(origin,yAxis),xAxis),
              add(add(add(origin,yAxis),xAxis),qAxis),3);
  drawSegment(add(add(origin,xAxis),yAxis),
              add(add(add(origin,xAxis),yAxis),qAxis),3);
  drawSegment(add(add(add(origin,yAxis),xAxis),qAxis),
              add(add(origin,yAxis),qAxis),3);
  drawSegment(add(add(add(origin,yAxis),xAxis),qAxis),
              add(add(origin,xAxis),qAxis),3);

  // CvB=0.5: front and right faces.
  drawSegment(add(origin,scale(yAxis,0.5)),
              add(add(origin,xAxis),scale(yAxis,0.5)),3);
  drawSegment(add(add(origin,xAxis),scale(yAxis,0.5)),
              add(add(add(origin,xAxis),scale(yAxis,0.5)),qAxis),3);
  // CvL=0.5 after the b veto: front and top faces.
  drawSegment(add(add(origin,scale(xAxis,0.5)),scale(yAxis,0.5)),
              add(add(origin,scale(xAxis,0.5)),yAxis),3);
  drawSegment(add(add(origin,yAxis),scale(xAxis,0.5)),
              add(add(add(origin,yAxis),scale(xAxis,0.5)),qAxis),3);
  // QvG=0.5 in the high-CvB, low-CvL region on the top face.
  drawSegment(add(add(origin,yAxis),scale(qAxis,0.5)),
              add(add(add(origin,yAxis),scale(xAxis,0.5)),scale(qAxis,0.5)),3);

  TLatex text;
  text.SetTextAlign(22);
  text.SetTextSize(0.034);
  text.DrawLatex(0.44,0.10,"UParT CvL");
  text.SetTextAngle(90.);
  text.DrawLatex(0.105,0.43,"UParT CvB");
  text.SetTextAngle(34.);
  text.DrawLatex(0.84,0.76,"UParT QvG");
  text.SetTextAngle(0.);
  text.SetTextSize(0.040);
  text.DrawLatex(0.44,0.28,"b");
  text.DrawLatex(0.58,0.56,"c");
  text.DrawLatex(0.31,0.56,"uds / g");
  text.SetTextSize(0.032);
  text.DrawLatex(0.34,0.76,"g");
  text.DrawLatex(0.40,0.82,"uds");

  TLatex annotation;
  annotation.SetTextAlign(32);
  annotation.SetTextSize(0.040);
  annotation.DrawLatex(0.96,0.94,
    (std::string("true ")+flavor.name).c_str());

  // Compact in-frame palette; this avoids ROOT's large 3D palette objects.
  const int colorBins = 30;
  std::vector<std::unique_ptr<TBox> > colorBoxes;
  for (int index=0; index<colorBins; ++index) {
    const double y0 = 0.18+0.45*index/colorBins;
    const double y1 = 0.18+0.45*(index+1)/colorBins;
    const double value = minimum+(maximum-minimum)*(index+0.5)/colorBins;
    std::unique_ptr<TBox> box(new TBox(1.00,y0,1.035,y1));
    box->SetFillColor(paletteColor(value,minimum,maximum));
    box->SetLineColor(box->GetFillColor());
    box->Draw();
    colorBoxes.push_back(std::move(box));
  }
  TLatex paletteLabel;
  paletteLabel.SetTextSize(0.026);
  paletteLabel.SetTextAlign(12);
  paletteLabel.DrawLatex(1.045,0.18,Form("%.2g",minimum));
  paletteLabel.DrawLatex(1.045,0.63,Form("%.2g",maximum));
  paletteLabel.SetTextAngle(90.);
  paletteLabel.SetTextAlign(22);
  paletteLabel.DrawLatex(1.105,0.405,"normalized density");

  gPad->RedrawAxis();
  saveBoth(canvas.get(),outputDirectory,
           std::string("tagger_cube_true")+flavor.fileName);
}

void drawFlavorCubes(TFile &mc, const std::string &outputDirectory) {
  std::map<int,CubeFaces> allFaces;
  double minimum = 0.;
  double maximum = 0.;
  for (const FlavorStyle &flavor : kTruthFlavors) {
    TH3D *cube = requireObject<TH3D>(&mc,
      (std::string("FlavorMatrix/controls/h3_cvb_cvl_qvg_true")+
       std::to_string(flavor.id)).c_str());
    allFaces[flavor.id] = cubeFaces(cube);
    for (const std::vector<double> *face : {
           &allFaces[flavor.id].front,
           &allFaces[flavor.id].top,
           &allFaces[flavor.id].right}) {
      for (double value : *face) {
        minimum = std::min(minimum,value);
        maximum = std::max(maximum,value);
      }
    }
  }
  if (minimum < -0.05*std::max(maximum,1.e-12)) {
    const double extent = std::max(std::fabs(minimum),std::fabs(maximum));
    minimum = -extent;
    maximum = extent;
  }
  else {
    minimum = 0.;
  }
  for (const FlavorStyle &flavor : kTruthFlavors)
    drawCube(allFaces.at(flavor.id),flavor,minimum,
             std::max(maximum,1.e-12),outputDirectory);
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
  std::unique_ptr<TH2D> histogram(dynamic_cast<TH2D*>(source->Clone(
    ("plot_"+stem).c_str())));
  histogram->SetDirectory(nullptr);
  histogram->SetTitle("");
  histogram->GetXaxis()->SetTitle("Reco UParT flavor");
  histogram->GetYaxis()->SetTitle("Generator parton flavor");
  // Worker files keep category axes numeric so hadd does not attempt to
  // extend labelled axes.  Labels belong only on these detached plot clones.
  labelFlavorAxis(histogram->GetXaxis(),true);
  labelFlavorAxis(histogram->GetYaxis(),false);
  histogram->GetZaxis()->SetTitle(zTitle);
  histogram->SetMinimum(minimum);
  histogram->SetMaximum(maximum);
  configureFlavorStyle("2024I + Summer24");
  gStyle->SetPalette(diverging ? kBlueRedYellow : kViridis);
  gStyle->SetPaintTextFormat(".3f");
  std::unique_ptr<TCanvas> canvas = heatmapCanvas(
    "c_"+stem,histogram.get(),"2024I + Summer24");
  gStyle->SetPalette(diverging ? kBlueRedYellow : kViridis);
  gStyle->SetPaintTextFormat(".3f");
  histogram->Draw("COLZ TEXT");
  drawHeatmapHeader(canvas.get());
  gPad->RedrawAxis();
  saveBoth(canvas.get(),outputDirectory,stem);
}

std::unique_ptr<TGraphErrors> graphFromHistogram(TH1 *histogram) {
  if (!histogram) return nullptr;
  std::unique_ptr<TGraphErrors> graph(new TGraphErrors());
  graph->SetName("g_response_residual_plot");
  for (int bin=1; bin<=histogram->GetNbinsX(); ++bin) {
    const double value = histogram->GetBinContent(bin);
    const double error = histogram->GetBinError(bin);
    if (!std::isfinite(value) || (value==0. && error==0.)) continue;
    const int point = graph->GetN();
    graph->SetPoint(point,histogram->GetBinCenter(bin),value);
    graph->SetPointError(point,0.5*histogram->GetBinWidth(bin),error);
  }
  return graph;
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
  TH1D frame("frame_response_residual","",7,-0.5,6.5);
  frame.SetDirectory(nullptr);
  frame.SetMinimum(1.-extent);
  frame.SetMaximum(1.+extent);
  frame.GetXaxis()->SetTitle("Generator parton flavor");
  frame.GetYaxis()->SetTitle("fitted data / MC flavor response");
  labelFlavorAxis(frame.GetXaxis(),false);
  std::unique_ptr<TCanvas> canvas(tdrCanvas(
    "c_flavor_response_residual",&frame,8,11,kSquare));
  TLine unity(-0.5,1.,6.5,1.);
  unity.SetLineStyle(kDashed);
  unity.SetLineColor(kGray+2);
  unity.DrawClone();
  graph->SetMarkerStyle(kFullCircle);
  graph->SetMarkerColor(kBlack);
  graph->SetLineColor(kBlack);
  graph->Draw("PZ SAME");
  TLatex annotation;
  annotation.SetNDC();
  annotation.SetTextSize(0.035);
  annotation.DrawLatex(0.19,0.78,
                       "matrix-inferred response, p_{T} > 30 GeV");
  annotation.SetTextSize(0.030);
  annotation.DrawLatex(0.19,0.73,
                       "one-file smoke; null modes are prior-controlled");
  gPad->RedrawAxis();
  saveBoth(canvas.get(),outputDirectory,"flavor_response_residual");
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
  drawResponseResidual(analysis,outputDirectory);
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
