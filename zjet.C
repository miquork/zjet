#define zjet_cxx
#include "zjet.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

#include "TLorentzVector.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TGraphErrors.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TProfile3D.h"
#include "TObjString.h"
#include "TSystem.h"

#include "ZJetLumi.h"
#include "ZJetJerResolution.h"
#include "ZJetMuonCorrections.h"
#include "FlavorMatrixTools.h"
#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <ctime>
#include <iomanip>
#include <limits>
#include <map>
#include <memory>
#include <random>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

namespace {

struct ResponseAxis1D {
  TH1D *statistics = nullptr;
  TProfile *rmpf = nullptr;
  TProfile *rmpfjet1 = nullptr;
  TProfile *rmpfjetn = nullptr;
  TProfile *rmpfuncl = nullptr;
  TProfile *rmpfjetnu = nullptr;
  TProfile *rbal = nullptr;
  TProfile *rgenjet1 = nullptr;
  TProfile *residual = nullptr;
  TProfile *chHEF = nullptr;
  TProfile *neEmEF = nullptr;
  TProfile *neHEF = nullptr;
  TProfile *chEmEF = nullptr;
  TProfile *muEF = nullptr;
  TProfile *rho = nullptr;
  TH2D *mass = nullptr;
};

struct ResponseProfiles1D {
  std::map<std::string,ResponseAxis1D> axes;
};

struct TruthHDMAxis1D {
  TH1D *statistics = nullptr;
  TH1D *matchedStatistics = nullptr;
  TH1D *generatorStatistics = nullptr;
  TProfile *matchFraction = nullptr;
  TProfile *genJetPt = nullptr;
  TProfile *recoOverGen = nullptr;
  TProfile *deltaR = nullptr;
  TProfile *genOverBin = nullptr;
  TProfile *recoOverBin = nullptr;
  TProfile *genZOverRecoZ = nullptr;
  TProfile *recoMpf1 = nullptr;
  TProfile *recoMpfn = nullptr;
  TProfile *recoMpfu = nullptr;
  TProfile *recoMpf1Matched = nullptr;
  TProfile *recoMpfnMatched = nullptr;
  TProfile *recoMpfuMatched = nullptr;
  TProfile *genMpf1RecoAxis = nullptr;
  TProfile *genMpfnRecoAxis = nullptr;
  TProfile *genMpfuRecoAxis = nullptr;
  TProfile *genMpf1GenAxis = nullptr;
  TProfile *genMpfnGenAxis = nullptr;
  TProfile *genMpfuGenAxis = nullptr;
  TProfile *mpf1RecoGenProduct = nullptr;
  TProfile *mpfnRecoGenProduct = nullptr;
  TProfile *mpfuRecoGenProduct = nullptr;
  TProfile *mpf1GenSquared = nullptr;
  TProfile *mpfnGenSquared = nullptr;
  TProfile *mpfuGenSquared = nullptr;
  TProfile *previousResidual = nullptr;
  TProfile *dbWithoutResidual = nullptr;
  TProfile *mpf1WithoutResidual = nullptr;
};

struct TruthHDMProfiles1D {
  std::map<std::string,std::map<std::string,TruthHDMAxis1D> > regions;
};

struct FlavorMatrixComponents {
  double m0 = 0.;
  double m2 = 0.;
  double mn = 0.;
  double mu = 0.;
  double mnu = 0.;
  double hdm = 0.;
};

struct FlavorMatrixVariant {
  double x = 0.;
  FlavorMatrixComponents response;
};

struct FlavorMatrixHistograms {
  TH3D *counts = nullptr;
  TH3D *parallelCounts = nullptr;
  TH3D *transverseCounts = nullptr;
  std::map<std::string,TProfile3D*> profiles;
  std::map<std::string,TH3D*> pairControls;
  std::map<int,TH3D*> cubeControls;
};

struct GeneratorRecoil {
  bool hasGeneratorZ = false;
  TLorentzVector z;
  TLorentzVector ht;
  TLorentzVector met;
  TLorentzVector metu;
};

struct GeneratorPairComponents {
  bool valid = false;
  double genJetPt = 0.;
  double deltaR = 0.;
  double genMpf1RecoAxis = 0.;
  double genMpfnRecoAxis = 0.;
  double genMpfuRecoAxis = 0.;
  double genMpf1GenAxis = 0.;
  double genMpfnGenAxis = 0.;
  double genMpfuGenAxis = 0.;
  TLorentzVector met1;
  TLorentzVector metn;
};

const double responseBins[] = {
  12.,15.,20.,25.,30.,35.,40.,45.,50.,60.,70.,85.,105.,130.,175.,
  230.,300.,400.,500.,700.,1000.,1500.
};
const int responseBinCount = sizeof(responseBins)/sizeof(responseBins[0])-1;

const double flavorMatrixPtBins[] = {
  10.,15.,20.,25.,30.,35.,40.,50.,60.,70.,85.,100.,125.,155.,180.,
  210.,250.,300.,350.,400.,500.,600.,800.,1000.,1200.,1500.,1800.,
  2100.,2400.,2700.,3000.,3500.,4000.
};
const int flavorMatrixPtBinCount =
  sizeof(flavorMatrixPtBins)/sizeof(flavorMatrixPtBins[0])-1;
const double flavorIdBins[] = {
  -0.5,0.5,1.5,2.5,3.5,4.5,5.5,6.5,
};

int reconstructedUParTFlavor(double cvb, double cvl, double qvg) {
  if (!std::isfinite(cvb) || !std::isfinite(cvl) || !std::isfinite(qvg) ||
      cvb<0. || cvl<0. || qvg<0.)
    return 0;
  if (cvb<0.5) return 5;
  if (cvl>=0.5) return 4;
  return (qvg>=0.5 ? 1 : 6);
}

int generatorFlavorId(int partonFlavor) {
  const int flavor = std::abs(partonFlavor);
  if (flavor==1 || flavor==2) return 1;
  if (flavor==3 || flavor==4 || flavor==5) return flavor;
  if (flavor==21) return 6;
  return 0;
}

FlavorMatrixHistograms bookFlavorMatrix(TDirectory *parent) {
  if (!parent)
    throw std::runtime_error("Cannot book FlavorMatrix without a directory");
  TDirectory *directory = parent->mkdir("FlavorMatrix");
  if (!directory)
    throw std::runtime_error("Failed to create FlavorMatrix directory");
  directory->cd();

  FlavorMatrixHistograms result;
  result.counts = new TH3D(
    "h3counts_flavormatrix",
    ";p_{T,Z} (GeV);Reco UParT flavor;Generator parton flavor",
    flavorMatrixPtBinCount,flavorMatrixPtBins,
    7,flavorIdBins,7,flavorIdBins);
  result.counts->Sumw2();
  result.parallelCounts = dynamic_cast<TH3D*>(result.counts->Clone(
    "h3counts_parallel_flavormatrix"));
  result.parallelCounts->SetTitle(
    ";p_{T,Z} (GeV);Reco UParT flavor;Generator parton flavor");
  result.transverseCounts = dynamic_cast<TH3D*>(result.counts->Clone(
    "h3counts_transverse_flavormatrix"));
  result.transverseCounts->SetTitle(
    ";p_{T,Z} (GeV);Reco UParT flavor;Generator parton flavor");

  const std::vector<std::string> observables = {
    "m0", "m2", "mn", "mu", "mnu", "hdm",
  };
  const std::vector<std::string> variants = {"ab","ad","tc","pf"};
  for (const std::string &observable : observables) {
    for (const std::string &variant : variants) {
      const std::string name =
        "p3"+observable+variant+"_flavormatrix";
      TProfile3D *profile = new TProfile3D(
        name.c_str(),
        ";Reference p_{T} (GeV);Reco UParT flavor;Generator parton flavor",
        flavorMatrixPtBinCount,flavorMatrixPtBins,
        7,flavorIdBins,7,flavorIdBins);
      result.profiles[name] = profile;
    }
  }

  TDirectory *controls = directory->mkdir("controls");
  if (!controls)
    throw std::runtime_error("Failed to create FlavorMatrix/controls");
  controls->cd();
  const std::vector<std::tuple<std::string,std::string,std::string> > pairs = {
    {"h3_cvb_cvl_trueflavor","UParT CvB","UParT CvL"},
    {"h3_cvb_qvg_trueflavor","UParT CvB","UParT QvG"},
    {"h3_cvl_qvg_trueflavor","UParT CvL","UParT QvG"},
  };
  for (const auto &pair : pairs) {
    const std::string &name = std::get<0>(pair);
    TH3D *histogram = new TH3D(
      name.c_str(),
      Form(";%s;%s;Generator parton flavor",
           std::get<1>(pair).c_str(),std::get<2>(pair).c_str()),
      50,0.,1.,50,0.,1.,7,-0.5,6.5);
    histogram->Sumw2();
    result.pairControls[name] = histogram;
  }
  for (int flavor=0; flavor<=6; ++flavor) {
    TH3D *histogram = new TH3D(
      Form("h3_cvb_cvl_qvg_true%d",flavor),
      ";UParT CvB;UParT CvL;UParT QvG",
      24,0.,1.,24,0.,1.,24,0.,1.);
    histogram->Sumw2();
    result.cubeControls[flavor] = histogram;
  }
  return result;
}

FlavorMatrixComponents flavorMatrixComponents(
  double m0, double m2, double mn, double mu, double mnu) {
  FlavorMatrixComponents result;
  result.m0 = m0;
  result.m2 = m2;
  result.mn = mn;
  result.mu = mu;
  result.mnu = mnu;
  const double denominator = 1.-result.mn/1.00-result.mu/0.92;
  result.hdm = (std::fabs(denominator)>1.e-9
                  ? (result.m0-result.mn-result.mu)/denominator
                  : std::numeric_limits<double>::quiet_NaN());
  return result;
}

std::map<std::string,FlavorMatrixVariant> flavorMatrixVariants(
  const TLorentzVector &z, const TLorentzVector &jet,
  double m0, double m2, double mn, double mu, double mnu,
  bool transverse) {
  std::map<std::string,FlavorMatrixVariant> result;
  if (z.Pt()<=0. || jet.Pt()<=0.) return result;

  TLorentzVector zAxis;
  zAxis.SetPtEtaPhiM(1.,0.,z.Phi(),0.);
  TLorentzVector probeAxis;
  probeAxis.SetPtEtaPhiM(1.,0.,jet.Phi()+TMath::Pi(),0.);
  const double ptave = 0.5*(z.Pt()+jet.Pt());

  TLorentzVector bisector = zAxis+probeAxis;
  if (bisector.Pt()>0.) bisector *= 1./bisector.Pt();
  const double ptavp = transverse
    ? ptave
    : 0.5*((z.Px()*bisector.Px()+z.Py()*bisector.Py())-
           (jet.Px()*bisector.Px()+jet.Py()*bisector.Py()));

  // Z+jet already evaluates every component on the Z response axis. The
  // ab/ad/tc/pf variants deliberately change only the reference-pT binning,
  // matching the existing p2m* convention and avoiding a hidden change of
  // response definition inside a flavor bookkeeping object.
  const std::vector<std::pair<std::string,double> > variants = {
    {"ab",ptavp}, {"ad",ptave}, {"tc",z.Pt()}, {"pf",jet.Pt()},
  };
  for (const auto &entry : variants) {
    FlavorMatrixVariant variant;
    variant.x = entry.second;
    variant.response = flavorMatrixComponents(m0,m2,mn,mu,mnu);
    result[entry.first] = variant;
  }
  return result;
}

void fillFlavorMatrix(
  FlavorMatrixHistograms &histograms, double ptz, double ptj,
  int recoFlavor, int trueFlavor, double cvb, double cvl, double qvg,
  const std::map<std::string,FlavorMatrixVariant> &variants,
  double weight, bool transverse) {
  histograms.counts->Fill(ptz,recoFlavor,trueFlavor,weight);
  (transverse ? histograms.transverseCounts : histograms.parallelCounts)
    ->Fill(ptz,recoFlavor,trueFlavor,transverse ? -weight : weight);
  for (const auto &variant : variants) {
    const FlavorMatrixComponents &value = variant.second.response;
    const std::map<std::string,double> observables = {
      {"m0",value.m0}, {"m2",value.m2}, {"mn",value.mn},
      {"mu",value.mu}, {"mnu",value.mnu},
    };
    for (const auto &observable : observables) {
      const std::string name =
        "p3"+observable.first+variant.first+"_flavormatrix";
      if (std::isfinite(observable.second) &&
          std::isfinite(variant.second.x))
        histograms.profiles.at(name)->Fill(
          variant.second.x,recoFlavor,trueFlavor,observable.second,weight);
    }
  }
  // Tagger-shape controls describe the un-subtracted signal-window sample.
  // Filling them with the signed transverse weight would make ordinary
  // discriminator densities negative and would obscure the tagger working
  // point study.  The signed estimator remains available in counts/profiles,
  // while parallelCounts records the matching raw signal-window population
  // (which can still contain negative generator weights in NLO samples).
  if (transverse || ptj<=30. || cvb<0. || cvl<0. || qvg<0. ||
      !std::isfinite(cvb) || !std::isfinite(cvl) || !std::isfinite(qvg))
    return;
  histograms.pairControls.at("h3_cvb_cvl_trueflavor")->Fill(
    cvb,cvl,trueFlavor,weight);
  histograms.pairControls.at("h3_cvb_qvg_trueflavor")->Fill(
    cvb,qvg,trueFlavor,weight);
  histograms.pairControls.at("h3_cvl_qvg_trueflavor")->Fill(
    cvl,qvg,trueFlavor,weight);
  histograms.cubeControls.at(trueFlavor)->Fill(cvb,cvl,qvg,weight);
}

ResponseProfiles1D bookResponseProfiles1D(TDirectory *parent,
                                          const char *directoryName) {
  if (!parent)
    throw std::runtime_error("Cannot book one-dimensional response profiles");
  TDirectory *directory = parent->mkdir(directoryName);
  if (!directory)
    throw std::runtime_error("Failed to create response directory " +
                             std::string(directoryName));

  ResponseProfiles1D result;
  for (const char *axisName : {"zmmjet","jetpt","ptave"}) {
    TDirectory *axisDirectory = directory->mkdir(axisName);
    if (!axisDirectory)
      throw std::runtime_error("Failed to create one-dimensional axis " +
                               std::string(axisName));
    axisDirectory->cd();
    ResponseAxis1D &axis = result.axes[axisName];
    axis.statistics = new TH1D("statistics_rmpf",";p_{T} (GeV);Pairs",
                               responseBinCount,responseBins);
    axis.statistics->Sumw2();
    axis.rmpf = new TProfile("rmpf",";p_{T} (GeV);MPF",
                             responseBinCount,responseBins);
    axis.rmpfjet1 = new TProfile("rmpfjet1",";p_{T} (GeV);MPF1",
                                 responseBinCount,responseBins);
    axis.rmpfjetn = new TProfile("rmpfjetn",";p_{T} (GeV);MPFn",
                                 responseBinCount,responseBins);
    axis.rmpfuncl = new TProfile("rmpfuncl",";p_{T} (GeV);MPFu",
                                 responseBinCount,responseBins);
    axis.rmpfjetnu = new TProfile("rmpfjetnu",";p_{T} (GeV);MPFnu",
                                  responseBinCount,responseBins);
    axis.rbal = new TProfile("rbal",";p_{T} (GeV);DB",
                             responseBinCount,responseBins);
    axis.rgenjet1 = new TProfile("rgenjet1",";p_{T} (GeV);Gen balance",
                                 responseBinCount,responseBins);
    axis.residual = new TProfile("residual",";p_{T} (GeV);Previous residual",
                                 responseBinCount,responseBins);
    axis.chHEF = new TProfile("chHEF",";p_{T} (GeV);chHEF",
                              responseBinCount,responseBins);
    axis.neEmEF = new TProfile("neEmEF",";p_{T} (GeV);neEmEF",
                               responseBinCount,responseBins);
    axis.neHEF = new TProfile("neHEF",";p_{T} (GeV);neHEF",
                              responseBinCount,responseBins);
    axis.chEmEF = new TProfile("chEmEF",";p_{T} (GeV);chEmEF",
                               responseBinCount,responseBins);
    axis.muEF = new TProfile("muEF",";p_{T} (GeV);muEF",
                             responseBinCount,responseBins);
    axis.rho = new TProfile("rho",";p_{T} (GeV);#rho",
                            responseBinCount,responseBins);
    axis.mass = new TH2D("mass",";p_{T} (GeV);m_{#mu#mu} (GeV)",
                         responseBinCount,responseBins,120,60.,120.);
  }
  return result;
}

TruthHDMProfiles1D bookTruthHDMProfiles1D(
  TDirectory *parent, const char *directoryName,
  const std::vector<std::string> &regions) {
  if (!parent)
    throw std::runtime_error("Cannot book generator-level HDM controls");
  TDirectory *directory = parent->mkdir(directoryName);
  if (!directory)
    throw std::runtime_error("Failed to create truth directory " +
                             std::string(directoryName));
  TruthHDMProfiles1D result;
  for (const std::string &regionName : regions) {
    TDirectory *region = directory->mkdir(regionName.c_str());
    if (!region)
      throw std::runtime_error("Failed to create truth region " + regionName);
    for (const char *axisName : {"zmmjet","jetpt","ptave"}) {
      TDirectory *axisDirectory = region->mkdir(axisName);
      if (!axisDirectory)
        throw std::runtime_error("Failed to create truth axis " +
                                 std::string(axisName));
      axisDirectory->cd();
      TruthHDMAxis1D &axis = result.regions[regionName][axisName];
      auto profile = [&](const char *name, const char *title) {
        return new TProfile(name,title,responseBinCount,responseBins);
      };
      axis.statistics = new TH1D(
        "statistics_all",";p_{T} (GeV);Selected pairs",
        responseBinCount,responseBins);
      axis.matchedStatistics = new TH1D(
        "statistics_matched",";p_{T} (GeV);Reco-gen matched pairs",
        responseBinCount,responseBins);
      axis.generatorStatistics = new TH1D(
        "statistics_generator_recoil",
        ";p_{T} (GeV);Pairs with matched jet and generator Z",
        responseBinCount,responseBins);
      for (TH1D *histogram : {axis.statistics,axis.matchedStatistics,
                              axis.generatorStatistics})
        histogram->Sumw2();
      axis.matchFraction = profile(
        "match_fraction",";p_{T} (GeV);Reco-gen match fraction");
      axis.genJetPt = profile(
        "genjet_pt",";p_{T} (GeV);p_{T,gen jet} (GeV)");
      axis.recoOverGen = profile(
        "reco_over_gen",";p_{T} (GeV);p_{T,reco}/p_{T,gen}");
      axis.deltaR = profile(
        "delta_r",";p_{T} (GeV);#DeltaR(reco jet,gen jet)");
      axis.genOverBin = profile(
        "gen_over_bin",";p_{T} (GeV);p_{T,gen jet}/p_{T,bin}");
      axis.recoOverBin = profile(
        "reco_over_bin",";p_{T} (GeV);p_{T,reco jet}/p_{T,bin}");
      axis.genZOverRecoZ = profile(
        "gen_z_over_reco_z",";p_{T} (GeV);p_{T,gen Z}/p_{T,reco Z}");
      axis.recoMpf1 = profile(
        "reco_mpf1",";p_{T} (GeV);Reco MPF1 component");
      axis.recoMpfn = profile(
        "reco_mpfn",";p_{T} (GeV);Reco MPFn component");
      axis.recoMpfu = profile(
        "reco_mpfu",";p_{T} (GeV);Reco MPFu component");
      axis.recoMpf1Matched = profile(
        "reco_mpf1_matched",
        ";p_{T} (GeV);Reco MPF1 for generator-recoil pairs");
      axis.recoMpfnMatched = profile(
        "reco_mpfn_matched",
        ";p_{T} (GeV);Reco MPFn for generator-recoil pairs");
      axis.recoMpfuMatched = profile(
        "reco_mpfu_matched",
        ";p_{T} (GeV);Reco MPFu for generator-recoil pairs");
      axis.genMpf1RecoAxis = profile(
        "gen_mpf1_reco_axis",
        ";p_{T} (GeV);Gen MPF1 on reco-Z axis");
      axis.genMpfnRecoAxis = profile(
        "gen_mpfn_reco_axis",
        ";p_{T} (GeV);Gen MPFn on reco-Z axis");
      axis.genMpfuRecoAxis = profile(
        "gen_mpfu_reco_axis",
        ";p_{T} (GeV);Gen MPFu on reco-Z axis");
      axis.genMpf1GenAxis = profile(
        "gen_mpf1_gen_axis",
        ";p_{T} (GeV);Gen MPF1 on gen-Z axis");
      axis.genMpfnGenAxis = profile(
        "gen_mpfn_gen_axis",
        ";p_{T} (GeV);Gen MPFn on gen-Z axis");
      axis.genMpfuGenAxis = profile(
        "gen_mpfu_gen_axis",
        ";p_{T} (GeV);Gen MPFu on gen-Z axis");
      axis.mpf1RecoGenProduct = profile(
        "mpf1_reco_gen_product",
        ";p_{T} (GeV);Reco MPF1 #times gen MPF1");
      axis.mpfnRecoGenProduct = profile(
        "mpfn_reco_gen_product",
        ";p_{T} (GeV);Reco MPFn #times gen MPFn");
      axis.mpfuRecoGenProduct = profile(
        "mpfu_reco_gen_product",
        ";p_{T} (GeV);Reco MPFu #times gen MPFu");
      axis.mpf1GenSquared = profile(
        "mpf1_gen_squared",";p_{T} (GeV);(Gen MPF1)^{2}");
      axis.mpfnGenSquared = profile(
        "mpfn_gen_squared",";p_{T} (GeV);(Gen MPFn)^{2}");
      axis.mpfuGenSquared = profile(
        "mpfu_gen_squared",";p_{T} (GeV);(Gen MPFu)^{2}");
      axis.previousResidual = profile(
        "previous_residual",
        ";p_{T} (GeV);Inverse previous residual correction");
      axis.dbWithoutResidual = profile(
        "db_without_previous_residual",
        ";p_{T} (GeV);DB with selected-jet residual removed");
      axis.mpf1WithoutResidual = profile(
        "mpf1_without_previous_residual",
        ";p_{T} (GeV);MPF1 with selected-jet residual removed");
    }
  }
  return result;
}

void fillTruthHDMProfiles1D(
  TruthHDMProfiles1D &profiles, const std::string &region,
  double ptz, double ptj, double ptave, bool matched,
  const GeneratorRecoil &generatorRecoil,
  const GeneratorPairComponents &generatorPair, double recoMpf1,
  double recoMpfn, double recoMpfu, double inverseResidual, double weight) {
  const std::map<std::string,double> x = {
    {"zmmjet",ptz}, {"jetpt",ptj}, {"ptave",ptave},
  };
  for (const auto &entry : x) {
    TruthHDMAxis1D &axis = profiles.regions.at(region).at(entry.first);
    const double binValue = entry.second;
    axis.statistics->Fill(binValue,weight);
    axis.matchFraction->Fill(binValue,matched ? 1. : 0.,weight);
    axis.recoOverBin->Fill(binValue,ptj/binValue,weight);
    axis.previousResidual->Fill(binValue,inverseResidual,weight);
    axis.dbWithoutResidual->Fill(
      binValue,(ptj/ptz)*inverseResidual,weight);
    axis.mpf1WithoutResidual->Fill(
      binValue,recoMpf1*inverseResidual,weight);
    axis.recoMpf1->Fill(binValue,recoMpf1,weight);
    axis.recoMpfn->Fill(binValue,recoMpfn,weight);
    axis.recoMpfu->Fill(binValue,recoMpfu,weight);
    if (!matched || !generatorPair.valid || generatorPair.genJetPt<=0.)
      continue;
    axis.matchedStatistics->Fill(binValue,weight);
    axis.genJetPt->Fill(binValue,generatorPair.genJetPt,weight);
    axis.recoOverGen->Fill(
      binValue,ptj/generatorPair.genJetPt,weight);
    axis.deltaR->Fill(binValue,generatorPair.deltaR,weight);
    axis.genOverBin->Fill(
      binValue,generatorPair.genJetPt/binValue,weight);
    if (!generatorRecoil.hasGeneratorZ || generatorRecoil.z.Pt()<=0.)
      continue;
    axis.generatorStatistics->Fill(binValue,weight);
    axis.recoMpf1Matched->Fill(binValue,recoMpf1,weight);
    axis.recoMpfnMatched->Fill(binValue,recoMpfn,weight);
    axis.recoMpfuMatched->Fill(binValue,recoMpfu,weight);
    axis.genZOverRecoZ->Fill(
      binValue,generatorRecoil.z.Pt()/ptz,weight);
    axis.genMpf1RecoAxis->Fill(
      binValue,generatorPair.genMpf1RecoAxis,weight);
    axis.genMpfnRecoAxis->Fill(
      binValue,generatorPair.genMpfnRecoAxis,weight);
    axis.genMpfuRecoAxis->Fill(
      binValue,generatorPair.genMpfuRecoAxis,weight);
    axis.genMpf1GenAxis->Fill(
      binValue,generatorPair.genMpf1GenAxis,weight);
    axis.genMpfnGenAxis->Fill(
      binValue,generatorPair.genMpfnGenAxis,weight);
    axis.genMpfuGenAxis->Fill(
      binValue,generatorPair.genMpfuGenAxis,weight);
    axis.mpf1RecoGenProduct->Fill(
      binValue,recoMpf1*generatorPair.genMpf1RecoAxis,weight);
    axis.mpfnRecoGenProduct->Fill(
      binValue,recoMpfn*generatorPair.genMpfnRecoAxis,weight);
    axis.mpfuRecoGenProduct->Fill(
      binValue,recoMpfu*generatorPair.genMpfuRecoAxis,weight);
    axis.mpf1GenSquared->Fill(
      binValue,std::pow(generatorPair.genMpf1RecoAxis,2),weight);
    axis.mpfnGenSquared->Fill(
      binValue,std::pow(generatorPair.genMpfnRecoAxis,2),weight);
    axis.mpfuGenSquared->Fill(
      binValue,std::pow(generatorPair.genMpfuRecoAxis,2),weight);
  }
}

void writeTruthHDMDerivedGraphs(TruthHDMProfiles1D &profiles) {
  for (auto &regionEntry : profiles.regions) {
    for (auto &axisEntry : regionEntry.second) {
      TruthHDMAxis1D &axis = axisEntry.second;
      TDirectory *directory = axis.recoMpf1->GetDirectory();
      if (!directory) continue;
      const std::vector<std::tuple<const char*,TProfile*,TProfile*> > ratios = {
        {"response_r1_reco_axis",axis.recoMpf1Matched,
         axis.genMpf1RecoAxis},
        {"response_rn_reco_axis",axis.recoMpfnMatched,
         axis.genMpfnRecoAxis},
        {"response_ru_reco_axis",axis.recoMpfuMatched,
         axis.genMpfuRecoAxis},
        {"closure_r1_gen_axis",axis.recoMpf1Matched,
         axis.genMpf1GenAxis},
        {"closure_rn_gen_axis",axis.recoMpfnMatched,
         axis.genMpfnGenAxis},
        {"closure_ru_gen_axis",axis.recoMpfuMatched,
         axis.genMpfuGenAxis},
      };
      for (const auto &ratio : ratios) {
        const char *name = std::get<0>(ratio);
        TProfile *numerator = std::get<1>(ratio);
        TProfile *denominator = std::get<2>(ratio);
        std::unique_ptr<TGraphErrors> graph(new TGraphErrors());
        graph->SetName(name);
        graph->SetTitle(
          ";p_{T} (GeV);Ratio of reco and generator component means");
        for (int bin=1; bin<=numerator->GetNbinsX(); ++bin) {
          const double n = numerator->GetBinContent(bin);
          const double d = denominator->GetBinContent(bin);
          if (numerator->GetBinEntries(bin)==0. ||
              denominator->GetBinEntries(bin)==0. ||
              !std::isfinite(n) || !std::isfinite(d) ||
              std::fabs(d)<1.e-9)
            continue;
          const double en = numerator->GetBinError(bin);
          const double ed = denominator->GetBinError(bin);
          const double value = n/d;
          const double error = std::hypot(en/d,n*ed/(d*d));
          const int point = graph->GetN();
          graph->SetPoint(point,numerator->GetBinCenter(bin),value);
          graph->SetPointError(
            point,0.5*numerator->GetBinWidth(bin),std::fabs(error));
        }
        directory->cd();
        graph->Write(name,TObject::kOverwrite);
      }
      const std::vector<std::tuple<const char*,TProfile*,TProfile*> > slopes = {
        {"slope_r1_reco_axis",axis.mpf1RecoGenProduct,
         axis.mpf1GenSquared},
        {"slope_rn_reco_axis",axis.mpfnRecoGenProduct,
         axis.mpfnGenSquared},
        {"slope_ru_reco_axis",axis.mpfuRecoGenProduct,
         axis.mpfuGenSquared},
      };
      for (const auto &slope : slopes) {
        const char *name = std::get<0>(slope);
        TProfile *product = std::get<1>(slope);
        TProfile *square = std::get<2>(slope);
        std::unique_ptr<TGraphErrors> graph(new TGraphErrors());
        graph->SetName(name);
        graph->SetTitle(
          ";p_{T} (GeV);Zero-intercept reco-versus-gen slope");
        for (int bin=1; bin<=product->GetNbinsX(); ++bin) {
          const double numerator = product->GetBinContent(bin);
          const double denominator = square->GetBinContent(bin);
          if (product->GetBinEntries(bin)==0. ||
              square->GetBinEntries(bin)==0. ||
              !std::isfinite(numerator) || !std::isfinite(denominator) ||
              std::fabs(denominator)<1.e-12)
            continue;
          const double value = numerator/denominator;
          const double error = std::hypot(
            product->GetBinError(bin)/denominator,
            numerator*square->GetBinError(bin)/(denominator*denominator));
          const int point = graph->GetN();
          graph->SetPoint(point,product->GetBinCenter(bin),value);
          graph->SetPointError(
            point,0.5*product->GetBinWidth(bin),std::fabs(error));
        }
        directory->cd();
        graph->Write(name,TObject::kOverwrite);
      }
    }
  }
}

void fillResponseAxis1D(ResponseAxis1D &axis, double x, double db,
                        double mpf, double mpf1, double mpfn, double mpfu,
                        double mpfnu, double residual, double chHEF,
                        double neEmEF, double neHEF, double chEmEF,
                        double muEF, double rho, bool hasGenResponse,
                        double genResponse, double weight) {
  axis.statistics->Fill(x,weight);
  axis.rmpf->Fill(x,mpf,weight);
  axis.rmpfjet1->Fill(x,mpf1,weight);
  axis.rmpfjetn->Fill(x,mpfn,weight);
  axis.rmpfuncl->Fill(x,mpfu,weight);
  axis.rmpfjetnu->Fill(x,mpfnu,weight);
  axis.rbal->Fill(x,db,weight);
  axis.residual->Fill(x,residual,weight);
  axis.chHEF->Fill(x,chHEF,weight);
  axis.neEmEF->Fill(x,neEmEF,weight);
  axis.neHEF->Fill(x,neHEF,weight);
  axis.chEmEF->Fill(x,chEmEF,weight);
  axis.muEF->Fill(x,muEF,weight);
  axis.rho->Fill(x,rho,weight);
  if (hasGenResponse) axis.rgenjet1->Fill(x,genResponse,weight);
}

void fillResponseProfiles1D(ResponseProfiles1D &profiles, double ptz,
                            double ptj, double ptave, double db, double mpf,
                            double mpf1, double mpfn, double mpfu,
                            double mpfnu, double residual, double chHEF,
                            double neEmEF, double neHEF, double chEmEF,
                            double muEF, double rho, bool hasGenResponse,
                            double genResponse, double weight) {
  fillResponseAxis1D(profiles.axes.at("zmmjet"),ptz,db,mpf,mpf1,mpfn,mpfu,
                     mpfnu,residual,chHEF,neEmEF,neHEF,chEmEF,muEF,rho,
                     hasGenResponse,genResponse,weight);
  fillResponseAxis1D(profiles.axes.at("jetpt"),ptj,db,mpf,mpf1,mpfn,mpfu,
                     mpfnu,residual,chHEF,neEmEF,neHEF,chEmEF,muEF,rho,
                     hasGenResponse,genResponse,weight);
  fillResponseAxis1D(profiles.axes.at("ptave"),ptave,db,mpf,mpf1,mpfn,mpfu,
                     mpfnu,residual,chHEF,neEmEF,neHEF,chEmEF,muEF,rho,
                     hasGenResponse,genResponse,weight);
}

} // namespace

void zjet::Loop()
{
//   In a ROOT session, you can do:
//      root> .L zjet.C
//      root> zjet t
//      root> t.GetEntry(12); // Fill t data members with entry number 12
//      root> t.Show();       // Show values of entry 12
//      root> t.Show(16);     // Read and show values of entry 16
//      root> t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
   if (fChain == 0) return;

   fChain->SetBranchStatus("*",0);  // disable all branches

   fChain->SetBranchStatus("run",1);
   fChain->SetBranchStatus("luminosityBlock",1);
   fChain->SetBranchStatus("event",1);

   fChain->SetBranchStatus("nMuon",1);
   fChain->SetBranchStatus("Muon_pt",1);
   fChain->SetBranchStatus("Muon_eta",1);
   fChain->SetBranchStatus("Muon_phi",1);
   fChain->SetBranchStatus("Muon_mass",1);
   fChain->SetBranchStatus("Muon_charge",1);
   fChain->SetBranchStatus("Muon_looseId",1);
   fChain->SetBranchStatus("Muon_mediumId",1);
   fChain->SetBranchStatus("Muon_tightId",1);
   fChain->SetBranchStatus("Muon_pfIsoId",1);
   fChain->SetBranchStatus("Muon_pfRelIso04_all",1);
   fChain->SetBranchStatus("Muon_nTrackerLayers",1);

   fChain->SetBranchStatus("nJet",1);
   fChain->SetBranchStatus("Jet_pt",1);
   fChain->SetBranchStatus("Jet_eta",1);
   fChain->SetBranchStatus("Jet_phi",1);
   fChain->SetBranchStatus("Jet_mass",1);
   fChain->SetBranchStatus("Jet_rawFactor",1);
   fChain->SetBranchStatus("Jet_area",1);
   fChain->SetBranchStatus("Jet_chMultiplicity",1);
   fChain->SetBranchStatus("Jet_neMultiplicity",1);
   fChain->SetBranchStatus("Jet_nConstituents",1);
   fChain->SetBranchStatus("Jet_chHEF",1);
   fChain->SetBranchStatus("Jet_chEmEF",1);
   fChain->SetBranchStatus("Jet_neHEF",1);
   fChain->SetBranchStatus("Jet_neEmEF",1);
   fChain->SetBranchStatus("Jet_muEF",1);
   fChain->SetBranchStatus("Jet_btagDeepFlavB",1);
   fChain->SetBranchStatus("Jet_btagDeepFlavCvB",1);
   fChain->SetBranchStatus("Jet_btagDeepFlavCvL",1);
   fChain->SetBranchStatus("Jet_btagDeepFlavQG",1);
   fChain->SetBranchStatus("Jet_btagUParTAK4CvB",1);
   fChain->SetBranchStatus("Jet_btagUParTAK4CvL",1);
   fChain->SetBranchStatus("Jet_btagUParTAK4QvG",1);

   fChain->SetBranchStatus("PV_npvs",1);
   fChain->SetBranchStatus("Rho_fixedGridRhoFastjetAll",1);

   fChain->SetBranchStatus("Flag_goodVertices",1);
   fChain->SetBranchStatus("Flag_globalSuperTightHalo2016Filter",1);
   fChain->SetBranchStatus("Flag_HBHENoiseFilter",1);
   fChain->SetBranchStatus("Flag_HBHENoiseIsoFilter",1);
   fChain->SetBranchStatus("Flag_EcalDeadCellTriggerPrimitiveFilter",1);
   fChain->SetBranchStatus("Flag_BadPFMuonFilter",1);
   fChain->SetBranchStatus("Flag_BadPFMuonDzFilter",1);
   fChain->SetBranchStatus("Flag_hfNoisyHitsFilter",1);
   fChain->SetBranchStatus("Flag_eeBadScFilter",1);
   fChain->SetBranchStatus("Flag_ecalBadCalibFilter",1);
   fChain->SetBranchStatus("HLT_IsoMu24",1);
   fChain->SetBranchStatus("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8",1);
   fChain->SetBranchStatus("nTrigObj",1);
   fChain->SetBranchStatus("TrigObj_id",1);
   fChain->SetBranchStatus("TrigObj_pt",1);
   fChain->SetBranchStatus("TrigObj_eta",1);
   fChain->SetBranchStatus("TrigObj_phi",1);

   if (isMC) {
     fChain->SetBranchStatus("genWeight",1);
     fChain->SetBranchStatus("Pileup_nTrueInt",1);
     fChain->SetBranchStatus("nGenJet",1);
     fChain->SetBranchStatus("GenJet_pt",1);
     fChain->SetBranchStatus("GenJet_eta",1);
     fChain->SetBranchStatus("GenJet_phi",1);
     fChain->SetBranchStatus("GenJet_mass",1);
     fChain->SetBranchStatus("GenJet_partonFlavour",1);
     fChain->SetBranchStatus("Jet_genJetIdx",1);
     fChain->SetBranchStatus("Jet_partonFlavour",1);
     fChain->SetBranchStatus("nGenPart",1);
     fChain->SetBranchStatus("GenPart_pdgId",1);
     fChain->SetBranchStatus("GenPart_status",1);
     fChain->SetBranchStatus("GenPart_pt",1);
     fChain->SetBranchStatus("GenPart_eta",1);
     fChain->SetBranchStatus("GenPart_phi",1);
     fChain->SetBranchStatus("GenPart_mass",1);
     fChain->SetBranchStatus("GenMET_pt",1);
     fChain->SetBranchStatus("GenMET_phi",1);
   }

   fChain->SetBranchStatus("PuppiMET_pt",1);
   fChain->SetBranchStatus("PuppiMET_phi",1);
   fChain->SetBranchStatus("RawPuppiMET_pt",1);
   fChain->SetBranchStatus("RawPuppiMET_phi",1);
   
   cout << "Opening input files and reading entry metadata. "
        << "The first remote access can take a while..." << endl << flush;
   const auto metadataStart = std::chrono::steady_clock::now();
   // GetEntries() opens the files. This is intentionally not GetEntriesFast(),
   // because an exact count is needed for reliable progress and ETA reports.
   const Long64_t nentries = fChain->GetEntries();
   const double metadataSeconds =
     std::chrono::duration<double>(std::chrono::steady_clock::now()-
                                  metadataStart).count();
   if (nentries<=0) {
     cout << "ERROR: the input chain has zero readable entries after "
          << Form("%.1f",metadataSeconds) << " s. No output file will be "
          << "created. Check the first URL in the input list separately."
          << endl;
     return;
   }
   cout << "Input opened in " << Form("%.1f",metadataSeconds)
        << " s. Starting loop over dataset with " << nentries
        << " entries." << endl;
   if (isMC)  cout << "Running over MC branches" << endl;
   if (!isMC) cout << "Running over DATA branches" << endl;

   ZJetLumiData lumiData;
   if (!isMC && !goldenJsonFile.empty()) {
     if (!lumiData.loadGoldenJson(goldenJsonFile)) {
       cout << "Failed to load golden JSON " << goldenJsonFile << endl;
       return;
     }
   }
   if (!isMC && !lumiPileupFile.empty()) {
     if (!lumiData.loadPileup(lumiPileupFile)) {
       cout << "Failed to load lumisection pileup file " << lumiPileupFile << endl;
       return;
     }
   }

   TH1 *pileupWeights(0);
   if (isMC && !pileupWeightFile.empty()) {
     TFile pileupFile(pileupWeightFile.c_str(), "READ");
     TH1 *source = (TH1*)pileupFile.Get("pileup_ratio");
     if (!source) source = (TH1*)pileupFile.Get("pileup");
     if (!source) {
       cout << "Could not find pileup_ratio or pileup in "
            << pileupWeightFile << endl;
       return;
     }
     pileupWeights = (TH1*)source->Clone("zjet_pileup_weights");
     pileupWeights->SetDirectory(0);
   }

   const bool applyMuonCorrections = !muonCorrectionFile.empty();
   if (applyMuonCorrections) {
     std::ifstream input(muonCorrectionFile);
     if (!input.good()) {
       cout << "Failed to open Summer24 muon correction source "
            << muonCorrectionFile << endl;
       return;
     }
     cout << "Applying Summer24 nominal muon scale corrections to "
          << (isMC ? "MC" : "data");
     if (isMC) cout << " and the deterministic MC resolution correction";
     cout << "." << endl;
   }
   zjetcorrections::Summer24MuonCorrections muonCorrections;

   FactorizedJetCorrector *jec = 0;
   if (!jecL2File.empty()) {
     try {
       std::vector<JetCorrectorParameters> corrections;
       corrections.push_back(JetCorrectorParameters(jecL2File));
       if (!isMC && !jecResidualFile.empty())
         corrections.push_back(JetCorrectorParameters(jecResidualFile));
       jec = new FactorizedJetCorrector(corrections);
     }
     catch (const std::exception &error) {
       cout << "Failed to initialize JEC: " << error.what() << endl;
       return;
     }
     cout << "Recomputing " << (isMC ? "MC" : "data")
          << " jets from raw pT with L2 correction " << jecL2File;
     if (!isMC && !jecResidualFile.empty())
       cout << " and residual correction " << jecResidualFile;
     cout << "." << endl;
   }

   std::unique_ptr<zjetcorrections::JetPtResolution> jerResolution;
   std::unique_ptr<zjetcorrections::JetResolutionScaleFactor> jerScaleFactor;
   if (isMC && (jerResolutionFile.empty()!=jerScaleFactorFile.empty())) {
     cout << "JER smearing requires both a resolution and a scale-factor file."
          << endl;
     return;
   }
   if (isMC && !jerResolutionFile.empty()) {
     try {
       cout << "Initializing MC JER resolution from " << jerResolutionFile
            << "." << endl << flush;
       jerResolution.reset(
         new zjetcorrections::JetPtResolution(jerResolutionFile));
       cout << "Initializing MC JER scale factors from "
            << jerScaleFactorFile << "." << endl << flush;
       jerScaleFactor.reset(
         new zjetcorrections::JetResolutionScaleFactor(jerScaleFactorFile));
     }
     catch (const std::exception &error) {
       cout << "Failed to initialize JER: " << error.what() << endl;
       return;
     }
     cout << "Applying MC JER smearing with resolution " << jerResolutionFile
          << " (" << jerResolution->size() << " bins) and scale factors "
          << jerScaleFactorFile << " (" << jerScaleFactor->size()
          << " bins)." << endl;
   }
   else if (isMC) {
     cout << "MC JER smearing disabled." << endl;
   }

   TH2 *jetVetoMap = 0;
   if (!isMC && !jetVetoMapFile.empty()) {
     std::unique_ptr<TFile> mapFile(TFile::Open(jetVetoMapFile.c_str(),"READ"));
     TH2 *source = mapFile && !mapFile->IsZombie()
       ? dynamic_cast<TH2*>(mapFile->Get("jetvetomap")) : nullptr;
     if (!source) {
       cout << "Failed to load jetvetomap from " << jetVetoMapFile << endl;
       return;
     }
     jetVetoMap = dynamic_cast<TH2*>(source->Clone("zjet_jetvetomap"));
     jetVetoMap->SetDirectory(0);
     cout << "Applying data-only jet veto map " << jetVetoMapFile
          << ":jetvetomap." << endl;
   }

   const bool recalculatePuppiMet =
     (jec!=nullptr) || (isMC && jerResolution && jerScaleFactor);
   std::mt19937 jerRandomNumberGenerator(92837465);
   if (!recalculatePuppiMet) {
     cout << "JEC recomputation disabled; using stored NanoAOD jet pT and "
          << "PuppiMET." << endl;
   }

   const std::string::size_type separator = outputFile.find_last_of("/\\");
   if (separator!=std::string::npos) {
     const std::string outputDirectory = outputFile.substr(0,separator);
     gSystem->mkdir(outputDirectory.c_str(),kTRUE);
     if (gSystem->AccessPathName(outputDirectory.c_str())) {
       cout << "Failed to create output directory " << outputDirectory << endl;
       return;
     }
   }

   TDirectory *curdir = gDirectory;
   TFile *fout = new TFile(outputFile.c_str(),"RECREATE");
   if (!fout || fout->IsZombie()) {
     cout << "Failed to create output file " << outputFile << endl;
     return;
   }
   fout->cd();
   TObjString jecMode(jec ? "raw-pT JEC recomputation" : "stored NanoAOD JEC");
   TObjString jecL2(jecL2File.c_str());
   TObjString jecResidual((!isMC ? jecResidualFile : "").c_str());
   TObjString jerResolutionMetadata(
     (isMC ? jerResolutionFile : "").c_str());
   TObjString jerScaleFactorMetadata(
     (isMC ? jerScaleFactorFile : "").c_str());
   TObjString muonCorrectionMetadata(muonCorrectionFile.c_str());
   TObjString muonCorrectionSha(
     applyMuonCorrections ? ZJetMuonCorrectionData::sourceSha256 : "");
   TObjString jetVetoMapMetadata(
     (!isMC ? jetVetoMapFile : "").c_str());
   TObjString type1MetMetadata(
     recalculatePuppiMet
       ? "RawPuppiMET plus raw-minus-JEC/JER-corrected lepton-cleaned jets with corrected pT>15 GeV"
       : "stored PuppiMET");
   TObjString flavorDefinition(
     "Bettina/Sami DeepJet: B>0.7527; C=0.5*(CvB+CvL)>0.3985 after B veto; QG split at 0.5 after B/C veto");
   TObjString flavorMatrixDefinition(
     "FlavorMatrix uses |eta(jet)|<1.3 and the signed all-pairs signal-minus-two-half-weight-sidebands estimator; UParTAK4 tag IDs: undefined=0, uds=1, c=4, b=5, g=6; true IDs: undefined=0, d+u=1, s=3, c=4, b=5, g=6; data uses true ID 0 because truth is unavailable; HDM is derived from component means and re-finalized after hadd with Rn=1.00 and Ru=0.92; category axes remain numerically unlabelled until plotting so ROOT merges them without extension");

   
   // Object pT plots
   fout->mkdir("control");
   fout->cd("control");
   TH1D *h_cutflow = new TH1D("h_cutflow","",7,0.5,7.5);
   h_cutflow->GetXaxis()->SetBinLabel(1,"all");
   h_cutflow->GetXaxis()->SetBinLabel(2,"golden JSON");
   h_cutflow->GetXaxis()->SetBinLabel(3,"MET filters");
   h_cutflow->GetXaxis()->SetBinLabel(4,"HLT dimuon Mass8");
   h_cutflow->GetXaxis()->SetBinLabel(5,"synchronized dimuon");
   h_cutflow->GetXaxis()->SetBinLabel(6,"Z pT and mass");
   h_cutflow->GetXaxis()->SetBinLabel(7,"paired probe veto");

   // Alternative Z selections evaluated on the same HLT+filter event sample.
   // Except for the first bin, every entry also includes the narrow Z window.
   TH1D *h_muon_selection = new TH1D("h_muon_selection","",12,0.5,12.5);
   const char *muonSelectionLabels[] = {
     "HLT + filters", "OS eta", "loose ID", "medium ID", "tight ID",
     "medium + loose iso", "medium + medium iso",
     "medium + tight iso", "tight + tight iso", "tag-probe 27/10",
     "medium loose iso 27/20", "tight tight iso 27/20"
   };
   for (int ibin = 1; ibin <= 12; ++ibin)
     h_muon_selection->GetXaxis()->SetBinLabel(ibin,
                                               muonSelectionLabels[ibin-1]);

   TH1D *h_probe_veto = new TH1D("h_probe_veto","",5,0.5,5.5);
   h_probe_veto->GetXaxis()->SetBinLabel(1,"Z mass");
   h_probe_veto->GetXaxis()->SetBinLabel(2,"signal clear");
   h_probe_veto->GetXaxis()->SetBinLabel(3,"+90 pair valid");
   h_probe_veto->GetXaxis()->SetBinLabel(4,"-90 pair valid");
   h_probe_veto->GetXaxis()->SetBinLabel(5,"effective signal");
   TH1D *h_probe_pair_state = new TH1D("h_probe_pair_state","",4,-0.5,3.5);
   h_probe_pair_state->GetXaxis()->SetBinLabel(1,"neither");
   h_probe_pair_state->GetXaxis()->SetBinLabel(2,"+90 only");
   h_probe_pair_state->GetXaxis()->SetBinLabel(3,"-90 only");
   h_probe_pair_state->GetXaxis()->SetBinLabel(4,"both");

   TH1D *h_muon_scale_factor = new TH1D(
     "h_muon_scale_factor",";Muon scale factor;Muons",200,0.98,1.02);
   TH1D *h_muon_resolution_factor = new TH1D(
     "h_muon_resolution_factor",";Muon resolution factor;Muons",
     300,0.85,1.15);
   TH1D *h_jer_smear_factor = new TH1D(
     "h_jer_smear_factor",";JER smearing factor;Jets",300,0.7,1.3);
   TH1D *h_jet_veto_map = new TH1D("h_jet_veto_map","",3,0.5,3.5);
   h_jet_veto_map->GetXaxis()->SetBinLabel(1,"checked");
   h_jet_veto_map->GetXaxis()->SetBinLabel(2,"passed");
   h_jet_veto_map->GetXaxis()->SetBinLabel(3,"vetoed");

   std::map<std::string, TProfile*> pileupControl;
   const char *observables[] = {"npvs", "rho", "mu"};
   const int observableBins[] = {100, 100, 100};
   const double observableMax[] = {100., 100., 100.};
   const char *regions[] = {"parallel", "transverse", "subtracted"};
   for (int io = 0; io != 3; ++io) {
     for (int ir = 0; ir != 3; ++ir) {
       for (const char *response : {"db", "mpf"}) {
         const string name = Form("p_%s_vs_%s_%s", response,
                                  observables[io], regions[ir]);
         pileupControl[name] = new TProfile(name.c_str(), "",
                                            observableBins[io], 0.,
                                            observableMax[io]);
       }
     }
   }

   std::map<std::string, TProfile*> truthControl;
   for (const char *region : regions) {
     for (const char *observable : {"ptz", "npvs", "rho", "mu"}) {
       const int bins = (string(observable)=="ptz" ? 100 : 100);
       const double xmax = (string(observable)=="ptz" ? 200. : 100.);
       const string name = Form("p_pujet_fraction_vs_%s_%s", observable,
                                region);
       truthControl[name] = new TProfile(name.c_str(), "", bins, 0., xmax);
     }
   }
   std::map<std::string, TProfile*> methodTruthControl;
   for (const char *region : regions) {
     for (const char *category : {"all","matched","pileup"}) {
       for (const char *response :
            {"db","mpf","mpf1","mpfn","mpfu","mpfnu"}) {
         const string ptName = Form("p_%s_vs_ptz_%s_%s_central",response,
                                    category,region);
         methodTruthControl[ptName] =
           new TProfile(ptName.c_str(),"",100,0.,200.);
         const string etaName = Form(
           "p_%s_vs_abseta_%s_%s_ptz15to30",response,category,region);
         methodTruthControl[etaName] =
           new TProfile(etaName.c_str(),"",50,0.,5.);
       }
     }
   }
   for (const char *region : {"parallel", "transverse"}) {
     for (const char *category : {"matched", "pileup"}) {
       for (const char *response : {"db", "mpf"}) {
         const string name = Form("p_%s_vs_ptz_%s_%s", response, category,
                                  region);
         truthControl[name] = new TProfile(name.c_str(), "",100,0.,200.);
       }
     }
   }
   for (const char *region : regions) {
     const string indexName = Form("p_no_gen_index_fraction_vs_ptz_%s",region);
     truthControl[indexName] = new TProfile(indexName.c_str(),"",100,0.,200.);
     const string extraCutsName =
       Form("p_extra_match_cuts_unmatched_fraction_vs_ptz_%s",region);
     truthControl[extraCutsName] =
       new TProfile(extraCutsName.c_str(),"",100,0.,200.);
     for (const char *etaRegion : {"central","endcap","forward"}) {
       const string name = Form("p_pujet_fraction_vs_ptz_%s_%s",region,
                                etaRegion);
       truthControl[name] = new TProfile(name.c_str(),"",100,0.,200.);
     }
   }
   std::map<std::string, TH1D*> truthYield;
   std::map<std::string, TH2D*> truthMatchQuality;
   for (const char *region : regions) {
     for (const char *category : {"all","unmatched","no_gen_index",
                                  "extra_cuts_unmatched"}) {
       const string name = Form("h_%s_jets_vs_ptz_%s",category,region);
       truthYield[name] = new TH1D(name.c_str(),"",100,0.,200.);
       truthYield[name]->Sumw2();
     }
     const string name = Form("h2_truth_match_quality_vs_ptz_%s",region);
     truthMatchQuality[name] = new TH2D(name.c_str(),"",100,0.,200.,
                                        4,-0.5,3.5);
     truthMatchQuality[name]->GetYaxis()->SetBinLabel(1,"no gen index");
     truthMatchQuality[name]->GetYaxis()->SetBinLabel(2,"index, gen pT <= 8");
     truthMatchQuality[name]->GetYaxis()->SetBinLabel(3,"index, DeltaR >= 0.4");
     truthMatchQuality[name]->GetYaxis()->SetBinLabel(4,"index, extra cuts pass");
   }
   TH1D *h_truth_parallel = new TH1D("h_truth_parallel","",2,-0.5,1.5);
   TH1D *h_truth_transverse = new TH1D("h_truth_transverse","",2,-0.5,1.5);
   TH1D *h_truth_subtracted = new TH1D("h_truth_subtracted","",2,-0.5,1.5);
   for (TH1D *h : {h_truth_parallel,h_truth_transverse,h_truth_subtracted}) {
     h->GetXaxis()->SetBinLabel(1,"pileup/unmatched");
     h->GetXaxis()->SetBinLabel(2,"truth matched");
   }
   TH1D *h_nlep = new TH1D("h_nlep","",20,0,20);
   TH1D *h_lep1pt = new TH1D("h_lep1pt","",200,0,200);
   TH1D *h_lep2pt = new TH1D("h_lep2pt","",200,0,200);
   TH1D *h_leppt = new TH1D("h_leppt","",200,0,200);
   TH1D *h_lepeta = new TH1D("h_lepeta","",100,-5,5);
   TH2D *h2_lepeta_vs_ptz = new TH2D("h2_lepeta_vs_ptz","",200,0,200,100,-5,5);
   TH1D *h_lepdphi = new TH1D("h_lepdphi","",120,-TMath::TwoPi(),TMath::TwoPi());
   TH2D *h2_lepdphi_vs_ptz = new TH2D("h2_lepdphi_vs_ptz","",200,0,200,120,-TMath::TwoPi(),TMath::TwoPi());

   TH1D *h_lepdphiabs = new TH1D("h_lepdphiabs","",120,0,TMath::Pi());
   TH2D *h2_lepdphiabs_vs_ptz = new TH2D("h2_lepdphiabs_vs_ptz","",200,0,200,120,0,TMath::Pi());
   TH1D *h_lepdphimax = new TH1D("h_lepdphimax","",120,0,TMath::Pi());
   TH2D *h2_lepdphimax_vs_ptz = new TH2D("h2_lepdphimax_vs_ptz","",200,0,200,120,0,TMath::Pi());
   TH1D *h_lepdphimin = new TH1D("h_lepdphimin","",120,0,TMath::Pi());
   TH2D *h2_lepdphimin_vs_ptz = new TH2D("h2_lepdphimin_vs_ptz","",200,0,200,120,0,TMath::Pi());

   TH1D *h_zpt_precut = new TH1D("h_zpt_precut","",200,0,200);
   TH1D *h_zmass_precut = new TH1D("h_zmass_precut","",300,60,120);
   TH1D *h_zeta_precut = new TH1D("h_zeta_precut","",100,-5,5);
   TH2D *h2_zeta_precut_vs_ptz = new TH2D("h2_zeta_precut_vs_ptz","",200,0,200,100,-5,5);
   TH2D *h2_zmass_precut_vs_pt = new TH2D("h2_zmass_precut_vs_pt","",200,0,200,300,60,120);
   TProfile *p_zmass_precut_vs_pt = new TProfile("p_zmass_precut_vs_pt","",200,0,200);
   
   TH1D *h_zpt = new TH1D("h_zpt","",200,0,200);
   TH1D *h_zeta = new TH1D("h_zeta","",100,-5,5);
   TH2D *h2_zeta_vs_ptz = new TH2D("h2_zeta_vs_ptz","",200,0,200,100,-5,5);
   TH1D *h_zmass = new TH1D("h_zmass","",300,60,120);
   TH2D *h2_zmass_vs_pt = new TH2D("h2_zmass_vs_pt","",200,0,200,300,60,120);
   TProfile *p_zmass_vs_pt = new TProfile("p_zmass_vs_pt","",200,0,200);

   TH1D *h_zpt_probeveto = new TH1D("h_zpt_probeveto","",200,0,200);
   
   TH1D *h_njet = new TH1D("h_njet","",20,0,20);
   TH1D *h_jet1pt = new TH1D("h_jet1pt","",200,0,200);
   TH1D *h_jetpt = new TH1D("h_jetpt","",200,0,200);
   TH1D *h_jet1eta = new TH1D("h_jet1eta","",100,-5,5);
   TH1D *h_jeteta = new TH1D("h_jeteta","",100,-5,5);

   TH1D *h_nsel = new TH1D("h_nsel","",40,0,20);
   TH1D *h_sel1pt = new TH1D("h_sel1pt","",200,0,200);
   TH1D *h_selpt = new TH1D("h_selpt","",200,0,200);
   TH1D *h_sel1eta = new TH1D("h_sel1eta","",100,-5,5);
   TH1D *h_seleta = new TH1D("h_seleta","",100,-5,5);
   TH1D *h_seldphi = new TH1D("h_seldphi","",120,-TMath::TwoPi(),TMath::TwoPi());

   TH2D *h2_parpteta = new TH2D("h2_parpteta","",200,0,200,100,-5,5);
   TH1D *h_parpt = new TH1D("h_parpt","",200,0,200);
   TH1D *h_pareta = new TH1D("h_pareta","",100,-5,5);
   TH1D *h_pardphi = new TH1D("h_pardphi","",120,-TMath::TwoPi(),TMath::TwoPi());
   
   TH2D *h2_tranpteta = new TH2D("h2_tranpteta","",200,0,200,100,-5,5);
   TH1D *h_ntran = new TH1D("h_tran","",40,0,20);
   TH1D *h_tran1pt = new TH1D("h_tran1pt","",200,0,200);
   TH1D *h_tranpt = new TH1D("h_tranpt","",200,0,200);
   TH1D *h_tran1eta = new TH1D("h_tran1eta","",100,-5,5);
   TH1D *h_traneta = new TH1D("h_traneta","",100,-5,5);
   TH1D *h_trandphi = new TH1D("h_trandphi","",120,-TMath::TwoPi(),TMath::TwoPi());

   TH2D *h2_mixpteta = new TH2D("h2_mixpteta","",200,0,200,100,-5,5);
   TH1D *h_mixpt = new TH1D("h_mixpt","",200,0,200);
   TH1D *h_mixeta = new TH1D("h_mixeta","",100,-5,5);
   TH1D *h_mixdphi = new TH1D("h_mixdphi","",120,-TMath::TwoPi(),TMath::TwoPi());
   
   TH1D *h_db = new TH1D("h_db","",200,0,2);
   TH1D *h_mpf = new TH1D("h_mpf","",700,-3,4);
   TProfile2D *p2_db = new TProfile2D("p2_db","",40,0,200,100,-5,5);
   TProfile2D *p2_mpf = new TProfile2D("p2_mpf","",40,0,200,100,-5,5);
   TProfile2D *p2_mpfn = new TProfile2D("p2_mpfn","",40,0,200,100,-5,5);
   TProfile2D *p2_mpfu = new TProfile2D("p2_mpfu","",40,0,200,100,-5,5);
   TProfile2D *p2_mpfnu = new TProfile2D("p2_mpfnu","",40,0,200,100,-5,5);
   TH2D *h2_db = new TH2D("h2_db","",200,0,200,200,0,200);
   TProfile *p_db_vsz = new TProfile("p_db_vsz","",200,0,200);
   TProfile *p_db_vsj = new TProfile("p_db_vsj","",200,0,200);
   TProfile *p_db_vsa = new TProfile("p_db_vsa","",200,0,200);

   TH1D *h_dbp = new TH1D("h_dbp","",200,0,2);
   TH1D *h_mpfp = new TH1D("h_mpfp","",700,-3,4);
   TProfile2D *p2_dbp = new TProfile2D("p2_dbp","",40,0,200,100,-5,5);
   TProfile2D *p2_mpfp = new TProfile2D("p2_mpfp","",40,0,200,100,-5,5);
   TProfile2D *p2_mpfnp = new TProfile2D("p2_mpfnp","",40,0,200,100,-5,5);
   TProfile2D *p2_mpfup = new TProfile2D("p2_mpfup","",40,0,200,100,-5,5);
   TProfile2D *p2_mpfnup = new TProfile2D("p2_mpfnup","",40,0,200,100,-5,5);
   TH2D *h2_dbp = new TH2D("h2_dbp","",200,0,200,200,0,200);
   TProfile *p_dbp_vsz = new TProfile("p_dbp_vsz","",200,0,200);
   TProfile *p_dbp_vsj = new TProfile("p_dbp_vsj","",200,0,200);
   TProfile *p_dbp_vsa = new TProfile("p_dbp_vsa","",200,0,200);

   TH1D *h_dbt = new TH1D("h_dbt","",200,0,2);
   TH1D *h_mpft = new TH1D("h_mpft","",700,-3,4);
   TProfile2D *p2_dbt = new TProfile2D("p2_dbt","",40,0,200,100,-5,5);
   TProfile2D *p2_mpft = new TProfile2D("p2_mpft","",40,0,200,100,-5,5);
   TProfile2D *p2_mpfnt = new TProfile2D("p2_mpfnt","",40,0,200,100,-5,5);
   TProfile2D *p2_mpfut = new TProfile2D("p2_mpfut","",40,0,200,100,-5,5);
   TProfile2D *p2_mpfnut = new TProfile2D("p2_mpfnut","",40,0,200,100,-5,5);
   TH2D *h2_dbt = new TH2D("h2_dbt","",200,0,200,200,0,200);
   TProfile *p_dbt_vsz = new TProfile("p_dbt_vsz","",200,0,200);
   TProfile *p_dbt_vsj = new TProfile("p_dbt_vsj","",200,0,200);
   TProfile *p_dbt_vsa = new TProfile("p_dbt_vsa","",200,0,200);

   
   // Actual JEC stuff
   fout->mkdir("l2res");
   fout->cd("l2res");
   
   double vs[] = {0, 0.087, 0.174, 0.261, 0.348, 0.435, 0.522, 0.609, 0.696, 0.783, 0.879, 0.957, 1.044, 1.131, 1.218, 1.305, 1.392, 1.479, 1.566, 1.653, 1.74, 1.83, 1.93, 2.043, 2.172, 2.322, 2.5, 2.65, 2.853, 2.964, 3.139, 3.314, 3.489, 3.664, 3.839, 4.013, 4.191, 4.363, 4.538, 4.716, 4.889, 5.191};
   const int ns = sizeof(vs)/sizeof(vs[0])-1;
   double vp[] = {15, 21, 28, 37, 49, 59, 86, 110, 132, 170, 204, 236, 279, 302, 373, 460, 575, 638, 737, 846, 967, 1101, 1248, 1410, 1588, 1784, 2000, 2238, 2500, 2787, 3103};
   const int np = sizeof(vp)/sizeof(vp[0])-1;

   TH2D *h2ptetapf_ = new TH2D("h2ptetapf",";eta;probe",ns,vs,np,vp);
   TH2D *h2pteta_   = new TH2D("h2pteta",";eta;avp",ns,vs,np,vp);
   TH2D *h2ptetatc_ = new TH2D("h2ptetatc",";eta;tag",ns,vs,np,vp);

   TProfile *pmzpf_ = new TProfile("pmzpf",";probe;mz",np,vp);
   TProfile *pmz_ = new TProfile("pmz",";avp;mz",np,vp);
   TProfile *pmztc_ = new TProfile("pmztc",";tag;mz",np,vp);
   
   TProfile2D *p2jespf_ = new TProfile2D("p2jespf",";eta;probe",ns,vs,np,vp);
   TProfile2D *p2jes_   = new TProfile2D("p2jes",";eta;avp",ns,vs,np,vp);
   TProfile2D *p2jestc_ = new TProfile2D("p2jestc",";eta;tag",ns,vs,np,vp);

   TProfile2D *p2jsfpf_ = new TProfile2D("p2jsfpf",";eta;probe",ns,vs,np,vp);
   TProfile2D *p2jsf_   = new TProfile2D("p2jsf",";eta;avp",ns,vs,np,vp);
   TProfile2D *p2jsftc_ = new TProfile2D("p2jsftc",";eta;tag",ns,vs,np,vp);

   TProfile2D *p2m0pf_ = new TProfile2D("p2m0pf",";eta;probe",ns,vs,np,vp);
   TProfile2D *p2m0_   = new TProfile2D("p2m0",";eta;avp",ns,vs,np,vp);
   TProfile2D *p2m0tc_ = new TProfile2D("p2m0tc",";eta;tag",ns,vs,np,vp);

   TProfile2D *p2m2pf_ = new TProfile2D("p2m2pf",";eta;probe",ns,vs,np,vp);
   TProfile2D *p2m2_   = new TProfile2D("p2m2",";eta;avp",ns,vs,np,vp);
   TProfile2D *p2m2tc_ = new TProfile2D("p2m2tc",";eta;tag",ns,vs,np,vp);
   
   TProfile2D *p2mnpf_ = new TProfile2D("p2mnpf",";eta;probe",ns,vs,np,vp);
   TProfile2D *p2mn_   = new TProfile2D("p2mn",";eta;avp",ns,vs,np,vp);
   TProfile2D *p2mntc_ = new TProfile2D("p2mntc",";eta;tag",ns,vs,np,vp);
   
   TProfile2D *p2mnupf_ = new TProfile2D("p2mnupf",";eta;probe",ns,vs,np,vp);
   TProfile2D *p2mnu_   = new TProfile2D("p2mnu",";eta;avp",ns,vs,np,vp);
   TProfile2D *p2mnutc_ = new TProfile2D("p2mnutc",";eta;tag",ns,vs,np,vp);

   TProfile2D *p2mupf_ = new TProfile2D("p2mupf",";eta;probe",ns,vs,np,vp);
   TProfile2D *p2mu_   = new TProfile2D("p2mu",";eta;avp",ns,vs,np,vp);
   TProfile2D *p2mutc_ = new TProfile2D("p2mutc",";eta;tag",ns,vs,np,vp);

   TProfile2D *p2respf_ = new TProfile2D("p2respf",";eta;probe",ns,vs,np,vp);
   TProfile2D *p2res_   = new TProfile2D("p2res",";eta;avp",ns,vs,np,vp);
   TProfile2D *p2restc_ = new TProfile2D("p2restc",";eta;tag",ns,vs,np,vp);

   // Inputs used by the reprocess.C -> softrad3.C -> globalFit.C chain.
   // They use the same signed signal-minus-sideband weights as the response.
   TProfile2D *p2chftc_ = new TProfile2D("p2chftc",";eta;tag",ns,vs,np,vp);
   TProfile2D *p2neftc_ = new TProfile2D("p2neftc",";eta;tag",ns,vs,np,vp);
   TProfile2D *p2nhftc_ = new TProfile2D("p2nhftc",";eta;tag",ns,vs,np,vp);
   TProfile2D *p2ceftc_ = new TProfile2D("p2ceftc",";eta;tag",ns,vs,np,vp);
   TProfile2D *p2muftc_ = new TProfile2D("p2muftc",";eta;tag",ns,vs,np,vp);
   TProfile2D *p2rhotc_ = new TProfile2D("p2rhotc",";eta;tag",ns,vs,np,vp);
   TProfile2D *p2rgentc_ = new TProfile2D("p2rgentc",";eta;tag",ns,vs,np,vp);
   TH2D *h2mztc_ = new TH2D("h2mztc",";p_{T,Z};m_{#mu#mu}",
                             np,vp,120,60.,120.);

   // Flavor subsets follow the gamma+jet implementation maintained by
   // Bettina and Sami.  The first index is the reconstructed DeepJet tag and
   // the second is the matched generator parton flavor.  "i" is inclusive
   // and "n" means no classified tag/flavor.  These profiles use the same
   // signed signal-minus-sideband weights as the inclusive response above.
   fout->mkdir("flavor");
   fout->cd("flavor");
   std::map<std::string,
            std::map<std::string,std::map<std::string,TH1*> > > flavorProfiles;
   const std::vector<std::string> flavorVariables = {
     "counts", "mpfchs1", "ptchs", "mpf1", "mpfn", "mpfu",
     "rjet", "gjet",
   };
   const std::vector<std::string> recoFlavorTags = {
     "i", "b", "c", "q", "g", "n",
   };
   const std::vector<std::string> genFlavorTags = {
     "i", "b", "c", "q", "g", "n",
   };
   for (const std::string &variable : flavorVariables) {
     for (const std::string &recoTag : recoFlavorTags) {
       for (const std::string &genTag : genFlavorTags) {
         const std::string name =
           variable + "_g" + recoTag + genTag;
         if (variable=="counts") {
           TH1D *histogram = new TH1D(name.c_str(),";p_{T,Z};Events",np,vp);
           histogram->Sumw2();
           flavorProfiles[variable][recoTag][genTag] = histogram;
         }
         else {
           flavorProfiles[variable][recoTag][genTag] =
             new TProfile(name.c_str(),";p_{T,Z};Response",np,vp);
         }
       }
     }
   }

   FlavorMatrixHistograms flavorMatrix = bookFlavorMatrix(fout);

   fout->mkdir("l2res1");
   fout->cd("l2res1");

   double vx[] = {-5.191, -4.889, -4.716, -4.538, -4.363, -4.191, -4.013, -3.839, -3.664, -3.489, -3.314, -3.139, -2.964, -2.853, -2.65, -2.5, -2.322, -2.172, -2.043, -1.93, -1.83, -1.74, -1.653, -1.566, -1.479, -1.392, -1.305, -1.218, -1.131, -1.044, -0.957, -0.879, -0.783, -0.696, -0.609, -0.522, -0.435, -0.348, -0.261, -0.174, -0.087, 0, 0.087, 0.174, 0.261, 0.348, 0.435, 0.522, 0.609, 0.696, 0.783, 0.879, 0.957, 1.044, 1.131, 1.218, 1.305, 1.392, 1.479, 1.566, 1.653, 1.74, 1.83, 1.93, 2.043, 2.172, 2.322, 2.5, 2.65, 2.853, 2.964, 3.139, 3.314, 3.489, 3.664, 3.839, 4.013, 4.191, 4.363, 4.538, 4.716, 4.889, 5.191};
   const int nx = sizeof(vx)/sizeof(vx[0])-1;
   double vy[] = {15, 21, 28, 37, 49, 59, 86, 110, 132, 170, 204, 236, 279, 302, 373, 460, 575, 638, 737, 846, 967, 1101, 1248, 1410, 1588, 1784, 2000, 2238, 2500, 2787, 3103};
   const int ny = sizeof(vy)/sizeof(vy[0])-1;

   TH2D *h2ptetapf = new TH2D("h2ptetapf",";eta;probe",nx,vx,ny,vy);
   TH2D *h2pteta   = new TH2D("h2pteta",";eta;avp",nx,vx,ny,vy);
   TH2D *h2ptetatc = new TH2D("h2ptetatc",";eta;tag",nx,vx,ny,vy);

   TProfile *pmzpf = new TProfile("pmzpf",";probe;mz",ny,vy);
   TProfile *pmz = new TProfile("pmz",";eta;avp",ny,vy);
   TProfile *pmztc = new TProfile("pmztc",";tag;mz",ny,vy);
      
   TProfile2D *p2jespf = new TProfile2D("p2jespf",";eta;probe",nx,vx,ny,vy);
   TProfile2D *p2jes   = new TProfile2D("p2jes",";eta;avp",nx,vx,ny,vy);
   TProfile2D *p2jestc = new TProfile2D("p2jestc",";eta;tag",nx,vx,ny,vy);

   TProfile2D *p2jsfpf = new TProfile2D("p2jsfpf",";eta;probe",nx,vx,ny,vy);
   TProfile2D *p2jsf   = new TProfile2D("p2jsf",";eta;avp",nx,vx,ny,vy);
   TProfile2D *p2jsftc = new TProfile2D("p2jsftc",";eta;tag",nx,vx,ny,vy);

   TProfile2D *p2m0pf = new TProfile2D("p2m0pf",";eta;probe",nx,vx,ny,vy);
   TProfile2D *p2m0   = new TProfile2D("p2m0",";eta;avp",nx,vx,ny,vy);
   TProfile2D *p2m0tc = new TProfile2D("p2m0tc",";eta;tag",nx,vx,ny,vy);

   TProfile2D *p2m2pf = new TProfile2D("p2m2pf",";eta;probe",nx,vx,ny,vy);
   TProfile2D *p2m2   = new TProfile2D("p2m2",";eta;avp",nx,vx,ny,vy);
   TProfile2D *p2m2tc = new TProfile2D("p2m2tc",";eta;tag",nx,vx,ny,vy);
   
   TProfile2D *p2mnpf = new TProfile2D("p2mnpf",";eta;probe",nx,vx,ny,vy);
   TProfile2D *p2mn   = new TProfile2D("p2mn",";eta;avp",nx,vx,ny,vy);
   TProfile2D *p2mntc = new TProfile2D("p2mntc",";eta;tag",nx,vx,ny,vy);
   
   TProfile2D *p2mnupf = new TProfile2D("p2mnupf",";eta;probe",nx,vx,ny,vy);
   TProfile2D *p2mnu   = new TProfile2D("p2mnu",";eta;avp",nx,vx,ny,vy);
   TProfile2D *p2mnutc = new TProfile2D("p2mnutc",";eta;tag",nx,vx,ny,vy);

   TProfile2D *p2mupf = new TProfile2D("p2mupf",";eta;probe",nx,vx,ny,vy);
   TProfile2D *p2mu   = new TProfile2D("p2mu",";eta;avp",nx,vx,ny,vy);
   TProfile2D *p2mutc = new TProfile2D("p2mutc",";eta;tag",nx,vx,ny,vy);

   TProfile2D *p2respf = new TProfile2D("p2respf",";eta;probe",nx,vx,ny,vy);
   TProfile2D *p2res   = new TProfile2D("p2res",";eta;avp",nx,vx,ny,vy);
   TProfile2D *p2restc = new TProfile2D("p2restc",";eta;tag",nx,vx,ny,vy);

   // Preserve the native ZbAnalysis pT binning in genuine one-dimensional
   // profiles.  writeJecsys3.C prefers these over projections of the coarse
   // L2Res two-dimensional profiles, while retaining the latter for the
   // relative-eta workflow.
   ResponseProfiles1D allPairsProfiles =
     bookResponseProfiles1D(fout,"profiles1d");
   TruthHDMProfiles1D allPairsTruthProfiles = bookTruthHDMProfiles1D(
     fout,"truth_hdm",{"parallel","transverse","subtracted"});
   TDirectory *legacyDirectory = fout->mkdir("legacy");
   ResponseProfiles1D legacyProfiles =
     bookResponseProfiles1D(legacyDirectory,"profiles1d");
   TruthHDMProfiles1D legacyTruthProfiles = bookTruthHDMProfiles1D(
     legacyDirectory,"truth_hdm",{"parallel"});

   // The legacy and all-pairs methods select different jet populations.
   // Preserve method-specific previous-JEC profiles instead of silently using
   // the sideband-subtracted all-pairs profile for the legacy compatibility
   // output. The unsuffixed legacy p2res follows ZbAnalysis exactly: despite
   // its historical "avp" title it is filled versus tag (Z) pT.
   TDirectory *legacyL2ResDirectory = legacyDirectory->mkdir("l2res");
   legacyL2ResDirectory->cd();
   TH2D *legacyH2PtEta = new TH2D(
     "h2pteta",";|#eta|;p_{T,Z} (GeV)",ns,vs,np,vp);
   TH2D *legacyH2PtEtaPf = new TH2D(
     "h2ptetapf",";|#eta|;p_{T,jet} (GeV)",ns,vs,np,vp);
   TH2D *legacyH2PtEtaTc = new TH2D(
     "h2ptetatc",";|#eta|;p_{T,Z} (GeV)",ns,vs,np,vp);
   TProfile2D *legacyP2Jes = new TProfile2D(
     "p2jes",";|#eta|;p_{T,Z} (GeV);Inverse full JEC",ns,vs,np,vp);
   TProfile2D *legacyP2JesPf = new TProfile2D(
     "p2jespf",";|#eta|;p_{T,jet} (GeV);Inverse full JEC",ns,vs,np,vp);
   TProfile2D *legacyP2JesTc = new TProfile2D(
     "p2jestc",";|#eta|;p_{T,Z} (GeV);Inverse full JEC",ns,vs,np,vp);
   TProfile2D *legacyP2Res = new TProfile2D(
     "p2res",";|#eta|;p_{T,Z} (GeV);Inverse L2L3Residual",ns,vs,np,vp);
   TProfile2D *legacyP2ResPf = new TProfile2D(
     "p2respf",";|#eta|;p_{T,jet} (GeV);Inverse L2L3Residual",
     ns,vs,np,vp);
   TProfile2D *legacyP2ResTc = new TProfile2D(
     "p2restc",";|#eta|;p_{T,Z} (GeV);Inverse L2L3Residual",
     ns,vs,np,vp);
   TDirectory *legacyControlDirectory = legacyDirectory->mkdir("control");
   legacyControlDirectory->cd();
   TH1D *h_legacy_cutflow = new TH1D("h_cutflow","",9,0.5,9.5);
   const char *legacyCutLabels[] = {
     "all", "golden JSON", "MET filters", "dimuon trigger",
     "trigger-matched tight muons", "Z pT and mass", "leading jet",
     "back-to-back", "alpha < 1"
   };
   for (int ibin=1; ibin<=9; ++ibin)
     h_legacy_cutflow->GetXaxis()->SetBinLabel(ibin,legacyCutLabels[ibin-1]);
   TH1D *h_legacy_zpt = new TH1D("h_zpt",";p_{T,Z} (GeV);Events",
                                 200,0.,200.);
   TH1D *h_legacy_jetpt = new TH1D("h_jetpt",";p_{T,jet} (GeV);Events",
                                   200,0.,200.);
   TH1D *h_legacy_dphi = new TH1D("h_dphi",";|#Delta#phi|-#pi;Events",
                                  120,-TMath::Pi(),TMath::Pi());
   TH1D *h_legacy_alpha = new TH1D(
     "h_alpha",";p_{T,jet2}/p_{T,Z};Events",150,0.,3.);
   TH1D *h_legacy_subleading_jetpt = new TH1D(
     "h_subleading_jetpt",";p_{T,jet2} (GeV);Events",200,0.,1000.);
   TH2D *h_legacy_alpha_vs_jetpt = new TH2D(
     "h_alpha_vs_jetpt",";p_{T,jet1} (GeV);p_{T,jet2}/p_{T,Z}",
     200,0.,1000.,150,0.,3.);
   TProfile *p_legacy_db_before_alpha = dynamic_cast<TProfile*>(
     legacyProfiles.axes.at("jetpt").rbal->Clone(
       "p_db_vs_jetpt_before_alpha"));
   TProfile *p_legacy_mpf1_before_alpha = dynamic_cast<TProfile*>(
     legacyProfiles.axes.at("jetpt").rmpfjet1->Clone(
       "p_mpf1_vs_jetpt_before_alpha"));
   TProfile *p_legacy_rho_before_alpha = dynamic_cast<TProfile*>(
     legacyProfiles.axes.at("zmmjet").rho->Clone(
       "p_rho_vs_zpt_before_alpha"));
   p_legacy_db_before_alpha->Reset();
   p_legacy_mpf1_before_alpha->Reset();
   p_legacy_rho_before_alpha->Reset();
   p_legacy_db_before_alpha->SetDirectory(legacyControlDirectory);
   p_legacy_mpf1_before_alpha->SetDirectory(legacyControlDirectory);
   p_legacy_rho_before_alpha->SetDirectory(legacyControlDirectory);
   std::map<int,TProfile*> legacyRhoAlpha;
   for (const int alphaPercent : {10,15,20,30}) {
     TProfile *profile = dynamic_cast<TProfile*>(
       legacyProfiles.axes.at("zmmjet").rho->Clone(
         Form("p_rho_vs_zpt_alpha%03d",alphaPercent)));
     profile->Reset();
     profile->SetDirectory(legacyControlDirectory);
     legacyRhoAlpha[alphaPercent] = profile;
   }
   fout->cd();
   TObjString synchronizedSelection(
     "ZbAnalysis master 46dbf340 with production JetID setting: HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8; both muons trigger matched within DeltaR<0.3; tight ID; pfRelIso04<0.15; pT>20/10 GeV; |eta|<2.3; pT(Z)>12 GeV; |m-90 GeV|<20 GeV; MC pileup<=100; leading lepton-cleaned jet pT>=12 GeV and |eta|<=5; alpha=pT(jet2)/pT(Z), set to zero below pT(jet2)=15 GeV, alpha<1; central profiles use 0<|eta(jet1)|<1.3");
   TObjString legacyJetId(
     "disabled to match the current ZbAnalysis synchronization reference");
   TObjString truthDefinition(
     "truth_hdm profiles use status-1 generator muons matched to the selected reco muons within DeltaR<0.1; GenJet AK4 with pT>15 GeV and DeltaR(mu)>0.3 for generator HT; GenMET for invisible momentum; matched selected jet from Jet_genJetIdx; reco-axis and gen-axis recoil projections are both stored; R1, Rn and Ru should be formed from ratios of profile means, not means of event ratios");
   TObjString residualDefinition(
     "p2res stores the inverse final L2L3Residual factor; all-pairs l2res uses signed parallel-minus-transverse weights; legacy/l2res uses the legacy leading-jet selection without transverse subtraction and reproduces ZbAnalysis p2res binning versus Z pT");
   
   curdir->cd();

   TLorentzVector p4lplus, p4lminus, p4z, p4jet1, p4jet, p4sel1, p4tran1;
   TLorentzVector p4p, p4pz, p4t1, p4t1z, p4t2, p4t2z;
   TLorentzVector met, ht, met1, metn, metu, metnu, meta;
   TLorentzVector rawjet, rawjets, corrjets;

   // JMENANOv15 does not store Jet_jetId. Reconstruct the Run-3 Tight
   // PF Jet ID used by the standard NanoAOD jetId tight bit.
   auto passTightJetId = [&](int ijet) {
     const double abseta = fabs(Jet_eta[ijet]);
     if (abseta <= 2.6)
       return (Jet_neHEF[ijet] < 0.99 && Jet_neEmEF[ijet] < 0.90 &&
               Jet_nConstituents[ijet] > 1 && Jet_chHEF[ijet] > 0.01 &&
               Jet_chMultiplicity[ijet] > 0);
     if (abseta <= 2.7)
       return (Jet_neHEF[ijet] < 0.90 && Jet_neEmEF[ijet] < 0.99);
     if (abseta <= 3.0)
       return (Jet_neHEF[ijet] < 0.99);
     if (abseta < 5.0)
       return (Jet_neEmEF[ijet] < 0.40 && Jet_neMultiplicity[ijet] >= 2);
     return false;
   };

   const double mz = 91.1880;
   const double dmz = 1.5*2.4955; // 1.5*Gamma,Z~3.7 GeV
   
   Long64_t nbytes = 0, nb = 0;
   bool readFailure = false;
   const auto loopStart = std::chrono::steady_clock::now();
   auto previousProgress = loopStart;
   auto reportProgress = [&](Long64_t processed) {
     const auto now = std::chrono::steady_clock::now();
     const double elapsed =
       std::chrono::duration<double>(now-loopStart).count();
     const double sincePrevious =
       std::chrono::duration<double>(now-previousProgress).count();
     const bool earlyReport =
       (processed==1000 || processed==10000 || processed==100000);
     const bool periodicReport = (sincePrevious>=60.);
     const bool finalReport = (processed==nentries);
     if (!earlyReport && !periodicReport && !finalReport) return;

     const double rate = (elapsed>0. ? processed/elapsed : 0.);
     const double remaining =
       (rate>0. ? (nentries-processed)/rate : 0.);
     const std::time_t completionTime =
       std::time(0)+static_cast<std::time_t>(std::llround(remaining));
     const std::tm *localCompletion = std::localtime(&completionTime);
     cout << "Processed " << processed << "/" << nentries << " ("
          << Form("%.1f",100.*processed/nentries) << "%) in "
          << Form("%.1f",elapsed/60.) << " min at "
          << Form("%.0f",rate) << " events/s; "
          << Form("%.1f",remaining/60.) << " min remaining";
     if (localCompletion)
       cout << ", estimated completion "
            << std::put_time(localCompletion,"%Y-%m-%d %H:%M:%S");
     cout << "." << endl << flush;
     previousProgress = now;
   };

   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) {
        cout << "ERROR: failed to load event " << jentry << "." << endl;
        readFailure = true;
        break;
      }
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      if (nb<=0) {
        cout << "ERROR: failed to read event " << jentry
             << ". Stopping to avoid writing a silently incomplete sample."
             << endl;
        readFailure = true;
        break;
      }
      if (jentry==0)
        cout << "First event read successfully; the analysis loop is running."
             << endl << flush;
      reportProgress(jentry+1);
      // if (Cut(ientry) < 0) continue;

      double eventWeight = 1.;
      if (isMC) {
        eventWeight = (genWeight >= 0. ? 1. : -1.);
        if (pileupWeights) {
          const int bin = pileupWeights->GetXaxis()->FindFixBin(Pileup_nTrueInt);
          eventWeight *= pileupWeights->GetBinContent(bin);
        }
      }
      // ZbAnalysis does not apply genWeight to its nominal response profiles.
      // Keep that convention in the synchronization control only; the new
      // all-pairs method retains its signed generator-event weighting.
      const double legacyEventWeight = (isMC ? 1. : eventWeight);
      const double mu = (isMC ? Pileup_nTrueInt
                              : lumiData.pileup(run, luminosityBlock));
      // Apply this before JER smearing, as in ZbAnalysis. Besides selecting
      // the same events, this preserves the random-number sequence used for
      // stochastic smearing of later accepted events.
      if (isMC && mu>100.) continue;

      h_cutflow->Fill(1., eventWeight);
      h_legacy_cutflow->Fill(1.,legacyEventWeight);
      if (!isMC && !lumiData.accept(run, luminosityBlock)) continue;
      h_cutflow->Fill(2., eventWeight);
      h_legacy_cutflow->Fill(2.,legacyEventWeight);

      const bool passMetFilters =
        (isMC ||
         (Flag_goodVertices && Flag_globalSuperTightHalo2016Filter &&
         Flag_HBHENoiseFilter && Flag_HBHENoiseIsoFilter &&
         Flag_EcalDeadCellTriggerPrimitiveFilter && Flag_BadPFMuonFilter &&
         Flag_BadPFMuonDzFilter && Flag_eeBadScFilter &&
         Flag_ecalBadCalibFilter));
      if (!passMetFilters) continue;
      h_cutflow->Fill(3., eventWeight);
      h_legacy_cutflow->Fill(3.,legacyEventWeight);

      if (!HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8) continue;
      h_cutflow->Fill(4., eventWeight);
      h_legacy_cutflow->Fill(4.,legacyEventWeight);

      if (nMuon>nMuonMax || nJet>nJetMax || nTrigObj>nMaxTrigObj ||
          (isMC && (nGenJet>nMaxGenJet || nGenPart>nMaxGenPart))) {
        cout << "ERROR: collection size exceeds fixed MakeClass buffer: nMuon="
             << nMuon << ", nJet=" << nJet << ", nTrigObj=" << nTrigObj
             << ", nGenJet=" << (isMC ? nGenJet : 0)
             << ", nGenPart=" << (isMC ? nGenPart : 0)
             << endl;
        continue;
      }

      double jetInverseResidual[nJetMax];
      double jetStoredInverseJec[nJetMax];
      double jetRawPt[nJetMax];
      double jetRawMass[nJetMax];
      for (int ijet=0; ijet<nJet; ++ijet) {
        jetInverseResidual[ijet] = 1.;
        // Preserve the NanoAOD correction before Jet_rawFactor is updated to
        // describe the correction recomputed below. ZbAnalysis stores this
        // original inverse JEC in its legacy p2jes profiles.
        jetStoredInverseJec[ijet] = 1.-Jet_rawFactor[ijet];
        jetRawPt[ijet] = Jet_pt[ijet]*(1.-Jet_rawFactor[ijet]);
        jetRawMass[ijet] = Jet_mass[ijet]*(1.-Jet_rawFactor[ijet]);
      }

      // Keep trigger matching on the uncorrected NanoAOD muon, as in the
      // reference analysis, but use corrected four-vectors for pair ranking,
      // kinematic cuts, the Z boson, and lepton-jet cleaning.
      std::vector<TLorentzVector> correctedMuonP4(nMuon);
      for (int ilep=0; ilep<nMuon; ++ilep) {
        const double scaleFactor = applyMuonCorrections
          ? muonCorrections.scaleFactor(!isMC,Muon_pt[ilep],Muon_eta[ilep],
                                        Muon_phi[ilep],Muon_charge[ilep])
          : 1.;
        const double scaledPt = Muon_pt[ilep]*scaleFactor;
        const double resolutionFactor =
          (isMC && applyMuonCorrections)
          ? muonCorrections.resolutionFactor(
              scaledPt,Muon_eta[ilep],Muon_phi[ilep],
              Muon_nTrackerLayers[ilep],static_cast<int>(event),
              static_cast<int>(luminosityBlock))
          : 1.;
        correctedMuonP4[ilep].SetPtEtaPhiM(
          scaledPt*resolutionFactor,Muon_eta[ilep],Muon_phi[ilep],
          Muon_mass[ilep]*scaleFactor*resolutionFactor);
        h_muon_scale_factor->Fill(scaleFactor,eventWeight);
        if (isMC)
          h_muon_resolution_factor->Fill(resolutionFactor,eventWeight);
      }
      h_muon_selection->Fill(1., eventWeight);

      p4lplus.SetPtEtaPhiM(0,0,0,0);
      p4lminus.SetPtEtaPhiM(0,0,0,0);
      p4z.SetPtEtaPhiM(0,0,0,0);
      p4t1.SetPtEtaPhiM(0,0,0,0);
      p4t2.SetPtEtaPhiM(0,0,0,0);
      p4jet1.SetPtEtaPhiM(0,0,0,0);
      p4jet.SetPtEtaPhiM(0,0,0,0);
      p4sel1.SetPtEtaPhiM(0,0,0,0);
      p4tran1.SetPtEtaPhiM(0,0,0,0);
      ht.SetPtEtaPhiM(0,0,0,0);
      met.SetPtEtaPhiM(0,0,0,0);
      met1.SetPtEtaPhiM(0,0,0,0);
      metn.SetPtEtaPhiM(0,0,0,0);
      metu.SetPtEtaPhiM(0,0,0,0);
      metnu.SetPtEtaPhiM(0,0,0,0);
      meta.SetPtEtaPhiM(0,0,0,0);
      rawjet.SetPtEtaPhiM(0,0,0,0);
      rawjets.SetPtEtaPhiM(0,0,0,0);
      corrjets.SetPtEtaPhiM(0,0,0,0);
      int nlep(0);
      double nsel(0.);
      double ntran(0.);

      // The configurable selector below is retained for the control cut-flow.
      // The nominal analysis uses selectSynchronizedMuonPair further below.
      auto triggerMatchedMuon = [&](int ilep) {
        TLorentzVector muon;
        muon.SetPtEtaPhiM(Muon_pt[ilep],Muon_eta[ilep],Muon_phi[ilep],
                          Muon_mass[ilep]);
        for (int iobj=0; iobj<nTrigObj; ++iobj) {
          if (std::abs(int(TrigObj_id[iobj]))!=13) continue;
          TLorentzVector triggerObject;
          triggerObject.SetPtEtaPhiM(TrigObj_pt[iobj],TrigObj_eta[iobj],
                                     TrigObj_phi[iobj],0.);
          if (muon.DeltaR(triggerObject)<0.3) return true;
        }
        return false;
      };

      auto selectMuonPair = [&](int idWorkingPoint, int isolationWorkingPoint,
                                double minimumPt, double leadingPt,
                                bool requireTag, TLorentzVector &plus,
                                TLorentzVector &minus) {
        plus.SetPtEtaPhiM(0,0,0,0);
        minus.SetPtEtaPhiM(0,0,0,0);
        double bestMassDistance = 1.e9;
        auto passMuon = [&](int ilep) {
          const bool passId =
            (idWorkingPoint==0 ||
             (idWorkingPoint==1 && Muon_looseId[ilep]) ||
             (idWorkingPoint==2 && Muon_mediumId[ilep]) ||
             (idWorkingPoint==3 && Muon_tightId[ilep]));
          return (passId && Muon_pfIsoId[ilep]>=isolationWorkingPoint &&
                  correctedMuonP4[ilep].Pt()>minimumPt &&
                  fabs(correctedMuonP4[ilep].Eta())<2.4);
        };
        auto isTag = [&](int ilep) {
          return (correctedMuonP4[ilep].Pt()>27. && Muon_tightId[ilep] &&
                  Muon_pfIsoId[ilep]>=4);
        };
        for (int iplus = 0; iplus != nMuon; ++iplus) {
          if (Muon_charge[iplus]<=0 || !passMuon(iplus)) continue;
          const TLorentzVector &plusCandidate = correctedMuonP4[iplus];
          for (int iminus = 0; iminus != nMuon; ++iminus) {
            if (Muon_charge[iminus]>=0 || !passMuon(iminus)) continue;
            const TLorentzVector &minusCandidate = correctedMuonP4[iminus];
            if (max(plusCandidate.Pt(),minusCandidate.Pt())<=leadingPt)
              continue;
            if (requireTag && !isTag(iplus) && !isTag(iminus)) continue;
            const double massDistance =
              fabs((plusCandidate+minusCandidate).M()-mz);
            if (massDistance<bestMassDistance) {
              bestMassDistance = massDistance;
              plus = plusCandidate;
              minus = minusCandidate;
            }
          }
        }
        return (bestMassDistance<1.e8);
      };

      // Synchronized with the 2024 dimuon selection in ZbAnalysis master.
      // The closest opposite-sign pair to 90 GeV is retained.
      std::vector<TLorentzVector> synchronizedMuons;
      auto selectSynchronizedMuonPair = [&](TLorentzVector &plus,
                                            TLorentzVector &minus) {
        plus.SetPtEtaPhiM(0,0,0,0);
        minus.SetPtEtaPhiM(0,0,0,0);
        synchronizedMuons.clear();
        double bestMassDistance = 1.e9;
        for (int ilep=0; ilep<nMuon; ++ilep) {
          if (Muon_pt[ilep]<=8. || !Muon_tightId[ilep] ||
              Muon_pfRelIso04_all[ilep]>=0.15 ||
              !triggerMatchedMuon(ilep)) continue;
          synchronizedMuons.push_back(correctedMuonP4[ilep]);
        }
        if (synchronizedMuons.size()<2 || synchronizedMuons.size()>3)
          return false;
        for (int iplus=0; iplus<nMuon; ++iplus) {
          if (Muon_charge[iplus]<=0 || Muon_pt[iplus]<=8. ||
              !Muon_tightId[iplus] ||
              Muon_pfRelIso04_all[iplus]>=0.15 ||
              !triggerMatchedMuon(iplus)) continue;
          const TLorentzVector &plusCandidate = correctedMuonP4[iplus];
          for (int iminus=0; iminus<nMuon; ++iminus) {
            if (Muon_charge[iminus]>=0 || Muon_pt[iminus]<=8. ||
                !Muon_tightId[iminus] ||
                Muon_pfRelIso04_all[iminus]>=0.15 ||
                !triggerMatchedMuon(iminus)) continue;
            const TLorentzVector &minusCandidate = correctedMuonP4[iminus];
            const double massDistance =
              fabs((plusCandidate+minusCandidate).M()-90.);
            if (massDistance<bestMassDistance) {
              bestMassDistance = massDistance;
              plus = plusCandidate;
              minus = minusCandidate;
            }
          }
        }
        if (bestMassDistance>=20. || fabs(plus.Eta())>2.3 ||
            fabs(minus.Eta())>2.3) return false;
        return (std::max(plus.Pt(),minus.Pt())>20. &&
                std::min(plus.Pt(),minus.Pt())>10.);
      };

      auto fillMuonSelection = [&](int bin, int idWorkingPoint,
                                   int isolationWorkingPoint,
                                   double minimumPt, double leadingPt,
                                   bool requireTag) {
        TLorentzVector plus, minus;
        if (!selectMuonPair(idWorkingPoint,isolationWorkingPoint,minimumPt,
                            leadingPt,requireTag,plus,minus)) return;
        const TLorentzVector candidate = plus+minus;
        if (fabs(candidate.M()-mz)<dmz)
          h_muon_selection->Fill(bin,eventWeight);
      };

      fillMuonSelection(2,0,0,0.,0.,false);
      fillMuonSelection(3,1,0,0.,0.,false);
      fillMuonSelection(4,2,0,0.,0.,false);
      fillMuonSelection(5,3,0,0.,0.,false);
      fillMuonSelection(6,2,2,0.,0.,false);
      fillMuonSelection(7,2,3,0.,0.,false);
      fillMuonSelection(8,2,4,0.,0.,false);
      fillMuonSelection(9,3,4,0.,0.,false);
      fillMuonSelection(10,2,2,10.,27.,true);
      fillMuonSelection(11,2,2,20.,27.,false);
      fillMuonSelection(12,3,4,20.,27.,false);

      // Select the nominal synchronized dimuon pair.
      h_nlep->Fill(nMuon, eventWeight);
      if (!selectSynchronizedMuonPair(p4lplus,p4lminus)) continue;
      h_legacy_cutflow->Fill(5.,legacyEventWeight);
      auto separatedFromSynchronizedMuons = [&](const TLorentzVector &jet) {
        for (const TLorentzVector &muon : synchronizedMuons)
          if (jet.DeltaR(muon)<0.3) return false;
        return true;
      };

      // Recompute JEC only for lepton-cleaned jets and smear the first three
      // such jets in their original NanoAOD order. This matches ZbAnalysis;
      // sorting by the smeared pT happens implicitly only when the leading jet
      // is selected below.
      int cleanedJetOrdinal = 0;
      for (int ijet=0; ijet<nJet; ++ijet) {
        TLorentzVector storedJet;
        storedJet.SetPtEtaPhiM(Jet_pt[ijet],Jet_eta[ijet],Jet_phi[ijet],
                               Jet_mass[ijet]);
        if (!separatedFromSynchronizedMuons(storedJet)) continue;

        double correctedPt = Jet_pt[ijet];
        double correctedMass = Jet_mass[ijet];
        if (jec) {
          jec->setJetPt(jetRawPt[ijet]);
          jec->setJetEta(Jet_eta[ijet]);
          jec->setJetPhi(Jet_phi[ijet]);
          jec->setJetA(Jet_area[ijet]);
          jec->setRho(Rho_fixedGridRhoFastjetAll);
          jec->setNPV(PV_npvs);
          const std::vector<float> subCorrections = jec->getSubCorrections();
          if (subCorrections.empty() ||
              !std::isfinite(subCorrections.back()) ||
              subCorrections.back()<=0.) {
            cout << "ERROR: invalid JEC for event " << jentry << ", jet "
                 << ijet << "." << endl;
            readFailure = true;
            break;
          }
          const double correction = subCorrections.back();
          if (subCorrections.size()>1)
            jetInverseResidual[ijet] =
              subCorrections[subCorrections.size()-2]/correction;
          correctedPt = correction*jetRawPt[ijet];
          correctedMass = correction*jetRawMass[ijet];
        }

        if (isMC && jerResolution && jerScaleFactor &&
            cleanedJetOrdinal<3) {
          TLorentzVector correctedJet;
          correctedJet.SetPtEtaPhiM(correctedPt,Jet_eta[ijet],Jet_phi[ijet],
                                    correctedMass);
          const double resolution = jerResolution->resolution(
            Jet_eta[ijet],Rho_fixedGridRhoFastjetAll,correctedPt);
          const double scaleFactor = jerScaleFactor->scaleFactor(
            Jet_eta[ijet],correctedPt);
          if (!std::isfinite(resolution) || resolution<0. ||
              !std::isfinite(scaleFactor) || scaleFactor<=0.) {
            cout << "ERROR: invalid JER for event " << jentry << ", jet "
                 << ijet << "." << endl;
            readFailure = true;
            break;
          }

          double smearFactor = 1.;
          bool matched = false;
          if (Jet_genJetIdx[ijet]>=0 && Jet_genJetIdx[ijet]<nGenJet) {
            const int igen = Jet_genJetIdx[ijet];
            TLorentzVector generatorJet;
            generatorJet.SetPtEtaPhiM(
              GenJet_pt[igen],GenJet_eta[igen],GenJet_phi[igen],
              GenJet_mass[igen]);
            matched = (generatorJet.Pt()>0.01 &&
                       correctedJet.DeltaR(generatorJet)<0.2 &&
                       fabs(correctedPt-generatorJet.Pt())/
                         generatorJet.Pt()<3.*resolution);
            if (matched)
              smearFactor += (scaleFactor-1.)*
                (correctedPt-generatorJet.Pt())/generatorJet.Pt();
          }
          if (!matched && scaleFactor>1.) {
            const double width =
              resolution*std::sqrt(scaleFactor*scaleFactor-1.);
            std::normal_distribution<double> gaussian(0.,width);
            smearFactor += gaussian(jerRandomNumberGenerator);
          }
          const double minimumSmear = 0.01/correctedJet.E();
          if (smearFactor<minimumSmear) smearFactor = minimumSmear;
          h_jer_smear_factor->Fill(smearFactor,eventWeight);
          correctedPt *= smearFactor;
          correctedMass *= smearFactor;
        }
        ++cleanedJetOrdinal;

        Jet_pt[ijet] = correctedPt;
        Jet_mass[ijet] = correctedMass;
        Jet_rawFactor[ijet] =
          (correctedPt>0. ? 1.-jetRawPt[ijet]/correctedPt : 0.);
      }
      if (readFailure) break;

      auto passesJetVetoMap = [&](const TLorentzVector &jet) {
        if (!jetVetoMap) return true;
        const int etaBin = jetVetoMap->GetXaxis()->FindBin(jet.Eta());
        const int phiBin = jetVetoMap->GetYaxis()->FindBin(jet.Phi());
        return !(jetVetoMap->GetBinContent(etaBin,phiBin)>0.);
      };

      // Reconstruct Z boson
      if (p4lplus.Pt()>0) ++nlep;
      if (p4lminus.Pt()>0) ++nlep;
      if (p4lplus.Pt()>0 && p4lminus.Pt()>0) {
	p4z += p4lplus; p4z += p4lminus;
      }
      if (nlep != 2) continue;
      h_cutflow->Fill(5., eventWeight);
      //if (p4z.Pt()>0 && p4z.M()>80 && p4z.M()<100) {
      if (p4z.Pt()>0) {
	h_zpt_precut->Fill(p4z.Pt());
	h_zeta_precut->Fill(p4z.Eta());
	h2_zeta_precut_vs_ptz->Fill(p4z.Pt(), p4z.Eta());
	h_zmass_precut->Fill(p4z.M());
	h2_zmass_precut_vs_pt->Fill(p4z.Pt(),p4z.M());
	p_zmass_precut_vs_pt->Fill(p4z.Pt(),p4z.M());
      }
      if (p4z.Pt()>12. && fabs(p4z.M()-90.)<20.) {
	h_lep1pt->Fill(max(p4lplus.Pt(),p4lminus.Pt()));
	h_lep2pt->Fill(min(p4lplus.Pt(),p4lminus.Pt()));
	h_leppt->Fill(p4lplus.Pt());
	h_leppt->Fill(p4lminus.Pt());
	h_lepeta->Fill(p4lplus.Eta());
	h_lepeta->Fill(p4lminus.Eta());
	h2_lepeta_vs_ptz->Fill(p4z.Pt(), p4lplus.Eta());
	h2_lepeta_vs_ptz->Fill(p4z.Pt(), p4lminus.Eta());

	h_lepdphi->Fill(p4z.DeltaPhi(p4lplus));
	h_lepdphi->Fill(p4z.DeltaPhi(p4lminus));
	h2_lepdphi_vs_ptz->Fill(p4z.Pt(), p4z.DeltaPhi(p4lplus));
	h2_lepdphi_vs_ptz->Fill(p4z.Pt(), p4z.DeltaPhi(p4lminus));

	h_lepdphiabs->Fill(fabs(p4z.DeltaPhi(p4lplus)));
	h_lepdphiabs->Fill(fabs(p4z.DeltaPhi(p4lminus)));
	h2_lepdphiabs_vs_ptz->Fill(p4z.Pt(), fabs(p4z.DeltaPhi(p4lplus)));
	h2_lepdphiabs_vs_ptz->Fill(p4z.Pt(), fabs(p4z.DeltaPhi(p4lminus)));

	double dphimax = max(fabs(p4z.DeltaPhi(p4lplus)),
			     fabs(p4z.DeltaPhi(p4lminus)));
	h_lepdphimax->Fill(dphimax);
	h2_lepdphimax_vs_ptz->Fill(p4z.Pt(),dphimax);
	double dphimin = min(fabs(p4z.DeltaPhi(p4lplus)),
			     fabs(p4z.DeltaPhi(p4lminus)));
	h_lepdphimin->Fill(dphimin);
	h2_lepdphimin_vs_ptz->Fill(p4z.Pt(),dphimin);
	
	h_zpt->Fill(p4z.Pt());
	h_zeta->Fill(p4z.Eta());
	h2_zeta_vs_ptz->Fill(p4z.Pt(), p4z.Eta());
	h_zmass->Fill(p4z.M());
	h2_zmass_vs_pt->Fill(p4z.Pt(),p4z.M());
	p_zmass_vs_pt->Fill(p4z.Pt(),p4z.M());
      }
      else
	continue;
      h_cutflow->Fill(6., eventWeight);
      h_legacy_cutflow->Fill(6.,legacyEventWeight);
      h_probe_veto->Fill(1., eventWeight);

      // Build an independent generator-level recoil from the two status-one
      // muons matched to the selected reconstructed pair, particle-level jets
      // and GenMET. Keeping this decomposition separate from reco MET makes it
      // possible to extract R1, Rn and Ru as ratios of mean projected recoil
      // components without unstable event-by-event divisions near zero.
      GeneratorRecoil generatorRecoil;
      TLorentzVector generatorMuonPlus, generatorMuonMinus;
      if (isMC) {
        auto matchedGeneratorMuon = [&](const TLorentzVector &recoMuon,
                                        int pdgId,
                                        TLorentzVector &generatorMuon) {
          int bestIndex = -1;
          double bestDeltaR = 0.1;
          for (int igen=0; igen<nGenPart; ++igen) {
            if (GenPart_pdgId[igen]!=pdgId || GenPart_status[igen]!=1)
              continue;
            TLorentzVector candidate;
            candidate.SetPtEtaPhiM(
              GenPart_pt[igen],GenPart_eta[igen],GenPart_phi[igen],
              GenPart_mass[igen]);
            const double deltaR = recoMuon.DeltaR(candidate);
            if (deltaR<bestDeltaR) {
              bestDeltaR = deltaR;
              bestIndex = igen;
              generatorMuon = candidate;
            }
          }
          return bestIndex;
        };
        const int generatorPlusIndex = matchedGeneratorMuon(
          p4lplus,-13,generatorMuonPlus);
        const int generatorMinusIndex = matchedGeneratorMuon(
          p4lminus,+13,generatorMuonMinus);
        if (generatorPlusIndex>=0 && generatorMinusIndex>=0 &&
            generatorPlusIndex!=generatorMinusIndex) {
          generatorRecoil.hasGeneratorZ = true;
          generatorRecoil.z = generatorMuonPlus+generatorMuonMinus;
          generatorRecoil.ht = generatorRecoil.z;
          for (int igen=0; igen<nGenJet; ++igen) {
            TLorentzVector generatorJet;
            generatorJet.SetPtEtaPhiM(
              GenJet_pt[igen],GenJet_eta[igen],GenJet_phi[igen],
              GenJet_mass[igen]);
            if (generatorJet.Pt()<=15. ||
                generatorJet.DeltaR(generatorMuonPlus)<0.3 ||
                generatorJet.DeltaR(generatorMuonMinus)<0.3)
              continue;
            generatorRecoil.ht += generatorJet;
          }
          generatorRecoil.ht.SetPtEtaPhiM(
            generatorRecoil.ht.Pt(),0.,generatorRecoil.ht.Phi(),0.);
          generatorRecoil.met.SetPtEtaPhiM(
            GenMET_pt,0.,GenMET_phi,0.);
          generatorRecoil.metu = generatorRecoil.met+generatorRecoil.ht;
        }
      }

      auto generatorPairComponents = [&]
        (int generatorJetIndex, const TLorentzVector &reconstructedJet,
         const TLorentzVector &recoProjectionAxis,
         const TLorentzVector &genProjectionAxis, bool transverse) {
        GeneratorPairComponents result;
        if (!isMC || generatorJetIndex<0 || generatorJetIndex>=nGenJet)
          return result;
        TLorentzVector generatorJet;
        generatorJet.SetPtEtaPhiM(
          GenJet_pt[generatorJetIndex],GenJet_eta[generatorJetIndex],
          GenJet_phi[generatorJetIndex],GenJet_mass[generatorJetIndex]);
        if (generatorJet.Pt()<=0.) return result;
        result.valid = true;
        result.genJetPt = generatorJet.Pt();
        result.deltaR = reconstructedJet.DeltaR(generatorJet);
        if (!generatorRecoil.hasGeneratorZ ||
            generatorRecoil.z.Pt()<=0. || p4z.Pt()<=0.)
          return result;

        result.met1 = -generatorRecoil.z-generatorJet;
        result.met1.SetPtEtaPhiM(
          result.met1.Pt(),0.,result.met1.Phi(),0.);
        result.metn = -generatorRecoil.ht+generatorRecoil.z+generatorJet;
        result.metn.SetPtEtaPhiM(
          result.metn.Pt(),0.,result.metn.Phi(),0.);

        auto project = [](const TLorentzVector &component,
                          const TLorentzVector &axis, double denominator,
                          double offset) {
          return offset+component.Vect().Dot(axis.Vect())/denominator;
        };
        const double recoDenominator = p4z.Pt()*p4z.Pt();
        const double genDenominator =
          generatorRecoil.z.Pt()*generatorRecoil.z.Pt();
        const double baseMpf1Reco = project(
          result.met1,p4z,recoDenominator,1.);
        const double baseMpfnReco = project(
          result.metn,p4z,recoDenominator,0.);
        const double baseMpfuReco = project(
          generatorRecoil.metu,p4z,recoDenominator,0.);
        const double baseMpf1Gen = project(
          result.met1,generatorRecoil.z,genDenominator,1.);
        const double baseMpfnGen = project(
          result.metn,generatorRecoil.z,genDenominator,0.);
        const double baseMpfuGen = project(
          generatorRecoil.metu,generatorRecoil.z,genDenominator,0.);
        if (!transverse) {
          result.genMpf1RecoAxis = baseMpf1Reco;
          result.genMpfnRecoAxis = baseMpfnReco;
          result.genMpfuRecoAxis = baseMpfuReco;
          result.genMpf1GenAxis = baseMpf1Gen;
          result.genMpfnGenAxis = baseMpfnGen;
          result.genMpfuGenAxis = baseMpfuGen;
        }
        else {
          // Mirror the transverse estimator used below: retain the parallel
          // component and add the projection on the rotated axis.
          result.genMpf1RecoAxis = project(
            result.met1,recoProjectionAxis,recoDenominator,1.)+
            (baseMpf1Reco-1.);
          result.genMpfnRecoAxis = project(
            result.metn,recoProjectionAxis,recoDenominator,0.)+baseMpfnReco;
          result.genMpfuRecoAxis = project(
            generatorRecoil.metu,recoProjectionAxis,recoDenominator,0.)+
            baseMpfuReco;
          result.genMpf1GenAxis = project(
            result.met1,genProjectionAxis,genDenominator,1.)+
            (baseMpf1Gen-1.);
          result.genMpfnGenAxis = project(
            result.metn,genProjectionAxis,genDenominator,0.)+baseMpfnGen;
          result.genMpfuGenAxis = project(
            generatorRecoil.metu,genProjectionAxis,genDenominator,0.)+
            baseMpfuGen;
        }
        return result;
      };

      // Legacy leading-jet reference. It shares the synchronized event and
      // dimuon selection above but retains the reference analysis choices,
      // including its alpha definition and JetID setting.
      {
        TLorentzVector legacyMet, legacyHt, legacyRawJets, legacyCorrJets;
        legacyMet.SetPtEtaPhiM(
          recalculatePuppiMet ? RawPuppiMET_pt : PuppiMET_pt,0.,
          recalculatePuppiMet ? RawPuppiMET_phi : PuppiMET_phi,0.);
        legacyHt = p4z;
        std::vector<int> legacyJetIndices;
        for (int ijet=0; ijet<nJet; ++ijet) {
          TLorentzVector candidate;
          candidate.SetPtEtaPhiM(Jet_pt[ijet],Jet_eta[ijet],Jet_phi[ijet],
                                 Jet_mass[ijet]);
          if (!separatedFromSynchronizedMuons(candidate)) continue;
          if (candidate.Pt()>=10.) legacyJetIndices.push_back(ijet);
          if (candidate.Pt()>15.) legacyHt += candidate;
          if (recalculatePuppiMet && candidate.Pt()>=15.) {
            TLorentzVector rawCandidate;
            rawCandidate.SetPtEtaPhiM(
              jetRawPt[ijet],Jet_eta[ijet],Jet_phi[ijet],jetRawMass[ijet]);
            legacyRawJets += rawCandidate;
            legacyCorrJets += candidate;
          }
        }
        std::sort(legacyJetIndices.begin(),legacyJetIndices.end(),
                  [&](int first, int second) {
                    return Jet_pt[first]>Jet_pt[second];
                  });
        legacyHt.SetPtEtaPhiM(legacyHt.Pt(),0.,legacyHt.Phi(),0.);
        if (recalculatePuppiMet) legacyMet += legacyRawJets-legacyCorrJets;
        legacyMet.SetPtEtaPhiM(legacyMet.Pt(),0.,legacyMet.Phi(),0.);
        const TLorentzVector legacyMetu = legacyMet+legacyHt;

        int legacyJetIndex = -1;
        for (int candidateIndex : legacyJetIndices) {
          if (Jet_pt[candidateIndex]>=12. &&
              fabs(Jet_eta[candidateIndex])<=5.) {
            legacyJetIndex = candidateIndex;
            break;
          }
        }
        if (legacyJetIndex>=0) {
          TLorentzVector legacyJet;
          legacyJet.SetPtEtaPhiM(
            Jet_pt[legacyJetIndex],Jet_eta[legacyJetIndex],
            Jet_phi[legacyJetIndex],Jet_mass[legacyJetIndex]);
          TLorentzVector legacySubleadingJet;
          legacySubleadingJet.SetPtEtaPhiM(0.,0.,0.,0.);
          const auto subleading = std::find_if(
            legacyJetIndices.begin(),legacyJetIndices.end(),
            [&](int candidateIndex) {
              return candidateIndex!=legacyJetIndex;
            });
          if (subleading!=legacyJetIndices.end()) {
            const int subleadingIndex = *subleading;
            legacySubleadingJet.SetPtEtaPhiM(
              Jet_pt[subleadingIndex],Jet_eta[subleadingIndex],
              Jet_phi[subleadingIndex],Jet_mass[subleadingIndex]);
          }

          h_legacy_cutflow->Fill(7.,legacyEventWeight);
          h_legacy_zpt->Fill(p4z.Pt(),legacyEventWeight);
          h_legacy_jetpt->Fill(legacyJet.Pt(),legacyEventWeight);
          // Reproduce the reference phiBB convention, including the
          // orientation-dependent forward spike veto used in that code.
          const double legacyPhiBB =
            fabs(legacyJet.DeltaPhi(p4z)-TMath::Pi());
          const double legacyDphiResidual =
            std::min(legacyPhiBB,2.*TMath::Pi()-legacyPhiBB);
          h_legacy_dphi->Fill(legacyDphiResidual,legacyEventWeight);
          const bool passesBackToBack =
            !(legacyPhiBB>0.44 &&
              legacyPhiBB<2.*TMath::Pi()-0.44);
          const bool rejectedAsSpike =
            (legacyJet.Pt()<70. && fabs(legacyJet.Eta())>2.65 &&
             fabs(legacyJet.Eta())<2.964 && legacyPhiBB<2.7);
          if (passesBackToBack && passesJetVetoMap(legacyJet) &&
              !rejectedAsSpike) {
            h_legacy_cutflow->Fill(8.,legacyEventWeight);
            const double ptz = p4z.Pt();
            const double ptj = legacyJet.Pt();
            const double ptave = 0.5*(ptz+ptj);
            const double alpha = (legacySubleadingJet.Pt()>=15.)
              ? legacySubleadingJet.Pt()/ptz : 0.;
            TLorentzVector legacyMet1 = -p4z-legacyJet;
            legacyMet1.SetPtEtaPhiM(legacyMet1.Pt(),0.,legacyMet1.Phi(),0.);
            TLorentzVector legacyMetn = -legacyHt+p4z+legacyJet;
            legacyMetn.SetPtEtaPhiM(legacyMetn.Pt(),0.,legacyMetn.Phi(),0.);
            const TLorentzVector legacyMetnu = legacyMetn+legacyMetu;
            const TLorentzVector legacyMeta =
              legacyMet1+legacyMetn+legacyMetu;
            const double denominator = ptz*ptz;
            const double db = ptj/ptz;
            const double mpf =
              1.+legacyMeta.Vect().Dot(p4z.Vect())/denominator;
            const double mpf1 =
              1.+legacyMet1.Vect().Dot(p4z.Vect())/denominator;
            const double mpfn =
              legacyMetn.Vect().Dot(p4z.Vect())/denominator;
            const double mpfu =
              legacyMetu.Vect().Dot(p4z.Vect())/denominator;
            const double mpfnu =
              legacyMetnu.Vect().Dot(p4z.Vect())/denominator;
            h_legacy_alpha->Fill(alpha,legacyEventWeight);
            h_legacy_subleading_jetpt->Fill(
              legacySubleadingJet.Pt(),legacyEventWeight);
            h_legacy_alpha_vs_jetpt->Fill(ptj,alpha,legacyEventWeight);
            p_legacy_db_before_alpha->Fill(ptj,db,legacyEventWeight);
            p_legacy_mpf1_before_alpha->Fill(ptj,mpf1,legacyEventWeight);
            p_legacy_rho_before_alpha->Fill(
              ptz,Rho_fixedGridRhoFastjetAll,legacyEventWeight);
            for (const auto &entry : legacyRhoAlpha)
              if (alpha<0.01*entry.first)
                entry.second->Fill(
                  ptz,Rho_fixedGridRhoFastjetAll,legacyEventWeight);

            if (alpha<1.) {
              h_legacy_cutflow->Fill(9.,legacyEventWeight);
              bool hasGenResponse = false;
              double genResponse = 0.;
              int legacyGeneratorJetIndex = -1;
              if (isMC && Jet_genJetIdx[legacyJetIndex]>=0 &&
                  Jet_genJetIdx[legacyJetIndex]<nGenJet) {
                const int igen = Jet_genJetIdx[legacyJetIndex];
                legacyGeneratorJetIndex = igen;
                TLorentzVector generatorJet;
                generatorJet.SetPtEtaPhiM(
                  GenJet_pt[igen],GenJet_eta[igen],GenJet_phi[igen],
                  GenJet_mass[igen]);
                genResponse =
                  -(generatorJet.Px()*p4z.Px()+
                    generatorJet.Py()*p4z.Py())/denominator;
                hasGenResponse = std::isfinite(genResponse);
              }
              const double absLegacyEta = fabs(legacyJet.Eta());
              const double inverseStoredJec =
                jetStoredInverseJec[legacyJetIndex];
              legacyH2PtEta->Fill(
                absLegacyEta,ptz,legacyEventWeight);
              legacyH2PtEtaPf->Fill(
                absLegacyEta,ptj,legacyEventWeight);
              legacyH2PtEtaTc->Fill(
                absLegacyEta,ptz,legacyEventWeight);
              legacyP2Jes->Fill(
                absLegacyEta,ptz,inverseStoredJec,legacyEventWeight);
              legacyP2JesPf->Fill(
                absLegacyEta,ptj,inverseStoredJec,legacyEventWeight);
              legacyP2JesTc->Fill(
                absLegacyEta,ptz,inverseStoredJec,legacyEventWeight);
              legacyP2Res->Fill(
                absLegacyEta,ptz,jetInverseResidual[legacyJetIndex],
                legacyEventWeight);
              legacyP2ResPf->Fill(
                absLegacyEta,ptj,jetInverseResidual[legacyJetIndex],
                legacyEventWeight);
              legacyP2ResTc->Fill(
                absLegacyEta,ptz,jetInverseResidual[legacyJetIndex],
                legacyEventWeight);
              if (absLegacyEta>0. && absLegacyEta<1.3) {
                const TLorentzVector &generatorAxis =
                  generatorRecoil.hasGeneratorZ ? generatorRecoil.z : p4z;
                const GeneratorPairComponents generatorPair =
                  generatorPairComponents(
                    legacyGeneratorJetIndex,legacyJet,p4z,generatorAxis,
                    false);
                fillTruthHDMProfiles1D(
                  legacyTruthProfiles,"parallel",ptz,ptj,ptave,
                  hasGenResponse,generatorRecoil,generatorPair,mpf1,mpfn,
                  mpfu,jetInverseResidual[legacyJetIndex],legacyEventWeight);
                fillResponseProfiles1D(
                  legacyProfiles,ptz,ptj,ptave,db,mpf,mpf1,mpfn,mpfu,mpfnu,
                  jetInverseResidual[legacyJetIndex],
                  Jet_chHEF[legacyJetIndex],Jet_neEmEF[legacyJetIndex],
                  Jet_neHEF[legacyJetIndex],Jet_chEmEF[legacyJetIndex],
                  Jet_muEF[legacyJetIndex],Rho_fixedGridRhoFastjetAll,
                  hasGenResponse,genResponse,legacyEventWeight);
                legacyProfiles.axes.at("zmmjet").mass->Fill(
                  ptz,p4z.M(),legacyEventWeight);
                legacyProfiles.axes.at("jetpt").mass->Fill(
                  ptj,p4z.M(),legacyEventWeight);
                legacyProfiles.axes.at("ptave").mass->Fill(
                  ptave,p4z.M(),legacyEventWeight);
              }
            }
          }
        }
      }

      // Set Z-parallel (probe) directions
      p4p.SetPtEtaPhiM(p4z.Pt(),p4z.Eta(),p4z.Phi()+TMath::Pi(),p4z.M());
      p4pz.SetPtEtaPhiM(p4z.Pt(),p4z.Eta(),p4z.Phi(),p4z.M());
	
      // Use both transverse sidebands. Their later weights are one half each,
      // so their average has the same azimuthal acceptance as the signal.
      p4t1.SetPtEtaPhiM(p4z.Pt(),p4z.Eta(),p4z.Phi()+TMath::Pi()*0.5,p4z.M());
      p4t1z.SetPtEtaPhiM(p4z.Pt(),p4z.Eta(),p4z.Phi()-TMath::Pi()*0.5,p4z.M());
      p4t2.SetPtEtaPhiM(p4z.Pt(),p4z.Eta(),p4z.Phi()-TMath::Pi()*0.5,p4z.M());
      p4t2z.SetPtEtaPhiM(p4z.Pt(),p4z.Eta(),p4z.Phi()+TMath::Pi()*0.5,p4z.M());

      // Treat (+90) and (-90) as independent signal-sideband pairs. A lepton
      // veto in one transverse window removes only that window and the
      // matching half of the parallel signal. This keeps the lepton-veto
      // acceptance identical in each subtraction and avoids penalising the
      // signal twice merely because two sidebands are evaluated.
      const double vetoWidth = TMath::Pi()/8.;
      auto leptonClear = [&](const TLorentzVector &probe) {
        return (fabs(probe.DeltaPhi(p4lplus))>=vetoWidth &&
                fabs(probe.DeltaPhi(p4lminus))>=vetoWidth);
      };
      const bool signalClear = leptonClear(p4p);
      const bool pairValid[] = {signalClear && leptonClear(p4t1),
                                signalClear && leptonClear(p4t2)};
      const int pairState = (pairValid[0] ? 1 : 0)+(pairValid[1] ? 2 : 0);
      const double signalAcceptance =
        0.5*((pairValid[0] ? 1. : 0.)+(pairValid[1] ? 1. : 0.));
      h_probe_pair_state->Fill(pairState,eventWeight);
      if (signalClear) h_probe_veto->Fill(2.,eventWeight);
      if (pairValid[0]) h_probe_veto->Fill(3.,eventWeight);
      if (pairValid[1]) h_probe_veto->Fill(4.,eventWeight);
      h_probe_veto->Fill(5.,eventWeight*signalAcceptance);
      if (signalAcceptance==0.) continue;
      h_cutflow->Fill(7., eventWeight*signalAcceptance);
      h_zpt_probeveto->Fill(p4z.Pt(), eventWeight*signalAcceptance);
      
      // Calculate MET and HT sum
      if (recalculatePuppiMet)
        met.SetPtEtaPhiM(RawPuppiMET_pt,0.,RawPuppiMET_phi,0.);
      else
        met.SetPtEtaPhiM(PuppiMET_pt,0.,PuppiMET_phi,0.);
      ht += p4z;
      for (int ijet = 0; ijet != nJet; ++ijet) {
	p4jet.SetPtEtaPhiM(Jet_pt[ijet], Jet_eta[ijet], Jet_phi[ijet],
			   Jet_mass[ijet]);
	//if (p4jet.DeltaR(p4lplus)>0.4 && p4jet.DeltaR(p4lminus)>0.4 &&
	if (separatedFromSynchronizedMuons(p4jet) && p4jet.Pt()>15.) {
	  ht += p4jet;
	  if (recalculatePuppiMet) {
	    rawjet.SetPtEtaPhiM(jetRawPt[ijet],Jet_eta[ijet],Jet_phi[ijet],
	                         jetRawMass[ijet]);
	    rawjets += rawjet;
	    corrjets += p4jet;
	  }
	}
      }
      ht.SetPtEtaPhiM(ht.Pt(),0,ht.Phi(),0);
      if (recalculatePuppiMet) met += rawjets-corrjets;
      met.SetPtEtaPhiM(met.Pt(),0,met.Phi(),0.);
      metu = met + ht;

      auto fillPileupResponse = [&](const char *region, double db,
                                    double mpfValue, double weight) {
        const double x[] = {double(PV_npvs),
                            double(Rho_fixedGridRhoFastjetAll), mu};
        for (int io = 0; io != 3; ++io) {
          if (x[io] < 0.) continue;
          pileupControl[Form("p_db_vs_%s_%s",observables[io],region)]
            ->Fill(x[io], db, weight);
          pileupControl[Form("p_mpf_vs_%s_%s",observables[io],region)]
            ->Fill(x[io], mpfValue, weight);
        }
      };

      auto fillTruth = [&](const char *region, bool matched,
                           bool passesExtraMatchCuts, bool hasGenIndex,
                           int matchCategory, double db, double mpfValue,
                           double weight) {
        if (!isMC) return;
        TH1D *hist = (string(region)=="parallel" ? h_truth_parallel :
                      string(region)=="transverse" ? h_truth_transverse :
                      h_truth_subtracted);
        hist->Fill(matched ? 1. : 0., weight);

        const double x[] = {p4z.Pt(), double(PV_npvs),
                            double(Rho_fixedGridRhoFastjetAll), mu};
        const char *names[] = {"ptz", "npvs", "rho", "mu"};
        for (int io = 0; io != 4; ++io) {
          if (x[io] < 0.) continue;
          truthControl[Form("p_pujet_fraction_vs_%s_%s",names[io],region)]
            ->Fill(x[io], matched ? 0. : 1., weight);
        }

        truthControl[Form("p_no_gen_index_fraction_vs_ptz_%s",region)]
          ->Fill(p4z.Pt(),hasGenIndex ? 0. : 1.,weight);
        truthControl[Form("p_extra_match_cuts_unmatched_fraction_vs_ptz_%s",
                          region)]
          ->Fill(p4z.Pt(),passesExtraMatchCuts ? 0. : 1.,weight);
        const double abseta = fabs(p4jet.Eta());
        const char *etaRegion = (abseta<1.3 ? "central" :
                                 abseta<2.5 ? "endcap" : "forward");
        truthControl[Form("p_pujet_fraction_vs_ptz_%s_%s",region,etaRegion)]
          ->Fill(p4z.Pt(),matched ? 0. : 1.,weight);
        truthYield[Form("h_all_jets_vs_ptz_%s",region)]
          ->Fill(p4z.Pt(),weight);
        if (!matched)
          truthYield[Form("h_unmatched_jets_vs_ptz_%s",region)]
            ->Fill(p4z.Pt(),weight);
        if (!hasGenIndex)
          truthYield[Form("h_no_gen_index_jets_vs_ptz_%s",region)]
            ->Fill(p4z.Pt(),weight);
        if (!passesExtraMatchCuts)
          truthYield[Form("h_extra_cuts_unmatched_jets_vs_ptz_%s",region)]
            ->Fill(p4z.Pt(),weight);
        truthMatchQuality[Form("h2_truth_match_quality_vs_ptz_%s",region)]
          ->Fill(p4z.Pt(),matchCategory,weight);

        if (string(region)!="subtracted") {
          const char *category = (matched ? "matched" : "pileup");
          truthControl[Form("p_db_vs_ptz_%s_%s",category,region)]
            ->Fill(p4z.Pt(), db, weight);
          truthControl[Form("p_mpf_vs_ptz_%s_%s",category,region)]
            ->Fill(p4z.Pt(), mpfValue, weight);
        }
      };

      auto fillMethodTruth = [&](const char *region, bool matched,
                                 double db, double mpfValue,
                                 double mpfJet, double mpfNeutral,
                                 double mpfUnclustered, double mpfNu,
                                 double weight) {
        if (!isMC) return;
        const double values[] = {db,mpfValue,mpfJet,mpfNeutral,
                                 mpfUnclustered,mpfNu};
        const char *responses[] = {"db","mpf","mpf1","mpfn","mpfu",
                                   "mpfnu"};
        const char *truthCategory = matched ? "matched" : "pileup";
        const double absEta = fabs(p4jet.Eta());
        for (int index = 0; index != 6; ++index) {
          if (absEta<1.305) {
            methodTruthControl[Form("p_%s_vs_ptz_all_%s_central",
                                    responses[index],region)]
              ->Fill(p4z.Pt(),values[index],weight);
            methodTruthControl[Form("p_%s_vs_ptz_%s_%s_central",
                                    responses[index],truthCategory,region)]
              ->Fill(p4z.Pt(),values[index],weight);
          }
          if (p4z.Pt()>15. && p4z.Pt()<30.) {
            methodTruthControl[Form("p_%s_vs_abseta_all_%s_ptz15to30",
                                    responses[index],region)]
              ->Fill(absEta,values[index],weight);
            methodTruthControl[Form("p_%s_vs_abseta_%s_%s_ptz15to30",
                                    responses[index],truthCategory,region)]
              ->Fill(absEta,values[index],weight);
          }
        }
      };

      auto fillFlavor = [&](int ijet, int generatorFlavor,
                            bool hasGeneratorResponse,
                            double generatorResponse, double db,
                            double mpfValue, double mpfJet,
                            double mpfNeutral, double mpfUnclustered,
                            double weight) {
        if (fabs(Jet_eta[ijet])>=1.305) return;

        // Keep these definitions synchronized with Bettina's gamma+jet
        // flavor analysis.  The b and c selections have priority over Q/G.
        const double bScore = Jet_btagDeepFlavB[ijet];
        const double cScore =
          0.5*(Jet_btagDeepFlavCvB[ijet]+Jet_btagDeepFlavCvL[ijet]);
        const double qgScore = Jet_btagDeepFlavQG[ijet];
        const bool isB = (bScore>0.7527);
        const bool isC = (cScore>0.3985 && !isB);
        const bool isQ = (qgScore>=0.5 && !isB && !isC);
        const bool isG = (qgScore>=0. && qgScore<0.5 && !isB && !isC);
        const std::string recoTag =
          (isB ? "b" : isC ? "c" : isQ ? "q" : isG ? "g" : "n");

        std::string genTag = "n";
        if (std::abs(generatorFlavor)==5) genTag = "b";
        else if (std::abs(generatorFlavor)==4) genTag = "c";
        else if (std::abs(generatorFlavor)>=1 &&
                 std::abs(generatorFlavor)<=3) genTag = "q";
        else if (generatorFlavor==21) genTag = "g";

        const std::vector<std::string> selectedRecoTags = {"i",recoTag};
        std::vector<std::string> selectedGenTags = {"i"};
        if (isMC) selectedGenTags.push_back(genTag);
        for (const std::string &selectedReco : selectedRecoTags) {
          for (const std::string &selectedGen : selectedGenTags) {
            dynamic_cast<TH1D*>(
              flavorProfiles["counts"][selectedReco][selectedGen])
              ->Fill(p4z.Pt(),weight);
            dynamic_cast<TProfile*>(
              flavorProfiles["mpfchs1"][selectedReco][selectedGen])
              ->Fill(p4z.Pt(),mpfValue,weight);
            dynamic_cast<TProfile*>(
              flavorProfiles["ptchs"][selectedReco][selectedGen])
              ->Fill(p4z.Pt(),db,weight);
            dynamic_cast<TProfile*>(
              flavorProfiles["mpf1"][selectedReco][selectedGen])
              ->Fill(p4z.Pt(),mpfJet,weight);
            dynamic_cast<TProfile*>(
              flavorProfiles["mpfn"][selectedReco][selectedGen])
              ->Fill(p4z.Pt(),mpfNeutral,weight);
            dynamic_cast<TProfile*>(
              flavorProfiles["mpfu"][selectedReco][selectedGen])
              ->Fill(p4z.Pt(),mpfUnclustered,weight);
            dynamic_cast<TProfile*>(
              flavorProfiles["rjet"][selectedReco][selectedGen])
              ->Fill(p4z.Pt(),db,weight);
            if (hasGeneratorResponse)
              dynamic_cast<TProfile*>(
                flavorProfiles["gjet"][selectedReco][selectedGen])
                ->Fill(p4z.Pt(),generatorResponse,weight);
          }
        }
      };
      
      // Select leading jet
      h_njet->Fill(nJet, eventWeight);
      for (int ijet = 0; ijet != nJet; ++ijet) {

	p4jet.SetPtEtaPhiM(Jet_pt[ijet], Jet_eta[ijet], Jet_phi[ijet],
			   Jet_mass[ijet]);
	if (!passTightJetId(ijet)) continue;
	if (!isMC) {
	  h_jet_veto_map->Fill(1.,eventWeight);
	  if (!passesJetVetoMap(p4jet)) {
	    h_jet_veto_map->Fill(3.,eventWeight);
	    continue;
	  }
	  h_jet_veto_map->Fill(2.,eventWeight);
	}

	double eta = p4jet.Eta();
	double ptz = p4z.Pt();
	double ptj = p4jet.Pt();
	double pta = 0.5*(ptz+ptj);

	double jes = (1-Jet_rawFactor[ijet]);
	bool hasGenIndex = false;
	bool truthMatched = false;
	bool passesExtraMatchCuts = false;
	int truthMatchCategory = 0;
	bool hasGenResponse = false;
	double genBalance = 0.;
	int genJetIndex = -1;
	if (isMC && Jet_genJetIdx[ijet]>=0 && Jet_genJetIdx[ijet]<nGenJet) {
	  hasGenIndex = true;
	  // Jet_genJetIdx is the NanoAOD reco-to-particle-level match. Do not
	  // impose a second generator-pT cut on the nominal pileup classification:
	  // it creates a strong migration bias precisely in the low-pT region.
	  const int igen = Jet_genJetIdx[ijet];
	  genJetIndex = igen;
	  TLorentzVector p4gen;
	  p4gen.SetPtEtaPhiM(GenJet_pt[igen],GenJet_eta[igen],GenJet_phi[igen],
			     GenJet_mass[igen]);
	  truthMatched = true;
	  genBalance = -(p4gen.Px()*p4z.Px()+p4gen.Py()*p4z.Py())/
	               (ptz*ptz);
	  hasGenResponse = std::isfinite(genBalance);
	  if (p4gen.Pt()<=8.) truthMatchCategory = 1;
	  else if (p4jet.DeltaR(p4gen)>=0.4) truthMatchCategory = 2;
	  else {
	    truthMatchCategory = 3;
	    passesExtraMatchCuts = true;
	  }
	}
	const double upartCvB = Jet_btagUParTAK4CvB[ijet];
	const double upartCvL = Jet_btagUParTAK4CvL[ijet];
	const double upartQvG = Jet_btagUParTAK4QvG[ijet];
	const int recoUParTFlavor = reconstructedUParTFlavor(
	  upartCvB,upartCvL,upartQvG);
	const int truePartonFlavor =
	  (isMC ? generatorFlavorId(Jet_partonFlavour[ijet]) : 0);

	met1 = -p4z - p4jet;
	met1.SetPtEtaPhiM(met1.Pt(),0,met1.Phi(),0.);
	metn = -ht + p4z + p4jet;
	metn.SetPtEtaPhiM(metn.Pt(),0,metn.Phi(),0.);
	metnu = metn + metu;
	meta = met1 + metn + metu;
	double mpf = 1 + meta.Vect().Dot(p4z.Vect()) / (ptz*ptz);
	double mpf1 = 1 + met1.Vect().Dot(p4z.Vect()) / (ptz*ptz);
	double mpfn = metn.Vect().Dot(p4z.Vect()) / (ptz*ptz);
	double mpfu = metu.Vect().Dot(p4z.Vect()) / (ptz*ptz);
	double mpfnu = metnu.Vect().Dot(p4z.Vect()) / (ptz*ptz);
	
	//if (p4jet.DeltaR(p4lplus)>0.4 && p4jet.DeltaR(p4lminus)>0.4) {
	if (separatedFromSynchronizedMuons(p4jet)) {

	  h_jetpt->Fill(p4jet.Pt());
	  h_jeteta->Fill(p4jet.Eta());
	  if (p4jet.Pt()>p4jet1.Pt()) p4jet1 = p4jet;
	  
	  // Parallel region
	  //if (fabs(p4jet.DeltaPhi(p4z))>3./4.*TMath::Pi() && // >2.36
	  //if (fabs(p4jet.DeltaPhi(p4z))>7./8.*TMath::Pi() && // >2.75
	  //if (fabs(p4jet.DeltaPhi(p4z))>15./16.*TMath::Pi() && // >2.945
	  if (fabs(p4jet.DeltaPhi(p4p))<1./16.*TMath::Pi() && // 0.1963
	      //p4jet.Pt()>0.5*p4z.Pt() && p4z.Pt()>0.5*p4jet.Pt()) {
	      //p4jet.Pt()>0.6*p4z.Pt() && p4z.Pt()>0.6*p4jet.Pt()) {
	      //p4jet.Pt()>0.5*p4z.Pt() && p4jet.Pt()<1.5*p4z.Pt()) {
	      //p4jet.Pt()>0.25*p4z.Pt() && p4jet.Pt()<2.0*p4z.Pt()) {
		      p4jet.Pt()>0.5*p4z.Pt() && p4jet.Pt()<2.0*p4z.Pt()) {

		    const double db = ptj/ptz;
		    const double abseta = fabs(eta);
		    const double wt = eventWeight*signalAcceptance;
		    nsel += signalAcceptance;
		    h2_mixpteta->Fill(ptj,eta,wt);
		    h_mixpt->Fill(ptj,wt);
		    h_mixeta->Fill(eta,wt);
		    h_mixdphi->Fill(p4jet.DeltaPhi(p4p),wt);
		    h2_parpteta->Fill(ptj,eta,wt);
		    h_parpt->Fill(ptj,wt);
		    h_pareta->Fill(eta,wt);
		    h_pardphi->Fill(p4jet.DeltaPhi(p4p),wt);
		    if (ptj>p4sel1.Pt()) p4sel1 = p4jet;

		    h_dbp->Fill(db,wt); h_mpfp->Fill(mpf,wt);
		    p2_dbp->Fill(ptz,eta,db,wt); p2_mpfp->Fill(ptz,eta,mpf,wt);
		    p2_mpfnp->Fill(ptz,eta,mpfn,wt); p2_mpfup->Fill(ptz,eta,mpfu,wt);
		    p2_mpfnup->Fill(ptz,eta,mpfnu,wt); h2_dbp->Fill(ptz,ptj,wt);
		    p_dbp_vsz->Fill(ptz,db,wt); p_dbp_vsj->Fill(ptj,db,wt);
		    p_dbp_vsa->Fill(pta,db,wt);

		    h_selpt->Fill(ptj,wt); h_seleta->Fill(eta,wt);
		    h_seldphi->Fill(p4jet.DeltaPhi(p4p),wt);
		    h_db->Fill(db,wt); h_mpf->Fill(mpf,wt);
		    p2_db->Fill(ptz,eta,db,wt); p2_mpf->Fill(ptz,eta,mpf,wt);
		    p2_mpfn->Fill(ptz,eta,mpfn,wt); p2_mpfu->Fill(ptz,eta,mpfu,wt);
		    p2_mpfnu->Fill(ptz,eta,mpfnu,wt); h2_db->Fill(ptz,ptj,wt);
		    p_db_vsz->Fill(ptz,db,wt); p_db_vsj->Fill(ptj,db,wt);
		    p_db_vsa->Fill(pta,db,wt);
		    fillPileupResponse("parallel",db,mpf,wt);
		    fillPileupResponse("subtracted",db,mpf,wt);
		    fillTruth("parallel",truthMatched,passesExtraMatchCuts,
		              hasGenIndex,truthMatchCategory,db,mpf,wt);
		    fillTruth("subtracted",truthMatched,passesExtraMatchCuts,
		              hasGenIndex,truthMatchCategory,db,mpf,wt);
		    fillMethodTruth("parallel",truthMatched,db,mpf,mpf1,mpfn,
		                    mpfu,mpfnu,wt);
		    fillMethodTruth("subtracted",truthMatched,db,mpf,mpf1,mpfn,
		                    mpfu,mpfnu,wt);
		    fillFlavor(ijet,
		               (isMC && genJetIndex>=0
		                  ? GenJet_partonFlavour[genJetIndex] : 0),
		               hasGenResponse,genBalance,db,mpf,mpf1,mpfn,mpfu,wt);
		    if (abseta<1.305) {
		      if (abseta<1.3)
		        fillFlavorMatrix(
		          flavorMatrix,ptz,ptj,recoUParTFlavor,truePartonFlavor,
		          upartCvB,upartCvL,upartQvG,
		          flavorMatrixVariants(
		            p4z,p4jet,mpf,mpf1,mpfn,mpfu,mpfnu,false),
		          wt,false);
		      const TLorentzVector &generatorAxis =
		        generatorRecoil.hasGeneratorZ ? generatorRecoil.z : p4z;
		      const GeneratorPairComponents generatorPair =
		        generatorPairComponents(
		          genJetIndex,p4jet,p4z,generatorAxis,false);
		      fillTruthHDMProfiles1D(
		        allPairsTruthProfiles,"parallel",ptz,ptj,pta,truthMatched,
		        generatorRecoil,generatorPair,mpf1,mpfn,mpfu,
		        jetInverseResidual[ijet],wt);
		      fillTruthHDMProfiles1D(
		        allPairsTruthProfiles,"subtracted",ptz,ptj,pta,truthMatched,
		        generatorRecoil,generatorPair,mpf1,mpfn,mpfu,
		        jetInverseResidual[ijet],wt);
		      fillResponseProfiles1D(
		        allPairsProfiles,ptz,ptj,pta,db,mpf,mpf1,mpfn,mpfu,mpfnu,
		        jetInverseResidual[ijet],Jet_chHEF[ijet],Jet_neEmEF[ijet],
		        Jet_neHEF[ijet],Jet_chEmEF[ijet],Jet_muEF[ijet],
		        Rho_fixedGridRhoFastjetAll,hasGenResponse,genBalance,wt);
		      allPairsProfiles.axes.at("zmmjet").mass->Fill(ptz,p4z.M(),wt);
		      allPairsProfiles.axes.at("jetpt").mass->Fill(ptj,p4z.M(),wt);
		      allPairsProfiles.axes.at("ptave").mass->Fill(pta,p4z.M(),wt);
		    }

		    h2ptetapf_->Fill(abseta,ptj,wt); h2pteta_->Fill(abseta,pta,wt);
		    h2ptetatc_->Fill(abseta,ptz,wt);
		    h2ptetapf->Fill(eta,ptj,wt); h2pteta->Fill(eta,pta,wt);
		    h2ptetatc->Fill(eta,ptz,wt);
		    pmzpf_->Fill(ptj,p4z.M(),wt); pmz_->Fill(pta,p4z.M(),wt);
		    pmztc_->Fill(ptz,p4z.M(),wt); pmzpf->Fill(ptj,p4z.M(),wt);
		    pmz->Fill(pta,p4z.M(),wt); pmztc->Fill(ptz,p4z.M(),wt);
		    h2mztc_->Fill(ptz,p4z.M(),wt);

		    p2jespf_->Fill(abseta,ptj,jes,wt); p2jes_->Fill(abseta,pta,jes,wt);
		    p2jestc_->Fill(abseta,ptz,jes,wt); p2jespf->Fill(eta,ptj,jes,wt);
		    p2jes->Fill(eta,pta,jes,wt); p2jestc->Fill(eta,ptz,jes,wt);
		    p2m0pf_->Fill(abseta,ptj,mpf,wt); p2m0_->Fill(abseta,pta,mpf,wt);
		    p2m0tc_->Fill(abseta,ptz,mpf,wt); p2m0pf->Fill(eta,ptj,mpf,wt);
		    p2m0->Fill(eta,pta,mpf,wt); p2m0tc->Fill(eta,ptz,mpf,wt);
		    p2m2pf_->Fill(abseta,ptj,mpf1,wt); p2m2_->Fill(abseta,pta,mpf1,wt);
		    p2m2tc_->Fill(abseta,ptz,mpf1,wt); p2m2pf->Fill(eta,ptj,mpf1,wt);
		    p2m2->Fill(eta,pta,mpf1,wt); p2m2tc->Fill(eta,ptz,mpf1,wt);
		    p2mnpf_->Fill(abseta,ptj,mpfn,wt); p2mn_->Fill(abseta,pta,mpfn,wt);
		    p2mntc_->Fill(abseta,ptz,mpfn,wt); p2mnpf->Fill(eta,ptj,mpfn,wt);
		    p2mn->Fill(eta,pta,mpfn,wt); p2mntc->Fill(eta,ptz,mpfn,wt);
		    p2mnupf_->Fill(abseta,ptj,mpfnu,wt); p2mnu_->Fill(abseta,pta,mpfnu,wt);
		    p2mnutc_->Fill(abseta,ptz,mpfnu,wt); p2mnupf->Fill(eta,ptj,mpfnu,wt);
		    p2mnu->Fill(eta,pta,mpfnu,wt); p2mnutc->Fill(eta,ptz,mpfnu,wt);
		    p2mupf_->Fill(abseta,ptj,mpfu,wt); p2mu_->Fill(abseta,pta,mpfu,wt);
		    p2mutc_->Fill(abseta,ptz,mpfu,wt); p2mupf->Fill(eta,ptj,mpfu,wt);
		    p2mu->Fill(eta,pta,mpfu,wt); p2mutc->Fill(eta,ptz,mpfu,wt);
		    p2respf_->Fill(abseta,ptj,jetInverseResidual[ijet],wt);
		    p2res_->Fill(abseta,pta,jetInverseResidual[ijet],wt);
		    p2restc_->Fill(abseta,ptz,jetInverseResidual[ijet],wt);
		    p2respf->Fill(eta,ptj,jetInverseResidual[ijet],wt);
		    p2res->Fill(eta,pta,jetInverseResidual[ijet],wt);
		    p2restc->Fill(eta,ptz,jetInverseResidual[ijet],wt);
		    p2chftc_->Fill(abseta,ptz,Jet_chHEF[ijet],wt);
		    p2neftc_->Fill(abseta,ptz,Jet_neEmEF[ijet],wt);
		    p2nhftc_->Fill(abseta,ptz,Jet_neHEF[ijet],wt);
		    p2ceftc_->Fill(abseta,ptz,Jet_chEmEF[ijet],wt);
		    p2muftc_->Fill(abseta,ptz,Jet_muEF[ijet],wt);
		    p2rhotc_->Fill(abseta,ptz,Rho_fixedGridRhoFastjetAll,wt);
		    if (hasGenResponse) p2rgentc_->Fill(abseta,ptz,genBalance,wt);
		  } // Parallel region

		  const TLorentzVector *transverseProbe[] = {&p4t1,&p4t2};
		  const TLorentzVector *transverseAxis[] = {&p4t1z,&p4t2z};
		  for (int idir = 0; idir != 2; ++idir) {
		    if (!pairValid[idir]) continue;
		    const TLorentzVector &probe = *transverseProbe[idir];
		    const TLorentzVector &axis = *transverseAxis[idir];
		    if (fabs(p4jet.DeltaPhi(probe))>=1./16.*TMath::Pi() ||
			ptj<=0.5*ptz || ptj>=2.0*ptz) continue;

		    TLorentzVector met1t = -p4z-p4jet;
		    met1t.SetPtEtaPhiM(met1t.Pt(),0,met1t.Phi(),0.);
		    TLorentzVector metnt = -ht+p4z+p4jet;
		    metnt.SetPtEtaPhiM(metnt.Pt(),0,metnt.Phi(),0.);
		    const TLorentzVector metnut = metnt+metu;
		    const TLorentzVector metat = met1t+metnt+metu;
		    double mpfT = 1+metat.Vect().Dot(axis.Vect())/(ptz*ptz)+(mpf-1.);
		    double mpf1T = 1+met1t.Vect().Dot(axis.Vect())/(ptz*ptz)+(mpf1-1.);
		    double mpfnT = metnt.Vect().Dot(axis.Vect())/(ptz*ptz)+mpfn;
		    double mpfuT = metu.Vect().Dot(axis.Vect())/(ptz*ptz)+mpfu;
		    double mpfnuT = metnut.Vect().Dot(axis.Vect())/(ptz*ptz)+mpfnu;

		    const double db = ptj/ptz;
		    const double abseta = fabs(eta);
		    const double wraw = 0.5*eventWeight;
		    const double wt = -wraw;
		    ntran += 0.5;

		    h2_tranpteta->Fill(ptj,eta,wraw); h_tranpt->Fill(ptj,wraw);
		    h_traneta->Fill(eta,wraw); h_trandphi->Fill(p4jet.DeltaPhi(probe),wraw);
		    if (ptj>p4tran1.Pt()) p4tran1 = p4jet;
		    h_dbt->Fill(db,wraw); h_mpft->Fill(mpfT,wraw);
		    p2_dbt->Fill(ptz,eta,db,wraw); p2_mpft->Fill(ptz,eta,mpfT,wraw);
		    p2_mpfnt->Fill(ptz,eta,mpfnT,wraw); p2_mpfut->Fill(ptz,eta,mpfuT,wraw);
		    p2_mpfnut->Fill(ptz,eta,mpfnuT,wraw); h2_dbt->Fill(ptz,ptj,wraw);
		    p_dbt_vsz->Fill(ptz,db,wraw); p_dbt_vsj->Fill(ptj,db,wraw);
		    p_dbt_vsa->Fill(pta,db,wraw);

		    h2_mixpteta->Fill(ptj,eta,wt); h_mixpt->Fill(ptj,wt);
		    h_mixeta->Fill(eta,wt); h_mixdphi->Fill(p4jet.DeltaPhi(probe),wt);
		    h_selpt->Fill(ptj,wt); h_seleta->Fill(eta,wt);
		    h_seldphi->Fill(p4jet.DeltaPhi(probe),wt);
		    h_db->Fill(db,wt); h_mpf->Fill(mpfT,wt);
		    p2_db->Fill(ptz,eta,db,wt); p2_mpf->Fill(ptz,eta,mpfT,wt);
		    p2_mpfn->Fill(ptz,eta,mpfnT,wt); p2_mpfu->Fill(ptz,eta,mpfuT,wt);
		    p2_mpfnu->Fill(ptz,eta,mpfnuT,wt); h2_db->Fill(ptz,ptj,wt);
		    p_db_vsz->Fill(ptz,db,wt); p_db_vsj->Fill(ptj,db,wt);
		    p_db_vsa->Fill(pta,db,wt);
		    fillPileupResponse("transverse",db,mpfT,wraw);
		    fillPileupResponse("subtracted",db,mpfT,wt);
		    fillTruth("transverse",truthMatched,passesExtraMatchCuts,
		              hasGenIndex,truthMatchCategory,db,mpfT,wraw);
		    fillTruth("subtracted",truthMatched,passesExtraMatchCuts,
		              hasGenIndex,truthMatchCategory,db,mpfT,wt);
		    fillMethodTruth("transverse",truthMatched,db,mpfT,mpf1T,
		                    mpfnT,mpfuT,mpfnuT,wraw);
		    fillMethodTruth("subtracted",truthMatched,db,mpfT,mpf1T,
		                    mpfnT,mpfuT,mpfnuT,wt);
		    fillFlavor(ijet,
		               (isMC && genJetIndex>=0
		                  ? GenJet_partonFlavour[genJetIndex] : 0),
		               hasGenResponse,genBalance,db,mpfT,mpf1T,mpfnT,
		               mpfuT,wt);
		    if (abseta<1.305) {
		      if (abseta<1.3)
		        fillFlavorMatrix(
		          flavorMatrix,ptz,ptj,recoUParTFlavor,truePartonFlavor,
		          upartCvB,upartCvL,upartQvG,
		          flavorMatrixVariants(
		            p4z,p4jet,mpfT,mpf1T,mpfnT,mpfuT,mpfnuT,true),
		          wt,true);
		      TLorentzVector generatorTransverseAxis = axis;
		      if (generatorRecoil.hasGeneratorZ)
		        generatorTransverseAxis.SetPtEtaPhiM(
		          generatorRecoil.z.Pt(),generatorRecoil.z.Eta(),
		          generatorRecoil.z.Phi()+
		            (idir==0 ? -0.5*TMath::Pi() : 0.5*TMath::Pi()),
		          generatorRecoil.z.M());
		      const GeneratorPairComponents generatorPair =
		        generatorPairComponents(
		          genJetIndex,p4jet,axis,generatorTransverseAxis,true);
		      fillTruthHDMProfiles1D(
		        allPairsTruthProfiles,"transverse",ptz,ptj,pta,truthMatched,
		        generatorRecoil,generatorPair,mpf1T,mpfnT,mpfuT,
		        jetInverseResidual[ijet],wraw);
		      fillTruthHDMProfiles1D(
		        allPairsTruthProfiles,"subtracted",ptz,ptj,pta,truthMatched,
		        generatorRecoil,generatorPair,mpf1T,mpfnT,mpfuT,
		        jetInverseResidual[ijet],wt);
		      fillResponseProfiles1D(
		        allPairsProfiles,ptz,ptj,pta,db,mpfT,mpf1T,mpfnT,mpfuT,
		        mpfnuT,jetInverseResidual[ijet],Jet_chHEF[ijet],
		        Jet_neEmEF[ijet],Jet_neHEF[ijet],Jet_chEmEF[ijet],
		        Jet_muEF[ijet],Rho_fixedGridRhoFastjetAll,hasGenResponse,
		        genBalance,wt);
		    }

		    h2ptetapf_->Fill(abseta,ptj,wt); h2pteta_->Fill(abseta,pta,wt);
		    h2ptetatc_->Fill(abseta,ptz,wt); h2ptetapf->Fill(eta,ptj,wt);
		    h2pteta->Fill(eta,pta,wt); h2ptetatc->Fill(eta,ptz,wt);
		    p2jespf_->Fill(abseta,ptj,jes,wt); p2jes_->Fill(abseta,pta,jes,wt);
		    p2jestc_->Fill(abseta,ptz,jes,wt); p2jespf->Fill(eta,ptj,jes,wt);
		    p2jes->Fill(eta,pta,jes,wt); p2jestc->Fill(eta,ptz,jes,wt);
		    p2m0pf_->Fill(abseta,ptj,mpfT,wt); p2m0_->Fill(abseta,pta,mpfT,wt);
		    p2m0tc_->Fill(abseta,ptz,mpfT,wt); p2m0pf->Fill(eta,ptj,mpfT,wt);
		    p2m0->Fill(eta,pta,mpfT,wt); p2m0tc->Fill(eta,ptz,mpfT,wt);
		    p2m2pf_->Fill(abseta,ptj,mpf1T,wt); p2m2_->Fill(abseta,pta,mpf1T,wt);
		    p2m2tc_->Fill(abseta,ptz,mpf1T,wt); p2m2pf->Fill(eta,ptj,mpf1T,wt);
		    p2m2->Fill(eta,pta,mpf1T,wt); p2m2tc->Fill(eta,ptz,mpf1T,wt);
		    p2mnpf_->Fill(abseta,ptj,mpfnT,wt); p2mn_->Fill(abseta,pta,mpfnT,wt);
		    p2mntc_->Fill(abseta,ptz,mpfnT,wt); p2mnpf->Fill(eta,ptj,mpfnT,wt);
		    p2mn->Fill(eta,pta,mpfnT,wt); p2mntc->Fill(eta,ptz,mpfnT,wt);
		    p2mnupf_->Fill(abseta,ptj,mpfnuT,wt); p2mnu_->Fill(abseta,pta,mpfnuT,wt);
		    p2mnutc_->Fill(abseta,ptz,mpfnuT,wt); p2mnupf->Fill(eta,ptj,mpfnuT,wt);
		    p2mnu->Fill(eta,pta,mpfnuT,wt); p2mnutc->Fill(eta,ptz,mpfnuT,wt);
		    p2mupf_->Fill(abseta,ptj,mpfuT,wt); p2mu_->Fill(abseta,pta,mpfuT,wt);
		    p2mutc_->Fill(abseta,ptz,mpfuT,wt); p2mupf->Fill(eta,ptj,mpfuT,wt);
		    p2mu->Fill(eta,pta,mpfuT,wt); p2mutc->Fill(eta,ptz,mpfuT,wt);
		    p2respf_->Fill(abseta,ptj,jetInverseResidual[ijet],wt);
		    p2res_->Fill(abseta,pta,jetInverseResidual[ijet],wt);
		    p2restc_->Fill(abseta,ptz,jetInverseResidual[ijet],wt);
		    p2respf->Fill(eta,ptj,jetInverseResidual[ijet],wt);
		    p2res->Fill(eta,pta,jetInverseResidual[ijet],wt);
		    p2restc->Fill(eta,ptz,jetInverseResidual[ijet],wt);
		    p2chftc_->Fill(abseta,ptz,Jet_chHEF[ijet],wt);
		    p2neftc_->Fill(abseta,ptz,Jet_neEmEF[ijet],wt);
		    p2nhftc_->Fill(abseta,ptz,Jet_neHEF[ijet],wt);
		    p2ceftc_->Fill(abseta,ptz,Jet_chEmEF[ijet],wt);
		    p2muftc_->Fill(abseta,ptz,Jet_muEF[ijet],wt);
		    p2rhotc_->Fill(abseta,ptz,Rho_fixedGridRhoFastjetAll,wt);
		    if (hasGenResponse) p2rgentc_->Fill(abseta,ptz,genBalance,wt);
		  } // transverse direction
	  
	}
      } // for ijet
      if (p4jet1.Pt()>0) {
	h_jet1pt->Fill(p4jet1.Pt());
	h_jet1eta->Fill(p4jet1.Eta());
      }

      h_nsel->Fill(nsel);
      if (p4sel1.Pt()>0) {
	h_sel1pt->Fill(p4sel1.Pt());
	h_sel1eta->Fill(p4sel1.Eta());
      }
      h_ntran->Fill(ntran);
      if (p4tran1.Pt()>0) {
	h_tran1pt->Fill(p4tran1.Pt());
	h_tran1eta->Fill(p4tran1.Eta());
      }
	   } // for jentry in nentries

   if (readFailure) {
     fout->Close();
     gSystem->Unlink(outputFile.c_str());
     cout << "Removed incomplete output " << outputFile << "." << endl;
     return;
   }
	   writeTruthHDMDerivedGraphs(allPairsTruthProfiles);
	   writeTruthHDMDerivedGraphs(legacyTruthProfiles);
	   for (const char *variant : {"ab","ad","tc","pf"}) {
	     ZJetFlavorMatrix::finalizeHDMProfile(
	       flavorMatrix.profiles.at(std::string("p3hdm")+variant+
	                                "_flavormatrix"),
	       flavorMatrix.profiles.at(std::string("p3m0")+variant+
	                                "_flavormatrix"),
	       flavorMatrix.profiles.at(std::string("p3mn")+variant+
	                                "_flavormatrix"),
	       flavorMatrix.profiles.at(std::string("p3mu")+variant+
	                                "_flavormatrix"));
	   }

	   cout << endl << "Finished loop, writing file " << outputFile << "." << endl << flush;
    cout << "Processed " << nentries << " events\n";
    //cout << "Skipped " << _nbadevents_json << " events due to JSON ("
    //	 << (100.*_nbadevents_json/_nevents) << "%) \n";
    //cout << "Skipped " << _nbadevents_trigger << " events due to trigger ("
    //	 << (100.*_nbadevents_trigger/_ntot) << "%) \n";
    //cout << "Skipped " << _nbadevents_veto << " events due to veto ("
    //	 << (100.*_nbadevents_veto/_nevents) << "%) \n";

   // Write the directory-owned histograms once, then write metadata only
   // after that pass. Writing TObjString objects before TFile::Write creates
   // another key cycle for every metadata item and multiplies those cycles
   // when job outputs are merged with hadd.
   fout->Write("",TObject::kOverwrite);
   fout->cd();
   jecMode.Write("zjet_jec_mode",TObject::kOverwrite);
   jecL2.Write("zjet_jec_l2_file",TObject::kOverwrite);
   jecResidual.Write("zjet_jec_residual_file",TObject::kOverwrite);
   jerResolutionMetadata.Write(
     "zjet_jer_resolution_file",TObject::kOverwrite);
   jerScaleFactorMetadata.Write(
     "zjet_jer_scale_factor_file",TObject::kOverwrite);
   muonCorrectionMetadata.Write(
     "zjet_muon_correction_file",TObject::kOverwrite);
   muonCorrectionSha.Write(
     "zjet_muon_correction_sha256",TObject::kOverwrite);
   jetVetoMapMetadata.Write(
     "zjet_jet_veto_map_file",TObject::kOverwrite);
   type1MetMetadata.Write(
     "zjet_type1_met_definition",TObject::kOverwrite);
   flavorDefinition.Write("zjet_flavor_definition",TObject::kOverwrite);
   flavorMatrixDefinition.Write(
     "zjet_flavor_matrix_definition",TObject::kOverwrite);
   synchronizedSelection.Write(
     "zjet_synchronized_selection",TObject::kOverwrite);
   legacyJetId.Write("zjet_legacy_jet_id",TObject::kOverwrite);
   truthDefinition.Write("zjet_truth_hdm_definition",TObject::kOverwrite);
   residualDefinition.Write(
     "zjet_previous_residual_definition",TObject::kOverwrite);
   fout->Purge();
   fout->Close();
   delete jec;
}
