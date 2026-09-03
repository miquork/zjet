#include "TArrayD.h"
#include "TAxis.h"
#include "TDirectory.h"
#include "TFile.h"
#include "TH3D.h"
#include "TObject.h"
#include "TProfile3D.h"

#include "FlavorMatrixTools.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <map>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

const double expectedPtEdges[] = {
  10.,15.,20.,25.,30.,35.,40.,50.,60.,70.,85.,100.,125.,155.,180.,
  210.,250.,300.,350.,400.,500.,600.,800.,1000.,1200.,1500.,1800.,
  2100.,2400.,2700.,3000.,3500.,4000.,
};
const int expectedPtBins =
  sizeof(expectedPtEdges)/sizeof(expectedPtEdges[0])-1;

bool closeEnough(double first, double second, double absolute=1.e-9,
                 double relative=1.e-10) {
  return std::fabs(first-second) <=
         absolute+relative*std::max(std::fabs(first),std::fabs(second));
}

void fail(std::vector<std::string> &failures, const std::string &message) {
  failures.push_back(message);
}

bool hasLabels(const TAxis *axis) {
  if (!axis) return false;
  for (int bin=1; bin<=axis->GetNbins(); ++bin)
    if (axis->GetBinLabel(bin) && axis->GetBinLabel(bin)[0]!='\0') return true;
  return false;
}

void checkPtAxis(const TAxis *axis, const std::string &name,
                 std::vector<std::string> &failures) {
  if (!axis || axis->GetNbins()!=expectedPtBins) {
    fail(failures,name+": pT axis does not have 32 bins");
    return;
  }
  for (int edge=0; edge<=expectedPtBins; ++edge) {
    const double actual = edge==expectedPtBins
      ? axis->GetBinUpEdge(expectedPtBins) : axis->GetBinLowEdge(edge+1);
    if (!closeEnough(actual,expectedPtEdges[edge],1.e-12,1.e-12)) {
      std::ostringstream message;
      message << name << ": pT edge " << edge << " is " << actual
              << ", expected " << expectedPtEdges[edge];
      fail(failures,message.str());
    }
  }
}

void checkFlavorAxis(const TAxis *axis, const std::string &name,
                     std::vector<std::string> &failures) {
  if (!axis || axis->GetNbins()!=7) {
    fail(failures,name+": flavor axis does not have 7 bins");
    return;
  }
  for (int edge=0; edge<=7; ++edge) {
    const double actual = edge==7
      ? axis->GetBinUpEdge(7) : axis->GetBinLowEdge(edge+1);
    const double expected = -0.5+edge;
    if (!closeEnough(actual,expected,1.e-12,1.e-12)) {
      std::ostringstream message;
      message << name << ": flavor edge " << edge << " is " << actual
              << ", expected " << expected;
      fail(failures,message.str());
    }
  }
  if (hasLabels(axis))
    fail(failures,name+": category labels are forbidden in mergeable files");
}

template <class Expected>
Expected *requireExact(TFile *file, const std::string &path,
                       std::vector<std::string> &failures) {
  TObject *object = file ? file->Get(path.c_str()) : nullptr;
  if (!object) {
    fail(failures,"missing "+path);
    return nullptr;
  }
  if (object->IsA()!=Expected::Class()) {
    fail(failures,path+": class is "+object->ClassName()+", expected "+
                  Expected::Class()->GetName());
    return nullptr;
  }
  return static_cast<Expected*>(object);
}

void checkMatrixGeometry(const TH3 *histogram, const std::string &name,
                         std::vector<std::string> &failures) {
  if (!histogram) return;
  if (histogram->GetNbinsX()!=32 || histogram->GetNbinsY()!=7 ||
      histogram->GetNbinsZ()!=7)
    fail(failures,name+": expected geometry 32x7x7");
  checkPtAxis(histogram->GetXaxis(),name+" x",failures);
  checkFlavorAxis(histogram->GetYaxis(),name+" y",failures);
  checkFlavorAxis(histogram->GetZaxis(),name+" z",failures);
}

double binSumW2(const TH1 *histogram, int bin) {
  const TArrayD *sumw2 = histogram ? histogram->GetSumw2() : nullptr;
  return sumw2 && bin<sumw2->GetSize() ? sumw2->At(bin) : 0.;
}

double profileNumerator(const TProfile3D *profile, int bin) {
  return profile && profile->GetArray() ? profile->GetArray()[bin] : 0.;
}

void checkEmptyCountSlice(const TH3D *histogram, bool recoAxis, int id,
                          const std::string &description,
                          std::vector<std::string> &failures) {
  if (!histogram) return;
  const int fixedBin = id+1;
  for (int ix=0; ix<=histogram->GetNbinsX()+1; ++ix) {
    for (int iy=0; iy<=histogram->GetNbinsY()+1; ++iy) {
      for (int iz=0; iz<=histogram->GetNbinsZ()+1; ++iz) {
        if ((recoAxis && iy!=fixedBin) || (!recoAxis && iz!=fixedBin))
          continue;
        const int bin = histogram->GetBin(ix,iy,iz);
        const double content = histogram->GetBinContent(bin);
        const double error2 = binSumW2(histogram,bin);
        if (!closeEnough(content,0.) || !closeEnough(error2,0.)) {
          std::ostringstream message;
          message << histogram->GetName() << ": " << description
                  << " has content at global bin " << bin
                  << " (sumw=" << content << ", sumw2=" << error2 << ")";
          fail(failures,message.str());
          return;
        }
      }
    }
  }
}

void checkUnitScoreAxis(const TAxis *axis, const std::string &name,
                        int bins, std::vector<std::string> &failures) {
  if (!axis || axis->GetNbins()!=bins ||
      !closeEnough(axis->GetBinLowEdge(1),0.,1.e-12,1.e-12) ||
      !closeEnough(axis->GetBinUpEdge(bins),1.,1.e-12,1.e-12))
    fail(failures,name+": score axis has wrong binning");
  if (hasLabels(axis)) fail(failures,name+": unexpected axis labels");
}

} // namespace

void validateFlavorMatrix(const char *fileName, bool isMC) {
  std::vector<std::string> failures;
  std::unique_ptr<TFile> file(TFile::Open(fileName,"READ"));
  if (!file || file->IsZombie())
    throw std::runtime_error(std::string("Cannot open ")+fileName);
  if (!file->GetDirectory("FlavorMatrix"))
    fail(failures,"missing FlavorMatrix directory");
  if (!file->GetDirectory("FlavorMatrix/controls"))
    fail(failures,"missing FlavorMatrix/controls directory");

  const std::vector<std::string> countNames = {
    "h3counts_flavormatrix", "h3counts_parallel_flavormatrix",
    "h3counts_transverse_flavormatrix",
  };
  std::vector<TH3D*> counts;
  for (const std::string &name : countNames) {
    TH3D *histogram = requireExact<TH3D>(
      file.get(),"FlavorMatrix/"+name,failures);
    counts.push_back(histogram);
    checkMatrixGeometry(histogram,"FlavorMatrix/"+name,failures);
  }

  if (counts.size()==3 && counts[0] && counts[1] && counts[2]) {
    for (int bin=0; bin<counts[0]->GetNcells(); ++bin) {
      const double signedCount = counts[0]->GetBinContent(bin);
      const double expected =
        counts[1]->GetBinContent(bin)-counts[2]->GetBinContent(bin);
      if (!closeEnough(signedCount,expected)) {
        std::ostringstream message;
        message << "counts != parallel-transverse at global bin " << bin
                << " (" << signedCount << " vs " << expected << ")";
        fail(failures,message.str());
        break;
      }
      const double signedSumw2 = binSumW2(counts[0],bin);
      const double expectedSumw2 =
        binSumW2(counts[1],bin)+binSumW2(counts[2],bin);
      if (!closeEnough(signedSumw2,expectedSumw2)) {
        std::ostringstream message;
        message << "count sumw2 != parallel+transverse at global bin " << bin
                << " (" << signedSumw2 << " vs " << expectedSumw2 << ")";
        fail(failures,message.str());
        break;
      }
    }
  }
  for (TH3D *histogram : counts) {
    checkEmptyCountSlice(histogram,true,2,"reserved reco ID 2",failures);
    checkEmptyCountSlice(histogram,true,3,"reserved reco ID 3",failures);
    if (isMC)
      checkEmptyCountSlice(histogram,false,2,
                           "reserved true ID 2",failures);
    if (!isMC)
      for (int truth=1; truth<=6; ++truth)
        checkEmptyCountSlice(histogram,false,truth,
                             "data truth ID "+std::to_string(truth),failures);
  }

  TH3D *heavyTopology = requireExact<TH3D>(
    file.get(),"FlavorMatrix/h3counts_heavytopology",failures);
  checkMatrixGeometry(
    heavyTopology,"FlavorMatrix/h3counts_heavytopology",failures);
  if (heavyTopology && !counts.empty() && counts[0]) {
    bool projectionMismatch = false;
    for (int ix=0;
         ix<=heavyTopology->GetNbinsX()+1 && !projectionMismatch; ++ix) {
      for (int iy=0; iy<=heavyTopology->GetNbinsY()+1; ++iy) {
        double truthSum = 0.;
        double truthSumw2 = 0.;
        double topologySum = 0.;
        double topologySumw2 = 0.;
        for (int iz=0; iz<=heavyTopology->GetNbinsZ()+1; ++iz) {
          const int truthBin = counts[0]->GetBin(ix,iy,iz);
          const int topologyBin = heavyTopology->GetBin(ix,iy,iz);
          truthSum += counts[0]->GetBinContent(truthBin);
          truthSumw2 += binSumW2(counts[0],truthBin);
          topologySum += heavyTopology->GetBinContent(topologyBin);
          topologySumw2 += binSumW2(heavyTopology,topologyBin);
        }
        if (!closeEnough(truthSum,topologySum) ||
            !closeEnough(truthSumw2,topologySumw2)) {
          std::ostringstream message;
          message << "heavy-topology projection differs from flavor counts"
                  << " at x/y bins " << ix << "/" << iy;
          fail(failures,message.str());
          projectionMismatch = true;
          break;
        }
      }
    }
    checkEmptyCountSlice(
      heavyTopology,true,2,"reserved reco ID 2",failures);
    checkEmptyCountSlice(
      heavyTopology,true,3,"reserved reco ID 3",failures);
    if (isMC)
      checkEmptyCountSlice(
        heavyTopology,false,3,"reserved heavy-topology ID 3",failures);
    else
      for (int topology=0; topology<=5; ++topology)
        checkEmptyCountSlice(
          heavyTopology,false,topology,
          "data heavy-topology ID "+std::to_string(topology),failures);
  }

  const std::vector<std::string> observables = {
    "m0", "m2", "mn", "mu", "mnu", "hdm",
  };
  const std::vector<std::string> variants = {"ab","ad","tc","pf"};
  for (const std::string &suffix : {
         std::string("_flavormatrix"),
         std::string("_parallel_flavormatrix")}) {
    for (const std::string &variant : variants) {
      std::map<std::string,TProfile3D*> profiles;
      for (const std::string &observable : observables) {
        const std::string name = "p3"+observable+variant+suffix;
        TProfile3D *profile = requireExact<TProfile3D>(
          file.get(),"FlavorMatrix/"+name,failures);
        profiles[observable] = profile;
        checkMatrixGeometry(profile,"FlavorMatrix/"+name,failures);
      }
      if (!profiles["m0"] || !profiles["mn"] || !profiles["mu"] ||
          !profiles["mnu"] || !profiles["hdm"])
        continue;

      for (int bin=0; bin<profiles["m0"]->GetNcells(); ++bin) {
        const double sumw = profiles["m0"]->GetBinEntries(bin);
        for (const char *component : {"mn","mu","mnu"}) {
          const double other = profiles[component]->GetBinEntries(bin);
          if (!closeEnough(other,sumw)) {
            std::ostringstream message;
            message << variant << suffix << ": " << component
                    << " entries differ from m0 at global bin " << bin;
            fail(failures,message.str());
          }
        }
        const double mnuNumerator = profileNumerator(profiles["mnu"],bin);
        const double componentNumerator =
          profileNumerator(profiles["mn"],bin)+
          profileNumerator(profiles["mu"],bin);
        if (!closeEnough(mnuNumerator,componentNumerator,1.e-8,2.e-10)) {
          std::ostringstream message;
          message << variant << suffix
                  << ": mnu != mn+mu at global bin " << bin
                  << " (weighted numerators " << mnuNumerator << " vs "
                  << componentNumerator << ")";
          fail(failures,message.str());
        }

        TProfile3D *hdm = profiles["hdm"];
        const double hdmSumw = hdm->GetBinEntries(bin);
        const double m0 = profiles["m0"]->GetBinContent(bin);
        const double mn = profiles["mn"]->GetBinContent(bin);
        const double mu = profiles["mu"]->GetBinContent(bin);
        const double denominator =
          1.-mn/ZJetFlavorMatrix::responseN-mu/ZJetFlavorMatrix::responseU;
        const bool expectedDefined =
          sumw>0. && std::isfinite(m0) && std::isfinite(mn) &&
          std::isfinite(mu) && std::isfinite(denominator) &&
          std::fabs(denominator)>=1.e-8;
        if (expectedDefined) {
          const double expected = (m0-mn-mu)/denominator;
          if (!closeEnough(hdmSumw,sumw) ||
              !closeEnough(hdm->GetBinContent(bin),expected,1.e-8,2.e-10)) {
            std::ostringstream message;
            message << variant << suffix
                    << ": final HDM mismatch at global bin " << bin
                    << " (value " << hdm->GetBinContent(bin)
                    << ", expected " << expected << ")";
            fail(failures,message.str());
          }
        }
        else if (!closeEnough(hdmSumw,0.) ||
                 !closeEnough(profileNumerator(hdm,bin),0.)) {
          std::ostringstream message;
          message << variant << suffix
                  << ": HDM is populated in undefined global bin " << bin;
          fail(failures,message.str());
        }
      }
    }
  }

  // Environmental and jet-area controls were added after the first
  // FlavorMatrix production. Keep old merged files readable, but once any
  // member of a generation is present require its complete variant/suffix
  // set and matching entries.
  const bool hasFlavorEnvironmentControls =
    file->Get("FlavorMatrix/p3rhotc_flavormatrix") ||
    file->Get("FlavorMatrix/p3mueftc_flavormatrix");
  if (hasFlavorEnvironmentControls) {
    for (const std::string &suffix : {
           std::string("_flavormatrix"),
           std::string("_parallel_flavormatrix")})
      for (const std::string &variant : variants) {
        TProfile3D *reference = dynamic_cast<TProfile3D*>(file->Get(
          ("FlavorMatrix/p3m0"+variant+suffix).c_str()));
        for (const std::string &observable : {"rho","muef"}) {
          const std::string name = "p3"+observable+variant+suffix;
          TProfile3D *profile = requireExact<TProfile3D>(
            file.get(),"FlavorMatrix/"+name,failures);
          checkMatrixGeometry(profile,"FlavorMatrix/"+name,failures);
          if (!reference || !profile) continue;
          for (int bin=0; bin<reference->GetNcells(); ++bin)
            if (!closeEnough(profile->GetBinEntries(bin),
                             reference->GetBinEntries(bin))) {
              fail(failures,name+" entries differ from p3m0 at global bin "+
                std::to_string(bin));
              break;
            }
        }
      }
  }

  const bool hasFlavorAreaControls =
    file->Get("FlavorMatrix/p3areasumtc_flavormatrix") ||
    file->Get("FlavorMatrix/p3areaprojtc_flavormatrix") ||
    file->Get("FlavorMatrix/p3ueholetc_flavormatrix") ||
    file->Get("FlavorMatrix/p3mnufsrtc_flavormatrix");
  if (hasFlavorAreaControls) {
    for (const std::string &suffix : {
           std::string("_flavormatrix"),
           std::string("_parallel_flavormatrix")})
      for (const std::string &variant : variants) {
        TProfile3D *reference = dynamic_cast<TProfile3D*>(file->Get(
          ("FlavorMatrix/p3m0"+variant+suffix).c_str()));
        for (const std::string &observable : {
               "areasum","areaproj","uehole","mnufsr"}) {
          const std::string name = "p3"+observable+variant+suffix;
          TProfile3D *profile = requireExact<TProfile3D>(
            file.get(),"FlavorMatrix/"+name,failures);
          checkMatrixGeometry(profile,"FlavorMatrix/"+name,failures);
          if (!reference || !profile) continue;
          for (int bin=0; bin<reference->GetNcells(); ++bin)
            if (!closeEnough(profile->GetBinEntries(bin),
                             reference->GetBinEntries(bin))) {
              fail(failures,name+" entries differ from p3m0 at global bin "+
                std::to_string(bin));
              break;
            }
        }
      }
  }

  const bool hasFlavorTruthRecoil =
    file->Get("FlavorMatrix/p3genmutc_flavormatrix") ||
    file->Get("FlavorMatrix/p3recogenmutc_flavormatrix");
  if (hasFlavorTruthRecoil) {
    for (const std::string &suffix : {
           std::string("_flavormatrix"),
           std::string("_parallel_flavormatrix")})
      for (const std::string &variant : variants)
        for (const std::string &observable : {
               "recomnmatched","recoumatched","genmn","genmu",
               "recogenmn","recogenmu",
               "genmn2","genmu2"}) {
          const std::string name = "p3"+observable+variant+suffix;
          TProfile3D *profile = requireExact<TProfile3D>(
            file.get(),"FlavorMatrix/"+name,failures);
          checkMatrixGeometry(profile,"FlavorMatrix/"+name,failures);
        }
  }

  const std::vector<std::string> pairControls = {
    "h3_cvb_cvl_trueflavor", "h3_cvb_qvg_trueflavor",
    "h3_cvl_qvg_trueflavor", "h3_upartqvg_pnetqvg_trueflavor",
  };
  for (const std::string &name : pairControls) {
    TH3D *histogram = requireExact<TH3D>(
      file.get(),"FlavorMatrix/controls/"+name,failures);
    if (!histogram) continue;
    checkUnitScoreAxis(histogram->GetXaxis(),name+" x",50,failures);
    checkUnitScoreAxis(histogram->GetYaxis(),name+" y",50,failures);
    checkFlavorAxis(histogram->GetZaxis(),name+" z",failures);
    if (isMC)
      checkEmptyCountSlice(histogram,false,2,
                           "reserved true ID 2",failures);
    if (!isMC)
      for (int truth=1; truth<=6; ++truth)
        checkEmptyCountSlice(histogram,false,truth,
                             "data truth ID "+std::to_string(truth),failures);
  }
  const std::vector<std::string> topologyControls = {
    "h3_cvb_cvl_heavytopology", "h3_cvb_qvg_heavytopology",
    "h3_cvl_qvg_heavytopology",
  };
  for (const std::string &name : topologyControls) {
    TH3D *histogram = requireExact<TH3D>(
      file.get(),"FlavorMatrix/controls/"+name,failures);
    if (!histogram) continue;
    checkUnitScoreAxis(histogram->GetXaxis(),name+" x",50,failures);
    checkUnitScoreAxis(histogram->GetYaxis(),name+" y",50,failures);
    checkFlavorAxis(histogram->GetZaxis(),name+" z",failures);
  }
  TH3D *heavyMultiplicity = requireExact<TH3D>(
    file.get(),"FlavorMatrix/controls/h3_genjet_nc_nb_trueflavor",failures);
  if (heavyMultiplicity) {
    if (heavyMultiplicity->GetNbinsX()!=5 ||
        heavyMultiplicity->GetNbinsY()!=5)
      fail(failures,
           "h3_genjet_nc_nb_trueflavor: expected 5x5 hadron-count axes");
    checkFlavorAxis(heavyMultiplicity->GetZaxis(),
                    "h3_genjet_nc_nb_trueflavor z",failures);
  }
  for (int truth=0; truth<=6; ++truth) {
    const std::string name = "h3_cvb_cvl_qvg_true"+std::to_string(truth);
    TH3D *histogram = requireExact<TH3D>(
      file.get(),"FlavorMatrix/controls/"+name,failures);
    if (!histogram) continue;
    checkUnitScoreAxis(histogram->GetXaxis(),name+" x",24,failures);
    checkUnitScoreAxis(histogram->GetYaxis(),name+" y",24,failures);
    checkUnitScoreAxis(histogram->GetZaxis(),name+" z",24,failures);
    if (isMC && truth==2) {
      for (int bin=0; bin<histogram->GetNcells(); ++bin)
        if (!closeEnough(histogram->GetBinContent(bin),0.) ||
            !closeEnough(binSumW2(histogram,bin),0.)) {
          fail(failures,"reserved "+name+" is not empty");
          break;
        }
    }
    if (!isMC && truth!=0) {
      for (int bin=0; bin<histogram->GetNcells(); ++bin)
        if (!closeEnough(histogram->GetBinContent(bin),0.) ||
            !closeEnough(binSumW2(histogram,bin),0.)) {
          fail(failures,"data "+name+" is not empty");
          break;
        }
    }
  }

  if (!failures.empty()) {
    std::cerr << "FlavorMatrix validation failed for " << fileName << ":\n";
    const std::size_t maximumPrintedFailures = 100;
    for (std::size_t index=0;
         index<std::min(failures.size(),maximumPrintedFailures); ++index)
      std::cerr << "  - " << failures[index] << '\n';
    if (failures.size()>maximumPrintedFailures)
      std::cerr << "  - ... " << failures.size()-maximumPrintedFailures
                << " additional problem(s) omitted\n";
    throw std::runtime_error(
      "FlavorMatrix validation found "+std::to_string(failures.size())+
      " problem(s)");
  }
  std::cout << "Validated FlavorMatrix in " << fileName
            << " (" << (isMC ? "MC" : "data")
            << "): exact schema, merge-safe axes, signed-weight identities, "
            << "component closure, HDM closure, and controls are valid."
            << std::endl;
}
