#include "TDirectory.h"
#include "TFile.h"
#include "TH1.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TKey.h"
#include "TNamed.h"
#include "TObjString.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TROOT.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <map>
#include <memory>
#include <string>
#include <utility>
#include <vector>

namespace {

template <typename ObjectType>
ObjectType *requireObject(TDirectory *source, const std::string &path) {
  if (!source)
    throw std::runtime_error("Cannot read " + path + " from a null directory");
  TObject *object = source->Get(path.c_str());
  if (!object)
    throw std::runtime_error(
      "Missing object " + std::string(source->GetName()) + ":" + path);
  ObjectType *typed = dynamic_cast<ObjectType*>(object);
  if (!typed)
    throw std::runtime_error(
      "Object " + std::string(source->GetName()) + ":" + path +
      " has class " + object->ClassName() + ", expected " +
      ObjectType::Class()->GetName());
  return typed;
}

TDirectory *requireDirectory(TDirectory *source, const std::string &path) {
  if (!source)
    throw std::runtime_error("Cannot read " + path + " from a null directory");
  TDirectory *directory = source->GetDirectory(path.c_str());
  if (!directory)
    throw std::runtime_error(
      "Missing directory " + std::string(source->GetName()) + ":" + path);
  return directory;
}

TDirectory *makeDirectory(TDirectory *parent, const char *name) {
  if (!parent)
    throw std::runtime_error(
      "Cannot create directory " + std::string(name) + " below null parent");
  TDirectory *directory = parent->mkdir(name);
  if (!directory)
    throw std::runtime_error(
      "Failed to create directory " + std::string(parent->GetName()) + "/" +
      name);
  return directory;
}

void copyDirectory(TDirectory *source, TDirectory *target) {
  if (!source || !target)
    throw std::runtime_error("copyDirectory received a null directory");

  // A ROOT directory may contain many cycles of an identically named object.
  // Copy only the newest cycle so compatibility files stay compact and do not
  // accumulate tmp_* or response-profile duplicates on regeneration.
  std::map<std::string,TKey*> newestKeys;
  TIter findKeys(source->GetListOfKeys());
  while (TKey *key = static_cast<TKey*>(findKeys())) {
    auto found = newestKeys.find(key->GetName());
    if (found==newestKeys.end() || key->GetCycle()>found->second->GetCycle())
      newestKeys[key->GetName()] = key;
  }
  for (const auto &entry : newestKeys) {
    TKey *key = entry.second;
    TClass *objectClass = gROOT->GetClass(key->GetClassName());
    if (!objectClass)
      throw std::runtime_error(
        "Unknown ROOT class " + std::string(key->GetClassName()) + " for " +
        source->GetPath() + "/" + key->GetName());
    if (objectClass->InheritsFrom(TDirectory::Class())) {
      TDirectory *sourceChild = requireDirectory(source,key->GetName());
      TDirectory *targetChild = makeDirectory(target,key->GetName());
      copyDirectory(sourceChild,targetChild);
    }
    else {
      TObject *object = key->ReadObj();
      if (!object)
        throw std::runtime_error(
          "Failed to read " + std::string(source->GetPath()) + "/" +
          key->GetName());
      target->cd();
      object->Write(key->GetName());
      delete object;
    }
  }
}

void copyObjectOverwrite(TDirectory *source, TDirectory *target,
                         const std::string &name) {
  if (!source || !target)
    throw std::runtime_error("copyObjectOverwrite received a null directory");
  TObject *object = source->Get(name.c_str());
  if (!object)
    throw std::runtime_error(
      "Missing method-specific object " + std::string(source->GetPath()) +
      "/" + name);
  target->cd();
  object->Write(name.c_str(),TObject::kOverwrite);
}

void writeStoredHistogram(TDirectory *source, TDirectory *target,
                          const std::string &sourceName,
                          const std::string &targetName) {
  TH1 *input = requireObject<TH1>(source,sourceName);
  TH1 *output = dynamic_cast<TH1*>(input->Clone(targetName.c_str()));
  if (!output)
    throw std::runtime_error("Failed to clone " + sourceName);
  output->SetDirectory(nullptr);
  target->cd();
  output->Write(targetName.c_str());
  delete output;
}

void writeProfile(TFile *source, TDirectory *target,
                  const char *sourceName, const char *targetName,
                  int firstEtaBin, int lastEtaBin) {
  TProfile2D *profile = requireObject<TProfile2D>(
    source,"l2res/" + std::string(sourceName));
  TDirectory::TContext directoryContext(gROOT);
  TProfile *projection =
    profile->ProfileY(Form("tmp_%s",targetName),firstEtaBin,lastEtaBin);
  if (!projection)
    throw std::runtime_error(
      "Failed to project " + std::string(source->GetName()) + ":l2res/" +
      sourceName);
  projection->SetName(targetName);
  projection->SetDirectory(nullptr);
  target->cd();
  projection->Write(targetName);
  delete projection;
}

void writeCounts(TFile *source, TDirectory *target,
                 const char *sourceName, const char *targetName,
                 int firstEtaBin, int lastEtaBin) {
  TH2D *counts = requireObject<TH2D>(
    source,"l2res/" + std::string(sourceName));
  TDirectory::TContext directoryContext(gROOT);
  TH1D *projection = counts->ProjectionY(
    targetName,firstEtaBin,lastEtaBin,"e");
  if (!projection)
    throw std::runtime_error(
      "Failed to project " + std::string(source->GetName()) +
      ":l2res/" + sourceName);
  projection->SetDirectory(nullptr);
  target->cd();
  projection->Write();
  delete projection;
}

void writeResponseBinning(TFile *source, TDirectory *target,
                          const char *sourceSuffix, const char *targetAxis,
                          const char *countsName,
                          int firstEtaBin, int lastEtaBin) {
  const std::vector<std::pair<std::string,std::string> > profiles = {
    {"m0", "rmpf"},
    {"m2", "rmpfjet1"},
    {"mn", "rmpfjetn"},
    {"mu", "rmpfuncl"},
  };
  for (const auto &entry : profiles) {
    const std::string sourceName =
      "p2" + entry.first + std::string(sourceSuffix);
    const std::string targetName =
      entry.second + "_" + targetAxis + "_a100";
    writeProfile(source,target,sourceName.c_str(),targetName.c_str(),
                 firstEtaBin,lastEtaBin);
  }
  const std::string targetCounts =
    "statistics_rmpf_" + std::string(targetAxis) + "_a100";
  writeCounts(source,target,countsName,targetCounts.c_str(),
              firstEtaBin,lastEtaBin);
}

void writeNativeResponseBinning(TDirectory *profiles, TDirectory *target,
                                const char *axisName) {
  TDirectory *axis = requireDirectory(profiles,axisName);
  const std::vector<std::string> names = {
    "statistics_rmpf", "rmpf", "rmpfjet1", "rmpfjetn", "rmpfuncl",
    "rmpfjetnu",
  };
  for (const std::string &name : names)
    writeStoredHistogram(axis,target,name,
                         name+"_"+axisName+"_a100");
}

void writeEmptyFlavorObject(TDirectory *models, TDirectory *target,
                            const std::string &baseName,
                            const std::string &tag) {
  const std::string modelName = baseName + "_zmmjet_a100";
  TH1 *model = requireObject<TH1>(models,modelName);
  const std::string targetName = baseName + "_zmmjet" + tag + "_a100";
  TH1 *placeholder = dynamic_cast<TH1*>(model->Clone(targetName.c_str()));
  if (!placeholder)
    throw std::runtime_error("Failed to clone flavor placeholder " + targetName);
  placeholder->Reset("ICES");
  placeholder->SetDirectory(nullptr);
  placeholder->SetTitle("PLACEHOLDER: Z+flavor analysis not implemented");
  target->cd();
  placeholder->Write();
  delete placeholder;
}

void writeStoredFlavorObject(TDirectory *source, TDirectory *target,
                             const std::string &sourceName,
                             const std::string &targetName) {
  TH1 *input = requireObject<TH1>(source,sourceName);
  TH1 *output = dynamic_cast<TH1*>(input->Clone(targetName.c_str()));
  if (!output)
    throw std::runtime_error("Failed to clone measured flavor input " +
                             sourceName);
  output->SetDirectory(nullptr);
  target->cd();
  output->Write();
  delete output;
}

bool writeStoredFlavorInputs(TFile *source, TDirectory *sampleDirectory,
                             bool isMC) {
  TDirectory *flavor = source->GetDirectory("flavor");
  if (!flavor)
    return false;

  const std::vector<std::pair<std::string,std::string> > objects = {
    {"statistics_rmpf", "counts"},
    {"rmpf", "mpfchs1"},
    {"rmpfjet1", "mpf1"},
    {"rmpfjetn", "mpfn"},
    {"rmpfuncl", "mpfu"},
    {"rbal", "rjet"},
    {"rgenjet1", "gjet"},
  };
  const std::vector<std::pair<std::string,std::string> > recoTags = {
    {"i", ""}, {"b", "_btag"}, {"c", "_ctag"},
    {"q", "_quarktag"}, {"g", "_gluontag"}, {"n", "_notag"},
  };
  const std::vector<std::pair<std::string,std::string> > genTags = {
    {"i", "eta_00_13"}, {"b", "eta_00_13_genb"},
    {"c", "eta_00_13_genc"}, {"q", "eta_00_13_genuds"},
    {"g", "eta_00_13_geng"}, {"n", "eta_00_13_unclassified"},
  };

  for (const auto &genTag : genTags) {
    if (!isMC && genTag.first!="i") continue;
    TDirectory *target =
      (genTag.first=="i"
         ? requireDirectory(sampleDirectory,"eta_00_13")
         : makeDirectory(sampleDirectory,genTag.second.c_str()));
    for (const auto &recoTag : recoTags) {
      // The inclusive object was already written from the main l2res profile.
      if (genTag.first=="i" && recoTag.first=="i") continue;
      for (const auto &object : objects) {
        if (!isMC && (object.first=="rbal" || object.first=="rgenjet1"))
          continue;
        const std::string sourceName =
          object.second + "_g" + recoTag.first + genTag.first;
        const std::string targetName =
          object.first + "_zmmjet" + recoTag.second + "_a100";
        writeStoredFlavorObject(flavor,target,sourceName,targetName);
      }
    }
  }
  return true;
}

void writeEmptyFlavorSet(TDirectory *models, TDirectory *target,
                         const std::string &tag, bool isMC) {
  const std::vector<std::string> common = {
    "statistics_rmpf", "rmpf", "rmpfjet1", "rmpfjetn", "rmpfuncl",
  };
  for (const std::string &baseName : common)
    writeEmptyFlavorObject(models,target,baseName,tag);
  if (isMC) {
    writeEmptyFlavorObject(models,target,"rbal",tag);
    writeEmptyFlavorObject(models,target,"rgenjet1",tag);
  }
}

void writeFlavorPlaceholders(TDirectory *sampleDirectory, bool isMC) {
  TDirectory *inclusive = requireDirectory(sampleDirectory,"eta_00_13");
  const std::vector<std::string> recoTags = {
    "_btag", "_ctag", "_quarktag", "_gluontag", "_notag",
  };
  for (const std::string &tag : recoTags)
    writeEmptyFlavorSet(inclusive,inclusive,tag,isMC);

  if (!isMC)
    return;
  const std::vector<std::string> genDirectories = {
    "eta_00_13_genb", "eta_00_13_genc", "eta_00_13_genuds",
    "eta_00_13_geng", "eta_00_13_unclassified",
  };
  const std::vector<std::string> allRecoTags = {
    "", "_btag", "_ctag", "_quarktag", "_gluontag", "_notag",
  };
  for (const std::string &directoryName : genDirectories) {
    TDirectory *directory = makeDirectory(sampleDirectory,directoryName.c_str());
    for (const std::string &tag : allRecoTags)
      writeEmptyFlavorSet(inclusive,directory,tag,true);
  }
}

bool writeGlobalFitInputs(TFile *source, TDirectory *sampleDirectory,
                          bool isMC, bool useLegacyMethod,
                          bool preferOneDimensional) {
  TProfile2D *etaReference =
    requireObject<TProfile2D>(source,"l2res/p2m0tc");
  const int firstEtaBin = etaReference->GetXaxis()->FindFixBin(0.+1.e-6);
  const int lastEtaBin = etaReference->GetXaxis()->FindFixBin(1.305-1.e-6);

  TDirectory *etaDirectory = makeDirectory(sampleDirectory,"eta_00_13");
  const std::string nativePath =
    (useLegacyMethod ? "legacy/profiles1d" : "profiles1d");
  TDirectory *nativeProfiles = source->GetDirectory(nativePath.c_str());
  const bool useNativeProfiles = preferOneDimensional && nativeProfiles;

  // reprocess.C uses three complementary reference-pT choices. The tc, pf
  // and unsuffixed raw profiles are binned in Z pT, jet pT and average pT,
  // respectively. All three are projections of the same accepted pairs.
  if (useNativeProfiles) {
    writeNativeResponseBinning(nativeProfiles,etaDirectory,"zmmjet");
    writeNativeResponseBinning(nativeProfiles,etaDirectory,"jetpt");
    writeNativeResponseBinning(nativeProfiles,etaDirectory,"ptave");
    // Keep the previous residual conditional on every reference-pT choice.
    // Current reprocess.C consumes residual_zmmjet_a100, while the two extra
    // profiles make an axis-consistent iterative treatment possible without
    // rerunning the event analysis.
    for (const char *axisName : {"zmmjet","jetpt","ptave"}) {
      TDirectory *axis = requireDirectory(nativeProfiles,axisName);
      writeStoredHistogram(axis,etaDirectory,"residual",
                           std::string("residual_")+axisName+"_a100");
    }
  }
  else {
    if (useLegacyMethod)
      throw std::runtime_error(
        "Legacy method requested but legacy/profiles1d is unavailable");
    writeResponseBinning(source,etaDirectory,"tc","zmmjet","h2ptetatc",
                         firstEtaBin,lastEtaBin);
    writeResponseBinning(source,etaDirectory,"pf","jetpt","h2ptetapf",
                         firstEtaBin,lastEtaBin);
    writeResponseBinning(source,etaDirectory,"","ptave","h2pteta",
                         firstEtaBin,lastEtaBin);
  }

  if (useNativeProfiles) {
    TDirectory *zmm = requireDirectory(nativeProfiles,"zmmjet");
    const std::vector<std::pair<std::string,std::string> > profiles = {
      {"chHEF", "h_Zpt_chHEF_alpha100"},
      {"neEmEF", "h_Zpt_neEmEF_alpha100"},
      {"neHEF", "h_Zpt_neHEF_alpha100"},
      {"chEmEF", "h_Zpt_chEmEF_alpha100"},
      {"muEF", "h_Zpt_muEF_alpha100"},
      {"rho", "h_Zpt_rho_alpha100"},
      {"rbal", "rbal_zmmjet_a100"},
      {"mass", "h_Zpt_mZ_alpha100"},
    };
    for (const auto &entry : profiles)
      writeStoredHistogram(zmm,etaDirectory,entry.first,entry.second);
    if (isMC)
      writeStoredHistogram(zmm,etaDirectory,"rgenjet1",
                           "rgenjet1_zmmjet_a100");
  }
  else {
    const std::vector<std::pair<std::string,std::string> > profiles = {
      {"p2chftc", "h_Zpt_chHEF_alpha100"},
      {"p2neftc", "h_Zpt_neEmEF_alpha100"},
      {"p2nhftc", "h_Zpt_neHEF_alpha100"},
      {"p2ceftc", "h_Zpt_chEmEF_alpha100"},
      {"p2muftc", "h_Zpt_muEF_alpha100"},
      {"p2rhotc", "h_Zpt_rho_alpha100"},
      {"p2m2tc", "rbal_zmmjet_a100"},
    };
    for (const auto &entry : profiles)
      writeProfile(source,etaDirectory,entry.first.c_str(),entry.second.c_str(),
                   firstEtaBin,lastEtaBin);
    if (isMC)
      writeProfile(source,etaDirectory,"p2rgentc","rgenjet1_zmmjet_a100",
                   firstEtaBin,lastEtaBin);
    TH2D *mass = requireObject<TH2D>(source,"l2res/h2mztc");
    etaDirectory->cd();
    mass->Write("h_Zpt_mZ_alpha100");
  }
  return useNativeProfiles;
}

TH1D *makeHistogramWithAxis(const TH1 *model, const std::string &name) {
  if (!model)
    throw std::runtime_error("Cannot create " + name + " from a null model");
  const TAxis *axis = model->GetXaxis();
  TH1D *histogram = nullptr;
  if (axis->GetXbins()->GetSize()>0)
    histogram = new TH1D(name.c_str(),"",axis->GetNbins(),
                         axis->GetXbins()->GetArray());
  else
    histogram = new TH1D(name.c_str(),"",axis->GetNbins(),
                         axis->GetXmin(),axis->GetXmax());
  histogram->SetDirectory(nullptr);
  histogram->Sumw2();
  histogram->GetXaxis()->SetTitle(axis->GetTitle());
  return histogram;
}

void requireSameAxis(const TH1 *reference, const TH1 *candidate,
                     const std::string &description) {
  if (!reference || !candidate ||
      reference->GetNbinsX()!=candidate->GetNbinsX())
    throw std::runtime_error("Incompatible HDM input binning for " +
                             description);
  for (int bin = 1; bin <= reference->GetNbinsX()+1; ++bin) {
    const double first = reference->GetXaxis()->GetBinLowEdge(bin);
    const double second = candidate->GetXaxis()->GetBinLowEdge(bin);
    const double scale = std::max(1.,std::max(std::fabs(first),
                                             std::fabs(second)));
    if (std::fabs(first-second)>1.e-9*scale)
      throw std::runtime_error("Incompatible HDM bin edges for " +
                               description);
  }
}

double solveZJetHDM(double r0, double rn, double ru,
                    double responseN, double responseU) {
  const double denominator = 1.-rn/responseN-ru/responseU;
  if (std::fabs(denominator)<1.e-12)
    return std::numeric_limits<double>::quiet_NaN();
  return (r0-rn-ru)/denominator;
}

void writeHDMCombination(TFile *output, double responseN, double responseU) {
  if (!output || responseN<=0. || responseU<=0.)
    throw std::runtime_error("Invalid output or HDM recoil response");
  TDirectory *dataDirectory = requireDirectory(output,"data/eta_00_13");
  TDirectory *mcDirectory = requireDirectory(output,"mc/eta_00_13");
  TDirectory *ratioRoot = makeDirectory(output,"ratio");
  TDirectory *ratioDirectory = makeDirectory(ratioRoot,"eta_00_13");
  const std::vector<std::pair<std::string,std::string> > axes = {
    {"zmmjet","zjet"}, {"jetpt","jetz"}, {"ptave","zjav"},
  };
  for (const auto &axis : axes) {
    const std::string suffix = "_"+axis.first+"_a100";
    TH1 *dataMpf = requireObject<TH1>(dataDirectory,"rmpf"+suffix);
    TH1 *dataNeutral = requireObject<TH1>(dataDirectory,"rmpfjetn"+suffix);
    TH1 *dataUnclustered = requireObject<TH1>(dataDirectory,"rmpfuncl"+suffix);
    TH1 *mcMpf = requireObject<TH1>(mcDirectory,"rmpf"+suffix);
    TH1 *mcNeutral = requireObject<TH1>(mcDirectory,"rmpfjetn"+suffix);
    TH1 *mcUnclustered = requireObject<TH1>(mcDirectory,"rmpfuncl"+suffix);
    for (TH1 *input : {dataNeutral,dataUnclustered,mcMpf,mcNeutral,
                       mcUnclustered})
      requireSameAxis(dataMpf,input,axis.first);

    const std::string name = "hdm_mpfchs1_"+axis.second;
    std::unique_ptr<TH1D> dataHDM(makeHistogramWithAxis(dataMpf,name));
    std::unique_ptr<TH1D> mcHDM(makeHistogramWithAxis(mcMpf,name));
    std::unique_ptr<TH1D> ratioHDM(makeHistogramWithAxis(dataMpf,name));
    for (int bin = 1; bin <= dataMpf->GetNbinsX(); ++bin) {
      const double r0 = dataMpf->GetBinContent(bin);
      const double rn = dataNeutral->GetBinContent(bin);
      const double ru = dataUnclustered->GetBinContent(bin);
      const double q0 = mcMpf->GetBinContent(bin);
      const double qn = mcNeutral->GetBinContent(bin);
      const double qu = mcUnclustered->GetBinContent(bin);
      if (r0==0. || q0==0.) continue;
      const double dataValue = solveZJetHDM(r0,rn,ru,responseN,responseU);
      const double mcValue = solveZJetHDM(q0,qn,qu,responseN,responseU);
      if (!std::isfinite(dataValue) || !std::isfinite(mcValue) || mcValue==0.)
        continue;
      dataHDM->SetBinContent(bin,dataValue);
      dataHDM->SetBinError(bin,dataMpf->GetBinError(bin));
      mcHDM->SetBinContent(bin,mcValue);
      mcHDM->SetBinError(bin,mcMpf->GetBinError(bin));
      ratioHDM->SetBinContent(bin,dataValue/mcValue);
      // Preserve softrad3.C's present statistical-error convention exactly.
      ratioHDM->SetBinError(
        bin,std::hypot(dataMpf->GetBinError(bin),mcMpf->GetBinError(bin)));
    }
    for (const auto &target :
         std::vector<std::pair<TDirectory*,TH1D*> >{
           {dataDirectory,dataHDM.get()}, {mcDirectory,mcHDM.get()},
           {ratioDirectory,ratioHDM.get()}}) {
      target.first->cd();
      target.second->Write(name.c_str(),TObject::kOverwrite);
    }
  }
}

bool copySample(TFile *source, TFile *target, const char *sample,
                bool isMC, bool addFlavorPlaceholders,
                bool useLegacyMethod, bool preferOneDimensional,
                bool &usedNativeProfiles, bool &usedMethodL2Res) {
  TDirectory *sampleDirectory = makeDirectory(target,sample);
  usedNativeProfiles = writeGlobalFitInputs(
    source,sampleDirectory,isMC,useLegacyMethod,preferOneDimensional);
  const bool storedFlavorInputs =
    writeStoredFlavorInputs(source,sampleDirectory,isMC);
  if (!storedFlavorInputs && addFlavorPlaceholders)
    writeFlavorPlaceholders(sampleDirectory,isMC);
  for (const char *name : {"l2res","l2res1"}) {
    TDirectory *sourceDirectory = requireDirectory(source,name);
    TDirectory *targetDirectory = makeDirectory(sampleDirectory,name);
    copyDirectory(sourceDirectory,targetDirectory);
  }
  usedMethodL2Res = false;
  if (useLegacyMethod) {
    TDirectory *legacyL2Res = source->GetDirectory("legacy/l2res");
    if (legacyL2Res) {
      TDirectory *targetL2Res = requireDirectory(sampleDirectory,"l2res");
      for (const char *name : {
             "h2pteta", "h2ptetapf", "h2ptetatc",
             "p2jes", "p2jespf", "p2jestc",
             "p2res", "p2respf", "p2restc",
           })
        copyObjectOverwrite(legacyL2Res,targetL2Res,name);
      usedMethodL2Res = true;
    }
  }
  return storedFlavorInputs;
}

} // namespace

// Create one compatibility file that can be assigned to both ZMM_*_DATA and
// ZMM_*_MC in jecsys3/Config.C. It contains both the raw L2/L3 profile layout
// and the eta_00_13 inputs consumed by reprocess.C, softrad3.C and globalFit.C.
void writeJecsys3(
  const char *dataFile="rootfiles/zjet_DATA.root",
  const char *mcFile="rootfiles/zjet_MC.root",
  const char *outputFile="rootfiles/zjet_JMENANO_compat.root",
  bool addFlavorPlaceholders=false,
  bool useLegacyMethod=false,
  bool preferOneDimensional=true,
  double hdmResponseN=1.00,
  double hdmResponseU=0.92) {
  if (hdmResponseN<=0. || hdmResponseU<=0.)
    throw std::runtime_error("HDM recoil responses Rn and Ru must be positive");
  std::unique_ptr<TFile> data(TFile::Open(dataFile,"READ"));
  std::unique_ptr<TFile> mc(TFile::Open(mcFile,"READ"));
  if (!data || data->IsZombie())
    throw std::runtime_error("Failed to open data input " +
                             std::string(dataFile));
  if (!mc || mc->IsZombie())
    throw std::runtime_error("Failed to open MC input " + std::string(mcFile));

  for (const auto &sample :
       std::vector<std::pair<TFile*,bool> >{{data.get(),false},
                                            {mc.get(),true}}) {
    TFile *source = sample.first;
    requireDirectory(source,"l2res");
    requireDirectory(source,"l2res1");
    const std::vector<std::string> requiredProfiles = {
      "p2m0tc", "p2m2tc", "p2mntc", "p2mutc", "p2restc", "p2jes",
      "p2m0pf", "p2m2pf", "p2mnpf", "p2mupf",
      "p2m0", "p2m2", "p2mn", "p2mu",
      "p2res", "p2chftc", "p2neftc", "p2nhftc", "p2ceftc", "p2muftc",
      "p2rhotc",
    };
    for (const std::string &name : requiredProfiles)
      requireObject<TProfile2D>(source,"l2res/" + name);
    if (sample.second)
      requireObject<TProfile2D>(source,"l2res/p2rgentc");
    requireObject<TH2D>(source,"l2res/h2ptetatc");
    requireObject<TH2D>(source,"l2res/h2ptetapf");
    requireObject<TH2D>(source,"l2res/h2pteta");
    requireObject<TH2D>(source,"l2res/h2mztc");
    if (requireObject<TProfile2D>(source,"l2res/p2m0tc")->GetEntries()<=0)
      throw std::runtime_error(
        "Response profile has no entries in " + std::string(source->GetName()));
    if (requireObject<TProfile2D>(source,"l2res/p2restc")->GetEntries()<=0)
      throw std::runtime_error(
        "Residual profile has no entries in " + std::string(source->GetName()));
  }

  TFile output(outputFile,"RECREATE");
  if (output.IsZombie())
    throw std::runtime_error("Failed to create output " +
                             std::string(outputFile));
  bool dataNativeProfiles = false;
  bool mcNativeProfiles = false;
  bool dataMethodL2Res = false;
  bool mcMethodL2Res = false;
  const bool dataFlavorInputs = copySample(
    data.get(),&output,"data",false,addFlavorPlaceholders,useLegacyMethod,
    preferOneDimensional,dataNativeProfiles,dataMethodL2Res);
  const bool mcFlavorInputs = copySample(
    mc.get(),&output,"mc",true,addFlavorPlaceholders,useLegacyMethod,
    preferOneDimensional,mcNativeProfiles,mcMethodL2Res);
  if (dataNativeProfiles!=mcNativeProfiles)
    throw std::runtime_error(
      "Native one-dimensional profiles exist in only one sample");
  if (dataFlavorInputs!=mcFlavorInputs)
    throw std::runtime_error(
      "Measured flavor inputs exist in only one of data and MC");
  if (dataMethodL2Res!=mcMethodL2Res)
    throw std::runtime_error(
      "Method-specific legacy L2Res inputs exist in only one sample");
  writeHDMCombination(&output,hdmResponseN,hdmResponseU);
  output.cd();
  TNamed method(
    "zjet_method",
    useLegacyMethod
      ? "legacy leading jet; synchronized dimuon selection; no transverse subtraction"
      : "all accepted Z-jet pairs; two transverse sidebands with half weight");
  method.Write();
  TNamed profileSource(
    "zjet_profile_source",
    dataNativeProfiles
      ? "native one-dimensional ZbAnalysis pT binning"
      : "central-eta projection of two-dimensional L2Res profiles");
  profileSource.Write();
  if (useLegacyMethod && dataFlavorInputs) {
    TNamed flavorMethod(
      "zjet_legacy_flavor_note",
      "inclusive inputs use legacy leading jet; flavor-tagged inputs remain from all-pairs method during synchronization");
    flavorMethod.Write();
  }
  if (useLegacyMethod) {
    TNamed l2resMethod(
      "zjet_legacy_l2res_note",
      dataMethodL2Res
        ? "legacy leading-jet p2jes and p2res profiles overlay l2res; l2res1 remains the all-pairs relative-eta control"
        : "WARNING: input predates method-specific legacy p2jes/p2res; l2res and l2res1 remain from all-pairs method");
    l2resMethod.Write();
  }
  TNamed compatibility(
    "jecsys3_compatibility",
    "reprocess.C, softrad3.C and globalFit.C central eta input contract");
  compatibility.Write();
  TNamed hdmDefinition(
    "zjet_hdm_definition",
    Form("softrad3 Z+jet master equation: R=(r0-rn-ru)/(1-rn/Rn-ru/Ru); Rn=%.6g; Ru=%.6g; ratio errors follow current softrad3 convention",
         hdmResponseN,hdmResponseU));
  hdmDefinition.Write();
  if (dataFlavorInputs) {
    TNamed flavorStatus(
      "zjet_flavor_status",
      "measured DeepJet reco-tag and matched generator-flavor inputs");
    flavorStatus.Write();
    if (TObjString *definition =
          dynamic_cast<TObjString*>(data->Get("zjet_flavor_definition")))
      definition->Write("zjet_flavor_definition",TObject::kOverwrite);
  }
  else if (addFlavorPlaceholders) {
    TNamed flavorStatus(
      "zjet_flavor_status",
      "PLACEHOLDER ONLY: empty reco-tag and generator-flavor inputs; do not use for physics");
    flavorStatus.Write();
    Warning("writeJecsys3",
            "Writing empty Z+flavor placeholders; they are not physics inputs");
  }
  if (TObjString *metadata =
        dynamic_cast<TObjString*>(data->Get("zjet_campaign_metadata"))) {
    metadata->Write("zjet_campaign_metadata",TObject::kOverwrite);
  }
  output.Close();
}
