#include "TDirectory.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TKey.h"
#include "TNamed.h"
#include "TObjString.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TROOT.h"

#include <stdexcept>
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

  TIter next(source->GetListOfKeys());
  while (TKey *key = static_cast<TKey*>(next())) {
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

void writeProfile(TFile *source, TDirectory *target,
                  const char *sourceName, const char *targetName,
                  int firstEtaBin, int lastEtaBin) {
  TProfile2D *profile = requireObject<TProfile2D>(
    source,"l2res/" + std::string(sourceName));
  TProfile *projection =
    profile->ProfileY(Form("tmp_%s",targetName),firstEtaBin,lastEtaBin);
  if (!projection)
    throw std::runtime_error(
      "Failed to project " + std::string(source->GetName()) + ":l2res/" +
      sourceName);
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

void writeGlobalFitInputs(TFile *source, TDirectory *sampleDirectory,
                          bool isMC) {
  TProfile2D *etaReference =
    requireObject<TProfile2D>(source,"l2res/p2m0tc");
  const int firstEtaBin = etaReference->GetXaxis()->FindFixBin(0.+1.e-6);
  const int lastEtaBin = etaReference->GetXaxis()->FindFixBin(1.305-1.e-6);

  TDirectory *etaDirectory = makeDirectory(sampleDirectory,"eta_00_13");
  // reprocess.C uses three complementary reference-pT choices. The tc, pf
  // and unsuffixed raw profiles are binned in Z pT, jet pT and average pT,
  // respectively. All three are projections of the same accepted pairs.
  writeResponseBinning(source,etaDirectory,"tc","zmmjet","h2ptetatc",
                       firstEtaBin,lastEtaBin);
  writeResponseBinning(source,etaDirectory,"pf","jetpt","h2ptetapf",
                       firstEtaBin,lastEtaBin);
  writeResponseBinning(source,etaDirectory,"","ptave","h2pteta",
                       firstEtaBin,lastEtaBin);

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

void copySample(TFile *source, TFile *target, const char *sample,
                bool isMC) {
  TDirectory *sampleDirectory = makeDirectory(target,sample);
  writeGlobalFitInputs(source,sampleDirectory,isMC);
  for (const char *name : {"l2res","l2res1"}) {
    TDirectory *sourceDirectory = requireDirectory(source,name);
    TDirectory *targetDirectory = makeDirectory(sampleDirectory,name);
    copyDirectory(sourceDirectory,targetDirectory);
  }
}

} // namespace

// Create one compatibility file that can be assigned to both ZMM_*_DATA and
// ZMM_*_MC in jecsys3/Config.C. It contains both the raw L2/L3 profile layout
// and the eta_00_13 inputs consumed by reprocess.C, softrad3.C and globalFit.C.
void writeJecsys3(
  const char *dataFile="rootfiles/zjet_DATA.root",
  const char *mcFile="rootfiles/zjet_MC.root",
  const char *outputFile="rootfiles/zjet_JMENANO_compat.root") {
  TFile data(dataFile,"READ");
  TFile mc(mcFile,"READ");
  if (data.IsZombie())
    throw std::runtime_error("Failed to open data input " +
                             std::string(dataFile));
  if (mc.IsZombie())
    throw std::runtime_error("Failed to open MC input " + std::string(mcFile));

  for (const auto &sample :
       std::vector<std::pair<TFile*,bool> >{{&data,false},{&mc,true}}) {
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
  copySample(&data,&output,"data",false);
  copySample(&mc,&output,"mc",true);
  output.cd();
  TNamed method("zjet_method",
                "all accepted Z-jet pairs; two transverse sidebands with half weight");
  method.Write();
  TNamed compatibility(
    "jecsys3_compatibility",
    "reprocess.C, softrad3.C and globalFit.C central eta input contract");
  compatibility.Write();
  if (TObjString *metadata =
        dynamic_cast<TObjString*>(data.Get("zjet_campaign_metadata"))) {
    metadata->Write("zjet_campaign_metadata",TObject::kOverwrite);
  }
  output.Write();
  output.Close();
}
