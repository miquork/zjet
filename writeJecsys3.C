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

#include <cassert>
#include <string>
#include <utility>
#include <vector>

namespace {

void copyDirectory(TDirectory *source, TDirectory *target) {
  assert(source);
  assert(target);

  TIter next(source->GetListOfKeys());
  while (TKey *key = static_cast<TKey*>(next())) {
    TClass *objectClass = gROOT->GetClass(key->GetClassName());
    assert(objectClass);
    if (objectClass->InheritsFrom(TDirectory::Class())) {
      TDirectory *sourceChild = source->GetDirectory(key->GetName());
      TDirectory *targetChild = target->mkdir(key->GetName());
      copyDirectory(sourceChild,targetChild);
    }
    else {
      TObject *object = key->ReadObj();
      target->cd();
      object->Write(key->GetName());
      delete object;
    }
  }
}

void writeProfile(TFile *source, TDirectory *target,
                  const char *sourceName, const char *targetName,
                  int firstEtaBin, int lastEtaBin) {
  TProfile2D *profile =
    dynamic_cast<TProfile2D*>(source->Get(Form("l2res/%s",sourceName)));
  assert(profile);
  TProfile *projection =
    profile->ProfileY(Form("tmp_%s",targetName),firstEtaBin,lastEtaBin);
  assert(projection);
  projection->SetDirectory(nullptr);
  target->cd();
  projection->Write(targetName);
  delete projection;
}

void writeCounts(TFile *source, TDirectory *target,
                 int firstEtaBin, int lastEtaBin) {
  TH2D *counts = dynamic_cast<TH2D*>(source->Get("l2res/h2ptetatc"));
  assert(counts);
  TH1D *projection = counts->ProjectionY(
    "statistics_rmpf_zmmjet_a100",firstEtaBin,lastEtaBin,"e");
  assert(projection);
  projection->SetDirectory(nullptr);
  target->cd();
  projection->Write();
  delete projection;
}

void writeGlobalFitInputs(TFile *source, TDirectory *sampleDirectory,
                          bool isMC) {
  TProfile2D *etaReference =
    dynamic_cast<TProfile2D*>(source->Get("l2res/p2m0tc"));
  assert(etaReference);
  const int firstEtaBin = etaReference->GetXaxis()->FindFixBin(0.+1.e-6);
  const int lastEtaBin = etaReference->GetXaxis()->FindFixBin(1.305-1.e-6);

  TDirectory *etaDirectory = sampleDirectory->mkdir("eta_00_13");
  assert(etaDirectory);
  const std::vector<std::pair<std::string,std::string> > profiles = {
    {"p2m0tc", "rmpf_zmmjet_a100"},
    {"p2m2tc", "rmpfjet1_zmmjet_a100"},
    {"p2mntc", "rmpfjetn_zmmjet_a100"},
    {"p2mutc", "rmpfuncl_zmmjet_a100"},
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
  writeCounts(source,etaDirectory,firstEtaBin,lastEtaBin);

  TH2D *mass = dynamic_cast<TH2D*>(source->Get("l2res/h2mztc"));
  assert(mass);
  etaDirectory->cd();
  mass->Write("h_Zpt_mZ_alpha100");
}

void copySample(TFile *source, TFile *target, const char *sample,
                bool isMC) {
  TDirectory *sampleDirectory = target->mkdir(sample);
  assert(sampleDirectory);
  for (const char *name : {"l2res","l2res1"}) {
    TDirectory *sourceDirectory = source->GetDirectory(name);
    assert(sourceDirectory);
    TDirectory *targetDirectory = sampleDirectory->mkdir(name);
    copyDirectory(sourceDirectory,targetDirectory);
  }
  writeGlobalFitInputs(source,sampleDirectory,isMC);
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
  assert(!data.IsZombie());
  assert(!mc.IsZombie());

  for (TFile *source : {&data,&mc}) {
    TProfile2D *response = dynamic_cast<TProfile2D*>(source->Get("l2res/p2m0tc"));
    TProfile2D *residual = dynamic_cast<TProfile2D*>(source->Get("l2res/p2restc"));
    assert(response && response->GetEntries()>0);
    assert(residual && residual->GetEntries()>0);
    assert(source->Get("l2res/p2jes"));
    assert(source->Get("l2res/p2res"));
    assert(source->Get("l2res/p2chftc"));
    assert(source->Get("l2res/h2mztc"));
  }

  TFile output(outputFile,"RECREATE");
  assert(!output.IsZombie());
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
