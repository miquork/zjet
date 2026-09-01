#include "TFile.h"
#include "TKey.h"
#include "TObjString.h"
#include "TProfile3D.h"

#include "FlavorMatrixTools.h"

#include <cassert>
#include <fstream>
#include <iterator>
#include <map>
#include <memory>
#include <stdexcept>
#include <string>

// Store the exact campaign provenance alongside the histograms in a ROOT file.
void embedCampaignMetadata(const char *rootFile, const char *metadataFile) {
  std::ifstream input(metadataFile);
  assert(input.good());
  const std::string contents((std::istreambuf_iterator<char>(input)),
                             std::istreambuf_iterator<char>());
  assert(!contents.empty());

  TFile output(rootFile,"UPDATE");
  assert(!output.IsZombie());
  output.cd();

  // Recompute nonlinear HDM only after hadd has merged the linear component
  // sums. Worker-level HDM profiles are diagnostics; this pass makes the
  // final campaign object equal to HDM(<m0>,<mn>,<mu>) bin by bin.
  for (const char *variant : {"ab","ad","tc","pf"}) {
    const std::string suffix = std::string(variant)+"_flavormatrix";
    TProfile3D *hdm = dynamic_cast<TProfile3D*>(
      output.Get(("FlavorMatrix/p3hdm"+suffix).c_str()));
    TProfile3D *m0 = dynamic_cast<TProfile3D*>(
      output.Get(("FlavorMatrix/p3m0"+suffix).c_str()));
    TProfile3D *mn = dynamic_cast<TProfile3D*>(
      output.Get(("FlavorMatrix/p3mn"+suffix).c_str()));
    TProfile3D *mu = dynamic_cast<TProfile3D*>(
      output.Get(("FlavorMatrix/p3mu"+suffix).c_str()));
    ZJetFlavorMatrix::finalizeHDMProfile(hdm,m0,mn,mu);
    output.cd("FlavorMatrix");
    hdm->Write(hdm->GetName(),TObject::kOverwrite);
    output.cd();
  }

  // hadd keeps one key cycle per worker for non-mergeable TObjString
  // metadata. Require every worker value to agree, remove every old cycle,
  // and rewrite each item once so TBrowser shows a clean ;1 key instead of
  // dozens of identical cycles.
  std::map<std::string,std::string> metadataStrings;
  TIter next(output.GetListOfKeys());
  while (TKey *key = dynamic_cast<TKey*>(next())) {
    if (std::string(key->GetClassName())!="TObjString") continue;
    const std::string name = key->GetName();
    std::unique_ptr<TObject> object(key->ReadObj());
    TObjString *value = dynamic_cast<TObjString*>(object.get());
    if (!value)
      throw std::runtime_error("Failed to read metadata key " + name);
    const std::string contents = value->GetString().Data();
    const auto existing = metadataStrings.find(name);
    if (existing!=metadataStrings.end() && existing->second!=contents)
      throw std::runtime_error(
        "Conflicting worker metadata values for " + name +
        "; refusing to merge inconsistent correction configurations");
    metadataStrings[name] = contents;
  }
  for (const auto &entry : metadataStrings) {
    output.Delete((entry.first+";*").c_str());
    TObjString value(entry.second.c_str());
    value.Write(entry.first.c_str(),TObject::kOverwrite);
  }

  TObjString metadata(contents.c_str());
  metadata.Write("zjet_campaign_metadata",TObject::kOverwrite);
  // The metadata key is already persistent after Write(). A following
  // TFile::Write() would create a second cycle, and hadd would multiply it by
  // the number of worker files. Purge any pre-existing duplicate cycles from
  // the merged file before it is uploaded.
  output.Purge();
  for (const auto &entry : metadataStrings) {
    int cycles = 0;
    TIter verify(output.GetListOfKeys());
    while (TKey *key = dynamic_cast<TKey*>(verify()))
      if (entry.first==key->GetName()) ++cycles;
    assert(cycles==1);
  }
  int campaignCycles = 0;
  TIter verifyCampaign(output.GetListOfKeys());
  while (TKey *key = dynamic_cast<TKey*>(verifyCampaign()))
    if (std::string(key->GetName())=="zjet_campaign_metadata")
      ++campaignCycles;
  assert(campaignCycles==1);
  output.Close();
}
