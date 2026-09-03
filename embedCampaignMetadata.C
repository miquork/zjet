#include "TFile.h"
#include "TGraphErrors.h"
#include "TKey.h"
#include "TObjString.h"
#include "TProfile.h"
#include "TProfile3D.h"

#include "FlavorMatrixTools.h"

#include <cassert>
#include <cmath>
#include <fstream>
#include <iterator>
#include <map>
#include <memory>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

namespace {

void rebuildTruthGraph(TFile &output, const std::string &directoryPath,
                       const char *name, const char *numeratorName,
                       const char *denominatorName, bool slope) {
  TDirectory *directory = output.GetDirectory(directoryPath.c_str());
  if (!directory) return;
  TProfile *numerator = dynamic_cast<TProfile*>(
    directory->Get(numeratorName));
  TProfile *denominator = dynamic_cast<TProfile*>(
    directory->Get(denominatorName));
  if (!numerator || !denominator) return;
  std::unique_ptr<TGraphErrors> graph(new TGraphErrors());
  graph->SetName(name);
  graph->SetTitle(slope
    ? ";p_{T} (GeV);Zero-intercept reco-versus-gen slope"
    : ";p_{T} (GeV);Ratio of reco and generator component means");
  for (int bin=1; bin<=numerator->GetNbinsX(); ++bin) {
    const double n = numerator->GetBinContent(bin);
    const double d = denominator->GetBinContent(bin);
    const double minimum = slope ? 1.e-12 : 1.e-9;
    if (numerator->GetBinEntries(bin)==0. ||
        denominator->GetBinEntries(bin)==0. ||
        !std::isfinite(n) || !std::isfinite(d) ||
        std::fabs(d)<minimum)
      continue;
    const double value = n/d;
    const double error = std::hypot(
      numerator->GetBinError(bin)/d,
      n*denominator->GetBinError(bin)/(d*d));
    const int point = graph->GetN();
    graph->SetPoint(point,numerator->GetBinCenter(bin),value);
    graph->SetPointError(
      point,0.5*numerator->GetBinWidth(bin),std::fabs(error));
  }
  directory->cd();
  directory->Delete((std::string(name)+";*").c_str());
  graph->Write(name,TObject::kOverwrite);
  output.cd();
}

void finalizeTruthGraphs(TFile &output) {
  const std::vector<std::string> bases = {
    "truth_hdm/parallel", "truth_hdm/transverse", "truth_hdm/subtracted",
    "legacy/truth_hdm/parallel",
  };
  for (const std::string &base : bases)
    for (const char *axis : {"zmmjet","jetpt","ptave"}) {
      const std::string directory = base+"/"+axis;
      for (const auto &definition : std::vector<
             std::tuple<const char*,const char*,const char*> >{
             {"response_r1_reco_axis","reco_mpf1_matched",
              "gen_mpf1_reco_axis"},
             {"response_rn_reco_axis","reco_mpfn_matched",
              "gen_mpfn_reco_axis"},
             {"response_ru_reco_axis","reco_mpfu_matched",
              "gen_mpfu_reco_axis"},
             {"closure_r1_gen_axis","reco_mpf1_matched",
              "gen_mpf1_gen_axis"},
             {"closure_rn_gen_axis","reco_mpfn_matched",
              "gen_mpfn_gen_axis"},
             {"closure_ru_gen_axis","reco_mpfu_matched",
              "gen_mpfu_gen_axis"}})
        rebuildTruthGraph(output,directory,std::get<0>(definition),
                          std::get<1>(definition),std::get<2>(definition),
                          false);
      for (const auto &definition : std::vector<
             std::tuple<const char*,const char*,const char*> >{
             {"slope_r1_reco_axis","mpf1_reco_gen_product",
              "mpf1_gen_squared"},
             {"slope_rn_reco_axis","mpfn_reco_gen_product",
              "mpfn_gen_squared"},
             {"slope_ru_reco_axis","mpfu_reco_gen_product",
              "mpfu_gen_squared"}})
        rebuildTruthGraph(output,directory,std::get<0>(definition),
                          std::get<1>(definition),std::get<2>(definition),
                          true);
    }
}

} // namespace

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
    for (const char *family : {
           "_flavormatrix", "_parallel_flavormatrix"}) {
      const std::string suffix = std::string(variant)+family;
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
  }

  // TGraphErrors are not hadd-mergeable: hadd concatenates the points from
  // every worker. Rebuild all truth-derived ratio and slope graphs from the
  // already merged TProfiles so every pT bin occurs exactly once.
  finalizeTruthGraphs(output);

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
