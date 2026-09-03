#include "zjet.h"

#include <TChain.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TObjString.h>
#include <TKey.h>
#include <TSystem.h>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <string>

R__LOAD_LIBRARY(CondFormats/JetMETObjects/src/Utilities_cc.so)
R__LOAD_LIBRARY(CondFormats/JetMETObjects/src/JetCorrectorParameters_cc.so)
R__LOAD_LIBRARY(CondFormats/JetMETObjects/src/SimpleJetCorrector_cc.so)
R__LOAD_LIBRARY(CondFormats/JetMETObjects/src/FactorizedJetCorrector_cc.so)
R__LOAD_LIBRARY(zjet_C.so)

namespace {

int addInputList(TChain *chain, const char *inputList) {
  std::ifstream stream(inputList);
  if (!stream.is_open()) {
    std::cerr << "Could not open input list " << inputList << std::endl;
    return 0;
  }

  int added = 0;
  std::string line;
  while (std::getline(stream,line)) {
    const std::string::size_type first = line.find_first_not_of(" \t\r");
    if (first==std::string::npos || line[first]=='#') continue;
    const std::string::size_type last = line.find_last_not_of(" \t\r");
    added += chain->AddFile(line.substr(first,last-first+1).c_str());
  }
  return added;
}

int keyCycles(TDirectory *directory, const char *name) {
  if (!directory || !directory->GetListOfKeys()) return 0;
  int cycles = 0;
  TIter next(directory->GetListOfKeys());
  while (TKey *key = dynamic_cast<TKey*>(next()))
    if (std::string(key->GetName())==name) ++cycles;
  return cycles;
}

} // namespace

void run_zjet_job(const char *inputList, bool isMC, const char *outputFile,
                  const char *goldenJson="", const char *lumiPileup="",
                  const char *pileupWeights="", const char *jecL2="",
                  const char *jecResidual="", const char *jerResolution="",
                  const char *jerScaleFactor="",
                  const char *muonCorrections="",
                  const char *jetVetoMap="") {
  TChain *chain = new TChain("Events","Events");
  const int added = addInputList(chain,inputList);
  if (added<=0) {
    std::cerr << "ERROR: no input files were added from " << inputList
              << std::endl;
    gSystem->Exit(2);
    return;
  }

  std::cout << "Added " << added << " file(s) for the "
            << (isMC ? "MC" : "data") << " batch job." << std::endl;
  zjet analysis(chain,isMC,outputFile,goldenJson,lumiPileup,pileupWeights,
                jecL2,jecResidual,jerResolution,jerScaleFactor,
                muonCorrections,jetVetoMap);
  analysis.Loop();

  TFile check(outputFile,"READ");
  TH2 *inclusiveCounts =
    dynamic_cast<TH2*>(check.Get("l2res/h2ptetatc"));
  TH1 *flavorCounts =
    dynamic_cast<TH1*>(check.Get("flavor/counts_gii"));
  TH1 *nativeCounts =
    dynamic_cast<TH1*>(check.Get("profiles1d/zmmjet/statistics_rmpf"));
  TH1 *nativeMpfNu =
    dynamic_cast<TH1*>(check.Get("profiles1d/zmmjet/rmpfjetnu"));
  TH1 *legacyCounts = dynamic_cast<TH1*>(
    check.Get("legacy/profiles1d/zmmjet/statistics_rmpf"));
  TH1 *legacyMpfNu =
    dynamic_cast<TH1*>(check.Get("legacy/profiles1d/zmmjet/rmpfjetnu"));
  TObjString *jerResolutionMetadata = dynamic_cast<TObjString*>(
    check.Get("zjet_jer_resolution_file"));
  TObjString *jerScaleFactorMetadata = dynamic_cast<TObjString*>(
    check.Get("zjet_jer_scale_factor_file"));
  TObjString *muonCorrectionMetadata = dynamic_cast<TObjString*>(
    check.Get("zjet_muon_correction_file"));
  TObjString *jetVetoMapMetadata = dynamic_cast<TObjString*>(
    check.Get("zjet_jet_veto_map_file"));
  const std::string expectedJerResolution = isMC ? jerResolution : "";
  const std::string expectedJerScaleFactor = isMC ? jerScaleFactor : "";
  const std::string expectedJetVetoMap = isMC ? "" : jetVetoMap;
  const bool correctionMetadataMatches =
    (jerResolutionMetadata && jerScaleFactorMetadata &&
     muonCorrectionMetadata && jetVetoMapMetadata &&
     jerResolutionMetadata->GetString()==expectedJerResolution.c_str() &&
     jerScaleFactorMetadata->GetString()==expectedJerScaleFactor.c_str() &&
     muonCorrectionMetadata->GetString()==muonCorrections &&
     jetVetoMapMetadata->GetString()==expectedJetVetoMap.c_str());
  bool metadataHasSingleCycles = true;
  for (const char *name : {
         "zjet_jec_mode", "zjet_jec_l2_file", "zjet_jec_residual_file",
         "zjet_jer_resolution_file", "zjet_jer_scale_factor_file",
         "zjet_muon_correction_file", "zjet_muon_correction_sha256",
         "zjet_jet_veto_map_file", "zjet_type1_met_definition",
         "zjet_flavor_definition", "zjet_flavor_matrix_definition",
         "zjet_synchronized_selection", "zjet_legacy_jet_id",
         "zjet_truth_hdm_definition", "zjet_previous_residual_definition",
       })
    metadataHasSingleCycles =
      metadataHasSingleCycles && keyCycles(&check,name)==1;
  bool flavorCountsClose = false;
  if (inclusiveCounts && flavorCounts) {
    const int firstEtaBin =
      inclusiveCounts->GetXaxis()->FindFixBin(0.+1.e-6);
    const int lastEtaBin =
      inclusiveCounts->GetXaxis()->FindFixBin(1.305-1.e-6);
    const double inclusiveIntegral = inclusiveCounts->Integral(
      firstEtaBin,lastEtaBin,1,inclusiveCounts->GetNbinsY());
    const double flavorIntegral = flavorCounts->Integral();
    flavorCountsClose =
      (std::fabs(inclusiveIntegral-flavorIntegral) <=
       1.e-9*std::max(1.,std::fabs(inclusiveIntegral)));
  }
  if (check.IsZombie() || !check.Get("control/h_cutflow") ||
      !flavorCounts || !check.Get("zjet_flavor_definition") ||
      !nativeCounts || !nativeMpfNu || !legacyCounts || !legacyMpfNu ||
      !check.Get("legacy/control/h_cutflow") ||
      !check.Get("legacy/control/h_alpha") ||
      !check.Get("legacy/control/h_subleading_jetpt") ||
      !check.Get("legacy/control/h_alpha_vs_jetpt") ||
      !check.Get("legacy/control/p_db_vs_jetpt_before_alpha") ||
      !check.Get("legacy/control/p_mpf1_vs_jetpt_before_alpha") ||
      !check.Get("legacy/control/p_rho_vs_zpt_before_alpha") ||
      !check.Get("legacy/control/p_rho_vs_zpt_alpha010") ||
      !check.Get("control/p_mpf_vs_ptz_all_subtracted_central") ||
      !check.Get("control/p_mpfnu_vs_abseta_matched_parallel_ptz15to30") ||
      !check.Get("truth_hdm/parallel/zmmjet/reco_over_gen") ||
      !check.Get("truth_hdm/parallel/zmmjet/reco_mpfn_matched") ||
      !check.Get("truth_hdm/parallel/zmmjet/reco_mpfu_matched") ||
      !check.Get("truth_hdm/subtracted/zmmjet/gen_mpfn_reco_axis") ||
      !check.Get("truth_hdm/subtracted/zmmjet/gen_mpfu_reco_axis") ||
      !check.Get("truth_hdm/subtracted/zmmjet/response_rn_reco_axis") ||
      !check.Get("truth_hdm/subtracted/zmmjet/response_ru_reco_axis") ||
      !check.Get("truth_hdm/subtracted/zmmjet/mpfn_reco_gen_product") ||
      !check.Get("truth_hdm/subtracted/zmmjet/mpfu_gen_squared") ||
      !check.Get("truth_hdm/subtracted/zmmjet/slope_rn_reco_axis") ||
      !check.Get("truth_hdm/subtracted/zmmjet/slope_ru_reco_axis") ||
      !check.Get("legacy/truth_hdm/parallel/zmmjet/reco_over_gen") ||
      !check.Get("legacy/l2res/p2res") ||
      !check.Get("legacy/l2res/p2restc") ||
      !check.Get("FlavorMatrix/h3counts_flavormatrix") ||
      !check.Get("FlavorMatrix/h3counts_parallel_flavormatrix") ||
      !check.Get("FlavorMatrix/h3counts_transverse_flavormatrix") ||
      !check.Get("FlavorMatrix/p3hdmtc_flavormatrix") ||
      !check.Get("FlavorMatrix/p3mnuab_flavormatrix") ||
      !check.Get("FlavorMatrix/p3hdmtc_parallel_flavormatrix") ||
      !check.Get("FlavorMatrix/p3mnuab_parallel_flavormatrix") ||
      !check.Get("FlavorMatrix/p3areasumtc_parallel_flavormatrix") ||
      !check.Get("FlavorMatrix/p3areaprojtc_parallel_flavormatrix") ||
      !check.Get("FlavorMatrix/p3ueholetc_parallel_flavormatrix") ||
      !check.Get("FlavorMatrix/p3mnufsrtc_parallel_flavormatrix") ||
      !check.Get("FlavorMatrix/p3genmutc_parallel_flavormatrix") ||
      !check.Get("FlavorMatrix/p3recoumatchedtc_parallel_flavormatrix") ||
      !check.Get("FlavorMatrix/p3recogenmutc_parallel_flavormatrix") ||
      !check.Get("FlavorMatrix/p3genmu2tc_parallel_flavormatrix") ||
      !check.Get("FlavorMatrix/h3counts_heavytopology") ||
      !check.Get("FlavorMatrix/controls/h3_cvb_cvl_trueflavor") ||
      !check.Get("FlavorMatrix/controls/h3_cvb_cvl_qvg_true0") ||
      !check.Get("FlavorMatrix/controls/h3_genjet_nc_nb_trueflavor") ||
      !check.Get("FlavorMatrix/controls/h3_upartqvg_pnetqvg_trueflavor") ||
      !check.Get("FlavorMatrix/controls/h3_cvb_cvl_heavytopology") ||
      !check.Get("zjet_truth_hdm_definition") ||
      !check.Get("zjet_previous_residual_definition") ||
      !check.Get("zjet_synchronized_selection") ||
      !check.Get("zjet_legacy_jet_id") ||
      !correctionMetadataMatches ||
      !check.Get("zjet_muon_correction_sha256") ||
      !check.Get("zjet_type1_met_definition") ||
      !metadataHasSingleCycles ||
      !flavorCountsClose) {
    std::cerr << "ERROR: output validation failed for " << outputFile
              << std::endl;
    gSystem->Exit(3);
    return;
  }
  std::cout << "Validated output " << outputFile
            << ": flavor, method-specific residual, generator-recoil, "
            << "native 1D all-pairs, and synchronized legacy profiles are "
            << "present." << std::endl;
}
