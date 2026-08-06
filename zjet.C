#define zjet_cxx
#include "zjet.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

#include "TLorentzVector.h"
#include "TH2D.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TSystem.h"

#include "ZJetLumi.h"

#include <chrono>
#include <cmath>
#include <ctime>
#include <iomanip>
#include <map>
#include <string>

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

   fChain->SetBranchStatus("nJet",1);
   fChain->SetBranchStatus("Jet_pt",1);
   fChain->SetBranchStatus("Jet_eta",1);
   fChain->SetBranchStatus("Jet_phi",1);
   fChain->SetBranchStatus("Jet_mass",1);
   fChain->SetBranchStatus("Jet_rawFactor",1);
   fChain->SetBranchStatus("Jet_chMultiplicity",1);
   fChain->SetBranchStatus("Jet_neMultiplicity",1);
   fChain->SetBranchStatus("Jet_nConstituents",1);
   fChain->SetBranchStatus("Jet_chHEF",1);
   fChain->SetBranchStatus("Jet_neHEF",1);
   fChain->SetBranchStatus("Jet_neEmEF",1);

   fChain->SetBranchStatus("PV_npvs",1);
   fChain->SetBranchStatus("Rho_fixedGridRhoFastjetAll",1);

   fChain->SetBranchStatus("Flag_goodVertices",1);
   fChain->SetBranchStatus("Flag_globalSuperTightHalo2016Filter",1);
   fChain->SetBranchStatus("Flag_EcalDeadCellTriggerPrimitiveFilter",1);
   fChain->SetBranchStatus("Flag_BadPFMuonFilter",1);
   fChain->SetBranchStatus("Flag_BadPFMuonDzFilter",1);
   fChain->SetBranchStatus("Flag_hfNoisyHitsFilter",1);
   fChain->SetBranchStatus("Flag_eeBadScFilter",1);
   fChain->SetBranchStatus("Flag_ecalBadCalibFilter",1);
   fChain->SetBranchStatus("HLT_IsoMu24",1);

   if (isMC) {
     fChain->SetBranchStatus("genWeight",1);
     fChain->SetBranchStatus("Pileup_nTrueInt",1);
     fChain->SetBranchStatus("nGenJet",1);
     fChain->SetBranchStatus("GenJet_pt",1);
     fChain->SetBranchStatus("GenJet_eta",1);
     fChain->SetBranchStatus("GenJet_phi",1);
     fChain->SetBranchStatus("GenJet_mass",1);
     fChain->SetBranchStatus("Jet_genJetIdx",1);
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

   
   // Object pT plots
   fout->mkdir("control");
   fout->cd("control");
   TH1D *h_cutflow = new TH1D("h_cutflow","",7,0.5,7.5);
   h_cutflow->GetXaxis()->SetBinLabel(1,"all");
   h_cutflow->GetXaxis()->SetBinLabel(2,"golden JSON");
   h_cutflow->GetXaxis()->SetBinLabel(3,"MET filters");
   h_cutflow->GetXaxis()->SetBinLabel(4,"HLT IsoMu24");
   h_cutflow->GetXaxis()->SetBinLabel(5,"tag-probe muons");
   h_cutflow->GetXaxis()->SetBinLabel(6,"Z mass");
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
   TH1D *h_zmass_precut = new TH1D("h_zmass_precut","",300,75,105);
   TH1D *h_zeta_precut = new TH1D("h_zeta_precut","",100,-5,5);
   TH2D *h2_zeta_precut_vs_ptz = new TH2D("h2_zeta_precut_vs_ptz","",200,0,200,100,-5,5);
   TH2D *h2_zmass_precut_vs_pt = new TH2D("h2_zmass_precut_vs_pt","",200,0,200,300,75,105);
   TProfile *p_zmass_precut_vs_pt = new TProfile("p_zmass_precut_vs_pt","",200,0,200);
   
   TH1D *h_zpt = new TH1D("h_zpt","",200,0,200);
   TH1D *h_zeta = new TH1D("h_zeta","",100,-5,5);
   TH2D *h2_zeta_vs_ptz = new TH2D("h2_zeta_vs_ptz","",200,0,200,100,-5,5);
   TH1D *h_zmass = new TH1D("h_zmass","",300,75,105);
   TH2D *h2_zmass_vs_pt = new TH2D("h2_zmass_vs_pt","",200,0,200,300,75,105);
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
   
   curdir->cd();

   TLorentzVector p4lplus, p4lminus, p4z, p4jet1, p4jet, p4sel1, p4tran1;
   TLorentzVector p4p, p4pz, p4t1, p4t1z, p4t2, p4t2z;
   TLorentzVector met, ht, met1, metn, metu, metnu, meta;

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
      const double mu = (isMC ? Pileup_nTrueInt
                              : lumiData.pileup(run, luminosityBlock));

      h_cutflow->Fill(1., eventWeight);
      if (!isMC && !lumiData.accept(run, luminosityBlock)) continue;
      h_cutflow->Fill(2., eventWeight);

      const bool passMetFilters =
        (Flag_goodVertices && Flag_globalSuperTightHalo2016Filter &&
         Flag_EcalDeadCellTriggerPrimitiveFilter && Flag_BadPFMuonFilter &&
         Flag_BadPFMuonDzFilter && Flag_hfNoisyHitsFilter &&
         Flag_eeBadScFilter && Flag_ecalBadCalibFilter);
      if (!passMetFilters) continue;
      h_cutflow->Fill(3., eventWeight);

      if (!HLT_IsoMu24) continue;
      h_cutflow->Fill(4., eventWeight);

      if (nMuon>nMuonMax || nJet>nJetMax) {
        cout << "ERROR: collection size exceeds fixed MakeClass buffer: nMuon="
             << nMuon << ", nJet=" << nJet << endl;
        continue;
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
      int nlep(0);
      double nsel(0.);
      double ntran(0.);

      // Scan all opposite-sign pairs passing a configurable working point and
      // keep the pair closest to the Z mass. The tag requirement models IsoMu24
      // path: one tight, tightly isolated muon above the offline plateau;
      // the other muon can be a lower-pT medium-ID, loose-isolation probe.
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
                  Muon_pt[ilep]>minimumPt && fabs(Muon_eta[ilep])<2.4);
        };
        auto isTag = [&](int ilep) {
          return (Muon_pt[ilep]>27. && Muon_tightId[ilep] &&
                  Muon_pfIsoId[ilep]>=4);
        };
        for (int iplus = 0; iplus != nMuon; ++iplus) {
          if (Muon_charge[iplus]<=0 || !passMuon(iplus)) continue;
          TLorentzVector plusCandidate;
          plusCandidate.SetPtEtaPhiM(Muon_pt[iplus],Muon_eta[iplus],
                                     Muon_phi[iplus],Muon_mass[iplus]);
          for (int iminus = 0; iminus != nMuon; ++iminus) {
            if (Muon_charge[iminus]>=0 || !passMuon(iminus)) continue;
            TLorentzVector minusCandidate;
            minusCandidate.SetPtEtaPhiM(Muon_pt[iminus],Muon_eta[iminus],
                                        Muon_phi[iminus],Muon_mass[iminus]);
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

      // Select the nominal tag-probe pair.
      h_nlep->Fill(nMuon, eventWeight);
      if (!selectMuonPair(2,2,10.,27.,true,p4lplus,p4lminus)) continue;

      // Reconstruct Z boson
      if (p4lplus.Pt()>0) ++nlep;
      if (p4lminus.Pt()>0) ++nlep;
      if (p4lplus.Pt()>0 && p4lminus.Pt()>0) {
	p4z += p4lplus; p4z += p4lminus;
      }
      if (nlep != 2) continue;
      h_cutflow->Fill(5., eventWeight);
      //if (p4z.Pt()>0 && p4z.M()>80 && p4z.M()<100) {
      if (p4z.Pt()>0 && p4z.M()>75 && p4z.M()<105) {
	h_zpt_precut->Fill(p4z.Pt());
	h_zeta_precut->Fill(p4z.Eta());
	h2_zeta_precut_vs_ptz->Fill(p4z.Pt(), p4z.Eta());
	h_zmass_precut->Fill(p4z.M());
	h2_zmass_precut_vs_pt->Fill(p4z.Pt(),p4z.M());
	p_zmass_precut_vs_pt->Fill(p4z.Pt(),p4z.M());
      }
      if (p4z.Pt()>0 && fabs(p4z.M()-mz)<dmz) {
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
      h_probe_veto->Fill(1., eventWeight);

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
      met.SetPtEtaPhiM(PuppiMET_pt, 0., PuppiMET_phi, 0.0);
      ht += p4z;
      for (int ijet = 0; ijet != nJet; ++ijet) {
	p4jet.SetPtEtaPhiM(Jet_pt[ijet], Jet_eta[ijet], Jet_phi[ijet],
			   Jet_mass[ijet]);
	//if (p4jet.DeltaR(p4lplus)>0.4 && p4jet.DeltaR(p4lminus)>0.4 &&
	if (passTightJetId(ijet) &&
            p4jet.DeltaR(p4lplus)>0.2 && p4jet.DeltaR(p4lminus)>0.2 &&
	    p4jet.Pt()>15.) {
	  ht += p4jet;
	}
      }
      ht.SetPtEtaPhiM(ht.Pt(),0,ht.Phi(),0);
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
      
      // Select leading jet
      h_njet->Fill(nJet, eventWeight);
      for (int ijet = 0; ijet != nJet; ++ijet) {

	p4jet.SetPtEtaPhiM(Jet_pt[ijet], Jet_eta[ijet], Jet_phi[ijet],
			   Jet_mass[ijet]);
	if (!passTightJetId(ijet)) continue;

	double eta = p4jet.Eta();
	double ptz = p4z.Pt();
	double ptj = p4jet.Pt();
	double pta = 0.5*(ptz+ptj);

	double jes = (1-Jet_rawFactor[ijet]);
	bool hasGenIndex = false;
	bool truthMatched = false;
	bool passesExtraMatchCuts = false;
	int truthMatchCategory = 0;
	if (isMC && Jet_genJetIdx[ijet]>=0 && Jet_genJetIdx[ijet]<nGenJet) {
	  hasGenIndex = true;
	  // Jet_genJetIdx is the NanoAOD reco-to-particle-level match. Do not
	  // impose a second generator-pT cut on the nominal pileup classification:
	  // it creates a strong migration bias precisely in the low-pT region.
	  const int igen = Jet_genJetIdx[ijet];
	  TLorentzVector p4gen;
	  p4gen.SetPtEtaPhiM(GenJet_pt[igen],GenJet_eta[igen],GenJet_phi[igen],
			     GenJet_mass[igen]);
	  truthMatched = true;
	  if (p4gen.Pt()<=8.) truthMatchCategory = 1;
	  else if (p4jet.DeltaR(p4gen)>=0.4) truthMatchCategory = 2;
	  else {
	    truthMatchCategory = 3;
	    passesExtraMatchCuts = true;
	  }
	}

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
	if (p4jet.DeltaR(p4lplus)>0.2 && p4jet.DeltaR(p4lminus)>0.2) {

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

		    h2ptetapf_->Fill(abseta,ptj,wt); h2pteta_->Fill(abseta,pta,wt);
		    h2ptetatc_->Fill(abseta,ptz,wt);
		    h2ptetapf->Fill(eta,ptj,wt); h2pteta->Fill(eta,pta,wt);
		    h2ptetatc->Fill(eta,ptz,wt);
		    pmzpf_->Fill(ptj,p4z.M(),wt); pmz_->Fill(pta,p4z.M(),wt);
		    pmztc_->Fill(ptz,p4z.M(),wt); pmzpf->Fill(ptj,p4z.M(),wt);
		    pmz->Fill(pta,p4z.M(),wt); pmztc->Fill(ptz,p4z.M(),wt);

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
		    p2respf_->Fill(abseta,ptj,1.,wt); p2res_->Fill(abseta,pta,1.,wt);
		    p2restc_->Fill(abseta,ptz,1.,wt); p2respf->Fill(eta,ptj,1.,wt);
		    p2res->Fill(eta,pta,1.,wt); p2restc->Fill(eta,ptz,1.,wt);
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
		    p2respf_->Fill(abseta,ptj,1.,wt); p2res_->Fill(abseta,pta,1.,wt);
		    p2restc_->Fill(abseta,ptz,1.,wt); p2respf->Fill(eta,ptj,1.,wt);
		    p2res->Fill(eta,pta,1.,wt); p2restc->Fill(eta,ptz,1.,wt);
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
   
   cout << endl << "Finished loop, writing file " << outputFile << "." << endl << flush;
    cout << "Processed " << nentries << " events\n";
    //cout << "Skipped " << _nbadevents_json << " events due to JSON ("
    //	 << (100.*_nbadevents_json/_nevents) << "%) \n";
    //cout << "Skipped " << _nbadevents_trigger << " events due to trigger ("
    //	 << (100.*_nbadevents_trigger/_ntot) << "%) \n";
    //cout << "Skipped " << _nbadevents_veto << " events due to veto ("
    //	 << (100.*_nbadevents_veto/_nevents) << "%) \n";

   fout->Write();
   fout->Close();
}
