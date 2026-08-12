//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Feb 11 15:15:55 2026 by ROOT version 6.34.06
// from TTree Events/Events
// found on file: ../data/zjet/DYto2Mu-4Jets_Bin-MLL-50_TuneCP5_13p6TeV_madgraphMLM-pythia8_Summer24.root
//////////////////////////////////////////////////////////

#ifndef zjet_h
#define zjet_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

#include <string>

// Header file for the classes stored in the TTree if any.

class zjet {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   UInt_t          run;
   UInt_t          luminosityBlock;
   ULong64_t       event;
   UInt_t          bunchCrossing;
   UInt_t          orbitNumber;
   Int_t           nGenJet;
  static const int nMaxGenJet = 200;
   Float_t         GenJet_eta[nMaxGenJet];   //[nGenJet]
   Float_t         GenJet_mass[nMaxGenJet];   //[nGenJet]
   Float_t         GenJet_phi[nMaxGenJet];   //[nGenJet]
   Float_t         GenJet_pt[nMaxGenJet];   //[nGenJet]
  static const int nMaxGenPart = 500;
   Int_t           nGenPart;
   Short_t         GenPart_genPartIdxMother[nMaxGenPart];   //[nGenPart]
   UShort_t        GenPart_statusFlags[nMaxGenPart];   //[nGenPart]
   Int_t           GenPart_pdgId[nMaxGenPart];   //[nGenPart]
   Int_t           GenPart_status[nMaxGenPart];   //[nGenPart]
   Float_t         GenPart_eta[nMaxGenPart];   //[nGenPart]
   Float_t         GenPart_mass[nMaxGenPart];   //[nGenPart]
   Float_t         GenPart_phi[nMaxGenPart];   //[nGenPart]
   Float_t         GenPart_pt[nMaxGenPart];   //[nGenPart]
   Float_t         GenPart_iso[nMaxGenPart];   //[nGenPart]
   Int_t           Generator_id1;
   Int_t           Generator_id2;
   Float_t         Generator_binvar;
   Float_t         Generator_scalePDF;
   Float_t         Generator_weight;
   Float_t         Generator_x1;
   Float_t         Generator_x2;
   Float_t         Generator_xpdf1;
   Float_t         Generator_xpdf2;
   Float_t         GenVtx_x;
   Float_t         GenVtx_y;
   Float_t         GenVtx_z;
   Float_t         genWeight;
   Float_t         LHEWeight_originalXWGTUP;
   Int_t           nLHEPdfWeight;
   Float_t         LHEPdfWeight[103];   //[nLHEPdfWeight]
   Int_t           nLHEReweightingWeight;
   Float_t         LHEReweightingWeight[1];   //[nLHEReweightingWeight]
   Int_t           nLHEScaleWeight;
   Float_t         LHEScaleWeight[9];   //[nLHEScaleWeight]
   Int_t           nPSWeight;
   Float_t         PSWeight[4];   //[nPSWeight]
   Int_t           nJet;
  static const int nJetMax = 200;
   UChar_t         Jet_chMultiplicity[nJetMax];   //[nJet]
   UChar_t         Jet_nConstituents[nJetMax];   //[nJet]
   UChar_t         Jet_nElectrons[nJetMax];   //[nJet]
   UChar_t         Jet_nMuons[nJetMax];   //[nJet]
   UChar_t         Jet_nSVs[nJetMax];   //[nJet]
   UChar_t         Jet_neMultiplicity[nJetMax];   //[nJet]
   Short_t         Jet_electronIdx1[nJetMax];   //[nJet]
   Short_t         Jet_electronIdx2[nJetMax];   //[nJet]
   Short_t         Jet_muonIdx1[nJetMax];   //[nJet]
   Short_t         Jet_muonIdx2[nJetMax];   //[nJet]
   Short_t         Jet_svIdx1[nJetMax];   //[nJet]
   Short_t         Jet_svIdx2[nJetMax];   //[nJet]
   Int_t           Jet_hfadjacentEtaStripsSize[nJetMax];   //[nJet]
   Int_t           Jet_hfcentralEtaStripSize[nJetMax];   //[nJet]
   Float_t         Jet_PNetRegPtRawCorr[nJetMax];   //[nJet]
   Float_t         Jet_PNetRegPtRawCorrNeutrino[nJetMax];   //[nJet]
   Float_t         Jet_PNetRegPtRawRes[nJetMax];   //[nJet]
   Float_t         Jet_UParTAK4RegPtRawCorr[nJetMax];   //[nJet]
   Float_t         Jet_UParTAK4RegPtRawCorrNeutrino[nJetMax];   //[nJet]
   Float_t         Jet_UParTAK4RegPtRawRes[nJetMax];   //[nJet]
   Float_t         Jet_UParTAK4V1RegPtRawCorr[nJetMax];   //[nJet]
   Float_t         Jet_UParTAK4V1RegPtRawCorrNeutrino[nJetMax];   //[nJet]
   Float_t         Jet_UParTAK4V1RegPtRawRes[nJetMax];   //[nJet]
   Float_t         Jet_area[nJetMax];   //[nJet]
   Float_t         Jet_btagDeepFlavB[nJetMax];   //[nJet]
   Float_t         Jet_btagDeepFlavCvB[nJetMax];   //[nJet]
   Float_t         Jet_btagDeepFlavCvL[nJetMax];   //[nJet]
   Float_t         Jet_btagDeepFlavQG[nJetMax];   //[nJet]
   Float_t         Jet_btagPNetB[nJetMax];   //[nJet]
   Float_t         Jet_btagPNetCvB[nJetMax];   //[nJet]
   Float_t         Jet_btagPNetCvL[nJetMax];   //[nJet]
   Float_t         Jet_btagPNetCvNotB[nJetMax];   //[nJet]
   Float_t         Jet_btagPNetQvG[nJetMax];   //[nJet]
   Float_t         Jet_btagPNetTauVJet[nJetMax];   //[nJet]
   Float_t         Jet_btagUParTAK4B[nJetMax];   //[nJet]
   Float_t         Jet_btagUParTAK4CvB[nJetMax];   //[nJet]
   Float_t         Jet_btagUParTAK4CvL[nJetMax];   //[nJet]
   Float_t         Jet_btagUParTAK4CvNotB[nJetMax];   //[nJet]
   Float_t         Jet_btagUParTAK4Ele[nJetMax];   //[nJet]
   Float_t         Jet_btagUParTAK4Mu[nJetMax];   //[nJet]
   Float_t         Jet_btagUParTAK4QvG[nJetMax];   //[nJet]
   Float_t         Jet_btagUParTAK4SvCB[nJetMax];   //[nJet]
   Float_t         Jet_btagUParTAK4SvUDG[nJetMax];   //[nJet]
   Float_t         Jet_btagUParTAK4TauVJet[nJetMax];   //[nJet]
   Float_t         Jet_btagUParTAK4UDG[nJetMax];   //[nJet]
   Float_t         Jet_btagUParTAK4probb[nJetMax];   //[nJet]
   Float_t         Jet_btagUParTAK4probbb[nJetMax];   //[nJet]
   Float_t         Jet_chEmEF[nJetMax];   //[nJet]
   Float_t         Jet_chHEF[nJetMax];   //[nJet]
   Float_t         Jet_eta[nJetMax];   //[nJet]
   Float_t         Jet_hfEmEF[nJetMax];   //[nJet]
   Float_t         Jet_hfHEF[nJetMax];   //[nJet]
   Float_t         Jet_hfsigmaEtaEta[nJetMax];   //[nJet]
   Float_t         Jet_hfsigmaPhiPhi[nJetMax];   //[nJet]
   Float_t         Jet_mass[nJetMax];   //[nJet]
   Float_t         Jet_muEF[nJetMax];   //[nJet]
   Float_t         Jet_muonSubtrDeltaEta[nJetMax];   //[nJet]
   Float_t         Jet_muonSubtrDeltaPhi[nJetMax];   //[nJet]
   Float_t         Jet_muonSubtrFactor[nJetMax];   //[nJet]
   Float_t         Jet_neEmEF[nJetMax];   //[nJet]
   Float_t         Jet_neHEF[nJetMax];   //[nJet]
   Float_t         Jet_phi[nJetMax];   //[nJet]
   Float_t         Jet_pt[nJetMax];   //[nJet]
   Float_t         Jet_puIdDisc[nJetMax];   //[nJet]
   Float_t         Jet_rawFactor[nJetMax];   //[nJet]
   UChar_t         LHE_Njets;
   UChar_t         LHE_Nb;
   UChar_t         LHE_Nc;
   UChar_t         LHE_Nuds;
   UChar_t         LHE_Nglu;
   UChar_t         LHE_NpNLO;
   UChar_t         LHE_NpLO;
   Float_t         LHE_HT;
   Float_t         LHE_HTIncoming;
   Float_t         LHE_Vpt;
   Float_t         LHE_AlphaS;
   Int_t           nLHEPart;
   Short_t         LHEPart_firstMotherIdx[9];   //[nLHEPart]
   Short_t         LHEPart_lastMotherIdx[9];   //[nLHEPart]
   Int_t           LHEPart_pdgId[9];   //[nLHEPart]
   Int_t           LHEPart_status[9];   //[nLHEPart]
   Int_t           LHEPart_spin[9];   //[nLHEPart]
   Float_t         LHEPart_pt[9];   //[nLHEPart]
   Float_t         LHEPart_eta[9];   //[nLHEPart]
   Float_t         LHEPart_phi[9];   //[nLHEPart]
   Float_t         LHEPart_mass[9];   //[nLHEPart]
   Float_t         LHEPart_incomingpz[9];   //[nLHEPart]
   Float_t         GenMET_phi;
   Float_t         GenMET_pt;
   Int_t           nMuon;
  static const int nMuonMax = 100;
   UChar_t         Muon_bestTrackType[nMuonMax];   //[nMuon]
   UChar_t         Muon_highPtId[nMuonMax];   //[nMuon]
   Bool_t          Muon_highPurity[nMuonMax];   //[nMuon]
   Bool_t          Muon_inTimeMuon[nMuonMax];   //[nMuon]
   Bool_t          Muon_isGlobal[nMuonMax];   //[nMuon]
   Bool_t          Muon_isPFcand[nMuonMax];   //[nMuon]
   Bool_t          Muon_isStandalone[nMuonMax];   //[nMuon]
   Bool_t          Muon_isTracker[nMuonMax];   //[nMuon]
   UChar_t         Muon_jetNDauCharged[nMuonMax];   //[nMuon]
   Bool_t          Muon_looseId[nMuonMax];   //[nMuon]
   Bool_t          Muon_mediumId[nMuonMax];   //[nMuon]
   Bool_t          Muon_mediumPromptId[nMuonMax];   //[nMuon]
   UChar_t         Muon_miniIsoId[nMuonMax];   //[nMuon]
   UChar_t         Muon_multiIsoId[nMuonMax];   //[nMuon]
   UChar_t         Muon_mvaMuID_WP[nMuonMax];   //[nMuon]
   UChar_t         Muon_nStations[nMuonMax];   //[nMuon]
   UChar_t         Muon_nTrackerLayers[nMuonMax];   //[nMuon]
   UChar_t         Muon_pfIsoId[nMuonMax];   //[nMuon]
   UChar_t         Muon_puppiIsoId[nMuonMax];   //[nMuon]
   Bool_t          Muon_softId[nMuonMax];   //[nMuon]
   Bool_t          Muon_softMvaId[nMuonMax];   //[nMuon]
   UChar_t         Muon_tightCharge[nMuonMax];   //[nMuon]
   Bool_t          Muon_tightId[nMuonMax];   //[nMuon]
   UChar_t         Muon_tkIsoId[nMuonMax];   //[nMuon]
   Bool_t          Muon_triggerIdLoose[nMuonMax];   //[nMuon]
   Short_t         Muon_jetIdx[nMuonMax];   //[nMuon]
   Short_t         Muon_svIdx[nMuonMax];   //[nMuon]
   Short_t         Muon_fsrPhotonIdx[nMuonMax];   //[nMuon]
   Int_t           Muon_charge[nMuonMax];   //[nMuon]
   Int_t           Muon_pdgId[nMuonMax];   //[nMuon]
   Float_t         Muon_VXBS_Cov00[nMuonMax];   //[nMuon]
   Float_t         Muon_VXBS_Cov03[nMuonMax];   //[nMuon]
   Float_t         Muon_VXBS_Cov33[nMuonMax];   //[nMuon]
   Float_t         Muon_dxy[nMuonMax];   //[nMuon]
   Float_t         Muon_dxyErr[nMuonMax];   //[nMuon]
   Float_t         Muon_dxybs[nMuonMax];   //[nMuon]
   Float_t         Muon_dxybsErr[nMuonMax];   //[nMuon]
   Float_t         Muon_dz[nMuonMax];   //[nMuon]
   Float_t         Muon_dzErr[nMuonMax];   //[nMuon]
   Float_t         Muon_eta[nMuonMax];   //[nMuon]
   Float_t         Muon_ip3d[nMuonMax];   //[nMuon]
   Float_t         Muon_jetDF[nMuonMax];   //[nMuon]
   Float_t         Muon_jetPtRelv2[nMuonMax];   //[nMuon]
   Float_t         Muon_jetRelIso[nMuonMax];   //[nMuon]
   Float_t         Muon_mass[nMuonMax];   //[nMuon]
   Float_t         Muon_miniPFRelIso_all[nMuonMax];   //[nMuon]
   Float_t         Muon_miniPFRelIso_chg[nMuonMax];   //[nMuon]
   Float_t         Muon_mvaMuID[nMuonMax];   //[nMuon]
   Float_t         Muon_pfRelIso03_all[nMuonMax];   //[nMuon]
   Float_t         Muon_pfRelIso03_chg[nMuonMax];   //[nMuon]
   Float_t         Muon_pfRelIso04_all[nMuonMax];   //[nMuon]
   Float_t         Muon_phi[nMuonMax];   //[nMuon]
   Float_t         Muon_pt[nMuonMax];   //[nMuon]
   Float_t         Muon_ptErr[nMuonMax];   //[nMuon]
   Float_t         Muon_segmentComp[nMuonMax];   //[nMuon]
   Float_t         Muon_sip3d[nMuonMax];   //[nMuon]
   Float_t         Muon_softMva[nMuonMax];   //[nMuon]
   Float_t         Muon_softMvaRun3[nMuonMax];   //[nMuon]
   Float_t         Muon_tkRelIso[nMuonMax];   //[nMuon]
   Float_t         Muon_tuneP_charge[nMuonMax];   //[nMuon]
   Float_t         Muon_tuneP_pterr[nMuonMax];   //[nMuon]
   Float_t         Muon_tunepRelPt[nMuonMax];   //[nMuon]
   Float_t         Muon_bsConstrainedChi2[nMuonMax];   //[nMuon]
   Float_t         Muon_bsConstrainedPt[nMuonMax];   //[nMuon]
   Float_t         Muon_bsConstrainedPtErr[nMuonMax];   //[nMuon]
   Float_t         Muon_mvaLowPt[nMuonMax];   //[nMuon]
   Float_t         Muon_pnScore_heavy[nMuonMax];   //[nMuon]
   Float_t         Muon_pnScore_light[nMuonMax];   //[nMuon]
   Float_t         Muon_pnScore_prompt[nMuonMax];   //[nMuon]
   Float_t         Muon_pnScore_tau[nMuonMax];   //[nMuon]
   Float_t         Muon_promptMVA[nMuonMax];   //[nMuon]
   Int_t           Pileup_nPU;
   Int_t           Pileup_sumEOOT;
   Int_t           Pileup_sumLOOT;
   Float_t         Pileup_nTrueInt;
   Float_t         Pileup_pudensity;
   Float_t         Pileup_gpudensity;
   Float_t         Pileup_pthatmax;
   Float_t         PuppiMET_covXX;
   Float_t         PuppiMET_covXY;
   Float_t         PuppiMET_covYY;
   Float_t         PuppiMET_phi;
   Float_t         PuppiMET_phiUnclusteredDown;
   Float_t         PuppiMET_phiUnclusteredUp;
   Float_t         PuppiMET_pt;
   Float_t         PuppiMET_ptUnclusteredDown;
   Float_t         PuppiMET_ptUnclusteredUp;
   Float_t         PuppiMET_significance;
   Float_t         PuppiMET_sumEt;
   Float_t         PuppiMET_sumPtUnclustered;
   Float_t         RawPuppiMET_phi;
   Float_t         RawPuppiMET_pt;
   Float_t         RawPuppiMET_sumEt;
   Float_t         Rho_fixedGridRhoAll;
   Float_t         Rho_fixedGridRhoFastjetAll;
   Float_t         Rho_fixedGridRhoFastjetCentral;
   Float_t         Rho_fixedGridRhoFastjetCentralCalo;
   Float_t         Rho_fixedGridRhoFastjetCentralChargedPileUp;
   Float_t         Rho_fixedGridRhoFastjetCentralNeutral;
   Float_t         FiducialMET_phi;
   Float_t         FiducialMET_pt;
  static const int nMaxTrigObj = 200;
   Int_t           nTrigObj;
   Short_t         TrigObj_l1charge[nMaxTrigObj];   //[nTrigObj]
   UShort_t        TrigObj_id[nMaxTrigObj];   //[nTrigObj]
   Int_t           TrigObj_l1iso[nMaxTrigObj];   //[nTrigObj]
   ULong64_t       TrigObj_filterBits[nMaxTrigObj];   //[nTrigObj]
   Float_t         TrigObj_pt[nMaxTrigObj];   //[nTrigObj]
   Float_t         TrigObj_eta[nMaxTrigObj];   //[nTrigObj]
   Float_t         TrigObj_phi[nMaxTrigObj];   //[nTrigObj]
   Float_t         TrigObj_l1pt[nMaxTrigObj];   //[nTrigObj]
   Float_t         TrigObj_l1pt_2[nMaxTrigObj];   //[nTrigObj]
   Float_t         TrigObj_l2pt[nMaxTrigObj];   //[nTrigObj]
   Int_t           nOtherPV;
   Float_t         OtherPV_z[3];   //[nOtherPV]
   Float_t         OtherPV_score[3];   //[nOtherPV]
   UChar_t         PV_npvs;
   UChar_t         PV_npvsGood;
   Float_t         PV_ndof;
   Float_t         PV_x;
   Float_t         PV_y;
   Float_t         PV_z;
   Float_t         PV_chi2;
   Float_t         PV_score;
   Float_t         PV_sumpt2;
   Float_t         PV_sumpx;
   Float_t         PV_sumpy;
   Int_t           nSV;
   Short_t         SV_charge[13];   //[nSV]
   Float_t         SV_dlen[13];   //[nSV]
   Float_t         SV_dlenSig[13];   //[nSV]
   Float_t         SV_dxy[13];   //[nSV]
   Float_t         SV_dxySig[13];   //[nSV]
   Float_t         SV_pAngle[13];   //[nSV]
   UChar_t         GenJet_hadronFlavour[nMaxGenJet];   //[nGenJet]
   UChar_t         GenJet_nBHadrons[nMaxGenJet];   //[nGenJet]
   UChar_t         GenJet_nCHadrons[nMaxGenJet];   //[nGenJet]
   Short_t         GenJet_partonFlavour[nMaxGenJet];   //[nGenJet]
   Float_t         GenVtx_t0;
   UChar_t         Jet_hadronFlavour[nJetMax];   //[nJet]
   Short_t         Jet_genJetIdx[nJetMax];   //[nJet]
   Short_t         Jet_partonFlavour[nJetMax];   //[nJet]
   UChar_t         Muon_genPartFlav[nMuonMax];   //[nMuon]
   Short_t         Muon_genPartIdx[nMuonMax];   //[nMuon]
   Float_t         Muon_IPx[nMuonMax];   //[nMuon]
   Float_t         Muon_IPy[nMuonMax];   //[nMuon]
   Float_t         Muon_IPz[nMuonMax];   //[nMuon]
   Float_t         Muon_ipLengthSig[nMuonMax];   //[nMuon]
   UChar_t         SV_ntracks[13];   //[nSV]
   Float_t         SV_chi2[13];   //[nSV]
   Float_t         SV_eta[13];   //[nSV]
   Float_t         SV_mass[13];   //[nSV]
   Float_t         SV_ndof[13];   //[nSV]
   Float_t         SV_phi[13];   //[nSV]
   Float_t         SV_pt[13];   //[nSV]
   Float_t         SV_x[13];   //[nSV]
   Float_t         SV_y[13];   //[nSV]
   Float_t         SV_z[13];   //[nSV]
   Bool_t          Flag_HBHENoiseFilter;
   Bool_t          Flag_HBHENoiseIsoFilter;
   Bool_t          Flag_CSCTightHaloFilter;
   Bool_t          Flag_CSCTightHaloTrkMuUnvetoFilter;
   Bool_t          Flag_CSCTightHalo2015Filter;
   Bool_t          Flag_globalTightHalo2016Filter;
   Bool_t          Flag_globalSuperTightHalo2016Filter;
   Bool_t          Flag_HcalStripHaloFilter;
   Bool_t          Flag_hcalLaserEventFilter;
   Bool_t          Flag_EcalDeadCellTriggerPrimitiveFilter;
   Bool_t          Flag_EcalDeadCellBoundaryEnergyFilter;
   Bool_t          Flag_ecalBadCalibFilter;
   Bool_t          Flag_goodVertices;
   Bool_t          Flag_eeBadScFilter;
   Bool_t          Flag_ecalLaserCorrFilter;
   Bool_t          Flag_trkPOGFilters;
   Bool_t          Flag_chargedHadronTrackResolutionFilter;
   Bool_t          Flag_muonBadTrackFilter;
   Bool_t          Flag_BadChargedCandidateFilter;
   Bool_t          Flag_BadPFMuonFilter;
   Bool_t          Flag_BadPFMuonDzFilter;
   Bool_t          Flag_hfNoisyHitsFilter;
   Bool_t          Flag_BadChargedCandidateSummer16Filter;
   Bool_t          Flag_BadPFMuonSummer16Filter;
   Bool_t          Flag_trkPOG_manystripclus53X;
   Bool_t          Flag_trkPOG_toomanystripclus53X;
   Bool_t          Flag_trkPOG_logErrorTooManyClusters;
   Bool_t          HLT_Mu37_TkMu27;
   Bool_t          HLT_IsoMu20;
   Bool_t          HLT_IsoMu24;
   Bool_t          HLT_IsoMu24_eta2p1;
   Bool_t          HLT_IsoMu27;
   Bool_t          HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL;
   Bool_t          HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL;
   Bool_t          HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ;
   Bool_t          HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ;
   Bool_t          HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8;
   Bool_t          HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_PFJet30;
   Bool_t          HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass8;
   Bool_t          HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8;
   Bool_t          HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_CaloJet30;
   Bool_t          HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_AK8PFJet30;
   Bool_t          HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_AK8CaloJet30;
   Bool_t          HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass3p8;
   Bool_t          HLT_Mu15;
   Bool_t          HLT_Mu20;
   Bool_t          HLT_Mu27;
   Bool_t          HLT_Mu50;
   Bool_t          HLT_Mu55;

   // List of branches
   TBranch        *b_run;   //!
   TBranch        *b_luminosityBlock;   //!
   TBranch        *b_event;   //!
   TBranch        *b_bunchCrossing;   //!
   TBranch        *b_orbitNumber;   //!
   TBranch        *b_nGenJet;   //!
   TBranch        *b_GenJet_eta;   //!
   TBranch        *b_GenJet_mass;   //!
   TBranch        *b_GenJet_phi;   //!
   TBranch        *b_GenJet_pt;   //!
   TBranch        *b_nGenPart;   //!
   TBranch        *b_GenPart_genPartIdxMother;   //!
   TBranch        *b_GenPart_statusFlags;   //!
   TBranch        *b_GenPart_pdgId;   //!
   TBranch        *b_GenPart_status;   //!
   TBranch        *b_GenPart_eta;   //!
   TBranch        *b_GenPart_mass;   //!
   TBranch        *b_GenPart_phi;   //!
   TBranch        *b_GenPart_pt;   //!
   TBranch        *b_GenPart_iso;   //!
   TBranch        *b_Generator_id1;   //!
   TBranch        *b_Generator_id2;   //!
   TBranch        *b_Generator_binvar;   //!
   TBranch        *b_Generator_scalePDF;   //!
   TBranch        *b_Generator_weight;   //!
   TBranch        *b_Generator_x1;   //!
   TBranch        *b_Generator_x2;   //!
   TBranch        *b_Generator_xpdf1;   //!
   TBranch        *b_Generator_xpdf2;   //!
   TBranch        *b_GenVtx_x;   //!
   TBranch        *b_GenVtx_y;   //!
   TBranch        *b_GenVtx_z;   //!
   TBranch        *b_genWeight;   //!
   TBranch        *b_LHEWeight_originalXWGTUP;   //!
   TBranch        *b_nLHEPdfWeight;   //!
   TBranch        *b_LHEPdfWeight;   //!
   TBranch        *b_nLHEReweightingWeight;   //!
   TBranch        *b_LHEReweightingWeight;   //!
   TBranch        *b_nLHEScaleWeight;   //!
   TBranch        *b_LHEScaleWeight;   //!
   TBranch        *b_nPSWeight;   //!
   TBranch        *b_PSWeight;   //!
   TBranch        *b_nJet;   //!
   TBranch        *b_Jet_chMultiplicity;   //!
   TBranch        *b_Jet_nConstituents;   //!
   TBranch        *b_Jet_nElectrons;   //!
   TBranch        *b_Jet_nMuons;   //!
   TBranch        *b_Jet_nSVs;   //!
   TBranch        *b_Jet_neMultiplicity;   //!
   TBranch        *b_Jet_electronIdx1;   //!
   TBranch        *b_Jet_electronIdx2;   //!
   TBranch        *b_Jet_muonIdx1;   //!
   TBranch        *b_Jet_muonIdx2;   //!
   TBranch        *b_Jet_svIdx1;   //!
   TBranch        *b_Jet_svIdx2;   //!
   TBranch        *b_Jet_hfadjacentEtaStripsSize;   //!
   TBranch        *b_Jet_hfcentralEtaStripSize;   //!
   TBranch        *b_Jet_PNetRegPtRawCorr;   //!
   TBranch        *b_Jet_PNetRegPtRawCorrNeutrino;   //!
   TBranch        *b_Jet_PNetRegPtRawRes;   //!
   TBranch        *b_Jet_UParTAK4RegPtRawCorr;   //!
   TBranch        *b_Jet_UParTAK4RegPtRawCorrNeutrino;   //!
   TBranch        *b_Jet_UParTAK4RegPtRawRes;   //!
   TBranch        *b_Jet_UParTAK4V1RegPtRawCorr;   //!
   TBranch        *b_Jet_UParTAK4V1RegPtRawCorrNeutrino;   //!
   TBranch        *b_Jet_UParTAK4V1RegPtRawRes;   //!
   TBranch        *b_Jet_area;   //!
   TBranch        *b_Jet_btagDeepFlavB;   //!
   TBranch        *b_Jet_btagDeepFlavCvB;   //!
   TBranch        *b_Jet_btagDeepFlavCvL;   //!
   TBranch        *b_Jet_btagDeepFlavQG;   //!
   TBranch        *b_Jet_btagPNetB;   //!
   TBranch        *b_Jet_btagPNetCvB;   //!
   TBranch        *b_Jet_btagPNetCvL;   //!
   TBranch        *b_Jet_btagPNetCvNotB;   //!
   TBranch        *b_Jet_btagPNetQvG;   //!
   TBranch        *b_Jet_btagPNetTauVJet;   //!
   TBranch        *b_Jet_btagUParTAK4B;   //!
   TBranch        *b_Jet_btagUParTAK4CvB;   //!
   TBranch        *b_Jet_btagUParTAK4CvL;   //!
   TBranch        *b_Jet_btagUParTAK4CvNotB;   //!
   TBranch        *b_Jet_btagUParTAK4Ele;   //!
   TBranch        *b_Jet_btagUParTAK4Mu;   //!
   TBranch        *b_Jet_btagUParTAK4QvG;   //!
   TBranch        *b_Jet_btagUParTAK4SvCB;   //!
   TBranch        *b_Jet_btagUParTAK4SvUDG;   //!
   TBranch        *b_Jet_btagUParTAK4TauVJet;   //!
   TBranch        *b_Jet_btagUParTAK4UDG;   //!
   TBranch        *b_Jet_btagUParTAK4probb;   //!
   TBranch        *b_Jet_btagUParTAK4probbb;   //!
   TBranch        *b_Jet_chEmEF;   //!
   TBranch        *b_Jet_chHEF;   //!
   TBranch        *b_Jet_eta;   //!
   TBranch        *b_Jet_hfEmEF;   //!
   TBranch        *b_Jet_hfHEF;   //!
   TBranch        *b_Jet_hfsigmaEtaEta;   //!
   TBranch        *b_Jet_hfsigmaPhiPhi;   //!
   TBranch        *b_Jet_mass;   //!
   TBranch        *b_Jet_muEF;   //!
   TBranch        *b_Jet_muonSubtrDeltaEta;   //!
   TBranch        *b_Jet_muonSubtrDeltaPhi;   //!
   TBranch        *b_Jet_muonSubtrFactor;   //!
   TBranch        *b_Jet_neEmEF;   //!
   TBranch        *b_Jet_neHEF;   //!
   TBranch        *b_Jet_phi;   //!
   TBranch        *b_Jet_pt;   //!
   TBranch        *b_Jet_puIdDisc;   //!
   TBranch        *b_Jet_rawFactor;   //!
   TBranch        *b_LHE_Njets;   //!
   TBranch        *b_LHE_Nb;   //!
   TBranch        *b_LHE_Nc;   //!
   TBranch        *b_LHE_Nuds;   //!
   TBranch        *b_LHE_Nglu;   //!
   TBranch        *b_LHE_NpNLO;   //!
   TBranch        *b_LHE_NpLO;   //!
   TBranch        *b_LHE_HT;   //!
   TBranch        *b_LHE_HTIncoming;   //!
   TBranch        *b_LHE_Vpt;   //!
   TBranch        *b_LHE_AlphaS;   //!
   TBranch        *b_nLHEPart;   //!
   TBranch        *b_LHEPart_firstMotherIdx;   //!
   TBranch        *b_LHEPart_lastMotherIdx;   //!
   TBranch        *b_LHEPart_pdgId;   //!
   TBranch        *b_LHEPart_status;   //!
   TBranch        *b_LHEPart_spin;   //!
   TBranch        *b_LHEPart_pt;   //!
   TBranch        *b_LHEPart_eta;   //!
   TBranch        *b_LHEPart_phi;   //!
   TBranch        *b_LHEPart_mass;   //!
   TBranch        *b_LHEPart_incomingpz;   //!
   TBranch        *b_GenMET_phi;   //!
   TBranch        *b_GenMET_pt;   //!
   TBranch        *b_nMuon;   //!
   TBranch        *b_Muon_bestTrackType;   //!
   TBranch        *b_Muon_highPtId;   //!
   TBranch        *b_Muon_highPurity;   //!
   TBranch        *b_Muon_inTimeMuon;   //!
   TBranch        *b_Muon_isGlobal;   //!
   TBranch        *b_Muon_isPFcand;   //!
   TBranch        *b_Muon_isStandalone;   //!
   TBranch        *b_Muon_isTracker;   //!
   TBranch        *b_Muon_jetNDauCharged;   //!
   TBranch        *b_Muon_looseId;   //!
   TBranch        *b_Muon_mediumId;   //!
   TBranch        *b_Muon_mediumPromptId;   //!
   TBranch        *b_Muon_miniIsoId;   //!
   TBranch        *b_Muon_multiIsoId;   //!
   TBranch        *b_Muon_mvaMuID_WP;   //!
   TBranch        *b_Muon_nStations;   //!
   TBranch        *b_Muon_nTrackerLayers;   //!
   TBranch        *b_Muon_pfIsoId;   //!
   TBranch        *b_Muon_puppiIsoId;   //!
   TBranch        *b_Muon_softId;   //!
   TBranch        *b_Muon_softMvaId;   //!
   TBranch        *b_Muon_tightCharge;   //!
   TBranch        *b_Muon_tightId;   //!
   TBranch        *b_Muon_tkIsoId;   //!
   TBranch        *b_Muon_triggerIdLoose;   //!
   TBranch        *b_Muon_jetIdx;   //!
   TBranch        *b_Muon_svIdx;   //!
   TBranch        *b_Muon_fsrPhotonIdx;   //!
   TBranch        *b_Muon_charge;   //!
   TBranch        *b_Muon_pdgId;   //!
   TBranch        *b_Muon_VXBS_Cov00;   //!
   TBranch        *b_Muon_VXBS_Cov03;   //!
   TBranch        *b_Muon_VXBS_Cov33;   //!
   TBranch        *b_Muon_dxy;   //!
   TBranch        *b_Muon_dxyErr;   //!
   TBranch        *b_Muon_dxybs;   //!
   TBranch        *b_Muon_dxybsErr;   //!
   TBranch        *b_Muon_dz;   //!
   TBranch        *b_Muon_dzErr;   //!
   TBranch        *b_Muon_eta;   //!
   TBranch        *b_Muon_ip3d;   //!
   TBranch        *b_Muon_jetDF;   //!
   TBranch        *b_Muon_jetPtRelv2;   //!
   TBranch        *b_Muon_jetRelIso;   //!
   TBranch        *b_Muon_mass;   //!
   TBranch        *b_Muon_miniPFRelIso_all;   //!
   TBranch        *b_Muon_miniPFRelIso_chg;   //!
   TBranch        *b_Muon_mvaMuID;   //!
   TBranch        *b_Muon_pfRelIso03_all;   //!
   TBranch        *b_Muon_pfRelIso03_chg;   //!
   TBranch        *b_Muon_pfRelIso04_all;   //!
   TBranch        *b_Muon_phi;   //!
   TBranch        *b_Muon_pt;   //!
   TBranch        *b_Muon_ptErr;   //!
   TBranch        *b_Muon_segmentComp;   //!
   TBranch        *b_Muon_sip3d;   //!
   TBranch        *b_Muon_softMva;   //!
   TBranch        *b_Muon_softMvaRun3;   //!
   TBranch        *b_Muon_tkRelIso;   //!
   TBranch        *b_Muon_tuneP_charge;   //!
   TBranch        *b_Muon_tuneP_pterr;   //!
   TBranch        *b_Muon_tunepRelPt;   //!
   TBranch        *b_Muon_bsConstrainedChi2;   //!
   TBranch        *b_Muon_bsConstrainedPt;   //!
   TBranch        *b_Muon_bsConstrainedPtErr;   //!
   TBranch        *b_Muon_mvaLowPt;   //!
   TBranch        *b_Muon_pnScore_heavy;   //!
   TBranch        *b_Muon_pnScore_light;   //!
   TBranch        *b_Muon_pnScore_prompt;   //!
   TBranch        *b_Muon_pnScore_tau;   //!
   TBranch        *b_Muon_promptMVA;   //!
   TBranch        *b_Pileup_nPU;   //!
   TBranch        *b_Pileup_sumEOOT;   //!
   TBranch        *b_Pileup_sumLOOT;   //!
   TBranch        *b_Pileup_nTrueInt;   //!
   TBranch        *b_Pileup_pudensity;   //!
   TBranch        *b_Pileup_gpudensity;   //!
   TBranch        *b_Pileup_pthatmax;   //!
   TBranch        *b_PuppiMET_covXX;   //!
   TBranch        *b_PuppiMET_covXY;   //!
   TBranch        *b_PuppiMET_covYY;   //!
   TBranch        *b_PuppiMET_phi;   //!
   TBranch        *b_PuppiMET_phiUnclusteredDown;   //!
   TBranch        *b_PuppiMET_phiUnclusteredUp;   //!
   TBranch        *b_PuppiMET_pt;   //!
   TBranch        *b_PuppiMET_ptUnclusteredDown;   //!
   TBranch        *b_PuppiMET_ptUnclusteredUp;   //!
   TBranch        *b_PuppiMET_significance;   //!
   TBranch        *b_PuppiMET_sumEt;   //!
   TBranch        *b_PuppiMET_sumPtUnclustered;   //!
   TBranch        *b_RawPuppiMET_phi;   //!
   TBranch        *b_RawPuppiMET_pt;   //!
   TBranch        *b_RawPuppiMET_sumEt;   //!
   TBranch        *b_Rho_fixedGridRhoAll;   //!
   TBranch        *b_Rho_fixedGridRhoFastjetAll;   //!
   TBranch        *b_Rho_fixedGridRhoFastjetCentral;   //!
   TBranch        *b_Rho_fixedGridRhoFastjetCentralCalo;   //!
   TBranch        *b_Rho_fixedGridRhoFastjetCentralChargedPileUp;   //!
   TBranch        *b_Rho_fixedGridRhoFastjetCentralNeutral;   //!
   TBranch        *b_FiducialMET_phi;   //!
   TBranch        *b_FiducialMET_pt;   //!
   TBranch        *b_nTrigObj;   //!
   TBranch        *b_TrigObj_l1charge;   //!
   TBranch        *b_TrigObj_id;   //!
   TBranch        *b_TrigObj_l1iso;   //!
   TBranch        *b_TrigObj_filterBits;   //!
   TBranch        *b_TrigObj_pt;   //!
   TBranch        *b_TrigObj_eta;   //!
   TBranch        *b_TrigObj_phi;   //!
   TBranch        *b_TrigObj_l1pt;   //!
   TBranch        *b_TrigObj_l1pt_2;   //!
   TBranch        *b_TrigObj_l2pt;   //!
   TBranch        *b_TrkMET_phi;   //!
   TBranch        *b_TrkMET_pt;   //!
   TBranch        *b_TrkMET_sumEt;   //!
   TBranch        *b_PV_npvs;   //!
   TBranch        *b_PV_npvsGood;   //!
   TBranch        *b_PV_ndof;   //!
   TBranch        *b_PV_x;   //!
   TBranch        *b_PV_y;   //!
   TBranch        *b_PV_z;   //!
   TBranch        *b_PV_chi2;   //!
   TBranch        *b_PV_score;   //!
   TBranch        *b_PV_sumpt2;   //!
   TBranch        *b_PV_sumpx;   //!
   TBranch        *b_PV_sumpy;   //!
   TBranch        *b_GenJet_hadronFlavour;   //!
   TBranch        *b_GenJet_nBHadrons;   //!
   TBranch        *b_GenJet_nCHadrons;   //!
   TBranch        *b_GenJet_partonFlavour;   //!
   TBranch        *b_GenVtx_t0;   //!
   TBranch        *b_Jet_hadronFlavour;   //!
   TBranch        *b_Jet_genJetIdx;   //!
   TBranch        *b_Jet_partonFlavour;   //!
   TBranch        *b_Muon_genPartFlav;   //!
   TBranch        *b_Muon_genPartIdx;   //!
   TBranch        *b_Muon_IPx;   //!
   TBranch        *b_Muon_IPy;   //!
   TBranch        *b_Muon_IPz;   //!
   TBranch        *b_Muon_ipLengthSig;   //!
   TBranch        *b_Flag_HBHENoiseFilter;   //!
   TBranch        *b_Flag_HBHENoiseIsoFilter;   //!
   TBranch        *b_Flag_CSCTightHaloFilter;   //!
   TBranch        *b_Flag_CSCTightHaloTrkMuUnvetoFilter;   //!
   TBranch        *b_Flag_CSCTightHalo2015Filter;   //!
   TBranch        *b_Flag_globalTightHalo2016Filter;   //!
   TBranch        *b_Flag_globalSuperTightHalo2016Filter;   //!
   TBranch        *b_Flag_HcalStripHaloFilter;   //!
   TBranch        *b_Flag_hcalLaserEventFilter;   //!
   TBranch        *b_Flag_EcalDeadCellTriggerPrimitiveFilter;   //!
   TBranch        *b_Flag_EcalDeadCellBoundaryEnergyFilter;   //!
   TBranch        *b_Flag_ecalBadCalibFilter;   //!
   TBranch        *b_Flag_goodVertices;   //!
   TBranch        *b_Flag_eeBadScFilter;   //!
   TBranch        *b_Flag_ecalLaserCorrFilter;   //!
   TBranch        *b_Flag_trkPOGFilters;   //!
   TBranch        *b_Flag_chargedHadronTrackResolutionFilter;   //!
   TBranch        *b_Flag_muonBadTrackFilter;   //!
   TBranch        *b_Flag_BadChargedCandidateFilter;   //!
   TBranch        *b_Flag_BadPFMuonFilter;   //!
   TBranch        *b_Flag_BadPFMuonDzFilter;   //!
   TBranch        *b_Flag_hfNoisyHitsFilter;   //!
   TBranch        *b_Flag_BadChargedCandidateSummer16Filter;   //!
   TBranch        *b_Flag_BadPFMuonSummer16Filter;   //!
   TBranch        *b_Flag_trkPOG_manystripclus53X;   //!
   TBranch        *b_Flag_trkPOG_toomanystripclus53X;   //!
   TBranch        *b_Flag_trkPOG_logErrorTooManyClusters;   //!
   TBranch        *b_HLT_IsoMu20;   //!
   TBranch        *b_HLT_IsoMu24;   //!
   TBranch        *b_HLT_IsoMu24_eta2p1;   //!
   TBranch        *b_HLT_IsoMu27;   //!
   TBranch        *b_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL;   //!
   TBranch        *b_HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL;   //!
   TBranch        *b_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ;   //!
   TBranch        *b_HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ;   //!
   TBranch        *b_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8;   //!
   TBranch        *b_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_PFJet30;   //!
   TBranch        *b_HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass8;   //!
   TBranch        *b_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8;   //!
   TBranch        *b_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_CaloJet30;   //!
   TBranch        *b_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_AK8PFJet30;   //!
   TBranch        *b_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_AK8CaloJet30;   //!
   TBranch        *b_HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass3p8;   //!
   TBranch        *b_HLT_Mu15;   //!
   TBranch        *b_HLT_Mu20;   //!
   TBranch        *b_HLT_Mu27;   //!
   TBranch        *b_HLT_Mu50;   //!
   TBranch        *b_HLT_Mu55;   //!
   TBranch        *b_HLT_Mu8;   //!
   TBranch        *b_HLT_Mu17;   //!
   TBranch        *b_HLT_Mu19;   //!

  bool isMC;
  std::string outputFile;
  std::string goldenJsonFile;
  std::string lumiPileupFile;
  std::string pileupWeightFile;
  std::string jecL2File;
  std::string jecResidualFile;
  std::string jerResolutionFile;
  std::string jerScaleFactorFile;
  std::string muonCorrectionFile;
  std::string jetVetoMapFile;
  zjet(TTree *tree=0, bool isMC=false,
       const char *outputFile="rootfiles/zjet.root",
       const char *goldenJsonFile="",
       const char *lumiPileupFile="",
       const char *pileupWeightFile="",
       const char *jecL2File="",
       const char *jecResidualFile="",
       const char *jerResolutionFile="",
       const char *jerScaleFactorFile="",
       const char *muonCorrectionFile="",
       const char *jetVetoMapFile="");
   virtual ~zjet();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual bool     Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef zjet_cxx
zjet::zjet(TTree *tree, bool _isMC, const char *_outputFile,
           const char *_goldenJsonFile, const char *_lumiPileupFile,
           const char *_pileupWeightFile, const char *_jecL2File,
           const char *_jecResidualFile, const char *_jerResolutionFile,
           const char *_jerScaleFactorFile, const char *_muonCorrectionFile,
           const char *_jetVetoMapFile)
  : fChain(0), isMC(_isMC), outputFile(_outputFile),
    goldenJsonFile(_goldenJsonFile), lumiPileupFile(_lumiPileupFile),
    pileupWeightFile(_pileupWeightFile), jecL2File(_jecL2File),
    jecResidualFile(_jecResidualFile), jerResolutionFile(_jerResolutionFile),
    jerScaleFactorFile(_jerScaleFactorFile),
    muonCorrectionFile(_muonCorrectionFile),
    jetVetoMapFile(_jetVetoMapFile)
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("../data/zjet/DYto2Mu-4Jets_Bin-MLL-50_TuneCP5_13p6TeV_madgraphMLM-pythia8_Summer24.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("../data/zjet/DYto2Mu-4Jets_Bin-MLL-50_TuneCP5_13p6TeV_madgraphMLM-pythia8_Summer24.root");
      }
      f->GetObject("Events",tree);

   }
   Init(tree);
}

zjet::~zjet()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t zjet::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t zjet::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void zjet::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("run", &run, &b_run);
   fChain->SetBranchAddress("luminosityBlock", &luminosityBlock, &b_luminosityBlock);
   fChain->SetBranchAddress("event", &event, &b_event);
   fChain->SetBranchAddress("bunchCrossing", &bunchCrossing, &b_bunchCrossing);
   fChain->SetBranchAddress("orbitNumber", &orbitNumber, &b_orbitNumber);
   if (isMC) {
   fChain->SetBranchAddress("nGenJet", &nGenJet, &b_nGenJet);
   fChain->SetBranchAddress("GenJet_eta", GenJet_eta, &b_GenJet_eta);
   fChain->SetBranchAddress("GenJet_mass", GenJet_mass, &b_GenJet_mass);
   fChain->SetBranchAddress("GenJet_phi", GenJet_phi, &b_GenJet_phi);
   fChain->SetBranchAddress("GenJet_pt", GenJet_pt, &b_GenJet_pt);
   fChain->SetBranchAddress("nGenPart", &nGenPart, &b_nGenPart);
   fChain->SetBranchAddress("GenPart_genPartIdxMother", GenPart_genPartIdxMother, &b_GenPart_genPartIdxMother);
   fChain->SetBranchAddress("GenPart_statusFlags", GenPart_statusFlags, &b_GenPart_statusFlags);
   fChain->SetBranchAddress("GenPart_pdgId", GenPart_pdgId, &b_GenPart_pdgId);
   fChain->SetBranchAddress("GenPart_status", GenPart_status, &b_GenPart_status);
   fChain->SetBranchAddress("GenPart_eta", GenPart_eta, &b_GenPart_eta);
   fChain->SetBranchAddress("GenPart_mass", GenPart_mass, &b_GenPart_mass);
   fChain->SetBranchAddress("GenPart_phi", GenPart_phi, &b_GenPart_phi);
   fChain->SetBranchAddress("GenPart_pt", GenPart_pt, &b_GenPart_pt);
   fChain->SetBranchAddress("GenPart_iso", GenPart_iso, &b_GenPart_iso);
   fChain->SetBranchAddress("Generator_id1", &Generator_id1, &b_Generator_id1);
   fChain->SetBranchAddress("Generator_id2", &Generator_id2, &b_Generator_id2);
   fChain->SetBranchAddress("Generator_binvar", &Generator_binvar, &b_Generator_binvar);
   fChain->SetBranchAddress("Generator_scalePDF", &Generator_scalePDF, &b_Generator_scalePDF);
   fChain->SetBranchAddress("Generator_weight", &Generator_weight, &b_Generator_weight);
   fChain->SetBranchAddress("Generator_x1", &Generator_x1, &b_Generator_x1);
   fChain->SetBranchAddress("Generator_x2", &Generator_x2, &b_Generator_x2);
   fChain->SetBranchAddress("Generator_xpdf1", &Generator_xpdf1, &b_Generator_xpdf1);
   fChain->SetBranchAddress("Generator_xpdf2", &Generator_xpdf2, &b_Generator_xpdf2);
   fChain->SetBranchAddress("GenVtx_x", &GenVtx_x, &b_GenVtx_x);
   fChain->SetBranchAddress("GenVtx_y", &GenVtx_y, &b_GenVtx_y);
   fChain->SetBranchAddress("GenVtx_z", &GenVtx_z, &b_GenVtx_z);
   fChain->SetBranchAddress("genWeight", &genWeight, &b_genWeight);
   fChain->SetBranchAddress("LHEWeight_originalXWGTUP", &LHEWeight_originalXWGTUP, &b_LHEWeight_originalXWGTUP);
   fChain->SetBranchAddress("nLHEPdfWeight", &nLHEPdfWeight, &b_nLHEPdfWeight);
   fChain->SetBranchAddress("LHEPdfWeight", LHEPdfWeight, &b_LHEPdfWeight);
   fChain->SetBranchAddress("nLHEReweightingWeight", &nLHEReweightingWeight, &b_nLHEReweightingWeight);
   fChain->SetBranchAddress("LHEReweightingWeight", &LHEReweightingWeight, &b_LHEReweightingWeight);
   fChain->SetBranchAddress("nLHEScaleWeight", &nLHEScaleWeight, &b_nLHEScaleWeight);
   fChain->SetBranchAddress("LHEScaleWeight", LHEScaleWeight, &b_LHEScaleWeight);
   fChain->SetBranchAddress("nPSWeight", &nPSWeight, &b_nPSWeight);
   fChain->SetBranchAddress("PSWeight", PSWeight, &b_PSWeight);
   } // isMC
   fChain->SetBranchAddress("nJet", &nJet, &b_nJet);
   fChain->SetBranchAddress("Jet_chMultiplicity", Jet_chMultiplicity, &b_Jet_chMultiplicity);
   fChain->SetBranchAddress("Jet_nConstituents", Jet_nConstituents, &b_Jet_nConstituents);
   fChain->SetBranchAddress("Jet_nElectrons", Jet_nElectrons, &b_Jet_nElectrons);
   fChain->SetBranchAddress("Jet_nMuons", Jet_nMuons, &b_Jet_nMuons);
   fChain->SetBranchAddress("Jet_nSVs", Jet_nSVs, &b_Jet_nSVs);
   fChain->SetBranchAddress("Jet_neMultiplicity", Jet_neMultiplicity, &b_Jet_neMultiplicity);
   fChain->SetBranchAddress("Jet_electronIdx1", Jet_electronIdx1, &b_Jet_electronIdx1);
   fChain->SetBranchAddress("Jet_electronIdx2", Jet_electronIdx2, &b_Jet_electronIdx2);
   fChain->SetBranchAddress("Jet_muonIdx1", Jet_muonIdx1, &b_Jet_muonIdx1);
   fChain->SetBranchAddress("Jet_muonIdx2", Jet_muonIdx2, &b_Jet_muonIdx2);
   fChain->SetBranchAddress("Jet_svIdx1", Jet_svIdx1, &b_Jet_svIdx1);
   fChain->SetBranchAddress("Jet_svIdx2", Jet_svIdx2, &b_Jet_svIdx2);
   fChain->SetBranchAddress("Jet_hfadjacentEtaStripsSize", Jet_hfadjacentEtaStripsSize, &b_Jet_hfadjacentEtaStripsSize);
   fChain->SetBranchAddress("Jet_hfcentralEtaStripSize", Jet_hfcentralEtaStripSize, &b_Jet_hfcentralEtaStripSize);
   fChain->SetBranchAddress("Jet_PNetRegPtRawCorr", Jet_PNetRegPtRawCorr, &b_Jet_PNetRegPtRawCorr);
   fChain->SetBranchAddress("Jet_PNetRegPtRawCorrNeutrino", Jet_PNetRegPtRawCorrNeutrino, &b_Jet_PNetRegPtRawCorrNeutrino);
   fChain->SetBranchAddress("Jet_PNetRegPtRawRes", Jet_PNetRegPtRawRes, &b_Jet_PNetRegPtRawRes);
   fChain->SetBranchAddress("Jet_UParTAK4RegPtRawCorr", Jet_UParTAK4RegPtRawCorr, &b_Jet_UParTAK4RegPtRawCorr);
   fChain->SetBranchAddress("Jet_UParTAK4RegPtRawCorrNeutrino", Jet_UParTAK4RegPtRawCorrNeutrino, &b_Jet_UParTAK4RegPtRawCorrNeutrino);
   fChain->SetBranchAddress("Jet_UParTAK4RegPtRawRes", Jet_UParTAK4RegPtRawRes, &b_Jet_UParTAK4RegPtRawRes);
   fChain->SetBranchAddress("Jet_UParTAK4V1RegPtRawCorr", Jet_UParTAK4V1RegPtRawCorr, &b_Jet_UParTAK4V1RegPtRawCorr);
   fChain->SetBranchAddress("Jet_UParTAK4V1RegPtRawCorrNeutrino", Jet_UParTAK4V1RegPtRawCorrNeutrino, &b_Jet_UParTAK4V1RegPtRawCorrNeutrino);
   fChain->SetBranchAddress("Jet_UParTAK4V1RegPtRawRes", Jet_UParTAK4V1RegPtRawRes, &b_Jet_UParTAK4V1RegPtRawRes);
   fChain->SetBranchAddress("Jet_area", Jet_area, &b_Jet_area);
   fChain->SetBranchAddress("Jet_btagDeepFlavB", Jet_btagDeepFlavB, &b_Jet_btagDeepFlavB);
   fChain->SetBranchAddress("Jet_btagDeepFlavCvB", Jet_btagDeepFlavCvB, &b_Jet_btagDeepFlavCvB);
   fChain->SetBranchAddress("Jet_btagDeepFlavCvL", Jet_btagDeepFlavCvL, &b_Jet_btagDeepFlavCvL);
   fChain->SetBranchAddress("Jet_btagDeepFlavQG", Jet_btagDeepFlavQG, &b_Jet_btagDeepFlavQG);
   fChain->SetBranchAddress("Jet_btagPNetB", Jet_btagPNetB, &b_Jet_btagPNetB);
   fChain->SetBranchAddress("Jet_btagPNetCvB", Jet_btagPNetCvB, &b_Jet_btagPNetCvB);
   fChain->SetBranchAddress("Jet_btagPNetCvL", Jet_btagPNetCvL, &b_Jet_btagPNetCvL);
   fChain->SetBranchAddress("Jet_btagPNetCvNotB", Jet_btagPNetCvNotB, &b_Jet_btagPNetCvNotB);
   fChain->SetBranchAddress("Jet_btagPNetQvG", Jet_btagPNetQvG, &b_Jet_btagPNetQvG);
   fChain->SetBranchAddress("Jet_btagPNetTauVJet", Jet_btagPNetTauVJet, &b_Jet_btagPNetTauVJet);
   fChain->SetBranchAddress("Jet_btagUParTAK4B", Jet_btagUParTAK4B, &b_Jet_btagUParTAK4B);
   fChain->SetBranchAddress("Jet_btagUParTAK4CvB", Jet_btagUParTAK4CvB, &b_Jet_btagUParTAK4CvB);
   fChain->SetBranchAddress("Jet_btagUParTAK4CvL", Jet_btagUParTAK4CvL, &b_Jet_btagUParTAK4CvL);
   fChain->SetBranchAddress("Jet_btagUParTAK4CvNotB", Jet_btagUParTAK4CvNotB, &b_Jet_btagUParTAK4CvNotB);
   fChain->SetBranchAddress("Jet_btagUParTAK4Ele", Jet_btagUParTAK4Ele, &b_Jet_btagUParTAK4Ele);
   fChain->SetBranchAddress("Jet_btagUParTAK4Mu", Jet_btagUParTAK4Mu, &b_Jet_btagUParTAK4Mu);
   fChain->SetBranchAddress("Jet_btagUParTAK4QvG", Jet_btagUParTAK4QvG, &b_Jet_btagUParTAK4QvG);
   fChain->SetBranchAddress("Jet_btagUParTAK4SvCB", Jet_btagUParTAK4SvCB, &b_Jet_btagUParTAK4SvCB);
   fChain->SetBranchAddress("Jet_btagUParTAK4SvUDG", Jet_btagUParTAK4SvUDG, &b_Jet_btagUParTAK4SvUDG);
   fChain->SetBranchAddress("Jet_btagUParTAK4TauVJet", Jet_btagUParTAK4TauVJet, &b_Jet_btagUParTAK4TauVJet);
   fChain->SetBranchAddress("Jet_btagUParTAK4UDG", Jet_btagUParTAK4UDG, &b_Jet_btagUParTAK4UDG);
   fChain->SetBranchAddress("Jet_btagUParTAK4probb", Jet_btagUParTAK4probb, &b_Jet_btagUParTAK4probb);
   fChain->SetBranchAddress("Jet_btagUParTAK4probbb", Jet_btagUParTAK4probbb, &b_Jet_btagUParTAK4probbb);
   fChain->SetBranchAddress("Jet_chEmEF", Jet_chEmEF, &b_Jet_chEmEF);
   fChain->SetBranchAddress("Jet_chHEF", Jet_chHEF, &b_Jet_chHEF);
   fChain->SetBranchAddress("Jet_eta", Jet_eta, &b_Jet_eta);
   fChain->SetBranchAddress("Jet_hfEmEF", Jet_hfEmEF, &b_Jet_hfEmEF);
   fChain->SetBranchAddress("Jet_hfHEF", Jet_hfHEF, &b_Jet_hfHEF);
   fChain->SetBranchAddress("Jet_hfsigmaEtaEta", Jet_hfsigmaEtaEta, &b_Jet_hfsigmaEtaEta);
   fChain->SetBranchAddress("Jet_hfsigmaPhiPhi", Jet_hfsigmaPhiPhi, &b_Jet_hfsigmaPhiPhi);
   fChain->SetBranchAddress("Jet_mass", Jet_mass, &b_Jet_mass);
   fChain->SetBranchAddress("Jet_muEF", Jet_muEF, &b_Jet_muEF);
   fChain->SetBranchAddress("Jet_muonSubtrDeltaEta", Jet_muonSubtrDeltaEta, &b_Jet_muonSubtrDeltaEta);
   fChain->SetBranchAddress("Jet_muonSubtrDeltaPhi", Jet_muonSubtrDeltaPhi, &b_Jet_muonSubtrDeltaPhi);
   fChain->SetBranchAddress("Jet_muonSubtrFactor", Jet_muonSubtrFactor, &b_Jet_muonSubtrFactor);
   fChain->SetBranchAddress("Jet_neEmEF", Jet_neEmEF, &b_Jet_neEmEF);
   fChain->SetBranchAddress("Jet_neHEF", Jet_neHEF, &b_Jet_neHEF);
   fChain->SetBranchAddress("Jet_phi", Jet_phi, &b_Jet_phi);
   fChain->SetBranchAddress("Jet_pt", Jet_pt, &b_Jet_pt);
   fChain->SetBranchAddress("Jet_puIdDisc", Jet_puIdDisc, &b_Jet_puIdDisc);
   fChain->SetBranchAddress("Jet_rawFactor", Jet_rawFactor, &b_Jet_rawFactor);
   if (isMC) {
   fChain->SetBranchAddress("LHE_Njets", &LHE_Njets, &b_LHE_Njets);
   fChain->SetBranchAddress("LHE_Nb", &LHE_Nb, &b_LHE_Nb);
   fChain->SetBranchAddress("LHE_Nc", &LHE_Nc, &b_LHE_Nc);
   fChain->SetBranchAddress("LHE_Nuds", &LHE_Nuds, &b_LHE_Nuds);
   fChain->SetBranchAddress("LHE_Nglu", &LHE_Nglu, &b_LHE_Nglu);
   fChain->SetBranchAddress("LHE_NpNLO", &LHE_NpNLO, &b_LHE_NpNLO);
   fChain->SetBranchAddress("LHE_NpLO", &LHE_NpLO, &b_LHE_NpLO);
   fChain->SetBranchAddress("LHE_HT", &LHE_HT, &b_LHE_HT);
   fChain->SetBranchAddress("LHE_HTIncoming", &LHE_HTIncoming, &b_LHE_HTIncoming);
   fChain->SetBranchAddress("LHE_Vpt", &LHE_Vpt, &b_LHE_Vpt);
   fChain->SetBranchAddress("LHE_AlphaS", &LHE_AlphaS, &b_LHE_AlphaS);
   fChain->SetBranchAddress("nLHEPart", &nLHEPart, &b_nLHEPart);
   fChain->SetBranchAddress("LHEPart_firstMotherIdx", LHEPart_firstMotherIdx, &b_LHEPart_firstMotherIdx);
   fChain->SetBranchAddress("LHEPart_lastMotherIdx", LHEPart_lastMotherIdx, &b_LHEPart_lastMotherIdx);
   fChain->SetBranchAddress("LHEPart_pdgId", LHEPart_pdgId, &b_LHEPart_pdgId);
   fChain->SetBranchAddress("LHEPart_status", LHEPart_status, &b_LHEPart_status);
   fChain->SetBranchAddress("LHEPart_spin", LHEPart_spin, &b_LHEPart_spin);
   fChain->SetBranchAddress("LHEPart_pt", LHEPart_pt, &b_LHEPart_pt);
   fChain->SetBranchAddress("LHEPart_eta", LHEPart_eta, &b_LHEPart_eta);
   fChain->SetBranchAddress("LHEPart_phi", LHEPart_phi, &b_LHEPart_phi);
   fChain->SetBranchAddress("LHEPart_mass", LHEPart_mass, &b_LHEPart_mass);
   fChain->SetBranchAddress("LHEPart_incomingpz", LHEPart_incomingpz, &b_LHEPart_incomingpz);
   fChain->SetBranchAddress("GenMET_phi", &GenMET_phi, &b_GenMET_phi);
   fChain->SetBranchAddress("GenMET_pt", &GenMET_pt, &b_GenMET_pt);
   } // isMC
   fChain->SetBranchAddress("nMuon", &nMuon, &b_nMuon);
   fChain->SetBranchAddress("Muon_bestTrackType", Muon_bestTrackType, &b_Muon_bestTrackType);
   fChain->SetBranchAddress("Muon_highPtId", Muon_highPtId, &b_Muon_highPtId);
   fChain->SetBranchAddress("Muon_highPurity", Muon_highPurity, &b_Muon_highPurity);
   fChain->SetBranchAddress("Muon_inTimeMuon", Muon_inTimeMuon, &b_Muon_inTimeMuon);
   fChain->SetBranchAddress("Muon_isGlobal", Muon_isGlobal, &b_Muon_isGlobal);
   fChain->SetBranchAddress("Muon_isPFcand", Muon_isPFcand, &b_Muon_isPFcand);
   fChain->SetBranchAddress("Muon_isStandalone", Muon_isStandalone, &b_Muon_isStandalone);
   fChain->SetBranchAddress("Muon_isTracker", Muon_isTracker, &b_Muon_isTracker);
   fChain->SetBranchAddress("Muon_jetNDauCharged", Muon_jetNDauCharged, &b_Muon_jetNDauCharged);
   fChain->SetBranchAddress("Muon_looseId", Muon_looseId, &b_Muon_looseId);
   fChain->SetBranchAddress("Muon_mediumId", Muon_mediumId, &b_Muon_mediumId);
   fChain->SetBranchAddress("Muon_mediumPromptId", Muon_mediumPromptId, &b_Muon_mediumPromptId);
   fChain->SetBranchAddress("Muon_miniIsoId", Muon_miniIsoId, &b_Muon_miniIsoId);
   fChain->SetBranchAddress("Muon_multiIsoId", Muon_multiIsoId, &b_Muon_multiIsoId);
   fChain->SetBranchAddress("Muon_mvaMuID_WP", Muon_mvaMuID_WP, &b_Muon_mvaMuID_WP);
   fChain->SetBranchAddress("Muon_nStations", Muon_nStations, &b_Muon_nStations);
   fChain->SetBranchAddress("Muon_nTrackerLayers", Muon_nTrackerLayers, &b_Muon_nTrackerLayers);
   fChain->SetBranchAddress("Muon_pfIsoId", Muon_pfIsoId, &b_Muon_pfIsoId);
   fChain->SetBranchAddress("Muon_puppiIsoId", Muon_puppiIsoId, &b_Muon_puppiIsoId);
   fChain->SetBranchAddress("Muon_softId", Muon_softId, &b_Muon_softId);
   fChain->SetBranchAddress("Muon_softMvaId", Muon_softMvaId, &b_Muon_softMvaId);
   fChain->SetBranchAddress("Muon_tightCharge", Muon_tightCharge, &b_Muon_tightCharge);
   fChain->SetBranchAddress("Muon_tightId", Muon_tightId, &b_Muon_tightId);
   fChain->SetBranchAddress("Muon_tkIsoId", Muon_tkIsoId, &b_Muon_tkIsoId);
   fChain->SetBranchAddress("Muon_triggerIdLoose", Muon_triggerIdLoose, &b_Muon_triggerIdLoose);
   fChain->SetBranchAddress("Muon_jetIdx", Muon_jetIdx, &b_Muon_jetIdx);
   fChain->SetBranchAddress("Muon_svIdx", Muon_svIdx, &b_Muon_svIdx);
   fChain->SetBranchAddress("Muon_fsrPhotonIdx", Muon_fsrPhotonIdx, &b_Muon_fsrPhotonIdx);
   fChain->SetBranchAddress("Muon_charge", Muon_charge, &b_Muon_charge);
   fChain->SetBranchAddress("Muon_pdgId", Muon_pdgId, &b_Muon_pdgId);
   fChain->SetBranchAddress("Muon_VXBS_Cov00", Muon_VXBS_Cov00, &b_Muon_VXBS_Cov00);
   fChain->SetBranchAddress("Muon_VXBS_Cov03", Muon_VXBS_Cov03, &b_Muon_VXBS_Cov03);
   fChain->SetBranchAddress("Muon_VXBS_Cov33", Muon_VXBS_Cov33, &b_Muon_VXBS_Cov33);
   fChain->SetBranchAddress("Muon_dxy", Muon_dxy, &b_Muon_dxy);
   fChain->SetBranchAddress("Muon_dxyErr", Muon_dxyErr, &b_Muon_dxyErr);
   fChain->SetBranchAddress("Muon_dxybs", Muon_dxybs, &b_Muon_dxybs);
   fChain->SetBranchAddress("Muon_dxybsErr", Muon_dxybsErr, &b_Muon_dxybsErr);
   fChain->SetBranchAddress("Muon_dz", Muon_dz, &b_Muon_dz);
   fChain->SetBranchAddress("Muon_dzErr", Muon_dzErr, &b_Muon_dzErr);
   fChain->SetBranchAddress("Muon_eta", Muon_eta, &b_Muon_eta);
   fChain->SetBranchAddress("Muon_ip3d", Muon_ip3d, &b_Muon_ip3d);
   fChain->SetBranchAddress("Muon_jetDF", Muon_jetDF, &b_Muon_jetDF);
   fChain->SetBranchAddress("Muon_jetPtRelv2", Muon_jetPtRelv2, &b_Muon_jetPtRelv2);
   fChain->SetBranchAddress("Muon_jetRelIso", Muon_jetRelIso, &b_Muon_jetRelIso);
   fChain->SetBranchAddress("Muon_mass", Muon_mass, &b_Muon_mass);
   fChain->SetBranchAddress("Muon_miniPFRelIso_all", Muon_miniPFRelIso_all, &b_Muon_miniPFRelIso_all);
   fChain->SetBranchAddress("Muon_miniPFRelIso_chg", Muon_miniPFRelIso_chg, &b_Muon_miniPFRelIso_chg);
   fChain->SetBranchAddress("Muon_mvaMuID", Muon_mvaMuID, &b_Muon_mvaMuID);
   fChain->SetBranchAddress("Muon_pfRelIso03_all", Muon_pfRelIso03_all, &b_Muon_pfRelIso03_all);
   fChain->SetBranchAddress("Muon_pfRelIso03_chg", Muon_pfRelIso03_chg, &b_Muon_pfRelIso03_chg);
   fChain->SetBranchAddress("Muon_pfRelIso04_all", Muon_pfRelIso04_all, &b_Muon_pfRelIso04_all);
   fChain->SetBranchAddress("Muon_phi", Muon_phi, &b_Muon_phi);
   fChain->SetBranchAddress("Muon_pt", Muon_pt, &b_Muon_pt);
   fChain->SetBranchAddress("Muon_ptErr", Muon_ptErr, &b_Muon_ptErr);
   fChain->SetBranchAddress("Muon_segmentComp", Muon_segmentComp, &b_Muon_segmentComp);
   fChain->SetBranchAddress("Muon_sip3d", Muon_sip3d, &b_Muon_sip3d);
   fChain->SetBranchAddress("Muon_softMva", Muon_softMva, &b_Muon_softMva);
   fChain->SetBranchAddress("Muon_softMvaRun3", Muon_softMvaRun3, &b_Muon_softMvaRun3);
   fChain->SetBranchAddress("Muon_tkRelIso", Muon_tkRelIso, &b_Muon_tkRelIso);
   fChain->SetBranchAddress("Muon_tuneP_charge", Muon_tuneP_charge, &b_Muon_tuneP_charge);
   fChain->SetBranchAddress("Muon_tuneP_pterr", Muon_tuneP_pterr, &b_Muon_tuneP_pterr);
   fChain->SetBranchAddress("Muon_tunepRelPt", Muon_tunepRelPt, &b_Muon_tunepRelPt);
   fChain->SetBranchAddress("Muon_bsConstrainedChi2", Muon_bsConstrainedChi2, &b_Muon_bsConstrainedChi2);
   fChain->SetBranchAddress("Muon_bsConstrainedPt", Muon_bsConstrainedPt, &b_Muon_bsConstrainedPt);
   fChain->SetBranchAddress("Muon_bsConstrainedPtErr", Muon_bsConstrainedPtErr, &b_Muon_bsConstrainedPtErr);
   fChain->SetBranchAddress("Muon_mvaLowPt", Muon_mvaLowPt, &b_Muon_mvaLowPt);
   fChain->SetBranchAddress("Muon_pnScore_heavy", Muon_pnScore_heavy, &b_Muon_pnScore_heavy);
   fChain->SetBranchAddress("Muon_pnScore_light", Muon_pnScore_light, &b_Muon_pnScore_light);
   fChain->SetBranchAddress("Muon_pnScore_prompt", Muon_pnScore_prompt, &b_Muon_pnScore_prompt);
   fChain->SetBranchAddress("Muon_pnScore_tau", Muon_pnScore_tau, &b_Muon_pnScore_tau);
   fChain->SetBranchAddress("Muon_promptMVA", Muon_promptMVA, &b_Muon_promptMVA);
   if (isMC) {
   fChain->SetBranchAddress("Pileup_nPU", &Pileup_nPU, &b_Pileup_nPU);
   fChain->SetBranchAddress("Pileup_sumEOOT", &Pileup_sumEOOT, &b_Pileup_sumEOOT);
   fChain->SetBranchAddress("Pileup_sumLOOT", &Pileup_sumLOOT, &b_Pileup_sumLOOT);
   fChain->SetBranchAddress("Pileup_nTrueInt", &Pileup_nTrueInt, &b_Pileup_nTrueInt);
   fChain->SetBranchAddress("Pileup_pudensity", &Pileup_pudensity, &b_Pileup_pudensity);
   fChain->SetBranchAddress("Pileup_gpudensity", &Pileup_gpudensity, &b_Pileup_gpudensity);
   fChain->SetBranchAddress("Pileup_pthatmax", &Pileup_pthatmax, &b_Pileup_pthatmax);
   } // isMC
   fChain->SetBranchAddress("PuppiMET_covXX", &PuppiMET_covXX, &b_PuppiMET_covXX);
   fChain->SetBranchAddress("PuppiMET_covXY", &PuppiMET_covXY, &b_PuppiMET_covXY);
   fChain->SetBranchAddress("PuppiMET_covYY", &PuppiMET_covYY, &b_PuppiMET_covYY);
   fChain->SetBranchAddress("PuppiMET_phi", &PuppiMET_phi, &b_PuppiMET_phi);
   fChain->SetBranchAddress("PuppiMET_phiUnclusteredDown", &PuppiMET_phiUnclusteredDown, &b_PuppiMET_phiUnclusteredDown);
   fChain->SetBranchAddress("PuppiMET_phiUnclusteredUp", &PuppiMET_phiUnclusteredUp, &b_PuppiMET_phiUnclusteredUp);
   fChain->SetBranchAddress("PuppiMET_pt", &PuppiMET_pt, &b_PuppiMET_pt);
   fChain->SetBranchAddress("PuppiMET_ptUnclusteredDown", &PuppiMET_ptUnclusteredDown, &b_PuppiMET_ptUnclusteredDown);
   fChain->SetBranchAddress("PuppiMET_ptUnclusteredUp", &PuppiMET_ptUnclusteredUp, &b_PuppiMET_ptUnclusteredUp);
   fChain->SetBranchAddress("PuppiMET_significance", &PuppiMET_significance, &b_PuppiMET_significance);
   fChain->SetBranchAddress("PuppiMET_sumEt", &PuppiMET_sumEt, &b_PuppiMET_sumEt);
   fChain->SetBranchAddress("PuppiMET_sumPtUnclustered", &PuppiMET_sumPtUnclustered, &b_PuppiMET_sumPtUnclustered);
   fChain->SetBranchAddress("RawPuppiMET_phi", &RawPuppiMET_phi, &b_RawPuppiMET_phi);
   fChain->SetBranchAddress("RawPuppiMET_pt", &RawPuppiMET_pt, &b_RawPuppiMET_pt);
   fChain->SetBranchAddress("RawPuppiMET_sumEt", &RawPuppiMET_sumEt, &b_RawPuppiMET_sumEt);
   fChain->SetBranchAddress("Rho_fixedGridRhoAll", &Rho_fixedGridRhoAll, &b_Rho_fixedGridRhoAll);
   fChain->SetBranchAddress("Rho_fixedGridRhoFastjetAll", &Rho_fixedGridRhoFastjetAll, &b_Rho_fixedGridRhoFastjetAll);
   fChain->SetBranchAddress("Rho_fixedGridRhoFastjetCentral", &Rho_fixedGridRhoFastjetCentral, &b_Rho_fixedGridRhoFastjetCentral);
   fChain->SetBranchAddress("Rho_fixedGridRhoFastjetCentralCalo", &Rho_fixedGridRhoFastjetCentralCalo, &b_Rho_fixedGridRhoFastjetCentralCalo);
   fChain->SetBranchAddress("Rho_fixedGridRhoFastjetCentralChargedPileUp", &Rho_fixedGridRhoFastjetCentralChargedPileUp, &b_Rho_fixedGridRhoFastjetCentralChargedPileUp);
   fChain->SetBranchAddress("Rho_fixedGridRhoFastjetCentralNeutral", &Rho_fixedGridRhoFastjetCentralNeutral, &b_Rho_fixedGridRhoFastjetCentralNeutral);
   if (isMC) {
   fChain->SetBranchAddress("FiducialMET_phi", &FiducialMET_phi, &b_FiducialMET_phi);
   fChain->SetBranchAddress("FiducialMET_pt", &FiducialMET_pt, &b_FiducialMET_pt);
   } // isMC
   fChain->SetBranchAddress("nTrigObj", &nTrigObj, &b_nTrigObj);
   fChain->SetBranchAddress("TrigObj_l1charge", TrigObj_l1charge, &b_TrigObj_l1charge);
   fChain->SetBranchAddress("TrigObj_id", TrigObj_id, &b_TrigObj_id);
   fChain->SetBranchAddress("TrigObj_l1iso", TrigObj_l1iso, &b_TrigObj_l1iso);
   fChain->SetBranchAddress("TrigObj_filterBits", TrigObj_filterBits, &b_TrigObj_filterBits);
   fChain->SetBranchAddress("TrigObj_pt", TrigObj_pt, &b_TrigObj_pt);
   fChain->SetBranchAddress("TrigObj_eta", TrigObj_eta, &b_TrigObj_eta);
   fChain->SetBranchAddress("TrigObj_phi", TrigObj_phi, &b_TrigObj_phi);
   fChain->SetBranchAddress("TrigObj_l1pt", TrigObj_l1pt, &b_TrigObj_l1pt);
   fChain->SetBranchAddress("TrigObj_l1pt_2", TrigObj_l1pt_2, &b_TrigObj_l1pt_2);
   fChain->SetBranchAddress("TrigObj_l2pt", TrigObj_l2pt, &b_TrigObj_l2pt);
   fChain->SetBranchAddress("PV_npvs", &PV_npvs, &b_PV_npvs);
   fChain->SetBranchAddress("PV_npvsGood", &PV_npvsGood, &b_PV_npvsGood);
   fChain->SetBranchAddress("PV_ndof", &PV_ndof, &b_PV_ndof);
   fChain->SetBranchAddress("PV_x", &PV_x, &b_PV_x);
   fChain->SetBranchAddress("PV_y", &PV_y, &b_PV_y);
   fChain->SetBranchAddress("PV_z", &PV_z, &b_PV_z);
   fChain->SetBranchAddress("PV_chi2", &PV_chi2, &b_PV_chi2);
   fChain->SetBranchAddress("PV_score", &PV_score, &b_PV_score);
   fChain->SetBranchAddress("PV_sumpt2", &PV_sumpt2, &b_PV_sumpt2);
   fChain->SetBranchAddress("PV_sumpx", &PV_sumpx, &b_PV_sumpx);
   fChain->SetBranchAddress("PV_sumpy", &PV_sumpy, &b_PV_sumpy);
   if (isMC) {
   fChain->SetBranchAddress("GenJet_hadronFlavour", GenJet_hadronFlavour, &b_GenJet_hadronFlavour);
   fChain->SetBranchAddress("GenJet_nBHadrons", GenJet_nBHadrons, &b_GenJet_nBHadrons);
   fChain->SetBranchAddress("GenJet_nCHadrons", GenJet_nCHadrons, &b_GenJet_nCHadrons);
   fChain->SetBranchAddress("GenJet_partonFlavour", GenJet_partonFlavour, &b_GenJet_partonFlavour);
   fChain->SetBranchAddress("GenVtx_t0", &GenVtx_t0, &b_GenVtx_t0);
   fChain->SetBranchAddress("Jet_hadronFlavour", Jet_hadronFlavour, &b_Jet_hadronFlavour);
   fChain->SetBranchAddress("Jet_genJetIdx", Jet_genJetIdx, &b_Jet_genJetIdx);
   fChain->SetBranchAddress("Jet_partonFlavour", Jet_partonFlavour, &b_Jet_partonFlavour);
   fChain->SetBranchAddress("Muon_genPartFlav", Muon_genPartFlav, &b_Muon_genPartFlav);
   fChain->SetBranchAddress("Muon_genPartIdx", Muon_genPartIdx, &b_Muon_genPartIdx);
   } // isMC
   fChain->SetBranchAddress("Muon_IPx", Muon_IPx, &b_Muon_IPx);
   fChain->SetBranchAddress("Muon_IPy", Muon_IPy, &b_Muon_IPy);
   fChain->SetBranchAddress("Muon_IPz", Muon_IPz, &b_Muon_IPz);
   fChain->SetBranchAddress("Muon_ipLengthSig", Muon_ipLengthSig, &b_Muon_ipLengthSig);
   fChain->SetBranchAddress("Flag_HBHENoiseFilter", &Flag_HBHENoiseFilter, &b_Flag_HBHENoiseFilter);
   fChain->SetBranchAddress("Flag_HBHENoiseIsoFilter", &Flag_HBHENoiseIsoFilter, &b_Flag_HBHENoiseIsoFilter);
   fChain->SetBranchAddress("Flag_CSCTightHaloFilter", &Flag_CSCTightHaloFilter, &b_Flag_CSCTightHaloFilter);
   fChain->SetBranchAddress("Flag_CSCTightHaloTrkMuUnvetoFilter", &Flag_CSCTightHaloTrkMuUnvetoFilter, &b_Flag_CSCTightHaloTrkMuUnvetoFilter);
   fChain->SetBranchAddress("Flag_CSCTightHalo2015Filter", &Flag_CSCTightHalo2015Filter, &b_Flag_CSCTightHalo2015Filter);
   fChain->SetBranchAddress("Flag_globalTightHalo2016Filter", &Flag_globalTightHalo2016Filter, &b_Flag_globalTightHalo2016Filter);
   fChain->SetBranchAddress("Flag_globalSuperTightHalo2016Filter", &Flag_globalSuperTightHalo2016Filter, &b_Flag_globalSuperTightHalo2016Filter);
   fChain->SetBranchAddress("Flag_HcalStripHaloFilter", &Flag_HcalStripHaloFilter, &b_Flag_HcalStripHaloFilter);
   fChain->SetBranchAddress("Flag_hcalLaserEventFilter", &Flag_hcalLaserEventFilter, &b_Flag_hcalLaserEventFilter);
   fChain->SetBranchAddress("Flag_EcalDeadCellTriggerPrimitiveFilter", &Flag_EcalDeadCellTriggerPrimitiveFilter, &b_Flag_EcalDeadCellTriggerPrimitiveFilter);
   fChain->SetBranchAddress("Flag_EcalDeadCellBoundaryEnergyFilter", &Flag_EcalDeadCellBoundaryEnergyFilter, &b_Flag_EcalDeadCellBoundaryEnergyFilter);
   fChain->SetBranchAddress("Flag_ecalBadCalibFilter", &Flag_ecalBadCalibFilter, &b_Flag_ecalBadCalibFilter);
   fChain->SetBranchAddress("Flag_goodVertices", &Flag_goodVertices, &b_Flag_goodVertices);
   fChain->SetBranchAddress("Flag_eeBadScFilter", &Flag_eeBadScFilter, &b_Flag_eeBadScFilter);
   fChain->SetBranchAddress("Flag_ecalLaserCorrFilter", &Flag_ecalLaserCorrFilter, &b_Flag_ecalLaserCorrFilter);
   fChain->SetBranchAddress("Flag_trkPOGFilters", &Flag_trkPOGFilters, &b_Flag_trkPOGFilters);
   fChain->SetBranchAddress("Flag_chargedHadronTrackResolutionFilter", &Flag_chargedHadronTrackResolutionFilter, &b_Flag_chargedHadronTrackResolutionFilter);
   fChain->SetBranchAddress("Flag_muonBadTrackFilter", &Flag_muonBadTrackFilter, &b_Flag_muonBadTrackFilter);
   fChain->SetBranchAddress("Flag_BadChargedCandidateFilter", &Flag_BadChargedCandidateFilter, &b_Flag_BadChargedCandidateFilter);
   fChain->SetBranchAddress("Flag_BadPFMuonFilter", &Flag_BadPFMuonFilter, &b_Flag_BadPFMuonFilter);
   fChain->SetBranchAddress("Flag_BadPFMuonDzFilter", &Flag_BadPFMuonDzFilter, &b_Flag_BadPFMuonDzFilter);
   fChain->SetBranchAddress("Flag_hfNoisyHitsFilter", &Flag_hfNoisyHitsFilter, &b_Flag_hfNoisyHitsFilter);
   fChain->SetBranchAddress("Flag_BadChargedCandidateSummer16Filter", &Flag_BadChargedCandidateSummer16Filter, &b_Flag_BadChargedCandidateSummer16Filter);
   fChain->SetBranchAddress("Flag_BadPFMuonSummer16Filter", &Flag_BadPFMuonSummer16Filter, &b_Flag_BadPFMuonSummer16Filter);
   fChain->SetBranchAddress("Flag_trkPOG_manystripclus53X", &Flag_trkPOG_manystripclus53X, &b_Flag_trkPOG_manystripclus53X);
   fChain->SetBranchAddress("Flag_trkPOG_toomanystripclus53X", &Flag_trkPOG_toomanystripclus53X, &b_Flag_trkPOG_toomanystripclus53X);
   fChain->SetBranchAddress("Flag_trkPOG_logErrorTooManyClusters", &Flag_trkPOG_logErrorTooManyClusters, &b_Flag_trkPOG_logErrorTooManyClusters);
   fChain->SetBranchAddress("HLT_IsoMu20", &HLT_IsoMu20, &b_HLT_IsoMu20);
   fChain->SetBranchAddress("HLT_IsoMu24", &HLT_IsoMu24, &b_HLT_IsoMu24);
   fChain->SetBranchAddress("HLT_IsoMu24_eta2p1", &HLT_IsoMu24_eta2p1, &b_HLT_IsoMu24_eta2p1);
   fChain->SetBranchAddress("HLT_IsoMu27", &HLT_IsoMu27, &b_HLT_IsoMu27);
   fChain->SetBranchAddress("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL", &HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL, &b_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL);
   fChain->SetBranchAddress("HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL", &HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL, &b_HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL);
   fChain->SetBranchAddress("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ", &HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ, &b_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ);
   fChain->SetBranchAddress("HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ", &HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ, &b_HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ);
   fChain->SetBranchAddress("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8", &HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8, &b_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8);
   fChain->SetBranchAddress("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_PFJet30", &HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_PFJet30, &b_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_PFJet30);
   fChain->SetBranchAddress("HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass8", &HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass8, &b_HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass8);
   fChain->SetBranchAddress("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8", &HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8, &b_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8);
   fChain->SetBranchAddress("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_CaloJet30", &HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_CaloJet30, &b_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_CaloJet30);
   fChain->SetBranchAddress("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_AK8PFJet30", &HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_AK8PFJet30, &b_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_AK8PFJet30);
   fChain->SetBranchAddress("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_AK8CaloJet30", &HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_AK8CaloJet30, &b_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_AK8CaloJet30);
   fChain->SetBranchAddress("HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass3p8", &HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass3p8, &b_HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass3p8);
   fChain->SetBranchAddress("HLT_Mu15", &HLT_Mu15, &b_HLT_Mu15);
   fChain->SetBranchAddress("HLT_Mu20", &HLT_Mu20, &b_HLT_Mu20);
   fChain->SetBranchAddress("HLT_Mu27", &HLT_Mu27, &b_HLT_Mu27);
   fChain->SetBranchAddress("HLT_Mu50", &HLT_Mu50, &b_HLT_Mu50);
   fChain->SetBranchAddress("HLT_Mu55", &HLT_Mu55, &b_HLT_Mu55);
   Notify();
}

bool zjet::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return true;
}

void zjet::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t zjet::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef zjet_cxx
