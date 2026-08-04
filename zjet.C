#define zjet_cxx
#include "zjet.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

#include "TLorentzVector.h"
#include "TH2D.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TStopwatch.h"

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

   TStopwatch fulltime, laptime;
   fulltime.Start();
   TDatime bgn;
   int nlap(0);
   
   fChain->SetBranchStatus("*",0);  // disable all branches

   fChain->SetBranchStatus("nMuon",1);
   fChain->SetBranchStatus("Muon_pt",1);
   fChain->SetBranchStatus("Muon_eta",1);
   fChain->SetBranchStatus("Muon_phi",1);
   fChain->SetBranchStatus("Muon_mass",1);
   fChain->SetBranchStatus("Muon_charge",1);

   fChain->SetBranchStatus("nJet",1);
   fChain->SetBranchStatus("Jet_pt",1);
   fChain->SetBranchStatus("Jet_eta",1);
   fChain->SetBranchStatus("Jet_phi",1);
   fChain->SetBranchStatus("Jet_mass",1);
   fChain->SetBranchStatus("Jet_rawFactor",1);

   fChain->SetBranchStatus("PuppiMET_pt",1);
   fChain->SetBranchStatus("PuppiMET_phi",1);
   fChain->SetBranchStatus("RawPuppiMET_pt",1);
   fChain->SetBranchStatus("RawPuppiMET_phi",1);
   
   //Long64_t nentries = fChain->GetEntriesFast();
   Long64_t nentries = fChain->GetEntries(); // Long startup time
   cout << "\nStarting loop over " << "dataset" << " with "
	<< nentries << " entries" << endl;
   if (isMC)  cout << "Running over MC branches" << endl;
   if (!isMC) cout << "Running over DATA branches" << endl;

   TDirectory *curdir = gDirectory;
   TFile *fout = new TFile("rootfiles/zjet.root","RECREATE");

   
   // Object pT plots
   fout->mkdir("control");
   fout->cd("control");
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

   TH1D *h_nsel = new TH1D("h_nsel","",20,0,20);
   TH1D *h_sel1pt = new TH1D("h_sel1pt","",200,0,200);
   TH1D *h_selpt = new TH1D("h_selpt","",200,0,200);
   TH1D *h_sel1eta = new TH1D("h_sel1eta","",100,-5,5);
   TH1D *h_seleta = new TH1D("h_seleta","",100,-5,5);
   TH1D *h_seldphi = new TH1D("h_seldphi","",120,-TMath::TwoPi(),TMath::TwoPi());

   TH2D *h2_parpteta = new TH2D("h_parpteta","",200,0,200,100,-5,5);
   TH1D *h_parpt = new TH1D("h_parpt","",200,0,200);
   TH1D *h_pareta = new TH1D("h_pareta","",100,-5,5);
   TH1D *h_pardphi = new TH1D("h_pardphi","",120,-TMath::TwoPi(),TMath::TwoPi());
   
   TH2D *h2_tranpteta = new TH2D("h_tranpteta","",200,0,200,100,-5,5);
   TH1D *h_ntran = new TH1D("h_tran","",20,0,20);
   TH1D *h_tran1pt = new TH1D("h_tran1pt","",200,0,200);
   TH1D *h_tranpt = new TH1D("h_tranpt","",200,0,200);
   TH1D *h_tran1eta = new TH1D("h_tran1eta","",100,-5,5);
   TH1D *h_traneta = new TH1D("h_traneta","",100,-5,5);
   TH1D *h_trandphi = new TH1D("h_trandphi","",120,-TMath::TwoPi(),TMath::TwoPi());

   TH2D *h2_mixpteta = new TH2D("h_mixpteta","",200,0,200,100,-5,5);
   TH1D *h_mixpt = new TH1D("h_mixpt","",200,0,200);
   TH1D *h_mixeta = new TH1D("h_mixeta","",100,-5,5);
   TH1D *h_mixdphi = new TH1D("h_mixdphi","",120,-TMath::TwoPi(),TMath::TwoPi());
   
   TH1D *h_db = new TH1D("h_db","",200,0,2);
   TH2D *h2_db = new TH2D("h2_db","",200,0,200,200,0,200);
   TProfile *p_db_vsz = new TProfile("p_db_vsz","",200,0,200);
   TProfile *p_db_vsj = new TProfile("p_db_vsj","",200,0,200);
   TProfile *p_db_vsa = new TProfile("p_db_vsa","",200,0,200);

   TH1D *h_dbp = new TH1D("h_dbp","",200,0,2);
   TH2D *h2_dbp = new TH2D("h2_dbp","",200,0,200,200,0,200);
   TProfile *p_dbp_vsz = new TProfile("p_dbp_vsz","",200,0,200);
   TProfile *p_dbp_vsj = new TProfile("p_dbp_vsj","",200,0,200);
   TProfile *p_dbp_vsa = new TProfile("p_dbp_vsa","",200,0,200);

   TH1D *h_dbt = new TH1D("h_dbt","",200,0,2);
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
   TLorentzVector p4p, p4t, p4t1, p4t2;
   TLorentzVector met, ht, met1, metn, metu, metnu;
   
   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;

      if (jentry==100000 || jentry==1000000 || jentry==1000000 ||
	  (jentry%1000000==0 && jentry<10000000) ||
	  (jentry%10000000==0 && jentry!=0) ||
	  jentry==nentries-1) {
	if (jentry==0) { laptime.Start(); }
	if (nentries!=0) {
	  cout << Form("\nProcessed %lld events (%1.1f%%) in %1.0f sec. "
		       "(%1.0f sec. for last %d)",
		       jentry, 100.*jentry/nentries, fulltime.RealTime(),
		       laptime.RealTime(), nlap);
	}
	if (jentry!=0 && nlap!=0) {
	cout << Form("\nEstimated runtime:  %1.0f sec. "
		     " (%1.0f sec. for last %d)",
		     1.*nentries/jentry*fulltime.RealTime(),
		     1.*nentries/nlap*laptime.RealTime(),nlap) << flush;
	laptime.Reset();
	nlap = 0;
	}
	if (jentry==0) fulltime.Reset(); // Leave out initialization time
	fulltime.Continue();
	laptime.Continue();
      }
      if (jentry%10000==0) cout << "." << flush;
      ++nlap;
      
      p4lplus.SetPtEtaPhiM(0,0,0,0);
      p4lminus.SetPtEtaPhiM(0,0,0,0);
      p4z.SetPtEtaPhiM(0,0,0,0);
      p4t.SetPtEtaPhiM(0,0,0,0);
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
      int nlep(0), nsel(0), ntran(0);

      // Select leading leptons
      h_nlep->Fill(nMuon);
      if (nMuon>nMuonMax) {
	cout << "ERROR: nMuon="<<nMuon<<" > nMuonMax="<<nMuonMax<<endl;
	continue;
      }
      for (int ilep = 0; ilep != nMuon; ++ilep) {

	//if (ilep==0) h_lep1pt->Fill(Muon_pt[ilep]);
	//if (ilep==1) h_lep2pt->Fill(Muon_pt[ilep]);
	//h_leppt->Fill(Muon_pt[ilep]);

	if (Muon_charge[ilep]>0 && fabs(Muon_eta[ilep])<2.5 &&
	    Muon_pt[ilep]>p4lplus.Pt())
	  p4lplus.SetPtEtaPhiM(Muon_pt[ilep], Muon_eta[ilep], Muon_phi[ilep],
			       Muon_mass[ilep]);
	if (Muon_charge[ilep]<0 && fabs(Muon_eta[ilep])<2.5 &&
	    Muon_pt[ilep]>p4lminus.Pt())
	  p4lminus.SetPtEtaPhiM(Muon_pt[ilep], Muon_eta[ilep], Muon_phi[ilep],
				Muon_mass[ilep]);
      } // for ilep

      // Reconstruct Z boson
      if (p4lplus.Pt()>0) ++nlep;
      if (p4lminus.Pt()>0) ++nlep;
      if (p4lplus.Pt()>0 && p4lminus.Pt()>0) {
	p4z += p4lplus; p4z += p4lminus;
      }
      const double mz = 91.1880;
      const double dmz = 1.5*2.4955; // 1.5*Gamma,Z~3.7; was 10-20 GeV
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

      // Set Z-parallel directions
      p4p.SetPtEtaPhiM(p4z.Pt(),p4z.Eta(),p4z.Phi()+TMath::Pi(),p4z.M());
	
      // Set Z-transverse direction(s)
      p4t.SetPtEtaPhiM(p4z.Pt(),p4z.Eta(),p4z.Phi()+(jentry%2==0 ? +1 : -1)*TMath::Pi()*0.5,p4z.M());
      p4t1.SetPtEtaPhiM(p4z.Pt(),p4z.Eta(),p4z.Phi()+TMath::Pi()*0.5,p4z.M());
      p4t2.SetPtEtaPhiM(p4z.Pt(),p4z.Eta(),p4z.Phi()-TMath::Pi()*0.5,p4z.M());

      // Keep leptons out of parallel and transverse probe directions to avoid
      // relative bias between parallel and transverse pileup jet counts
      // (to be supersafe would need also margin of deltaPhi=0.2)
      // (using DeltaPhi<1./16.*pi for probes effectively provides this)
      // (outside tracker coverage could still have paired anti-parallel PU jet
      //  inside trackser coverage overlap with lepton from Z boson?)
      //if (fabs(p4p.DeltaPhi(p4z))<1./8.*TMath::Pi() ||
      //fabs(p4t.DeltaPhi(p4z))<1./8.*TMath::Pi())
      if (fabs(p4p.DeltaPhi(p4lplus))<1./8.*TMath::Pi() ||
	  fabs(p4p.DeltaPhi(p4lminus))<1./8.*TMath::Pi() ||
	  fabs(p4t.DeltaPhi(p4lplus))<1./8.*TMath::Pi() ||
	  fabs(p4t.DeltaPhi(p4lminus))<1./8.*TMath::Pi())
	continue;
      h_zpt_probeveto->Fill(p4z.Pt());
      
      // Calculate MET and HT sum
      met.SetPtEtaPhiM(PuppiMET_pt, 0., PuppiMET_phi, 0.0);
      ht += p4z;
      for (int ijet = 0; ijet != nJet; ++ijet) {
	p4jet.SetPtEtaPhiM(Jet_pt[ijet], Jet_eta[ijet], Jet_phi[ijet],
			   Jet_mass[ijet]);
	//if (p4jet.DeltaR(p4lplus)>0.4 && p4jet.DeltaR(p4lminus)>0.4 &&
	if (p4jet.DeltaR(p4lplus)>0.2 && p4jet.DeltaR(p4lminus)>0.2 &&
	    p4jet.Pt()>15.) {
	  ht += p4jet;
	}
      }
      ht.SetPtEtaPhiM(ht.Pt(),0,ht.Phi(),0);
      metu = met + ht;
      
      // Select leading jet
      h_njet->Fill(nJet-nlep);
      if (nJet>nJetMax) {
	cout << "ERROR: nJet="<<nJet<<" > nJetMax="<<nJetMax<<endl;
	continue;
      }
      for (int ijet = 0; ijet != nJet; ++ijet) {

	p4jet.SetPtEtaPhiM(Jet_pt[ijet], Jet_eta[ijet], Jet_phi[ijet],
			   Jet_mass[ijet]);
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
	      p4jet.Pt()>0.6*p4z.Pt() && p4z.Pt()>0.6*p4jet.Pt()) {


	    //if (fabs(p4jet.Eta())<1.305) {
	    if (true) {
	      ++nsel;
	      h2_mixpteta->Fill(p4jet.Pt(),p4jet.Eta());
	      h_mixpt->Fill(p4jet.Pt());
	      h_mixeta->Fill(p4jet.Eta());
	      h_mixdphi->Fill(p4jet.DeltaPhi(p4p));
	      
	      h2_parpteta->Fill(p4jet.Pt(),p4jet.Eta());
	      h_parpt->Fill(p4jet.Pt());
	      h_pareta->Fill(p4jet.Eta());
	      h_pardphi->Fill(p4jet.DeltaPhi(p4p));
	      if (p4jet.Pt()>p4sel1.Pt()) p4sel1 = p4jet;
	      
	      h_dbp->Fill(p4jet.Pt() / p4z.Pt());
	      h2_dbp->Fill(p4z.Pt(), p4jet.Pt());
	      p_dbp_vsz->Fill(p4z.Pt(), p4jet.Pt() / p4z.Pt());
	      p_dbp_vsj->Fill(p4jet.Pt(), p4jet.Pt() / p4z.Pt());
	      p_dbp_vsa->Fill(0.5*(p4z.Pt()+p4jet.Pt()), p4jet.Pt() / p4z.Pt());
	      
	      h_selpt->Fill(p4jet.Pt(), +1);
	      h_seleta->Fill(p4jet.Eta(), +1);
	      h_seldphi->Fill(p4jet.DeltaPhi(p4p), +1);
	      
	      h_db->Fill(p4jet.Pt() / p4z.Pt(), +1);
	      h2_db->Fill(p4z.Pt(), p4jet.Pt(), +1);
	      p_db_vsz->Fill(p4z.Pt(), p4jet.Pt() / p4z.Pt(), +1);
	      p_db_vsj->Fill(p4jet.Pt(), p4jet.Pt() / p4z.Pt(), +1);
	      p_db_vsa->Fill(0.5*(p4z.Pt()+p4jet.Pt()), p4jet.Pt() / p4z.Pt());
	    } // barrel

	    double eta = p4jet.Eta();
	    double ptz = p4z.Pt();
	    double ptj = p4jet.Pt();
	    double pta = 0.5*(ptz+ptj);
	    h2ptetapf_->Fill(eta, ptj);
	    h2pteta_->Fill(eta, pta);
	    h2ptetatc_->Fill(eta, ptz);
	    h2ptetapf->Fill(eta, ptj);
	    h2pteta->Fill(eta, pta);
	    h2ptetatc->Fill(eta, ptz);

	    pmzpf_->Fill(ptj, p4z.M());
	    pmz_->Fill(pta, p4z.M());
	    pmztc_->Fill(ptz, p4z.M());
	    pmzpf->Fill(ptj, p4z.M());
	    pmz->Fill(pta, p4z.M());
	    pmztc->Fill(ptz, p4z.M());
	    
	    double jes = (1-Jet_rawFactor[ijet]);
	    p2jespf_->Fill(eta, ptj, jes);
	    p2jes_->Fill(eta, pta, jes);
	    p2jestc_->Fill(eta, ptz, jes);
	    p2jespf->Fill(eta, ptj, jes);
	    p2jes->Fill(eta, pta, jes);
	    p2jestc->Fill(eta, ptz, jes);

	    met1 = -p4z - p4jet;
	    met1.SetPtEtaPhiM(met1.Pt(),0,met1.Phi(),0.);
	    metn = -ht + p4z + p4jet;
	    metn.SetPtEtaPhiM(metn.Pt(),0,metn.Phi(),0.);
	    metnu = metn + metu;
	    double mpf = 1 + met.Vect().Dot(p4z.Vect()) / (ptz*ptz);
	    double mpf1 = 1 + met1.Vect().Dot(p4z.Vect()) / (ptz*ptz);
	    double mpfn = metn.Vect().Dot(p4z.Vect()) / (ptz*ptz);
	    double mpfu = metu.Vect().Dot(p4z.Vect()) / (ptz*ptz);
	    double mpfnu = metnu.Vect().Dot(p4z.Vect()) / (ptz*ptz);

	    p2m0pf_->Fill(eta, ptj, mpf);
	    p2m0_->Fill(eta, pta, mpf);
	    p2m0tc_->Fill(eta, ptz, mpf);
	    p2m0pf->Fill(eta, ptj, mpf);
	    p2m0->Fill(eta, pta, mpf);
	    p2m0tc->Fill(eta, ptz, mpf);

	    p2m2pf_->Fill(eta, ptj, mpf1);
	    p2m2_->Fill(eta, pta, mpf1);
	    p2m2tc_->Fill(eta, ptz, mpf1);
	    p2m2pf->Fill(eta, ptj, mpf1);
	    p2m2->Fill(eta, pta, mpf1);
	    p2m2tc->Fill(eta, ptz, mpf1);

	    p2mnpf_->Fill(eta, ptj, mpfn);
	    p2mn_->Fill(eta, pta, mpfn);
	    p2mntc_->Fill(eta, ptz, mpfn);
	    p2mnpf->Fill(eta, ptj, mpfn);
	    p2mn->Fill(eta, pta, mpfn);
	    p2mntc->Fill(eta, ptz, mpfn);

	    p2mnupf_->Fill(eta, ptj, mpfnu);
	    p2mnu_->Fill(eta, pta, mpfnu);
	    p2mnutc_->Fill(eta, ptz, mpfnu);
	    p2mnupf->Fill(eta, ptj, mpfnu);
	    p2mnu->Fill(eta, pta, mpfnu);
	    p2mnutc->Fill(eta, ptz, mpfnu);

	    p2mupf_->Fill(eta, ptj, mpfu);
	    p2mu_->Fill(eta, pta, mpfu);
	    p2mutc_->Fill(eta, ptz, mpfu);
	    p2mupf->Fill(eta, ptj, mpfu);
	    p2mu->Fill(eta, pta, mpfu);
	    p2mutc->Fill(eta, ptz, mpfu);
	  } // Parallel region

	  // Transverse region(s)
	  //if ((fabs(p4jet.DeltaPhi(p4t1))>3./4.*TMath::Pi() || // >2.36
	  //   fabs(p4jet.DeltaPhi(p4t2))>3./4.*TMath::Pi()) && // >2.36
	  //if (fabs(p4jet.DeltaPhi(p4t1))>3./4.*TMath::Pi() && // >2.36
	  //if (fabs(p4jet.DeltaPhi(p4t1))>7./8.*TMath::Pi() && // >2.75
	  //if (fabs(p4jet.DeltaPhi(p4t))>15./16.*TMath::Pi() && // >2.945
	  if (fabs(p4jet.DeltaPhi(p4t))<1./16.*TMath::Pi() && // >0.1963
	      //p4jet.Pt()>0.5*p4z.Pt() && p4z.Pt()>0.5*p4jet.Pt()) {
	      p4jet.Pt()>0.6*p4z.Pt() && p4z.Pt()>0.6*p4jet.Pt()) {
	    
	    double wt = -1;//-0.5;
	    //if (fabs(p4jet.Eta())<1.305) {
	    if (true) {

	      ++ntran;
	      h2_mixpteta->Fill(p4jet.Pt(),p4jet.Eta(),wt);
	      h_mixpt->Fill(p4jet.Pt(),wt);
	      h_mixeta->Fill(p4jet.Eta(),wt);
	      h_mixdphi->Fill(p4jet.DeltaPhi(p4t),wt);
	      
	      h2_tranpteta->Fill(p4jet.Pt(),p4jet.Eta());
	      h_tranpt->Fill(p4jet.Pt());
	      h_traneta->Fill(p4jet.Eta());
	      h_trandphi->Fill(p4jet.DeltaPhi(p4t));
	      //h_trandphi->Fill(p4jet.DeltaPhi(p4t1));
	      //h_trandphi->Fill(p4jet.DeltaPhi(p4t2));
	      if (p4jet.Pt()>p4tran1.Pt()) p4tran1 = p4jet;
	      
	      h_dbt->Fill(p4jet.Pt() / p4z.Pt());
	      h2_dbt->Fill(p4z.Pt(), p4jet.Pt());
	      p_dbt_vsz->Fill(p4z.Pt(), p4jet.Pt() / p4z.Pt());
	      p_dbt_vsj->Fill(p4jet.Pt(), p4jet.Pt() / p4z.Pt());
	      p_dbt_vsa->Fill(0.5*(p4z.Pt()+p4jet.Pt()), p4jet.Pt() / p4z.Pt());
	      
	      // Negative addition to normal region to cancel out pileup
	      h_selpt->Fill(p4jet.Pt(), wt);
	      h_seleta->Fill(p4jet.Eta(), wt);
	      h_seldphi->Fill(p4jet.DeltaPhi(p4t), wt);
	      //h_seldphi->Fill(p4jet.DeltaPhi(p4t1), 0.5*wt);
	      //h_seldphi->Fill(p4jet.DeltaPhi(p4t2), 0.5*wt);
	      
	      h_db->Fill(p4jet.Pt() / p4z.Pt(), wt);
	      h2_db->Fill(p4z.Pt(), p4jet.Pt(), wt);
	      p_db_vsz->Fill(p4z.Pt(), p4jet.Pt() / p4z.Pt(), wt);
	      p_db_vsj->Fill(p4jet.Pt(), p4jet.Pt() / p4z.Pt(), wt);
	      p_db_vsa->Fill(0.5*(p4z.Pt()+p4jet.Pt()),p4jet.Pt()/p4z.Pt(), wt);
	    } // barrel

	    double eta = p4jet.Eta();
	    double ptz = p4z.Pt();
	    double ptj = p4jet.Pt();
	    double pta = 0.5*(ptz+ptj);
	    h2ptetapf_->Fill(eta, ptj, wt);
	    h2pteta_->Fill(eta, pta, wt);
	    h2ptetatc_->Fill(eta, ptz, wt);
	    h2ptetapf->Fill(eta, ptj, wt);
	    h2pteta->Fill(eta, pta, wt);
	    h2ptetatc->Fill(eta, ptz, wt);

	    double jes = (1-Jet_rawFactor[ijet]);
	    p2jespf_->Fill(eta, ptj, jes, wt);
	    p2jes_->Fill(eta, pta, jes, wt);
	    p2jestc_->Fill(eta, ptz, jes, wt);
	    p2jespf->Fill(eta, ptj, jes, wt);
	    p2jes->Fill(eta, pta, jes, wt);
	    p2jestc->Fill(eta, ptz, jes, wt);

	    met1 = -p4z - p4jet;
	    met1.SetPtEtaPhiM(met1.Pt(),0,met1.Phi(),0);
	    metn = -ht + p4z + p4jet;
	    metn.SetPtEtaPhiM(metn.Pt(),0,metn.Phi(),0);
	    metnu = metn + metu;
	    double mpf = 1 + met.Vect().Dot(p4t.Vect()) / (ptz*ptz);
	    double mpf1 = 1 + met1.Vect().Dot(p4t.Vect()) / (ptz*ptz);
	    double mpfn = metn.Vect().Dot(p4t.Vect()) / (ptz*ptz);
	    double mpfu = metu.Vect().Dot(p4t.Vect()) / (ptz*ptz);
	    double mpfnu = metnu.Vect().Dot(p4t.Vect()) / (ptz*ptz);

	    p2m0pf_->Fill(eta, ptj, mpf, wt);
	    p2m0_->Fill(eta, pta, mpf, wt);
	    p2m0tc_->Fill(eta, ptz, mpf, wt);
	    p2m0pf->Fill(eta, ptj, mpf, wt);
	    p2m0->Fill(eta, pta, mpf, wt);
	    p2m0tc->Fill(eta, ptz, mpf, wt);

	    p2m2pf_->Fill(eta, ptj, mpf1, wt);
	    p2m2_->Fill(eta, pta, mpf1, wt);
	    p2m2tc_->Fill(eta, ptz, mpf1, wt);
	    p2m2pf->Fill(eta, ptj, mpf1, wt);
	    p2m2->Fill(eta, pta, mpf1, wt);
	    p2m2tc->Fill(eta, ptz, mpf1, wt);

	    p2mnpf_->Fill(eta, ptj, mpfn, wt);
	    p2mn_->Fill(eta, pta, mpfn, wt);
	    p2mntc_->Fill(eta, ptz, mpfn, wt);
	    p2mnpf->Fill(eta, ptj, mpfn, wt);
	    p2mn->Fill(eta, pta, mpfn, wt);
	    p2mntc->Fill(eta, ptz, mpfn, wt);

	    p2mnupf_->Fill(eta, ptj, mpfnu, wt);
	    p2mnu_->Fill(eta, pta, mpfnu, wt);
	    p2mnutc_->Fill(eta, ptz, mpfnu, wt);
	    p2mnupf->Fill(eta, ptj, mpfnu, wt);
	    p2mnu->Fill(eta, pta, mpfnu, wt);
	    p2mnutc->Fill(eta, ptz, mpfnu, wt);

	    p2mupf_->Fill(eta, ptj, mpfu, wt);
	    p2mu_->Fill(eta, pta, mpfu, wt);
	    p2mutc_->Fill(eta, ptz, mpfu, wt);
	    p2mupf->Fill(eta, ptj, mpfu, wt);
	    p2mu->Fill(eta, pta, mpfu, wt);
	    p2mutc->Fill(eta, ptz, mpfu, wt);
	  } // Transverse region
	  
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
   
   cout << endl << "Finished loop, writing file rootfiles/zjet.root." << endl << flush;
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
