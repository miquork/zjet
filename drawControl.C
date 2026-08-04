// Purpose: draw control plots of the metho
#include "TH2D.h"
#include "TProfile2D.h"
#include "TF1.h"
#include "TFile.h"

#include "tdrstyle_mod22.C"

void drawControl() {

  gROOT->ProcessLine(".! mkdir pdf");
  gROOT->ProcessLine(".! touch pdf");
  gROOT->ProcessLine(".! mkdir pdf/drawControl");
  gROOT->ProcessLine(".! touch pdf/drawControl");
  
  setTDRStyle();

  TFile *f = new TFile("rootfiles/zjet_MC.root","READ");
  assert(f && !f->IsZombie());

  TH1D *h_parpt = (TH1D*)f->Get("control/h_parpt");   assert(h_parpt);
  TH1D *h_tranpt = (TH1D*)f->Get("control/h_tranpt"); assert(h_tranpt);
  TH1D *h_mixpt = (TH1D*)f->Get("control/h_mixpt");   assert(h_mixpt);

  TH1D *h1_pt = tdrHist("h1_pt","N_{events} / GeV",1,5e3,"p_{T,Z} (GeV)",0,200);
  lumi_136TeV = "JMENano Summer24 DY";
  TCanvas *c1_pt = tdrCanvas("c1_pt",h1_pt,8,11,kSquare);
  gPad->SetLogy();

  tdrDraw(h_parpt,"HIST",kNone,kGray+1,kSolid,-1,1001,kGray);
  tdrDraw(h_tranpt,"HIST",kNone,kRed-9+1,kSolid,-1,1001,kRed-9);
  tdrDraw(h_mixpt,"Pz",kFullCircle,kGreen+2,kSolid,-1,kNone,0,0.5);

  TF1 *f1_pt = new TF1("f1_pt","[0]*exp([1]*x)",10,200);
  f1_pt->SetParameters(610,-0.0276);
  h_mixpt->Fit(f1_pt,"RN");
  f1_pt->SetRange(0,200);
  f1_pt->SetLineColor(kGreen+2);
  f1_pt->Draw("SAME");
  
  TLegend *leg_pt = tdrLeg(0.52,0.87-4*0.05,0.82,0.87);
  leg_pt->SetHeader("|#Delta#phi|<#pi/16, |#Deltap_{T}/p_{T}|<0.4");
  leg_pt->AddEntry(h_parpt,"Parallel");
  leg_pt->AddEntry(h_tranpt,"Transverse");
  leg_pt->AddEntry(h_mixpt,"Par - tran");

  TLatex *tex = new TLatex();
  tex->SetNDC();
  tex->SetTextSize(0.045);
  tex->DrawLatex(0.52,0.62,"Veto leptons from");
  tex->DrawLatex(0.52,0.57,"par+tran |#Delta#phi|<#pi/8");
  
  gPad->RedrawAxis();

  c1_pt->SaveAs("pdf/drawControl/drawControl_c1_1_pt.pdf");
  

  TH1D *h_pareta = (TH1D*)f->Get("control/h_pareta");   assert(h_pareta);
  TH1D *h_traneta = (TH1D*)f->Get("control/h_traneta"); assert(h_traneta);
  TH1D *h_mixeta = (TH1D*)f->Get("control/h_mixeta");   assert(h_mixeta);

  TH1D *h1_eta = tdrHist("h1_eta","N_{events} / bin",0,1400,"#eta_{jet}",-5,5);
  lumi_136TeV = "JMENano Summer24 DY";
  TCanvas *c1_eta = tdrCanvas("c1_eta",h1_eta,8,11,kSquare);

  tdrDraw(h_pareta,"HIST",kNone,kGray+1,kSolid,-1,1001,kGray);
  tdrDraw(h_traneta,"HIST",kNone,kRed-9+1,kSolid,-1,1001,kRed-9);
  tdrDraw(h_mixeta,"Pz",kFullCircle,kGreen+2,kSolid,-1,kNone,0,0.5);

  TF1 *f1_eta = new TF1("f1_eta","[0]*exp([1]*cosh(x))",-5,5);
  f1_eta->SetParameters(610,-0.0276);
  h_mixeta->Fit(f1_eta,"RN");
  f1_eta->SetRange(-5,5);
  f1_eta->SetLineColor(kGreen+2);
  f1_eta->Draw("SAME");
  
  gPad->RedrawAxis();
  
  TLegend *leg_eta = tdrLeg(0.35,0.80-4*0.05,0.65,0.80);
  leg_eta->SetHeader("|#Delta#phi|<#pi/16, |#Deltap_{T}/p_{T}|<0.4");
  leg_eta->AddEntry(h_pareta,"Parallel");
  leg_eta->AddEntry(h_traneta,"Transverse");
  leg_eta->AddEntry(h_mixeta,"Par-tran");

  tex->DrawLatex(0.35,0.55,"Veto leptons from");
  tex->DrawLatex(0.35,0.50,"par+tran |#Delta#phi|<#pi/8");

  c1_eta->SaveAs("pdf/drawControl/drawControl_c1_2_eta.pdf");

  
  TH1D *h_dbp = (TH1D*)f->Get("control/h_dbp"); assert(h_dbp);
  TH1D *h_dbt = (TH1D*)f->Get("control/h_dbt"); assert(h_dbt);
  TH1D *h_db = (TH1D*)f->Get("control/h_db");   assert(h_db);

  TH1D *h1_db = tdrHist("h1_db","N_{events} / bin",0,900,"DB",0.,2.0);
  lumi_136TeV = "JMENano Summer24 DY";
  TCanvas *c1_db = tdrCanvas("c1_db",h1_db,8,11,kSquare);

  tdrDraw(h_dbp,"HIST",kNone,kGray+1,kSolid,-1,1001,kGray);
  tdrDraw(h_dbt,"HIST",kNone,kRed-9+1,kSolid,-1,1001,kRed-9);
  tdrDraw(h_db,"Pz",kFullCircle,kGreen+2,kSolid,-1,kNone,0,0.5);

  TF1 *f1_db = new TF1("f1_db","gaus",0.6,1./0.6);
  f1_eta->SetParameters(200,0.9,0.3);
  h_db->Fit(f1_db,"RN");
  f1_db->SetRange(0.5,1.8);
  f1_db->SetLineColor(kGreen+2);
  f1_db->Draw("SAME");
  
  gPad->RedrawAxis();
  
  TLegend *leg_db = tdrLeg(0.50,0.87-4*0.05,0.80,0.87);
  leg_db->SetHeader("|#Delta#phi|<#pi/16, |#Deltap_{T}/p_{T}|<0.4");
  leg_db->AddEntry(h_dbp,"Parallel");
  leg_db->AddEntry(h_dbt,"Transverse");
  leg_db->AddEntry(h_db,"Par-tran");

  tex->DrawLatex(0.50,0.62,"Veto leptons from");
  tex->DrawLatex(0.50,0.57,"par+tran |#Delta#phi|<#pi/8");

  c1_db->SaveAs("pdf/drawControl/drawControl_c1_3_db.pdf");


  TH1D *h_mpfp = (TH1D*)f->Get("control/h_mpfp"); assert(h_dbp);
  TH1D *h_mpft = (TH1D*)f->Get("control/h_mpft"); assert(h_mpft);
  TH1D *h_mpf = (TH1D*)f->Get("control/h_mpf");   assert(h_mpf);

  TH1D *h1_mpf = tdrHist("h1_mpf","N_{events} / bin",0,0.5*900,"MPF",-3.0,4.0);
  lumi_136TeV = "JMENano Summer24 DY";
  TCanvas *c1_mpf = tdrCanvas("c1_mpf",h1_mpf,8,11,kSquare);

  tdrDraw(h_mpfp,"HIST",kNone,kGray+1,kSolid,-1,1001,kGray);
  tdrDraw(h_mpft,"HIST",kNone,kRed-9+1,kSolid,-1,1001,kRed-9);
  tdrDraw(h_mpf,"Pz",kFullCircle,kGreen+2,kSolid,-1,kNone,0,0.5);

  TF1 *f1_mpf = new TF1("f1_mpf","gaus",0.6,1./0.6);
  f1_eta->SetParameters(200,0.9,0.3);
  h_mpf->Fit(f1_mpf,"RN");
  f1_mpf->SetRange(0.5,1.8);
  f1_mpf->SetLineColor(kGreen+2);
  f1_mpf->Draw("SAME");
  
  gPad->RedrawAxis();
  
  TLegend *leg_mpf = tdrLeg(0.50,0.87-4*0.05,0.80,0.87);
  leg_mpf->SetHeader("|#Delta#phi|<#pi/16, |#Deltap_{T}/p_{T}|<0.4");
  leg_mpf->AddEntry(h_mpfp,"Parallel");
  leg_mpf->AddEntry(h_mpft,"Transverse");
  leg_mpf->AddEntry(h_mpf,"Par-tran");

  tex->DrawLatex(0.50,0.62,"Veto leptons from");
  tex->DrawLatex(0.50,0.57,"par+tran |#Delta#phi|<#pi/8");

  c1_mpf->SaveAs("pdf/drawControl/drawControl_c1_4_mpf.pdf");


  TH2D *h2_parpteta = (TH2D*)f->Get("control/h2_parpteta");
  assert(h2_parpteta);
  TH2D *h2_tranpteta = (TH2D*)f->Get("control/h2_tranpteta");
  assert(h2_tranpteta);
  TH2D *h2_mixpteta = (TH2D*)f->Get("control/h2_mixpteta");
  assert(h2_mixpteta);

  //TH1D *h2_pteta = tdrHist("h2_pteta","p_{T,Z} (GeV)",0,200,"#eta_{jet}",-5,5);
  TCanvas *c2_pteta = new TCanvas("c2_pteta","",1800,600);
  c2_pteta->Divide(3,1,0,0);

  c2_pteta->cd(1);
  h2_parpteta->Draw("COLZ");
  h2_parpteta->UseCurrentStyle();
  h2_parpteta->GetZaxis()->SetRangeUser(1,150);
  h2_parpteta->GetXaxis()->SetTitle("p_{T,Z} (GeV)");
  h2_parpteta->GetYaxis()->SetTitle("#eta_{jet}");

  tex->DrawLatex(0.65,0.90,"Parallel");
  gPad->RedrawAxis();
  
  c2_pteta->cd(2);
  h2_tranpteta->Draw("COLZ");
  h2_tranpteta->UseCurrentStyle();
  h2_tranpteta->GetZaxis()->SetRangeUser(1,150);
  h2_tranpteta->GetXaxis()->SetTitle("p_{T,Z} (GeV)");
  h2_tranpteta->GetYaxis()->SetTitle("#eta_{jet}");
  
  tex->DrawLatex(0.65,0.90,"Transverse");
  gPad->RedrawAxis();
  
  c2_pteta->cd(3);
  h2_mixpteta->Draw("COLZ");
  h2_mixpteta->UseCurrentStyle();
  h2_mixpteta->GetZaxis()->SetRangeUser(1,150);
  h2_mixpteta->GetXaxis()->SetTitle("p_{T,Z} (GeV)");
  h2_mixpteta->GetYaxis()->SetTitle("#eta_{jet}");
  
  tex->DrawLatex(0.65,0.90,"Par - tran");
  gPad->RedrawAxis();
  
  c2_pteta->SaveAs("pdf/drawControl/drawControl_c2_1_pteta.pdf");


  TProfile2D *p2_db = (TProfile2D*)f->Get("control/p2_db");
  assert(p2_db);
  TProfile2D *p2_dbp = (TProfile2D*)f->Get("control/p2_dbp");
  assert(p2_db;);
  TProfile2D *p2_dbt = (TProfile2D*)f->Get("control/p2_dbt");
  assert(p2_dbt);

  TCanvas *c2_db = new TCanvas("c2_db","",1800,600);
  c2_db->Divide(3,1,0,0);

  c2_db->cd(1);
  p2_dbp->Draw("COLZ");
  p2_dbp->UseCurrentStyle();
  p2_dbp->GetZaxis()->SetRangeUser(0.5,1.5);
  p2_dbp->GetXaxis()->SetTitle("p_{T,Z} (GeV)");
  p2_dbp->GetYaxis()->SetTitle("#eta_{jet}");

  tex->DrawLatex(0.65,0.90,"Parallel");
  gPad->RedrawAxis();
  
  c2_db->cd(2);
  p2_dbt->Draw("COLZ");
  p2_dbt->UseCurrentStyle();
  p2_dbt->GetZaxis()->SetRangeUser(0.5,1.5);
  p2_dbt->GetXaxis()->SetTitle("p_{T,Z} (GeV)");
  p2_dbt->GetYaxis()->SetTitle("#eta_{jet}");
  
  tex->DrawLatex(0.65,0.90,"Transverse");
  gPad->RedrawAxis();
  
  c2_db->cd(3);
  p2_db->Draw("COLZ");
  p2_db->UseCurrentStyle();
  p2_db->GetZaxis()->SetRangeUser(0.5,1.5);
  p2_db->GetXaxis()->SetTitle("p_{T,Z} (GeV)");
  p2_db->GetYaxis()->SetTitle("#eta_{jet}");
  
  tex->DrawLatex(0.65,0.90,"Par - tran");
  gPad->RedrawAxis();
  
  c2_db->SaveAs("pdf/drawControl/drawControl_c2_2_db.pdf");
  

  TProfile2D *p2_mpf = (TProfile2D*)f->Get("control/p2_mpf");
  assert(p2_mpf);
  TProfile2D *p2_mpfp = (TProfile2D*)f->Get("control/p2_mpfp");
  assert(p2_mpf;);
  TProfile2D *p2_mpft = (TProfile2D*)f->Get("control/p2_mpft");
  assert(p2_mpft);

  TCanvas *c2_mpf = new TCanvas("c2_mpf","",1800,600);
  c2_mpf->Divide(3,1,0,0);

  c2_mpf->cd(1);
  p2_mpfp->Draw("COLZ");
  p2_mpfp->UseCurrentStyle();
  p2_mpfp->GetZaxis()->SetRangeUser(0.5,1.5);
  p2_mpfp->GetXaxis()->SetTitle("p_{T,Z} (GeV)");
  p2_mpfp->GetYaxis()->SetTitle("#eta_{jet}");

  tex->DrawLatex(0.65,0.90,"Parallel");
  gPad->RedrawAxis();
  
  c2_mpf->cd(2);
  p2_mpft->Draw("COLZ");
  p2_mpft->UseCurrentStyle();
  p2_mpft->GetZaxis()->SetRangeUser(0.5,1.5);
  p2_mpft->GetXaxis()->SetTitle("p_{T,Z} (GeV)");
  p2_mpft->GetYaxis()->SetTitle("#eta_{jet}");
  
  tex->DrawLatex(0.65,0.90,"Transverse");
  gPad->RedrawAxis();
  
  c2_mpf->cd(3);
  p2_mpf->Draw("COLZ");
  p2_mpf->UseCurrentStyle();
  p2_mpf->GetZaxis()->SetRangeUser(0.5,1.5);
  p2_mpf->GetXaxis()->SetTitle("p_{T,Z} (GeV)");
  p2_mpf->GetYaxis()->SetTitle("#eta_{jet}");
  
  tex->DrawLatex(0.65,0.90,"Par - tran");
  gPad->RedrawAxis();
  
  c2_mpf->SaveAs("pdf/drawControl/drawControl_c2_3_mpf.pdf");


  TProfile2D *p2_mpfn = (TProfile2D*)f->Get("control/p2_mpfn");
  assert(p2_mpfn);
  TProfile2D *p2_mpfnp = (TProfile2D*)f->Get("control/p2_mpfnp");
  assert(p2_mpfn;);
  TProfile2D *p2_mpfnt = (TProfile2D*)f->Get("control/p2_mpfnt");
  assert(p2_mpfnt);

  TCanvas *c2_mpfn = new TCanvas("c2_mpfn","",1800,600);
  c2_mpfn->Divide(3,1,0,0);

  c2_mpfn->cd(1);
  p2_mpfnp->Draw("COLZ");
  p2_mpfnp->UseCurrentStyle();
  p2_mpfnp->GetZaxis()->SetRangeUser(-0.5,0.5);
  p2_mpfnp->GetXaxis()->SetTitle("p_{T,Z} (GeV)");
  p2_mpfnp->GetYaxis()->SetTitle("#eta_{jet}");

  tex->DrawLatex(0.65,0.90,"Parallel");
  gPad->RedrawAxis();
  
  c2_mpfn->cd(2);
  p2_mpfnt->Draw("COLZ");
  p2_mpfnt->UseCurrentStyle();
  p2_mpfnt->GetZaxis()->SetRangeUser(-0.5,0.5);
  p2_mpfnt->GetXaxis()->SetTitle("p_{T,Z} (GeV)");
  p2_mpfnt->GetYaxis()->SetTitle("#eta_{jet}");
  
  tex->DrawLatex(0.65,0.90,"Transverse");
  gPad->RedrawAxis();
  
  c2_mpfn->cd(3);
  p2_mpfn->Draw("COLZ");
  p2_mpfn->UseCurrentStyle();
  p2_mpfn->GetZaxis()->SetRangeUser(-0.5,0.5);
  p2_mpfn->GetXaxis()->SetTitle("p_{T,Z} (GeV)");
  p2_mpfn->GetYaxis()->SetTitle("#eta_{jet}");
  
  tex->DrawLatex(0.65,0.90,"Par - tran");
  gPad->RedrawAxis();
  
  c2_mpfn->SaveAs("pdf/drawControl/drawControl_c2_4_mpfn.pdf");


  TProfile2D *p2_mpfu = (TProfile2D*)f->Get("control/p2_mpfu");
  assert(p2_mpfu);
  TProfile2D *p2_mpfup = (TProfile2D*)f->Get("control/p2_mpfup");
  assert(p2_mpfu;);
  TProfile2D *p2_mpfut = (TProfile2D*)f->Get("control/p2_mpfut");
  assert(p2_mpfut);

  TCanvas *c2_mpfu = new TCanvas("c2_mpfu","",1800,600);
  c2_mpfu->Divide(3,1,0,0);

  c2_mpfu->cd(1);
  p2_mpfup->Draw("COLZ");
  p2_mpfup->UseCurrentStyle();
  p2_mpfup->GetZaxis()->SetRangeUser(-0.5,0.5);
  p2_mpfup->GetXaxis()->SetTitle("p_{T,Z} (GeV)");
  p2_mpfup->GetYaxis()->SetTitle("#eta_{jet}");

  tex->DrawLatex(0.65,0.90,"Parallel");
  gPad->RedrawAxis();
  
  c2_mpfu->cd(2);
  p2_mpfut->Draw("COLZ");
  p2_mpfut->UseCurrentStyle();
  p2_mpfut->GetZaxis()->SetRangeUser(-0.5,0.5);
  p2_mpfut->GetXaxis()->SetTitle("p_{T,Z} (GeV)");
  p2_mpfut->GetYaxis()->SetTitle("#eta_{jet}");
  
  tex->DrawLatex(0.65,0.90,"Transverse");
  gPad->RedrawAxis();
  
  c2_mpfu->cd(3);
  p2_mpfu->Draw("COLZ");
  p2_mpfu->UseCurrentStyle();
  p2_mpfu->GetZaxis()->SetRangeUser(-0.5,0.5);
  p2_mpfu->GetXaxis()->SetTitle("p_{T,Z} (GeV)");
  p2_mpfu->GetYaxis()->SetTitle("#eta_{jet}");
  
  tex->DrawLatex(0.65,0.90,"Par - tran");
  gPad->RedrawAxis();
  
  c2_mpfu->SaveAs("pdf/drawControl/drawControl_c2_5_mpfu.pdf");


  TProfile2D *p2_mpfnu = (TProfile2D*)f->Get("control/p2_mpfnu");
  assert(p2_mpfnu);
  TProfile2D *p2_mpfnup = (TProfile2D*)f->Get("control/p2_mpfnup");
  assert(p2_mpfnu;);
  TProfile2D *p2_mpfnut = (TProfile2D*)f->Get("control/p2_mpfnut");
  assert(p2_mpfnut);

  TCanvas *c2_mpfnu = new TCanvas("c2_mpfnu","",1800,600);
  c2_mpfnu->Divide(3,1,0,0);

  c2_mpfnu->cd(1);
  p2_mpfnup->Draw("COLZ");
  p2_mpfnup->UseCurrentStyle();
  p2_mpfnup->GetZaxis()->SetRangeUser(-0.5,0.5);
  p2_mpfnup->GetXaxis()->SetTitle("p_{T,Z} (GeV)");
  p2_mpfnup->GetYaxis()->SetTitle("#eta_{jet}");

  tex->DrawLatex(0.65,0.90,"Parallel");
  gPad->RedrawAxis();
  
  c2_mpfnu->cd(2);
  p2_mpfnut->Draw("COLZ");
  p2_mpfnut->UseCurrentStyle();
  p2_mpfnut->GetZaxis()->SetRangeUser(-0.5,0.5);
  p2_mpfnut->GetXaxis()->SetTitle("p_{T,Z} (GeV)");
  p2_mpfnut->GetYaxis()->SetTitle("#eta_{jet}");
  
  tex->DrawLatex(0.65,0.90,"Transverse");
  gPad->RedrawAxis();
  
  c2_mpfnu->cd(3);
  p2_mpfnu->Draw("COLZ");
  p2_mpfnu->UseCurrentStyle();
  p2_mpfnu->GetZaxis()->SetRangeUser(-0.5,0.5);
  p2_mpfnu->GetXaxis()->SetTitle("p_{T,Z} (GeV)");
  p2_mpfnu->GetYaxis()->SetTitle("#eta_{jet}");
  
  tex->DrawLatex(0.65,0.90,"Par - tran");
  gPad->RedrawAxis();
  
  c2_mpfnu->SaveAs("pdf/drawControl/drawControl_c2_6_mpfnu.pdf");

} // drawControl
