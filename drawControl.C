// Purpose: draw control plots of the metho
#include "TH2D.h"
#include "TF1.h"

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

  c1_pt->SaveAs("pdf/drawControl/drawControl_c1_pt.pdf");
  

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

  c1_eta->SaveAs("pdf/drawControl/drawControl_c1_eta.pdf");
  
} // drawControl
