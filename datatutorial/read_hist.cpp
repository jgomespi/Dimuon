#include "TH1.h"
#include "TTree.h"
#include "TInterpreter.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TH2.h"
#include "TNtuple.h"
#include "TPaveLabel.h"
#include "TPaveText.h"
#include "TFrame.h"
#include "TSystem.h"
#include "TInterpreter.h"

void read_hist() {
  gROOT->ForceStyle();
  gStyle->SetOptStat(0); // for the sumperimposed plots
  //gStyle->SetOptStat("n");
  //gStyle->SetOptFit(1);	// for the fit plot

  TFile *file = new TFile("plots/histos.root","read");
  file->cd();

  TH1F * mypeak_pTcut0 = (TH1F*)file->Get("hDimuonMass_mypeak_pTcut0");
  TH1F * mypeak_pTcut10 = (TH1F*)file->Get("hDimuonMass_mypeak_pTcut10");
  TH1F * mypeak_pTcut20 = (TH1F*)file->Get("hDimuonMass_mypeak_pTcut20");
  TH1F * mypeak_pTcut30 = (TH1F*)file->Get("hDimuonMass_mypeak_pTcut30");


  mypeak_pTcut0->SetLineColor(2);  
  mypeak_pTcut10->SetLineColor(3);  
  mypeak_pTcut20->SetLineColor(6);
  mypeak_pTcut30->SetLineColor(9);
 
 /*
  Double_t mean_multH50 = multH50->GetMean();
  TString mean_multH50string ;
  mean_multH50string = Form("%.3f",mean_multH50);

  Double_t mean_multH100 = multH100->GetMean();
  TString mean_multH100string ;
  mean_multH100string = Form("%.3f",mean_multH100);

  Double_t mean_multH150 = multH150->GetMean();
  TString mean_multH150string ;
  mean_multH150string = Form("%.3f",mean_multH150);

  Double_t mean_multH200 = multH200->GetMean();
  TString mean_multH200string ;
  mean_multH200string = Form("%.3f",mean_multH200);

  Double_t mean_multH250 = multH250->GetMean();
  TString mean_multH250string ;
  mean_multH250string = Form("%.3f",mean_multH250);  

  Double_t mean_multH300 = multH300->GetMean();
  TString mean_multH300string ;
  mean_multH300string = Form("%.3f",mean_multH300);  

  Double_t mean_multH350 = multH350->GetMean();
  TString mean_multH350string ;
  mean_multH350string = Form("%.3f",mean_multH350);  

  Double_t mean_multH400 = multH400->GetMean();
  TString mean_multH400string ;
  mean_multH400string = Form("%.3f",mean_multH400);  

  TF1 pois("pois",poissonf,0,10,2);
  pois.SetParameter(0,1); 
  pois.SetParameter(1,1);  
*/
  auto legend = new TLegend(0.1,0.6,0.4,0.9);
  legend->SetHeader("p_{T} each muon","C"); // option "C" allows to center the header  
  legend->AddEntry(mypeak_pTcut0, "p_T > 0 GeV","f");
  legend->AddEntry(mypeak_pTcut10, "p_T > 10 GeV","f");
  legend->AddEntry(mypeak_pTcut20, "p_T > 20 GeV","f");
  legend->AddEntry(mypeak_pTcut30, "p_T > 30 GeV","f");

  TCanvas* c1 = new TCanvas("c1","#mu_{-}#mu_{+} invariant mass",800,800);  
  mypeak_pTcut0->SetTitle("#mu_{-}#mu_{+} invariant mass");
  mypeak_pTcut0->GetYaxis()->SetTitleOffset(1.5);
  mypeak_pTcut0->GetXaxis()->SetTitle("M_{#mu#mu} (GeV)");
  mypeak_pTcut0->GetYaxis()->SetTitle("Entries");
  mypeak_pTcut0->Draw();
  legend->Draw("SAME");
  mypeak_pTcut10->Draw("SAME");
  mypeak_pTcut20->Draw("SAME");
  mypeak_pTcut30->Draw("SAME");
  c1->SaveAs("inv_mass_pt.pdf","pdf");
  c1->SaveAs("inv_mass_pt.png","png");
  c1->Close();
}