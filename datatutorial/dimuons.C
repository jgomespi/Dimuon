#define dimuon_cxx
#include "dimuon.h"
#include "TGraph.h"
#include "TPaveLabel.h"
#include  <math.h>

void dimuons() {

  TChain * chain = new TChain("oniaTree","");
  chain->Add("Skim4.root");

  //chain->Show();
  //chain->Scan("*");

  dimuon a(chain);
  a._nev = 100000;
  //Jpsi:2.9-3.3GeV
  //Z:75-110 GeV

  a.GetSpectrum();

  a.SelectPeak(20.);

  a.FitPeak("hDimuonMass_mypeak_pTcut20");

  //a.PtEta(20.);

}

// this is the main processing function
void dimuon::GetSpectrum() {

  // check tree
  if (fChain == 0) return;

  // create and fill a simple mass histogram
  TH1F *hDimuonMass_normal = new TH1F("hDimuonMass_normal","hDimuonMass_normal",10000,0.2,300);
  FillHisto(hDimuonMass_normal);
  SaveHisto(hDimuonMass_normal);

  // now set log scales
  SaveHisto(hDimuonMass_normal,kTRUE);

  //define another (special) histogram: with variable (!) bin widths 

  double xbins[100000];
  xbins[0] = .1; 
  int nbins = 0;
  double binWidth=0.005; 
  for (int i=1; xbins[i-1]<500; i++) {
    xbins[i] = xbins[i-1]*(1+binWidth);
    nbins++;
  }
  TH1F *hDimuonMass = new TH1F("hDimuonMass","hDimuonMass",nbins,xbins);
  FillHisto(hDimuonMass);
  SaveHisto(hDimuonMass,kTRUE);

  // now: normalize yields (to adapt to variable binning!)
  for (int i=1; i<=hDimuonMass->GetNbinsX(); i++) {
    hDimuonMass->SetBinContent(i, hDimuonMass->GetBinContent(i)/hDimuonMass->GetBinWidth(i));
  }
  SaveHisto(hDimuonMass,kTRUE);

}

void dimuon::SaveHisto(TH1F* hist, Int_t log) {

  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);

  hist->GetXaxis()->SetTitle("#mu^{+}#mu^{-} invariant mass [GeV]");
  hist->GetYaxis()->SetTitle("Events / GeV");

  TCanvas *c = new TCanvas("c","c",800,600);

  if(log) {
    c->SetLogx();
    c->SetLogy();
  }

  hist->Draw("HIST");

  TString hn ("");
  hn += "plots/";
  hn += hist->GetName();
  if(log) hn += "_log";
  hn += ".png";
  c->SaveAs(hn);
  delete c;

  //TH1F* h2 = (TH1F*)hist->Clone();
  //h2->SetName(hn);
  _outFile->cd();
  hist->Draw("HIST");
  hist->Write();
  _outFile->Write();

}
// add one more par to the FillHisto function: the cut on each muon pT
void dimuon::FillHisto(TH1F* hist, Double_t pT_cut) {

  // loop over the tree, and fill the histograms
  Long64_t nentries = fChain->GetEntriesFast();
  nentries = _nev>0 ? _nev : nentries;
  //cout << "nentries: " << nentries << endl;
  TH2F *pTHist = new TH2F("pTHist","pT Hist",100,0.,300,100,0,0.3);
  TH2F *etaHist = new TH2F("etaHist","pT Hist",100,-5.,5.,100,0,0.3);

  Long64_t nbytes = 0, nb = 0;
  Double_t eta[nentries], pT[nentries], err[nentries];

  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    
    if ( Cut(ientry) < 0) continue;

    double mass = dimuon_p4->M();
    double pT = dimuon_p4->Pt();
    double eta = dimuon_p4->Eta();
    double err = abs(mass*mass-90.731*90.731)/90.731;

    if (muonN_p4->Pt() >= pT_cut && muonP_p4->Pt() >= pT_cut) hist->Fill(mass);  
    if (muonN_p4->Pt() >= pT_cut && muonP_p4->Pt() >= pT_cut && mass <= 110. && mass >= 65.){
       pTHist->Fill(pT,err);
       etaHist->Fill(eta,err);
    }
    
  }

  TCanvas *c10 = new TCanvas("pt_eta","pT_eta",800,600);
  //TPaveLabel* title = new TPaveLabel(0.1,0.96,0.9,0.99,"Espaco de fase na regiao de sinal (65 < M < 110 GeV)");
  //title->Draw();  
  pTHist->GetXaxis()->SetTitle("pT do par #mu^{+}#mu^{-} (GeV)");
  pTHist->GetYaxis()->SetTitle("#sigma_{rel}");
  pTHist->GetYaxis()->SetTitleOffset(1.2);
  etaHist->GetXaxis()->SetTitle("#eta do par #mu^{+}#mu^{-}");
  etaHist->GetYaxis()->SetTitle("#sigma_{rel}");  
  etaHist->GetYaxis()->SetTitleOffset(1.2);  
  c10->Divide(2,1);
  c10->cd(1);
  c10->SetLogy();
  pTHist->Draw("colz");
  c10->cd(2);
  etaHist->Draw("colz");
  c10->SaveAs("resultadocarai.png");
  c10->Close();
  c10->Delete();
}
/*
void dimuon::PtEta(Double_t pT_cut) {
  Long64_t nentries = fChain->GetEntriesFast();
  nentries = _nev>0 ? _nev : nentries;
  cout << "nentries: " << nentries << endl;

  Long64_t nbytes = 0, nb = 0;
  Double_t eta[33371], pT[33371], err[33371];
  int i=0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    
    Long64_t ientry = LoadTree(jentry);
    int entry = int(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    
    if ( Cut(ientry) < 0) continue;
    
    double mass = dimuon_p4->M();

    if (muonN_p4->Pt() >= pT_cut && muonP_p4->Pt() >= pT_cut && dimuon_p4->M() >= 65. && dimuon_p4->M() <= 110.) {
    	//std::cout << dimuon_p4->Eta() << " " << endl;
    	eta[i] = dimuon_p4->Eta();
    	pT[i] = dimuon_p4->Pt();
    	err[i] = abs((mass-90.731)/90.731);
    	i++; 
    }
  }
  //TCanvas *c = new TCanvas("c","A Simple Graph Example",800,600);
  //TGraph* Eta = new TGraph(nentries,eta,err);
  //TGraph* PT = new TGraph(nentries,pT,err);
  //Eta->Draw("AC*");
  //PT->Draw();
  //c->SaveAs("yields.png");
  //for (Long64_t jentry=0; jentry<nentries;jentry++)  std::cout << eta[jentry] << " " << endl;
  std::cout << std::vector::size(eta) << std::endl;
  //for (int i=0; i<33371;i++) std::cout << eta[i] << " " << endl;    

}
*/
void dimuon::SelectPeak(Double_t pT_cut) {
  
  Double_t mmin(65.), mmax(110.);
  _mmin = mmin;  _mmax = mmax;  
  if (pT_cut!=0) SpT_cut = Form("%i",int(pT_cut));

  // create an histogram around a peak
  TH1F *hDimuonMass_mypeak = new TH1F("hDimuonMass_mypeak_pTcut"+SpT_cut,"hDimuonMass_mypeak",100,mmin,mmax);
  FillHisto(hDimuonMass_mypeak,pT_cut);
  SaveHisto(hDimuonMass_mypeak);

}

void dimuon::FitPeak(TString hist_name) {
 
  if(!_outFile) {cout << "Check input file." << endl; return;}

  // retrive histogram with selected peak
  TH1F* hpeak= 0;
  TString hname(hist_name);
  _outFile->GetObject(hname,hpeak);
  if (!hpeak) {
    cout << "Check input histogram:" << hname <<  endl;
    return;
  }

  // define fit function and fit the histogram
  const Int_t nfitpar(5);
  TF1* f = new TF1("f",fitfun2,_mmin,_mmax,nfitpar);
  // change the mean parameter to 95.
  f->SetParameters(100,95.,5.,0,0);
  hpeak->Fit("f");

  // write fit results into array
  Double_t par[nfitpar];
  f->GetParameters(par);

  printf("\nFitResults:\n\tResonance mass: %5.3f +/- %5.3f GeV/c^2.\n",
	 par[1],f->GetParErrors()[1]);

  //return;   // comment for continuing

  // what follows is aesthetics, mostly ...

  gROOT->LoadMacro("tdrstyle.C");

  TCanvas *c0 = new TCanvas("peak","peak",800,600);
  //c0->SetFillColor(0);
  //c0->SetFrameFillColor(0);
  //c0->SetGrid();

  hpeak->GetXaxis()->SetTitle("#mu^{+}#mu^{-} invariant mass [GeV/c^{2}]");
  hpeak->GetYaxis()->SetTitle(Form("Events / %3.1f MeV/c^{2}",hpeak->GetBinWidth(1)*1000));
  hpeak->SetStats(0);
  hpeak->SetTitle("");
  hpeak->SetMarkerStyle(21);
  hpeak->SetMarkerSize(0.8);

  hpeak->Fit("f","V+","ep");

  // get the individual functions for separate representation 
  TF1 *signalFcn = new TF1("signalFcn",signal,_mmin,_mmax,3);
  signalFcn->SetLineColor(kBlue);
  signalFcn->SetNpx(500);
  TF1 *backFcn = new TF1("backFcn",backgr,_mmin,_mmax,2);
  backFcn->SetLineColor(kGray);
  backFcn->SetLineStyle(2);
  TF1 *signalFcn2 = new TF1("signalFcn2",signal2,_mmin,_mmax,3);
  signalFcn2->SetLineColor(kGreen);

  signalFcn2->SetParameters(par);
  signalFcn2->Draw("same");
  
  backFcn->SetParameters(&par[3]);
  backFcn->Draw("same");
    
  // draw the legend
  TLegend *legend=new TLegend(0.7,0.65,0.88,0.85);
  legend->SetBorderSize(0);
  legend->SetTextFont(40);
  legend->SetTextSize(0.03);
  legend->AddEntry(hpeak,"Data","lpe");
  legend->AddEntry(backFcn,"Background fit","l");
  legend->AddEntry(signalFcn2,"Signal fit","l");
  legend->AddEntry(f,"Global Fit","l");
  legend->Draw("same");

  // display info + fit results  
  TLatex L;
  L.SetNDC();
  L.SetTextSize(0.04);
  L.DrawLatex(0.15,0.8,"Dimuon Spectrum");
  L.SetTextSize(0.03);
  L.DrawLatex(0.15,0.75,"resonance: Z^{0} boson");
  L.DrawLatex(0.15,0.70,Form("mass: %5.3f #pm %5.3f GeV/c^{2}",
			     par[1], f->GetParErrors()[1]));
  L.DrawLatex(0.15,0.65,Form("with: %5.3f #pm %5.3f MeV/c^{2}", 
			     par[2]*1000, f->GetParErrors()[2]*1000));
  L.DrawLatex(0.15,0.6, Form("Entries: %5.0f ", hpeak->GetEntries()));
  // save the fitted histogram 
  c0->SaveAs("plots/mypeak_"+hist_name+".png");

}


Double_t signal(Double_t *x, Double_t *par) {
  //a simple gaussian
  return par[0]*exp(-0.5*TMath::Power(((x[0]-par[1])/(par[2])),2)); 
}

Double_t signal2(Double_t *x, Double_t *par) {
  //a Breit-Wigner
  Double_t arg1 = 14.0/22.0; // 2 over pi
  Double_t arg2 = par[2]*par[2]*par[1]*par[1]; //Gamma=par[1]  M=par[2]
  Double_t arg3 = ((x[0]*x[0]) - (par[1]*par[1]))*((x[0]*x[0]) - (par[1]*par[1]));
  Double_t arg4 = x[0]*x[0]*x[0]*x[0]*((par[2]*par[2])/(par[1]*par[1]));
  return par[0]*arg1*arg2/(arg3 + arg4);
  }

Double_t backgr(Double_t *x, Double_t *par) {
  //a simple polynomial
  return par[0]+par[2]*x[0];
}

Double_t fitfun(Double_t *x, Double_t *par) {
  //the total PDF function, sum of the above
  return signal(x,par) + backgr(x,&par[3]); 
}

Double_t fitfun2(Double_t *x, Double_t *par) {
  //the total PDF function, sum of the above
  return signal2(x,par) + backgr(x,&par[3]); 
}
