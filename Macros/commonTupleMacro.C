#include <TH2.h>
#include <TStyle.h>
#include <TFile.h>
#include "TH2F.h"
#include "TCanvas.h"
#include "TRandom.h"
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <strstream>
#include <iomanip.h>
//#include "TStopwatch.h"
#include "TROOT.h"
#include "TF2.h"
#include "TMath.h"
#include "TVirtualPad.h"
#include "TStyle.h"
#include "Riostream.h"
#include "TColor.h"
#include "TClass.h"
//#include <vector>

using namespace std;


void commonTupleMacro(){

  /* Macro to take first look at Zbb minitupes
   * hbb variable produced by Zbb_minituples.cc
   * Author:  Rachel Wilken - UNL
   */

    TFile *H115_Zbb = TFile::Open("/Users/rachelwilken/Documents/CMS/Zbb/CommonNtuples/output.root");
	//TFile *H115_Zbb = TFile::Open("/afs/cern.ch/user/w/wilken/scratch0/CMSSW_4_2_5/src/UserCode/wilken/output.root");


TString suffixps = ".gif";
TString directory = "Halloween/Signal/";
TString plot = directory+"somereallylongnameherehopethisarrayisbiggenough"+suffixps;
// Files for histogram output --> set suffixps to desired file type:  e.g. .eps, .jpg, ...

 
 	gStyle->SetOptStat(kFALSE);
	TCanvas *c1 = new TCanvas("c1","");
	gStyle->SetOptStat("emruo");
	c1->SetFillColor(10);
	c1->SetFillColor(10);
	double normalization;

// ********************************************************************
// Invariant Mass of Dijet pair
// ********************************************************************

    TH1F* hnJets_H115_Zbb			= (TH1F*) H115_Zbb->Get("hnJets");
    TH1F* hnMuons_H115_Zbb		= (TH1F*) H115_Zbb->Get("hnMuons");
	TH1F* hPtb1_H115_Zbb			= (TH1F*) H115_Zbb->Get("hPtb1");
    TH1F* hPtb2_H115_Zbb			= (TH1F*) H115_Zbb->Get("hPtb2");
    TH1F* hPtjj_H115_Zbb			= (TH1F*) H115_Zbb->Get("hPtjj");
    TH1F* hPtmumu_H115_Zbb		= (TH1F*) H115_Zbb->Get("hPtmumu");
    TH1F* hPtbalZH_H115_Zbb		= (TH1F*) H115_Zbb->Get("hPtbalZH");
    TH1F* hEtamu1_H115_Zbb			= (TH1F*) H115_Zbb->Get("hPtmu1");
    TH1F* hEtamu2_H115_Zbb			= (TH1F*) H115_Zbb->Get("hPtmu2");
    TH1F* hPtmu1_H115_Zbb			= (TH1F*) H115_Zbb->Get("hEtamu1");
    TH1F* hPtmu2_H115_Zbb			= (TH1F*) H115_Zbb->Get("hEtamu2");
    TH1F* hCSV1_H115_Zbb			= (TH1F*) H115_Zbb->Get("hCSV1");
    TH1F* hCSV2_H115_Zbb			= (TH1F*) H115_Zbb->Get("hCSV2");
    TH1F* hdphiVH_H115_Zbb		= (TH1F*) H115_Zbb->Get("hdphiVH");
    TH1F* hdetaJJ_H115_Zbb		= (TH1F*) H115_Zbb->Get("hdetaJJ");
    TH1F* hdphiJJ_H115_Zbb		= (TH1F*) H115_Zbb->Get("hdphiJJ");
    TH1F* hNaj_H115_Zbb			= (TH1F*) H115_Zbb->Get("hNaj");
    TH1F* hMmumu_H115_Zbb			= (TH1F*) H115_Zbb->Get("hMmumu");
    TH1F* hMjj_H115_Zbb			= (TH1F*) H115_Zbb->Get("hMjj");
    TH1F* hCutFlow_H115_Zbb		= (TH1F*) H115_Zbb->Get("hCutFlow");
    TH2F* hDphiDetajj_H115_Zbb		= (TH2F*) H115_Zbb->Get("hDphiDetajj");
	TH1F* hMjj_aftercuts_H115_Zbb		= (TH1F*) H115_Zbb->Get("hMjj_aftercuts");

//Histograms
	hnJets_H115_Zbb->SetTitle("Number of Good Jets");
	hnJets_H115_Zbb->GetXaxis()->SetTitle("Njets");
	hnJets_H115_Zbb->GetYaxis()->SetTitle(" ");
	hnJets_H115_Zbb->SetLineColor(kBlue);
	hnJets_H115_Zbb->Draw("hist");
	//gStyle->SetOptStat(kFALSE);
	plot = directory+"nJets"+suffixps;
	c1->Print(plot);
	c1->Clear();

	hnMuons_H115_Zbb->SetTitle("Number of Good Muonns");
	hnMuons_H115_Zbb->GetXaxis()->SetTitle("Nmuons");
	hnMuons_H115_Zbb->GetYaxis()->SetTitle(" ");
	hnMuons_H115_Zbb->SetLineColor(kBlue);
	hnMuons_H115_Zbb->Draw("hist");
	//gStyle->SetOptStat(kFALSE);
	plot = directory+"Nmuons"+suffixps;
	c1->Print(plot);
	c1->Clear();
	
	hPtb1_H115_Zbb->SetTitle("Pt of Jet with highest CSV");
	hPtb1_H115_Zbb->GetXaxis()->SetTitle("p_{T} (GeV)");
	hPtb1_H115_Zbb->GetYaxis()->SetTitle(" ");
	hPtb1_H115_Zbb->SetLineColor(kBlue);
	hPtb1_H115_Zbb->Draw("hist");
	//gStyle->SetOptStat(kFALSE);
	plot = directory+"Ptb1"+suffixps;
	c1->Print(plot);
	c1->Clear();
	
	hPtb2_H115_Zbb->SetTitle("Pt of Jet with second highest CSV");
	hPtb2_H115_Zbb->GetXaxis()->SetTitle("p_{T} (GeV)");
	hPtb2_H115_Zbb->GetYaxis()->SetTitle(" ");
	hPtb2_H115_Zbb->SetLineColor(kBlue);
	hPtb2_H115_Zbb->Draw("hist");
	//gStyle->SetOptStat(kFALSE);
	plot = directory+"Ptb2"+suffixps;
	c1->Print(plot);
	c1->Clear();	

	hPtjj_H115_Zbb->SetTitle("Pt of dijet pair");
	hPtjj_H115_Zbb->GetXaxis()->SetTitle("p_{T} (GeV)");
	hPtjj_H115_Zbb->GetYaxis()->SetTitle(" ");
	hPtjj_H115_Zbb->SetLineColor(kBlue);
	hPtjj_H115_Zbb->Draw("hist");
	//gStyle->SetOptStat(kFALSE);
	plot = directory+"Ptjj"+suffixps;
	c1->Print(plot);
	c1->Clear();

	hPtmumu_H115_Zbb->SetTitle("Pt of dimuon pair");
	hPtmumu_H115_Zbb->GetXaxis()->SetTitle("p_{T} (GeV)");
	hPtmumu_H115_Zbb->GetYaxis()->SetTitle(" ");
	hPtmumu_H115_Zbb->SetLineColor(kBlue);
	hPtmumu_H115_Zbb->Draw("hist");
	//gStyle->SetOptStat(kFALSE);
	plot = directory+"Ptmumu"+suffixps;
	c1->Print(plot);
	c1->Clear();

	hPtmu1_H115_Zbb->SetTitle("Muon with highest Pt");
	hPtmu1_H115_Zbb->GetXaxis()->SetTitle("p_{T} (GeV)");
	hPtmu1_H115_Zbb->GetYaxis()->SetTitle(" ");
	hPtmu1_H115_Zbb->SetLineColor(kBlue);
	hPtmu1_H115_Zbb->Draw("hist");
	//gStyle->SetOptStat(kFALSE);
	plot = directory+"Ptmu1"+suffixps;
	c1->Print(plot);
	c1->Clear();

	hPtmu2_H115_Zbb->SetTitle("Muon with second highest P");
	hPtmu2_H115_Zbb->GetXaxis()->SetTitle("p_{T} (GeV)");
	hPtmu2_H115_Zbb->GetYaxis()->SetTitle(" ");
	hPtmu2_H115_Zbb->SetLineColor(kBlue);
	hPtmu2_H115_Zbb->Draw("hist");
	//gStyle->SetOptStat(kFALSE);
	plot = directory+"Ptmu2"+suffixps;
	c1->Print(plot);
	c1->Clear();

	hEtamu1_H115_Zbb->SetTitle("Eta of Muon with highest Pt");
	hEtamu1_H115_Zbb->GetXaxis()->SetTitle("#eta");
	hEtamu1_H115_Zbb->GetYaxis()->SetTitle(" ");
	hEtamu1_H115_Zbb->SetLineColor(kBlue);
	hEtamu1_H115_Zbb->Draw("hist");
	//gStyle->SetOptStat(kFALSE);
	plot = directory+"Etamu1"+suffixps;
	c1->Print(plot);
	c1->Clear();
	
	hEtamu2_H115_Zbb->SetTitle("Eta of Muon with second highest P");
	hEtamu2_H115_Zbb->GetXaxis()->SetTitle("#eta");
	hEtamu2_H115_Zbb->GetYaxis()->SetTitle(" ");
	hEtamu2_H115_Zbb->SetLineColor(kBlue);
	hEtamu2_H115_Zbb->Draw("hist");
	//gStyle->SetOptStat(kFALSE);
	plot = directory+"Etamu2"+suffixps;
	c1->Print(plot);
	c1->Clear();
	

	hCSV1_H115_Zbb->SetTitle("First Jet ordered in CSV");
	hCSV1_H115_Zbb->GetXaxis()->SetRangeUser(0,1.2);
	hCSV1_H115_Zbb->GetXaxis()->SetTitle("CSV");
	hCSV1_H115_Zbb->GetYaxis()->SetTitle(" ");
	hCSV1_H115_Zbb->SetLineColor(kBlue);
	hCSV1_H115_Zbb->Draw("hist");
	//gStyle->SetOptStat(kFALSE);
	plot = directory+"CSV1"+suffixps;
	c1->Print(plot);
	c1->Clear();

	hCSV2_H115_Zbb->SetTitle("Second Jet ordered in CSV");
	hCSV2_H115_Zbb->GetXaxis()->SetRangeUser(0,1.2);
	hCSV2_H115_Zbb->GetXaxis()->SetTitle("CSV");
	hCSV2_H115_Zbb->GetYaxis()->SetTitle(" ");
	hCSV2_H115_Zbb->SetLineColor(kBlue);
	hCSV2_H115_Zbb->Draw("hist");
	//gStyle->SetOptStat(kFALSE);
	plot = directory+"CSV2"+suffixps;
	c1->Print(plot);
	c1->Clear();
	

	hdphiVH_H115_Zbb->SetTitle("Delta Phi of Vector Boson and Higgs");
	hdphiVH_H115_Zbb->GetXaxis()->SetTitle("#Delta#phi");
	hdphiVH_H115_Zbb->GetYaxis()->SetTitle(" ");
	hdphiVH_H115_Zbb->SetLineColor(kBlue);
	hdphiVH_H115_Zbb->Draw("hist");
	//gStyle->SetOptStat(kFALSE);
	plot = directory+"dphiVH"+suffixps;
	c1->Print(plot);
	c1->Clear();

	hdetaJJ_H115_Zbb->SetTitle("Delta Eta of dijet pair");
	hdetaJJ_H115_Zbb->GetXaxis()->SetTitle("#Delta#eta");
	hdetaJJ_H115_Zbb->GetYaxis()->SetTitle(" ");
	hdetaJJ_H115_Zbb->SetLineColor(kBlue);
	hdetaJJ_H115_Zbb->Draw("hist");
	//gStyle->SetOptStat(kFALSE);
	plot = directory+"DeltaEtaJJ"+suffixps;
	c1->Print(plot);
	c1->Clear();

	hdphiJJ_H115_Zbb->SetTitle("Delta Phi of dijet pair");
	hdphiJJ_H115_Zbb->GetXaxis()->SetTitle("#Delta#phi");
	hdphiJJ_H115_Zbb->GetYaxis()->SetTitle(" ");
	hdphiJJ_H115_Zbb->SetLineColor(kBlue);
	hdphiJJ_H115_Zbb->Draw("hist");
	//gStyle->SetOptStat(kFALSE);
	plot = directory+"DeltaPhiJJ"+suffixps;
	c1->Print(plot);
	c1->Clear();	

	hNaj_H115_Zbb->SetTitle("Number of Jets in Addition to two b jets");
	hNaj_H115_Zbb->GetXaxis()->SetTitle("Njets");
	hNaj_H115_Zbb->GetYaxis()->SetTitle(" ");
	hNaj_H115_Zbb->SetLineColor(kBlue);
	hNaj_H115_Zbb->Draw("hist");
	//gStyle->SetOptStat(kFALSE);
	plot = directory+"Naj"+suffixps;
	c1->Print(plot);
	c1->Clear();

	hMmumu_H115_Zbb->SetTitle("Mass of dimuon pair");
	hMmumu_H115_Zbb->GetXaxis()->SetRangeUser(30,130);
	hMmumu_H115_Zbb->GetXaxis()->SetTitle("Mass (GeV)");
	hMmumu_H115_Zbb->GetYaxis()->SetTitle(" ");
	hMmumu_H115_Zbb->SetLineColor(kBlue);
	hMmumu_H115_Zbb->Draw("hist");
	//gStyle->SetOptStat(kFALSE);
	plot = directory+"Mmumu"+suffixps;
	c1->Print(plot);
	c1->Clear();
	
	hPtbalZH_H115_Zbb->SetTitle("Pt balance of Z and H");
	//hPtbalZH_H115_Zbb->GetXaxis()->SetRangeUser(30,130);
	hPtbalZH_H115_Zbb->GetXaxis()->SetTitle("#Delta p_{T}^{HZ}");
	hPtbalZH_H115_Zbb->GetYaxis()->SetTitle(" ");
	hPtbalZH_H115_Zbb->SetLineColor(kBlue);
	hPtbalZH_H115_Zbb->Draw("hist");
	//gStyle->SetOptStat(kFALSE);
	plot = directory+"PtBalZH"+suffixps;
	c1->Print(plot);
	c1->Clear();
	
	hMjj_aftercuts_H115_Zbb->SetTitle("Mass of dijet pair");
	hMjj_aftercuts_H115_Zbb->GetXaxis()->SetTitle("Mass (GeV)");
	hMjj_aftercuts_H115_Zbb->GetYaxis()->SetTitle(" ");
	hMjj_aftercuts_H115_Zbb->SetLineColor(kBlue);
	hMjj_aftercuts_H115_Zbb->Draw("hist");
	//gStyle->SetOptStat(kFALSE);
	plot = directory+"Mjj_afterCuts"+suffixps;
	c1->Print(plot);
	c1->Clear();
	nsigmas = 1.5;    // Number of sigmas around mean to fit gaussian.  It uses 2 iterations 
	// i.e. range is = [mu - nsigmas * sigma, mu + nsigmas * sigma]
	hMjj_aftercuts_H115_Zbb->Fit("gaus");
	TF1 *myfunc = hMjj_aftercuts_H115_Zbb->GetFunction("gaus");
	float par0 = gaus->GetParameter(0); 
	float par1 = gaus->GetParameter(1); 
	float par2 = gaus->GetParameter(2); 
	cout << "Parameters are: " << "P0: " << par0  
	<<  " P1: " << par1 << " P2: " <<par2 << endl;
	float low = par1 -nsigmas * par2;
	float hi = par1 + nsigmas * par2;
	hMjj_aftercuts_H115_Zbb->Fit("gaus","R","",low,hi);
	par0 = gaus->GetParameter(0); 
	par1 = gaus->GetParameter(1); 
	par2 = gaus->GetParameter(2); 
	cout << "********* Second fit *********" << endl;
	cout << "Parameters are: " << "P0: " << par0  
	<<  " P1: " << par1 << " P2: " << par2 << endl;
	gStyle->SetOptStat(kTRUE);
	gStyle->SetOptFit(0111);
	low = par1 -nsigmas * par2;
	hi = par1 + nsigmas * par2;
	hMjj_aftercuts_H115_Zbb->Fit("gaus","R","",low,hi);	
	plot = directory+"Mjj_aftercuts_GausFit"+suffixps;
	c1->Print(plot);
	c1->Clear();
	
	hMjj_H115_Zbb->SetTitle("Mass of dijet pair Afer Cuts");
	hMjj_H115_Zbb->GetXaxis()->SetTitle("Mass (GeV)");
	hMjj_H115_Zbb->GetYaxis()->SetTitle(" ");
	hMjj_H115_Zbb->SetLineColor(kBlue);
	hMjj_H115_Zbb->Draw("hist");
	//gStyle->SetOptStat(kFALSE);
	plot = directory+"Mjj"+suffixps;
	c1->Print(plot);
	c1->Clear();
	nsigmas = 1.5;    // Number of sigmas around mean to fit gaussian.  It uses 2 iterations 
	// i.e. range is = [mu - nsigmas * sigma, mu + nsigmas * sigma]
	hMjj_H115_Zbb->Fit("gaus");
	TF1 *myfunc = hMjj_H115_Zbb->GetFunction("gaus");
	float par0 = gaus->GetParameter(0); 
	float par1 = gaus->GetParameter(1); 
	float par2 = gaus->GetParameter(2); 
	 cout << "Parameters are: " << "P0: " << par0  
	     <<  " P1: " << par1 << " P2: " <<par2 << endl;
	 float low = par1 -nsigmas * par2;
	 float hi = par1 + nsigmas * par2;
	 hMjj_H115_Zbb->Fit("gaus","R","",low,hi);
	 par0 = gaus->GetParameter(0); 
	 par1 = gaus->GetParameter(1); 
	 par2 = gaus->GetParameter(2); 
	 cout << "********* Second fit *********" << endl;
	 cout << "Parameters are: " << "P0: " << par0  
	      <<  " P1: " << par1 << " P2: " << par2 << endl;
	 gStyle->SetOptStat(kTRUE);
	 gStyle->SetOptFit(0111);
	 low = par1 -nsigmas * par2;
	 hi = par1 + nsigmas * par2;
	 hMjj_H115_Zbb->Fit("gaus","R","",low,hi);	
	plot = directory+"MjjGausFit"+suffixps;
	c1->Print(plot);
	c1->Clear();
	
	
	hDphiDetajj_H115_Zbb->SetTitle(" #Delta#eta vs #Delta#phi of jets");
	hDphiDetajj_H115_Zbb->GetXaxis()->SetTitle("#Delta#eta");
	hDphiDetajj_H115_Zbb->GetYaxis()->SetTitle("#Delta#phi");
	gStyle->SetPalette(1);
	hDphiDetajj_H115_Zbb->SetContour(5);
	hDphiDetajj_H115_Zbb->Draw("colz");
	plot = directory+"DphiDetajj"+suffixps;
	c1->Print(plot);
	c1->Clear();	
	
	
	 hCutFlow_H115_Zbb->GetXaxis()->SetBinLabel(1,"PreSelection");
	 hCutFlow_H115_Zbb->GetXaxis()->SetBinLabel(2,"p_{T}(jj)");
	 hCutFlow_H115_Zbb->GetXaxis()->SetBinLabel(3,"p_{T}(Z)");
	 hCutFlow_H115_Zbb->GetXaxis()->SetBinLabel(4,"CSV1");
	 hCutFlow_H115_Zbb->GetXaxis()->SetBinLabel(5,"CSV2");
	 hCutFlow_H115_Zbb->GetXaxis()->SetBinLabel(6,"#Delta#phi");
	 hCutFlow_H115_Zbb->GetXaxis()->SetBinLabel(7,"N_{aj}");
	 hCutFlow_H115_Zbb->GetXaxis()->SetBinLabel(8,"M_{jj}");	
	hCutFlow_H115_Zbb->SetTitle("Cut Flow HZ");
	hCutFlow_H115_Zbb->GetYaxis()->SetTitle("Number of Events");
	hCutFlow_H115_Zbb->SetLineColor(kBlack);
	hCutFlow_H115_Zbb->Draw("hist");
	gStyle->SetOptStat(kFALSE);
	plot = directory+"CutFlow"+suffixps;
	c1->Print(plot);
	c1->Clear();
	
	
	
}
