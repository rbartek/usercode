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


void SpacialVariables2d(){

  /* Macro to take first look at Zbb minitupes
   * hbb variable produced by Zbb_minituples.cc
   * Author:  Rachel Wilken - UNL
   */

    TFile *ZH_H115_Zbb = TFile::Open("/Users/rachelwilken/Documents/CMS/Zbb/AliceVariables/Nov15/ZH_M115.root");
	TFile *ZZ_Zbb = TFile::Open("/Users/rachelwilken/Documents/CMS/Zbb/AliceVariables/Nov15/ZZ.root");
	TFile *DYJetsToLL_M50_Zbb = TFile::Open("/Users/rachelwilken/Documents/CMS/Zbb/AliceVariables/Nov15/DYJetsToLL_M50.root");

	//TFile *ZH_H115_Zbb = TFile::Open("/afs/cern.ch/user/w/wilken/scratch0/CMSSW_4_2_5/src/UserCode/wilken/output.root");


TString suffixps = ".gif";
TString directory = "Nov15/Signal/";
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


    TH1F* hMjj_ZH_H115_Zbb			= (TH1F*) ZH_H115_Zbb->Get("hMjj");
    TH1F* hCutFlow_ZH_H115_Zbb		= (TH1F*) ZH_H115_Zbb->Get("hCutFlow");
    TH2F* hDphiDetajj_ZH_H115_Zbb		= (TH2F*) ZH_H115_Zbb->Get("hDphiDetajj");
	TH1F* hMjj_aftercuts_ZH_H115_Zbb		= (TH1F*) ZH_H115_Zbb->Get("hMjj_aftercuts");
	TH2F* hqtvsalphaJJ_ZH_H115_Zbb		= (TH2F*) ZH_H115_Zbb->Get("hqtvsalphaJJ");
	TH2F* hqtvsalphaZ_ZH_H115_Zbb		= (TH2F*) ZH_H115_Zbb->Get("hqtvsalphaZ");
	
	TH2F* hqtvsalphaJJ_ZZ_Zbb	= (TH2F*) ZZ_Zbb->Get("hqtvsalphaJJ");
    TH2F* hqtvsalphaZ_ZZ_Zbb	= (TH2F*) ZZ_Zbb->Get("hqtvsalphaZ");
    TH2F* hqtvsalphaJJ_DYJetsToLL_M50_Zbb	= (TH2F*) DYJetsToLL_M50_Zbb->Get("hqtvsalphaJJ");
    TH2F* hqtvsalphaZ_DYJetsToLL_M50_Zbb	= (TH2F*) DYJetsToLL_M50_Zbb->Get("hqtvsalphaZ");
	
	TH2F* hDphiDetaJJ_aftercuts_ZH_H115_Zbb		= (TH2F*) ZH_H115_Zbb->Get("hDphiDetajj_aftercuts");
	TH2F* hqtvsalphaJJ_aftercuts_ZH_H115_Zbb		= (TH2F*) ZH_H115_Zbb->Get("hqtvsalphaJJ_aftercuts");
	TH2F* hqtvsalphaZ_aftercuts_ZH_H115_Zbb		= (TH2F*) ZH_H115_Zbb->Get("hqtvsalphaZ_aftercuts");
	TH2F* hqtvsalphaJJ_aftercuts_ZZ_Zbb	= (TH2F*) ZZ_Zbb->Get("hqtvsalphaJJ_aftercuts");
    TH2F* hqtvsalphaZ_aftercuts_ZZ_Zbb	= (TH2F*) ZZ_Zbb->Get("hqtvsalphaZ_aftercuts");
    TH2F* hqtvsalphaJJ_aftercuts_DYJetsToLL_M50_Zbb	= (TH2F*) DYJetsToLL_M50_Zbb->Get("hqtvsalphaJJ_aftercuts");
    TH2F* hqtvsalphaZ_aftercuts_DYJetsToLL_M50_Zbb	= (TH2F*) DYJetsToLL_M50_Zbb->Get("hqtvsalphaZ_aftercuts");
	

//Histograms	
	hMjj_aftercuts_ZH_H115_Zbb->SetTitle("Mass of dijet pair");
	hMjj_aftercuts_ZH_H115_Zbb->GetXaxis()->SetTitle("Mass (GeV)");
	hMjj_aftercuts_ZH_H115_Zbb->GetYaxis()->SetTitle(" ");
	hMjj_aftercuts_ZH_H115_Zbb->SetLineColor(kBlue);
	hMjj_aftercuts_ZH_H115_Zbb->Draw("hist");
	//gStyle->SetOptStat(kFALSE);
	plot = directory+"Mjj_afterCuts"+suffixps;
	c1->Print(plot);
	c1->Clear();
	nsigmas = 1.5;    // Number of sigmas around mean to fit gaussian.  It uses 2 iterations 
	// i.e. range is = [mu - nsigmas * sigma, mu + nsigmas * sigma]
	hMjj_aftercuts_ZH_H115_Zbb->Fit("gaus");
	TF1 *myfunc = hMjj_aftercuts_ZH_H115_Zbb->GetFunction("gaus");
	float par0 = gaus->GetParameter(0); 
	float par1 = gaus->GetParameter(1); 
	float par2 = gaus->GetParameter(2); 
	cout << "Parameters are: " << "P0: " << par0  
	<<  " P1: " << par1 << " P2: " <<par2 << endl;
	float low = par1 -nsigmas * par2;
	float hi = par1 + nsigmas * par2;
	hMjj_aftercuts_ZH_H115_Zbb->Fit("gaus","R","",low,hi);
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
	hMjj_aftercuts_ZH_H115_Zbb->Fit("gaus","R","",low,hi);	
	plot = directory+"Mjj_aftercuts_GausFit"+suffixps;
	c1->Print(plot);
	c1->Clear();
	
	hMjj_ZH_H115_Zbb->SetTitle("Mass of dijet pair Afer Cuts");
	hMjj_ZH_H115_Zbb->GetXaxis()->SetTitle("Mass (GeV)");
	hMjj_ZH_H115_Zbb->GetYaxis()->SetTitle(" ");
	hMjj_ZH_H115_Zbb->SetLineColor(kBlue);
	hMjj_ZH_H115_Zbb->Draw("hist");
	//gStyle->SetOptStat(kFALSE);
	plot = directory+"Mjj"+suffixps;
	c1->Print(plot);
	c1->Clear();
	nsigmas = 1.5;    // Number of sigmas around mean to fit gaussian.  It uses 2 iterations 
	// i.e. range is = [mu - nsigmas * sigma, mu + nsigmas * sigma]
	hMjj_ZH_H115_Zbb->Fit("gaus");
	TF1 *myfunc = hMjj_ZH_H115_Zbb->GetFunction("gaus");
	float par0 = gaus->GetParameter(0); 
	float par1 = gaus->GetParameter(1); 
	float par2 = gaus->GetParameter(2); 
	 cout << "Parameters are: " << "P0: " << par0  
	     <<  " P1: " << par1 << " P2: " <<par2 << endl;
	 float low = par1 -nsigmas * par2;
	 float hi = par1 + nsigmas * par2;
	 hMjj_ZH_H115_Zbb->Fit("gaus","R","",low,hi);
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
	 hMjj_ZH_H115_Zbb->Fit("gaus","R","",low,hi);	
	plot = directory+"MjjGausFit"+suffixps;
	c1->Print(plot);
	c1->Clear();
	
	
	hDphiDetajj_ZH_H115_Zbb->SetTitle("#Delta#phi vs #Delta#eta of jets");
	hDphiDetajj_ZH_H115_Zbb->GetXaxis()->SetTitle("#Delta#eta");
	hDphiDetajj_ZH_H115_Zbb->GetYaxis()->SetTitle("#Delta#phi");
	gStyle->SetPalette(1);
	hDphiDetajj_ZH_H115_Zbb->SetContour(5);
	hDphiDetajj_ZH_H115_Zbb->Draw("colz");
	plot = directory+"DphiDetajj"+suffixps;
	c1->Print(plot);
	c1->Clear();	
	
	hqtvsalphaJJ_ZH_H115_Zbb->SetTitle(" Amenteros-Poldolansky of jets");
	hqtvsalphaJJ_ZH_H115_Zbb->GetXaxis()->SetTitle("#alpha");
	hqtvsalphaJJ_ZH_H115_Zbb->GetYaxis()->SetTitle("q_{T} (GeV)");
	gStyle->SetPalette(1);
	hqtvsalphaJJ_ZH_H115_Zbb->SetContour(5);
	hqtvsalphaJJ_ZH_H115_Zbb->Draw("colz");
	plot = directory+"qtvsalphaJJ"+suffixps;
	c1->Print(plot);
	c1->Clear();	
	
	hqtvsalphaZ_ZH_H115_Zbb->SetTitle(" Amenteros-Poldolansky of muons");
	hqtvsalphaZ_ZH_H115_Zbb->GetXaxis()->SetTitle("#alpha");
	hqtvsalphaZ_ZH_H115_Zbb->GetYaxis()->SetTitle("q_{T} (GeV)");
	gStyle->SetPalette(1);
	hqtvsalphaZ_ZH_H115_Zbb->SetContour(5);
	hqtvsalphaZ_ZH_H115_Zbb->Draw("colz");
	plot = directory+"qtvsalphaZ"+suffixps;
	c1->Print(plot);
	c1->Clear();	
	
	
	hqtvsalphaJJ_ZZ_Zbb->SetTitle(" Amenteros-Poldolansky of jets");
	hqtvsalphaJJ_ZZ_Zbb->GetXaxis()->SetTitle("#alpha");
	hqtvsalphaJJ_ZZ_Zbb->GetYaxis()->SetTitle("q_{T} (GeV)");
	gStyle->SetPalette(1);
	hqtvsalphaJJ_ZZ_Zbb->SetContour(5);
	hqtvsalphaJJ_ZZ_Zbb->Draw("colz");
	plot = "Nov15/qtvsalphaJJ_ZZ"+suffixps;
	c1->Print(plot);
	c1->Clear();	
	
	hqtvsalphaZ_ZZ_Zbb->SetTitle(" Amenteros-Poldolansky of muons");
	hqtvsalphaZ_ZZ_Zbb->GetXaxis()->SetTitle("#alpha");
	hqtvsalphaZ_ZZ_Zbb->GetYaxis()->SetTitle("q_{T} (GeV)");
	gStyle->SetPalette(1);
	hqtvsalphaZ_ZZ_Zbb->SetContour(5);
	hqtvsalphaZ_ZZ_Zbb->Draw("colz");
	plot = "Nov15/qtvsalphaZ_ZZ"+suffixps;
	c1->Print(plot);
	c1->Clear();	
	
	hqtvsalphaJJ_DYJetsToLL_M50_Zbb->SetTitle(" Amenteros-Poldolansky of jets");
	hqtvsalphaJJ_DYJetsToLL_M50_Zbb->GetXaxis()->SetTitle("#alpha");
	hqtvsalphaJJ_DYJetsToLL_M50_Zbb->GetYaxis()->SetTitle("q_{T} (GeV)");
	gStyle->SetPalette(1);
	hqtvsalphaJJ_DYJetsToLL_M50_Zbb->SetContour(5);
	hqtvsalphaJJ_DYJetsToLL_M50_Zbb->Draw("colz");
	plot = "Nov15/qtvsalphaJJ_DY"+suffixps;
	c1->Print(plot);
	c1->Clear();	
	
	hqtvsalphaZ_DYJetsToLL_M50_Zbb->SetTitle(" Amenteros-Poldolansky of muons");
	hqtvsalphaZ_DYJetsToLL_M50_Zbb->GetXaxis()->SetTitle("#alpha");
	hqtvsalphaZ_DYJetsToLL_M50_Zbb->GetYaxis()->SetTitle("q_{T} (GeV)");
	gStyle->SetPalette(1);
	hqtvsalphaZ_DYJetsToLL_M50_Zbb->SetContour(5);
	hqtvsalphaZ_DYJetsToLL_M50_Zbb->Draw("colz");
	plot = "Nov15/qtvsalphaZ_DY"+suffixps;
	c1->Print(plot);
	c1->Clear();
	
	 hCutFlow_ZH_H115_Zbb->GetXaxis()->SetBinLabel(1,"PreSelection");
	 hCutFlow_ZH_H115_Zbb->GetXaxis()->SetBinLabel(2,"p_{T}(jj)");
	 hCutFlow_ZH_H115_Zbb->GetXaxis()->SetBinLabel(3,"p_{T}(Z)");
	 hCutFlow_ZH_H115_Zbb->GetXaxis()->SetBinLabel(4,"CSV1");
	 hCutFlow_ZH_H115_Zbb->GetXaxis()->SetBinLabel(5,"CSV2");
	 hCutFlow_ZH_H115_Zbb->GetXaxis()->SetBinLabel(6,"#Delta#phi");
	 hCutFlow_ZH_H115_Zbb->GetXaxis()->SetBinLabel(7,"N_{aj}");
	 hCutFlow_ZH_H115_Zbb->GetXaxis()->SetBinLabel(8,"M_{jj}");	
	hCutFlow_ZH_H115_Zbb->SetTitle("Cut Flow HZ");
	hCutFlow_ZH_H115_Zbb->GetYaxis()->SetTitle("Number of Events");
	hCutFlow_ZH_H115_Zbb->SetLineColor(kBlack);
	hCutFlow_ZH_H115_Zbb->Draw("hist");
	gStyle->SetOptStat(kFALSE);
	plot = directory+"CutFlow"+suffixps;
	c1->Print(plot);
	c1->Clear();
	
	hDphiDetaJJ_aftercuts_ZH_H115_Zbb->SetTitle(" #Delta#eta vs #Delta#phi of jets After Cuts");
	hDphiDetaJJ_aftercuts_ZH_H115_Zbb->GetXaxis()->SetTitle("#Delta#eta");
	hDphiDetaJJ_aftercuts_ZH_H115_Zbb->GetYaxis()->SetTitle("#Delta#phi");
	gStyle->SetPalette(1);
	hDphiDetaJJ_aftercuts_ZH_H115_Zbb->SetContour(5);
	hDphiDetaJJ_aftercuts_ZH_H115_Zbb->Draw("colz");
	plot = directory+"DphiDetajj_aftercuts"+suffixps;
	c1->Print(plot);
	c1->Clear();	
	
	hqtvsalphaJJ_aftercuts_ZH_H115_Zbb->SetTitle(" Amenteros-Poldolansky of jets After Cuts H115");
	hqtvsalphaJJ_aftercuts_ZH_H115_Zbb->GetXaxis()->SetTitle("#alpha");
	hqtvsalphaJJ_aftercuts_ZH_H115_Zbb->GetYaxis()->SetTitle("q_{T} (GeV)");
	gStyle->SetPalette(1);
	hqtvsalphaJJ_aftercuts_ZH_H115_Zbb->SetContour(5);
	hqtvsalphaJJ_aftercuts_ZH_H115_Zbb->Draw("colz");
	plot = directory+"qtvsalphaJJ_aftercuts"+suffixps;
	c1->Print(plot);
	c1->Clear();	
	
	hqtvsalphaZ_aftercuts_ZH_H115_Zbb->SetTitle(" Amenteros-Poldolansky of muons After Cuts H115");
	hqtvsalphaZ_aftercuts_ZH_H115_Zbb->GetXaxis()->SetTitle("#alpha");
	hqtvsalphaZ_aftercuts_ZH_H115_Zbb->GetYaxis()->SetTitle("q_{T} (GeV)");
	gStyle->SetPalette(1);
	hqtvsalphaZ_aftercuts_ZH_H115_Zbb->SetContour(5);
	hqtvsalphaZ_aftercuts_ZH_H115_Zbb->Draw("colz");
	plot = directory+"qtvsalphaZ_aftercuts"+suffixps;
	c1->Print(plot);
	c1->Clear();	
	
	
	hqtvsalphaJJ_aftercuts_ZZ_Zbb->SetTitle(" Amenteros-Poldolansky of jets After Cuts ZZ");
	hqtvsalphaJJ_aftercuts_ZZ_Zbb->GetXaxis()->SetTitle("#alpha");
	hqtvsalphaJJ_aftercuts_ZZ_Zbb->GetYaxis()->SetTitle("q_{T} (GeV)");
	gStyle->SetPalette(1);
	hqtvsalphaJJ_aftercuts_ZZ_Zbb->SetContour(5);
	hqtvsalphaJJ_aftercuts_ZZ_Zbb->Draw("colz");
	plot = "Nov15/qtvsalphaJJ_aftercuts_ZZ"+suffixps;
	c1->Print(plot);
	c1->Clear();	
	
	hqtvsalphaZ_aftercuts_ZZ_Zbb->SetTitle(" Amenteros-Poldolansky of muons After Cuts ZZ");
	hqtvsalphaZ_aftercuts_ZZ_Zbb->GetXaxis()->SetTitle("#alpha");
	hqtvsalphaZ_aftercuts_ZZ_Zbb->GetYaxis()->SetTitle("q_{T} (GeV)");
	gStyle->SetPalette(1);
	hqtvsalphaZ_aftercuts_ZZ_Zbb->SetContour(5);
	hqtvsalphaZ_aftercuts_ZZ_Zbb->Draw("colz");
	plot = "Nov15/qtvsalphaZ_aftercuts_ZZ"+suffixps;
	c1->Print(plot);
	c1->Clear();	
	
	hqtvsalphaJJ_aftercuts_DYJetsToLL_M50_Zbb->SetTitle(" Amenteros-Poldolansky of jets After Cuts DY");
	hqtvsalphaJJ_aftercuts_DYJetsToLL_M50_Zbb->GetXaxis()->SetTitle("#alpha");
	hqtvsalphaJJ_aftercuts_DYJetsToLL_M50_Zbb->GetYaxis()->SetTitle("q_{T} (GeV)");
	gStyle->SetPalette(1);
	hqtvsalphaJJ_aftercuts_DYJetsToLL_M50_Zbb->SetContour(5);
	hqtvsalphaJJ_aftercuts_DYJetsToLL_M50_Zbb->Draw("colz");
	plot = "Nov15/qtvsalphaJJ_aftercuts_DY"+suffixps;
	c1->Print(plot);
	c1->Clear();	
	
	hqtvsalphaZ_aftercuts_DYJetsToLL_M50_Zbb->SetTitle(" Amenteros-Poldolansky of muons After Cuts DY");
	hqtvsalphaZ_aftercuts_DYJetsToLL_M50_Zbb->GetXaxis()->SetTitle("#alpha");
	hqtvsalphaZ_aftercuts_DYJetsToLL_M50_Zbb->GetYaxis()->SetTitle("q_{T} (GeV)");
	gStyle->SetPalette(1);
	hqtvsalphaZ_aftercuts_DYJetsToLL_M50_Zbb->SetContour(5);
	hqtvsalphaZ_aftercuts_DYJetsToLL_M50_Zbb->Draw("colz");
	plot = "Nov15/qtvsalphaZ_aftercuts_DY"+suffixps;
	c1->Print(plot);
	c1->Clear();	
	
}
