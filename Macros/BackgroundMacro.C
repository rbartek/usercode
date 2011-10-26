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


void BackgroundMacro(){

  /* Macro to take first look at Zbb minitupes
   * hbb variable produced by Zbb_minituples.cc
   * Author:  Rachel Wilken - UNL
   */

    TFile *ZZ_Zbb = TFile::Open("/Users/rachelwilken/Documents/CMS/Zbb/RecreateAnalysis/AddVar_ZZ.root");
    TFile *DY_Zbb = TFile::Open("/Users/rachelwilken/Documents/CMS/Zbb/RecreateAnalysis/AddVar_DY.root");
	//TFile *ZZ_Zbb = TFile::Open("/afs/cern.ch/user/w/wilken/scratch0/CMSSW_4_2_5/src/UserCode/wilken/output.root");


TString suffixps = "_Bkg.gif";
TString directory = "MacroHistos/";
TString plot = directory+"somereallylongnameherehopethisarrayisbiggenough"+suffixps;
// Files for histogram output --> set suffixps to desired file type:  e.g. .eps, .jpg, ...

 
 	gStyle->SetOptStat(kFALSE);
	TCanvas *c1 = new TCanvas("c1","");
	gStyle->SetOptStat("kTRUE");
	c1->SetFillColor(10);
	c1->SetFillColor(10);
	double normalization;

// ********************************************************************
// Invariant Mass of Dijet pair
// ********************************************************************

    TH1F* hjetCSV3_ZZ_Zbb		= (TH1F*) ZZ_Zbb->Get("jetCSV3");
    TH1F* hjetPt3_ZZ_Zbb		= (TH1F*) ZZ_Zbb->Get("jetPt3");
    TH1F* hnJets_ZZ_Zbb			= (TH1F*) ZZ_Zbb->Get("hnJets");
    TH1F* hnMuons_ZZ_Zbb		= (TH1F*) ZZ_Zbb->Get("hnMuons");
    TH1F* hjetPt1_ZZ_Zbb		= (TH1F*) ZZ_Zbb->Get("jetPt1");
    TH1F* hjetPt2_ZZ_Zbb		= (TH1F*) ZZ_Zbb->Get("jetPt2");
    TH1F* hPtb1_ZZ_Zbb			= (TH1F*) ZZ_Zbb->Get("hPtb1");
    TH1F* hPtb2_ZZ_Zbb			= (TH1F*) ZZ_Zbb->Get("hPtb2");
    TH1F* hPtjj_ZZ_Zbb			= (TH1F*) ZZ_Zbb->Get("hPtjj");
    TH1F* hPtmumu_ZZ_Zbb		= (TH1F*) ZZ_Zbb->Get("hPtmumu");
    TH1F* hPtbalZH_ZZ_Zbb		= (TH1F*) ZZ_Zbb->Get("hPtbalZH");
    TH1F* hPtmu1_ZZ_Zbb			= (TH1F*) ZZ_Zbb->Get("hPtmu1");
    TH1F* hEtamu2_ZZ_Zbb			= (TH1F*) ZZ_Zbb->Get("hEtamu1");
    TH1F* hEtamu1_ZZ_Zbb			= (TH1F*) ZZ_Zbb->Get("hEtamu2");
    TH1F* hPtmu2_ZZ_Zbb			= (TH1F*) ZZ_Zbb->Get("hPtmu2");
    TH1F* hCSV1_ZZ_Zbb			= (TH1F*) ZZ_Zbb->Get("hCSV1");
    TH1F* hCSV2_ZZ_Zbb			= (TH1F*) ZZ_Zbb->Get("hCSV2");
    TH1F* hdphiVH_ZZ_Zbb		= (TH1F*) ZZ_Zbb->Get("hdphiVH");
    TH1F* hdetaJJ_ZZ_Zbb		= (TH1F*) DY_Zbb->Get("hdetaJJ");
    TH1F* hdphiJJ_ZZ_Zbb		= (TH1F*) DY_Zbb->Get("hdphiJJ");
    TH1F* hNaj_ZZ_Zbb			= (TH1F*) ZZ_Zbb->Get("hNaj");
    TH1F* hMmumu_ZZ_Zbb			= (TH1F*) ZZ_Zbb->Get("hMmumu");
    TH1F* hMjj_ZZ_Zbb			= (TH1F*) ZZ_Zbb->Get("hMjj");
	TH1F* hMjj_aftercuts_ZZ_Zbb			= (TH1F*) ZZ_Zbb->Get("hMjj_aftercuts");
    TH1F* hCutFlow_ZZ_Zbb		= (TH1F*) ZZ_Zbb->Get("hCutFlow");
    TH2F* hDphiDetajj_ZZ_Zbb		= (TH2F*) ZZ_Zbb->Get("hDphiDetajj");
	
	TH1F* hjetCSV3_DY_Zbb		= (TH1F*) DY_Zbb->Get("jetCSV3");
    TH1F* hjetPt3_DY_Zbb		= (TH1F*) DY_Zbb->Get("jetPt3");
    TH1F* hnJets_DY_Zbb			= (TH1F*) DY_Zbb->Get("hnJets");
    TH1F* hnMuons_DY_Zbb		= (TH1F*) DY_Zbb->Get("hnMuons");
    TH1F* hjetPt1_DY_Zbb		= (TH1F*) DY_Zbb->Get("jetPt1");
    TH1F* hjetPt2_DY_Zbb		= (TH1F*) DY_Zbb->Get("jetPt2");
    TH1F* hPtb1_DY_Zbb			= (TH1F*) DY_Zbb->Get("hPtb1");
    TH1F* hPtb2_DY_Zbb			= (TH1F*) DY_Zbb->Get("hPtb2");
    TH1F* hPtjj_DY_Zbb			= (TH1F*) DY_Zbb->Get("hPtjj");
    TH1F* hPtmumu_DY_Zbb		= (TH1F*) DY_Zbb->Get("hPtmumu");
    TH1F* hPtbalZH_DY_Zbb		= (TH1F*) DY_Zbb->Get("hPtbalZH");
    TH1F* hPtmu1_DY_Zbb			= (TH1F*) DY_Zbb->Get("hPtmu1");
    TH1F* hPtmu2_DY_Zbb			= (TH1F*) DY_Zbb->Get("hPtmu2");
    TH1F* hEtamu1_DY_Zbb		= (TH1F*) DY_Zbb->Get("hEtamu1");
    TH1F* hEtamu2_DY_Zbb		= (TH1F*) DY_Zbb->Get("hEtamu2");
    TH1F* hCSV1_DY_Zbb			= (TH1F*) DY_Zbb->Get("hCSV1");
    TH1F* hCSV2_DY_Zbb			= (TH1F*) DY_Zbb->Get("hCSV2");
    TH1F* hdphiVH_DY_Zbb		= (TH1F*) DY_Zbb->Get("hdphiVH");
    TH1F* hdetaJJ_DY_Zbb		= (TH1F*) DY_Zbb->Get("hdetaJJ");
    TH1F* hdphiJJ_DY_Zbb		= (TH1F*) DY_Zbb->Get("hdphiJJ");
    TH1F* hNaj_DY_Zbb			= (TH1F*) DY_Zbb->Get("hNaj");
    TH1F* hMmumu_DY_Zbb			= (TH1F*) DY_Zbb->Get("hMmumu");
    TH1F* hMjj_DY_Zbb			= (TH1F*) DY_Zbb->Get("hMjj");
	TH1F* hMjj_aftercuts_DY_Zbb			= (TH1F*) DY_Zbb->Get("hMjj_aftercuts");
    TH1F* hCutFlow_DY_Zbb		= (TH1F*) DY_Zbb->Get("hCutFlow");
	TH2F* hDphiDetajj_DY_Zbb	= (TH2F*) DY_Zbb->Get("hDphiDetajj");

	//H115 signal weight for 10fb-1: 0.0013272727, 0.0013272727, 0.0013272704 = 0.00132727
	//ZZ weight for 10fb-1 0.014088593,0.014088595, 0.014088571 = 0.0140886
	//DY weight for 10fb-1 0.84017999, 0.84017778, 0.8401768 = 0.84018
	double H115_weight = 0.00132727;
	double ZZ_weight = 0.0140886;
	double DY_weight = 0.84018;
	
	//scale histograms
     hjetCSV3_ZZ_Zbb	->Scale(ZZ_weight);
     hjetPt3_ZZ_Zbb		->Scale(ZZ_weight);
     hnJets_ZZ_Zbb		->Scale(ZZ_weight);
     hnMuons_ZZ_Zbb		->Scale(ZZ_weight);
     hjetPt1_ZZ_Zbb		->Scale(ZZ_weight);
     hjetPt2_ZZ_Zbb		->Scale(ZZ_weight);
     hPtb1_ZZ_Zbb		->Scale(ZZ_weight);
     hPtb2_ZZ_Zbb		->Scale(ZZ_weight);
     hPtjj_ZZ_Zbb		->Scale(ZZ_weight);
     hPtmumu_ZZ_Zbb		->Scale(ZZ_weight);
	 hPtbalZH_ZZ_Zbb ->Scale(ZZ_weight);
     hPtmu1_ZZ_Zbb		->Scale(ZZ_weight);
     hPtmu2_ZZ_Zbb		->Scale(ZZ_weight);
	 hEtamu1_ZZ_Zbb->Scale(ZZ_weight);
	 hEtamu2_ZZ_Zbb->Scale(ZZ_weight);
     hCSV1_ZZ_Zbb		->Scale(ZZ_weight);
     hCSV2_ZZ_Zbb		->Scale(ZZ_weight);
     hdphiVH_ZZ_Zbb		->Scale(ZZ_weight);
	hdetaJJ_ZZ_Zbb		->Scale(ZZ_weight);
	hdphiJJ_ZZ_Zbb		->Scale(ZZ_weight);
     hNaj_ZZ_Zbb			->Scale(ZZ_weight);
     hMmumu_ZZ_Zbb		->Scale(ZZ_weight);
     hMjj_ZZ_Zbb			->Scale(ZZ_weight);
     hCutFlow_ZZ_Zbb		->Scale(ZZ_weight);
	 hDphiDetajj_ZZ_Zbb->Scale(DY_weight);
	
     hjetCSV3_DY_Zbb	->Scale(DY_weight);
     hjetPt3_DY_Zbb		->Scale(DY_weight);
     hnJets_DY_Zbb		->Scale(DY_weight);
     hnMuons_DY_Zbb		->Scale(DY_weight);
     hjetPt1_DY_Zbb		->Scale(DY_weight);
     hjetPt2_DY_Zbb		->Scale(DY_weight);
     hPtb1_DY_Zbb		->Scale(DY_weight);
     hPtb2_DY_Zbb		->Scale(DY_weight);
     hPtjj_DY_Zbb		->Scale(DY_weight);
     hPtmumu_DY_Zbb		->Scale(DY_weight);
	 hPtbalZH_DY_Zbb->Scale(DY_weight);
     hPtmu1_DY_Zbb		->Scale(DY_weight);
     hPtmu2_DY_Zbb		->Scale(DY_weight);
	hEtamu1_DY_Zbb->Scale(ZZ_weight);
	 hEtamu2_DY_Zbb->Scale(ZZ_weight);
     hCSV1_DY_Zbb		->Scale(DY_weight);
     hCSV2_DY_Zbb		->Scale(DY_weight);
     hdphiVH_DY_Zbb		->Scale(DY_weight);
	hdetaJJ_DY_Zbb		->Scale(DY_weight);
	hdphiJJ_DY_Zbb		->Scale(DY_weight);
     hNaj_DY_Zbb		->Scale(DY_weight);
     hMmumu_DY_Zbb		->Scale(DY_weight);
     hMjj_DY_Zbb		->Scale(DY_weight);
	 hMjj_aftercuts_DY_Zbb->Scale(DY_weight);
     hCutFlow_DY_Zbb		->Scale(DY_weight);
	 hDphiDetajj_DY_Zbb->Scale(DY_weight);
	
	
	// ********************************************************************
	// Have to create output file here
	// ********************************************************************	
  	// Create the root file for the histograms
	TFile *theHistogramFile;
    theHistogramFile = new TFile("Backgrounds.root", "RECREATE");
    theHistogramFile->cd();
	
	//Stack plot
	THStack *hjetCSV3_BkgStack = new THStack("hjetCSV3_BkgStack","Stacked Background Third Jet ordered in CSV");
	hjetCSV3_ZZ_Zbb->GetXaxis()->SetTitle("CSV");
	hjetCSV3_DY_Zbb->SetFillColor(kRed);
    hjetCSV3_BkgStack->Add(hjetCSV3_DY_Zbb);
	hjetCSV3_ZZ_Zbb->SetFillColor(kGreen);
    hjetCSV3_BkgStack->Add(hjetCSV3_ZZ_Zbb);
/*	hjetCSV3_tt_Zbb->SetFillColor(kBlue); 
	hjetCSV3_BkgStack->Add(hjetCSV3_tt_Zbb);	
	hjetCSV3_WZ_Zbb->SetFillColor(kOrange);
    hjetCSV3_BkgStack->Add(hjetCSV3_WZ_Zbb);
	hjetCSV3_QCD_Zbb->SetFillColor(kMagenta); 
    hjetCSV3_BkgStack->Add(hjetCSV3_QCD_Zbb);
  hjetCSV3_SingleTop_Zbb->SetFillColor(kCyan); 
  hjetCSV3_BkgStack->Add(hjetCSV3_SingleTop_Zbb);
	hjetCSV3_ZJets_Zbb->SetFillColor(kCyan); 
    hjetCSV3_BkgStack->Add(hjetCSV3_ZJets_Zbb);*/
	hjetCSV3_BkgStack->SetTitle("Third Jet ordered in CSV");
	hjetCSV3_BkgStack->Draw("hist");
	TLegend myLegend(0.7, 0.5, 0.89, 0.8);
	myLegend.SetTextSize(0.03);
	myLegend.AddEntry(hjetCSV3_DY_Zbb, "DY", "f");	
	myLegend.AddEntry(hjetCSV3_ZZ_Zbb, "ZZ", "f");	
	myLegend.Draw();		
	plot = directory+"CSV3"+suffixps;
	c1->Print(plot);
	c1->Clear();	
	
	
	//Stack plot
	THStack *hJetPt3_BkgStack = new THStack("hJetPt3_BkgStack","Stacked Background Third Jet ordered in p_{T}");
	hjetPt3_ZZ_Zbb->GetXaxis()->SetTitle("p_{T} (GeV)");
	hjetPt3_ZZ_Zbb->SetFillColor(kGreen);
    hJetPt3_BkgStack->Add(hjetPt3_ZZ_Zbb);
	hjetPt3_DY_Zbb->SetFillColor(kRed);
    hJetPt3_BkgStack->Add(hjetPt3_DY_Zbb);
	hJetPt3_BkgStack->SetTitle("Jet with third highest Pt");
	hJetPt3_BkgStack->Draw("hist");
	TLegend myLegend(0.7, 0.5, 0.89, 0.8);
	myLegend.SetTextSize(0.03);
	myLegend.AddEntry(hjetPt3_DY_Zbb, "DY", "f");	
	myLegend.AddEntry(hjetPt3_ZZ_Zbb, "ZZ", "f");	
	myLegend.Draw();		
	plot = directory+"JetPt3"+suffixps;
	c1->Print(plot);
	c1->Clear();	
		
	//Stack plot
	THStack *hnJets_BkgStack = new THStack("hnJets_BkgStack","Stacked Background Njets");
	hnJets_ZZ_Zbb->GetXaxis()->SetTitle("Njets");
	hnJets_ZZ_Zbb->SetFillColor(kGreen);
    hnJets_BkgStack->Add(hnJets_ZZ_Zbb);
	hnJets_DY_Zbb->SetFillColor(kRed);
    hnJets_BkgStack->Add(hnJets_DY_Zbb);
	hnJets_BkgStack->SetTitle("Number of Good Jets");
	hnJets_BkgStack->Draw("hist");
	TLegend myLegend(0.7, 0.5, 0.89, 0.8);
	myLegend.SetTextSize(0.03);
	myLegend.AddEntry(hnJets_DY_Zbb, "DY", "f");	
	myLegend.AddEntry(hnJets_ZZ_Zbb, "ZZ", "f");	
	myLegend.Draw();		
	plot = directory+"nJets"+suffixps;
	c1->Print(plot);
	c1->Clear();	
	
	THStack *hnMuons_BkgStack = new THStack("hnMuons_BkgStack","Stacked Background Nmuons");
	hnMuons_ZZ_Zbb->GetXaxis()->SetTitle("Nmuons");
	hnMuons_ZZ_Zbb->SetFillColor(kGreen);
    hnMuons_BkgStack->Add(hnMuons_ZZ_Zbb);
	hnMuons_DY_Zbb->SetFillColor(kRed);
    hnMuons_BkgStack->Add(hnMuons_DY_Zbb);
	hnMuons_BkgStack->SetTitle("Number of Good Muons");
	hnMuons_BkgStack->Draw("hist");
	TLegend myLegend(0.7, 0.5, 0.89, 0.8);
	myLegend.SetTextSize(0.03);
	myLegend.AddEntry(hnMuons_DY_Zbb, "DY", "f");	
	myLegend.AddEntry(hnMuons_ZZ_Zbb, "ZZ", "f");	
	myLegend.Draw();		
	plot = directory+"nMuons"+suffixps;
	c1->Print(plot);
	c1->Clear();	
	
	THStack *hJetPt2_BkgStack = new THStack("hJetPt2_BkgStack","Stacked Background Second Jet ordered in p_{T}");
	hjetPt2_ZZ_Zbb->GetXaxis()->SetTitle("p_{T} (GeV)");
	hjetPt2_ZZ_Zbb->SetFillColor(kGreen);
    hJetPt2_BkgStack->Add(hjetPt2_ZZ_Zbb);
	hjetPt2_DY_Zbb->SetFillColor(kRed);
    hJetPt2_BkgStack->Add(hjetPt2_DY_Zbb);
	hJetPt2_BkgStack->SetTitle("Jet with second highest Pt");
	hJetPt2_BkgStack->Draw("hist");
	TLegend myLegend(0.7, 0.5, 0.89, 0.8);
	myLegend.SetTextSize(0.03);
	myLegend.AddEntry(hjetPt2_DY_Zbb, "DY", "f");	
	myLegend.AddEntry(hjetPt2_ZZ_Zbb, "ZZ", "f");	
	myLegend.Draw();		
	plot = directory+"JetPt2"+suffixps;
	c1->Print(plot);
	c1->Clear();	
	
	THStack *hJetPt1_BkgStack = new THStack("hJetPt1_BkgStack","Stacked Background Highest p_{T} Jet");
	hjetPt1_ZZ_Zbb->GetXaxis()->SetTitle("p_{T} (GeV)");
	hjetPt1_ZZ_Zbb->SetFillColor(kGreen);
    hJetPt1_BkgStack->Add(hjetPt1_ZZ_Zbb);
	hjetPt1_DY_Zbb->SetFillColor(kRed);
    hJetPt1_BkgStack->Add(hjetPt1_DY_Zbb);
	hJetPt1_BkgStack->SetTitle("Jet withhighest Pt");
	hJetPt1_BkgStack->Draw("hist");
	TLegend myLegend(0.7, 0.5, 0.89, 0.8);
	myLegend.SetTextSize(0.03);
	myLegend.AddEntry(hjetPt1_DY_Zbb, "DY", "f");	
	myLegend.AddEntry(hjetPt1_ZZ_Zbb, "ZZ", "f");	
	myLegend.Draw();		
	plot = directory+"JetPt1"+suffixps;
	c1->Print(plot);
	c1->Clear();	

	THStack *hPtb1_BkgStack = new THStack("hPtb1_BkgStack","Stacked Background Pt of Jet with highest CSV");
	hPtb1_ZZ_Zbb->GetXaxis()->SetTitle("p_{T} (GeV)");
	hPtb1_ZZ_Zbb->SetFillColor(kGreen);
    hPtb1_BkgStack->Add(hPtb1_ZZ_Zbb);
	hPtb1_DY_Zbb->SetFillColor(kRed);
    hPtb1_BkgStack->Add(hPtb1_DY_Zbb);
	hPtb1_BkgStack->SetTitle("Pt of Jet with highest CSV");
	hPtb1_BkgStack->Draw("hist");
	TLegend myLegend(0.7, 0.5, 0.89, 0.8);
	myLegend.SetTextSize(0.03);
	myLegend.AddEntry(hPtb1_DY_Zbb, "DY", "f");	
	myLegend.AddEntry(hPtb1_ZZ_Zbb, "ZZ", "f");	
	myLegend.Draw();		
	plot = directory+"Ptb1"+suffixps;
	c1->Print(plot);
	c1->Clear();	
	
	THStack *hPtb2_BkgStack = new THStack("hPtb2_BkgStack","Stacked Background Pt of Jet with second highest CSV");
	hPtb2_ZZ_Zbb->GetXaxis()->SetTitle("p_{T} (GeV)");
	hPtb2_ZZ_Zbb->SetFillColor(kGreen);
    hPtb2_BkgStack->Add(hPtb2_ZZ_Zbb);
	hPtb2_DY_Zbb->SetFillColor(kRed);
    hPtb2_BkgStack->Add(hPtb2_DY_Zbb);
	hPtb2_BkgStack->SetTitle("Pt of Jet with second highest CSV");
	hPtb2_BkgStack->Draw("hist");
	TLegend myLegend(0.7, 0.5, 0.89, 0.8);
	myLegend.SetTextSize(0.03);
	myLegend.AddEntry(hPtb2_DY_Zbb, "DY", "f");	
	myLegend.AddEntry(hPtb2_ZZ_Zbb, "ZZ", "f");	
	myLegend.Draw();		
	plot = directory+"Ptb2"+suffixps;
	c1->Print(plot);
	c1->Clear();	
	
	THStack *hPtjj_BkgStack = new THStack("hPtjj_BkgStack","Stacked Background Pt of dijet pair");
	hPtjj_ZZ_Zbb->GetXaxis()->SetTitle("p_{T} (GeV)");
	hPtjj_ZZ_Zbb->SetFillColor(kGreen);
    hPtjj_BkgStack->Add(hPtjj_ZZ_Zbb);
	hPtjj_DY_Zbb->SetFillColor(kRed);
    hPtjj_BkgStack->Add(hPtjj_DY_Zbb);
	hPtjj_BkgStack->SetTitle("Pt of dijet pair");
	hPtjj_BkgStack->Draw("hist");
	TLegend myLegend(0.7, 0.5, 0.89, 0.8);
	myLegend.SetTextSize(0.03);
	myLegend.AddEntry(hPtjj_DY_Zbb, "DY", "f");	
	myLegend.AddEntry(hPtjj_ZZ_Zbb, "ZZ", "f");	
	myLegend.Draw();		
	plot = directory+"Ptjj"+suffixps;
	c1->Print(plot);
	c1->Clear();	
	
	THStack *hPtmumu_BkgStack = new THStack("hPtmumu_BkgStack","Stacked Background Pt of dimuon pair");
	hPtmumu_ZZ_Zbb->GetXaxis()->SetTitle("p_{T} (GeV)");
	hPtmumu_ZZ_Zbb->SetFillColor(kGreen);
    hPtmumu_BkgStack->Add(hPtmumu_ZZ_Zbb);
	hPtmumu_DY_Zbb->SetFillColor(kRed);
    hPtmumu_BkgStack->Add(hPtmumu_DY_Zbb);
	hPtmumu_BkgStack->SetTitle("Pt of dimuon pair");
	hPtmumu_BkgStack->Draw("hist");
	TLegend myLegend(0.7, 0.5, 0.89, 0.8);
	myLegend.SetTextSize(0.03);
	myLegend.AddEntry(hPtmumu_DY_Zbb, "DY", "f");	
	myLegend.AddEntry(hPtmumu_ZZ_Zbb, "ZZ", "f");	
	myLegend.Draw();		
	plot = directory+"Ptmumu"+suffixps;
	c1->Print(plot);
	c1->Clear();	
	
	THStack *hPtmu1_BkgStack = new THStack("hPtmu1_BkgStack","Stacked Background Pt of Muon with highest CSV");
	hPtmu1_ZZ_Zbb->GetXaxis()->SetTitle("p_{T} (GeV)");
	hPtmu1_ZZ_Zbb->SetFillColor(kGreen);
    hPtmu1_BkgStack->Add(hPtmu1_ZZ_Zbb);
	hPtmu1_DY_Zbb->SetFillColor(kRed);
    hPtmu1_BkgStack->Add(hPtmu1_DY_Zbb);
	hPtmu1_BkgStack->SetTitle("Pt of Muon with highest CSV");
	hPtmu1_BkgStack->Draw("hist");
	TLegend myLegend(0.7, 0.5, 0.89, 0.8);
	myLegend.SetTextSize(0.03);
	myLegend.AddEntry(hPtmu1_DY_Zbb, "DY", "f");	
	myLegend.AddEntry(hPtmu1_ZZ_Zbb, "ZZ", "f");	
	myLegend.Draw();		
	plot = directory+"Ptmu1"+suffixps;
	c1->Print(plot);
	c1->Clear();	
	
	THStack *hPtmu2_BkgStack = new THStack("hPtmu2_BkgStack","Stacked Background Pt of Muon with seccond highest CSV");
	hPtmu2_ZZ_Zbb->GetXaxis()->SetTitle("p_{T} (GeV)");
	hPtmu2_ZZ_Zbb->SetFillColor(kGreen);
    hPtmu2_BkgStack->Add(hPtmu2_ZZ_Zbb);
	hPtmu2_DY_Zbb->SetFillColor(kRed);
    hPtmu2_BkgStack->Add(hPtmu2_DY_Zbb);
	hPtmu2_BkgStack->SetTitle("Pt of Muon with second highest Pt");
	hPtmu2_BkgStack->Draw("hist");
	TLegend myLegend(0.7, 0.5, 0.89, 0.8);
	myLegend.SetTextSize(0.03);
	myLegend.AddEntry(hPtmu2_DY_Zbb, "DY", "f");	
	myLegend.AddEntry(hPtmu2_ZZ_Zbb, "ZZ", "f");	
	myLegend.Draw();		
	plot = directory+"Ptmu2"+suffixps;
	c1->Print(plot);
	c1->Clear();	
	
	THStack *hEtamu1_BkgStack = new THStack("hEtamu1_BkgStack","Stacked Background Eta of Muon with highest CSV");
	hEtamu1_ZZ_Zbb->GetXaxis()->SetTitle("#eta");
	hEtamu1_ZZ_Zbb->SetFillColor(kGreen);
    hEtamu1_BkgStack->Add(hEtamu1_ZZ_Zbb);
	hEtamu1_DY_Zbb->SetFillColor(kRed);
    hEtamu1_BkgStack->Add(hEtamu1_DY_Zbb);
	hEtamu1_BkgStack->SetTitle("Eta of Muon with highest Pt");
	hEtamu1_BkgStack->Draw("hist");
	TLegend myLegend(0.7, 0.5, 0.89, 0.8);
	myLegend.SetTextSize(0.03);
	myLegend.AddEntry(hEtamu1_DY_Zbb, "DY", "f");	
	myLegend.AddEntry(hEtamu1_ZZ_Zbb, "ZZ", "f");	
	myLegend.Draw();		
	plot = directory+"Etamu1"+suffixps;
	c1->Print(plot);
	c1->Clear();	
	
	THStack *hEtamu2_BkgStack = new THStack("hEtamu2_BkgStack","Stacked Background Eta of Muon with seccond highest CSV");
	hEtamu2_ZZ_Zbb->GetXaxis()->SetTitle("#eta");
	hEtamu2_ZZ_Zbb->SetFillColor(kGreen);
    hEtamu2_BkgStack->Add(hEtamu2_ZZ_Zbb);
	hEtamu2_DY_Zbb->SetFillColor(kRed);
    hEtamu2_BkgStack->Add(hEtamu2_DY_Zbb);
	hEtamu2_BkgStack->SetTitle("Eta of Muon with second highest Pt");
	hEtamu2_BkgStack->Draw("hist");
	TLegend myLegend(0.7, 0.5, 0.89, 0.8);
	myLegend.SetTextSize(0.03);
	myLegend.AddEntry(hEtamu2_DY_Zbb, "DY", "f");	
	myLegend.AddEntry(hEtamu2_ZZ_Zbb, "ZZ", "f");	
	myLegend.Draw();		
	plot = directory+"Etamu2"+suffixps;
	c1->Print(plot);
	c1->Clear();	
	
	THStack *hCSV2_BkgStack = new THStack("hCSV2_BkgStack","Stacked Background Second Jet ordered in CSV");
	hCSV2_ZZ_Zbb->GetXaxis()->SetTitle("CSV");
	hCSV2_DY_Zbb->SetFillColor(kRed);
    hCSV2_BkgStack->Add(hCSV2_DY_Zbb);
	hCSV2_ZZ_Zbb->SetFillColor(kGreen);
    hCSV2_BkgStack->Add(hCSV2_ZZ_Zbb);
	hCSV2_BkgStack->SetTitle("Second Jet ordered in CSV");
	hCSV2_BkgStack->Draw("hist");
	TLegend myLegend(0.7, 0.5, 0.89, 0.8);
	myLegend.SetTextSize(0.03);
	myLegend.AddEntry(hCSV2_DY_Zbb, "DY", "f");	
	myLegend.AddEntry(hCSV2_ZZ_Zbb, "ZZ", "f");	
	myLegend.Draw();		
	plot = directory+"CSV2"+suffixps;
	c1->Print(plot);
	c1->Clear();	
	
	THStack *hCSV1_BkgStack = new THStack("hCSV1_BkgStack","Stacked Background First Jet ordered in CSV");
	hCSV1_ZZ_Zbb->GetXaxis()->SetTitle("CSV");
	hCSV1_DY_Zbb->SetFillColor(kRed);
    hCSV1_BkgStack->Add(hCSV1_DY_Zbb);
	hCSV1_ZZ_Zbb->SetFillColor(kGreen);
    hCSV1_BkgStack->Add(hCSV1_ZZ_Zbb);
	hCSV1_BkgStack->SetTitle("Jet with largest CSV");
	hCSV1_BkgStack->Draw("hist");
	TLegend myLegend(0.7, 0.5, 0.89, 0.8);
	myLegend.SetTextSize(0.03);
	myLegend.AddEntry(hCSV1_DY_Zbb, "DY", "f");	
	myLegend.AddEntry(hCSV1_ZZ_Zbb, "ZZ", "f");	
	myLegend.Draw();		
	plot = directory+"CSV1"+suffixps;
	c1->Print(plot);
	c1->Clear();	
	
	THStack *hdphiVH_BkgStack = new THStack("hdphiVH_BkgStack","Stacked Background Delta Phi of Vector Boson and Higgs");
	hdphiVH_ZZ_Zbb->GetXaxis()->SetTitle("#Delta#phi");
	hdphiVH_ZZ_Zbb->SetFillColor(kGreen);
    hdphiVH_BkgStack->Add(hdphiVH_ZZ_Zbb);
	hdphiVH_DY_Zbb->SetFillColor(kRed);
    hdphiVH_BkgStack->Add(hdphiVH_DY_Zbb);
	hdphiVH_BkgStack->SetTitle("Delta Phi of Vector Boson and Higgs");
	hdphiVH_BkgStack->Draw("hist");
	TLegend myLegend(0.3, 0.6, 0.49, 0.9);
	myLegend.SetTextSize(0.03);
	myLegend.AddEntry(hdphiVH_DY_Zbb, "DY", "f");	
	myLegend.AddEntry(hdphiVH_ZZ_Zbb, "ZZ", "f");	
	myLegend.Draw();		
	plot = directory+"dphiVH"+suffixps;
	c1->Print(plot);
	c1->Clear();

	THStack *hdEtaJJ_BkgStack = new THStack("hdEtaJJ_BkgStack","Stacked Background #Delta#eta of the dijet pair");
	hdetaJJ_ZZ_Zbb->GetXaxis()->SetTitle("#Delta#eta");
	hdetaJJ_ZZ_Zbb->SetFillColor(kGreen);
    hdEtaJJ_BkgStack->Add(hdetaJJ_ZZ_Zbb);
	hdetaJJ_DY_Zbb->SetFillColor(kRed);
    hdEtaJJ_BkgStack->Add(hdetaJJ_DY_Zbb);
	hdEtaJJ_BkgStack->SetTitle("#Delta#eta of the dijet pair");
	hdEtaJJ_BkgStack->Draw("hist");
	TLegend myLegend(0.15, 0.7, 0.28, 0.93);
	myLegend.SetTextSize(0.03);
	myLegend.AddEntry(hdetaJJ_DY_Zbb, "DY", "f");	
	myLegend.AddEntry(hdetaJJ_ZZ_Zbb, "ZZ", "f");	
	myLegend.Draw();		
	plot = directory+"DeltaEtaJJ"+suffixps;
	c1->Print(plot);
	c1->Clear();	
	
	THStack *hNaj_BkgStack = new THStack("hNaj_BkgStack","Stacked Background Additional Jets");
	hNaj_ZZ_Zbb->GetXaxis()->SetTitle("Njets");
	hNaj_ZZ_Zbb->SetFillColor(kGreen);
    hNaj_BkgStack->Add(hNaj_ZZ_Zbb);
	hNaj_DY_Zbb->SetFillColor(kRed);
    hNaj_BkgStack->Add(hNaj_DY_Zbb);
	hNaj_BkgStack->SetTitle("Number of Jets in Addition to two b jets");
	hNaj_BkgStack->Draw("hist");
	TLegend myLegend(0.7, 0.5, 0.89, 0.8);
	myLegend.SetTextSize(0.03);
	myLegend.AddEntry(hNaj_DY_Zbb, "DY", "f");	
	myLegend.AddEntry(hNaj_ZZ_Zbb, "ZZ", "f");	
	myLegend.Draw();		
	plot = directory+"Naj"+suffixps;
	c1->Print(plot);
	c1->Clear();	
	
	THStack *hMmumu_BkgStack = new THStack("hMmumu_BkgStack","Stacked Background Mass of dimuon pair");
	hMmumu_ZZ_Zbb->GetXaxis()->SetTitle("Mass (GeV)");
	hMmumu_ZZ_Zbb->SetFillColor(kGreen);
    hMmumu_BkgStack->Add(hMmumu_ZZ_Zbb);
	hMmumu_DY_Zbb->SetFillColor(kRed);
    hMmumu_BkgStack->Add(hMmumu_DY_Zbb);
	hMmumu_BkgStack->SetTitle("Mass of dimuon pair");
	hMmumu_BkgStack->Draw("hist");
	TLegend myLegend(0.7, 0.5, 0.89, 0.8);
	myLegend.SetTextSize(0.03);
	myLegend.AddEntry(hMmumu_DY_Zbb, "DY", "f");	
	myLegend.AddEntry(hMmumu_ZZ_Zbb, "ZZ", "f");	
	myLegend.Draw();		
	plot = directory+"Mmumu"+suffixps;
	c1->Print(plot);
	c1->Clear();	
	
	THStack *hMjj_BkgStack = new THStack("hMjj_BkgStack","Stacked Background Mass of dijet pair");
	hMjj_ZZ_Zbb->GetXaxis()->SetTitle("Mass (GeV)");
	hMjj_ZZ_Zbb->SetFillColor(kGreen);
    hMjj_BkgStack->Add(hMjj_ZZ_Zbb);
	hMjj_DY_Zbb->SetFillColor(kRed);
    hMjj_BkgStack->Add(hMjj_DY_Zbb);
	hMjj_BkgStack->SetTitle("Mass of dijet pair");
	hMjj_BkgStack->Draw("hist");
	TLegend myLegend(0.7, 0.5, 0.89, 0.8);
	myLegend.SetTextSize(0.03);
	myLegend.AddEntry(hMjj_DY_Zbb, "DY", "f");	
	myLegend.AddEntry(hMjj_ZZ_Zbb, "ZZ", "f");	
	myLegend.Draw();		
	plot = directory+"Mjj"+suffixps;
	c1->Print(plot);
	c1->Clear();	

	THStack *hPtBalZH_BkgStack = new THStack("hPtBalZH_BkgStack","Stacked Background Pt Balance ZH");
	hPtbalZH_ZZ_Zbb->GetXaxis()->SetTitle("#Delta p_{T} (GeV)");
	hPtbalZH_ZZ_Zbb->SetFillColor(kGreen);
    hPtBalZH_BkgStack->Add(hPtbalZH_ZZ_Zbb);
	hPtbalZH_DY_Zbb->SetFillColor(kRed);
    hPtBalZH_BkgStack->Add(hPtbalZH_DY_Zbb);
	hPtBalZH_BkgStack->SetTitle("Pt balance of Z and H");
	hPtBalZH_BkgStack->Draw("hist");
	TLegend myLegend(0.7, 0.5, 0.89, 0.8);
	myLegend.SetTextSize(0.03);
	myLegend.AddEntry(hPtbalZH_DY_Zbb, "DY", "f");	
	myLegend.AddEntry(hPtbalZH_ZZ_Zbb, "ZZ", "f");	
	myLegend.Draw();		
	plot = directory+"PtBalZH"+suffixps;
	c1->Print(plot);
	c1->Clear();	
	
	THStack *hdphiJJ_BkgStack = new THStack("hdphiJJ_BkgStack","Stacked Background #Delta#phi JJ");
	hdphiJJ_ZZ_Zbb->GetXaxis()->SetTitle("#Delta#phi");
	hdphiJJ_ZZ_Zbb->SetFillColor(kGreen);
    hdphiJJ_BkgStack->Add(hdphiJJ_ZZ_Zbb);
	hdphiJJ_DY_Zbb->SetFillColor(kRed);
    hdphiJJ_BkgStack->Add(hdphiJJ_DY_Zbb);
	hdphiJJ_BkgStack->SetTitle("#Delta#phi JJ");
	hdphiJJ_BkgStack->Draw("hist");
	TLegend myLegend(0.7, 0.5, 0.89, 0.8);
	myLegend.SetTextSize(0.03);
	myLegend.AddEntry(hdphiJJ_DY_Zbb, "DY", "f");	
	myLegend.AddEntry(hdphiJJ_ZZ_Zbb, "ZZ", "f");	
	myLegend.Draw();		
	plot = directory+"dphiJJ"+suffixps;
	c1->Print(plot);
	c1->Clear();
	
	THStack *hMjj_aftercuts_BkgStack = new THStack("hMjj_aftercuts_BkgStack","Stacked Background JJ Mass after Cuts");
	hMjj_aftercuts_ZZ_Zbb->GetXaxis()->SetTitle("Mass (GeV)");
	hMjj_aftercuts_ZZ_Zbb->SetFillColor(kGreen);
    hMjj_aftercuts_BkgStack->Add(hMjj_aftercuts_ZZ_Zbb);
	hMjj_aftercuts_DY_Zbb->SetFillColor(kRed);
    hMjj_aftercuts_BkgStack->Add(hMjj_aftercuts_DY_Zbb);
	hMjj_aftercuts_BkgStack->SetTitle("JJ Mass After Cuts");
	hMjj_aftercuts_BkgStack->Draw("hist");
	TLegend myLegend(0.7, 0.5, 0.89, 0.8);
	myLegend.SetTextSize(0.03);
	myLegend.AddEntry(hMjj_aftercuts_DY_Zbb, "DY", "f");	
	myLegend.AddEntry(hMjj_aftercuts_ZZ_Zbb, "ZZ", "f");	
	myLegend.Draw();		
	plot = directory+"Mjj_aftercuts"+suffixps;
	c1->Print(plot);
	c1->Clear();
	
	THStack *hDphiDetajj_BkgStack = new THStack("hDphiDetajj_BkgStack","Stacked Background #Delta#eta vs #Delta#phi of jets");
	hDphiDetajj_ZZ_Zbb->GetXaxis()->SetTitle("#Delta#eta");
	hDphiDetajj_ZZ_Zbb->GetYaxis()->SetTitle("#Delta#phi");
	hDphiDetajj_ZZ_Zbb->SetMarkerColor(kGreen);
	hDphiDetajj_ZZ_Zbb->SetMarkerStyle(20);
	hDphiDetajj_ZZ_Zbb->SetMarkerSize(.1);
    hDphiDetajj_BkgStack->Add(hDphiDetajj_ZZ_Zbb);
	hDphiDetajj_DY_Zbb->SetMarkerColor(kRed);
	hDphiDetajj_DY_Zbb->SetMarkerStyle(20);
	hDphiDetajj_DY_Zbb->SetMarkerSize(.4);
    hDphiDetajj_BkgStack->Add(hDphiDetajj_DY_Zbb);
	hDphiDetajj_BkgStack->SetTitle(" #Delta#eta vs #Delta#phi of jets");
	hDphiDetajj_BkgStack->Draw("scat");
	TLegend myLegend(0.7, 0.5, 0.89, 0.8);
	myLegend.SetTextSize(0.03);
	myLegend.AddEntry(hDphiDetajj_DY_Zbb, "DY", "p");	
	myLegend.AddEntry(hDphiDetajj_ZZ_Zbb, "ZZ", "p");	
	myLegend.Draw();		
	plot = directory+"DphiDetajj"+suffixps;
	c1->Print(plot);
	c1->Clear();	
	hDphiDetajj_ZZ_Zbb->Add(hDphiDetajj_DY_Zbb);
	gStyle->SetPalette(1);
	hDphiDetajj_ZZ_Zbb->Draw("colz");
	plot = directory+"DphiDetajj_Colz"+suffixps;
	c1->Print(plot);
	c1->Clear();	
	
	
	
	hCutFlow_ZZ_Zbb->GetXaxis()->SetBinLabel(1,"HLT");
	hCutFlow_ZZ_Zbb->GetXaxis()->SetBinLabel(2,"PreSelection");
	hCutFlow_ZZ_Zbb->GetXaxis()->SetBinLabel(3,"p_{T}(jj)");
	hCutFlow_ZZ_Zbb->GetXaxis()->SetBinLabel(4,"p_{T}(Z)");
	hCutFlow_ZZ_Zbb->GetXaxis()->SetBinLabel(5,"CSV1");
	hCutFlow_ZZ_Zbb->GetXaxis()->SetBinLabel(6,"CSV2");
	hCutFlow_ZZ_Zbb->GetXaxis()->SetBinLabel(7,"#Delta#phi");
	hCutFlow_ZZ_Zbb->GetXaxis()->SetBinLabel(8,"N_{aj}");
	hCutFlow_ZZ_Zbb->GetXaxis()->SetBinLabel(9,"M_{jj}");		
	THStack *hMjj_BkgStack = new THStack("hMjj_BkgStack","Stacked Background Cut Flow");
	hCutFlow_ZZ_Zbb->GetYaxis()->SetTitle("Number of Events");
	hCutFlow_ZZ_Zbb->GetXaxis()->SetRangeUser(1,9);
	hCutFlow_DY_Zbb->GetXaxis()->SetRangeUser(1,9);
	hCutFlow_ZZ_Zbb->SetFillColor(kGreen);
    hMjj_BkgStack->Add(hCutFlow_ZZ_Zbb);
	hCutFlow_DY_Zbb->SetFillColor(kRed);
    hMjj_BkgStack->Add(hCutFlow_DY_Zbb);
	hMjj_BkgStack->SetTitle("Cut Flow HZ");
	hMjj_BkgStack->Draw("hist");
	TLegend myLegend(0.7, 0.5, 0.89, 0.8);
	myLegend.SetTextSize(0.03);
	myLegend.AddEntry(hCutFlow_DY_Zbb, "DY", "f");	
	myLegend.AddEntry(hCutFlow_ZZ_Zbb, "ZZ", "f");	
	myLegend.Draw();		
	plot = directory+"CutFlow"+suffixps;
	c1->Print(plot);
	c1->Clear();	

	theHistogramFile->Write();
	
//	hHighestPt_TTTo2L2Nu2B->Delete();
		
}
