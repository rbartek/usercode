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


void minitupleMacro(){

  /* Macro to take first look at Zbb minitupes
   * hbb variable produced by Zbb_minituples.cc
   * Author:  Rachel Wilken - UNL
   */

    TFile *ggH115_Zbb = TFile::Open("/Users/rachelwilken/Documents/CMS/Zbb/RecreateAnalysis/CutFlowHisto.root");
	//TFile *ggH115_Zbb = TFile::Open("/afs/cern.ch/user/w/wilken/scratch0/CMSSW_4_2_5/src/UserCode/wilken/output.root");


TString suffixps = ".gif";
TString directory = "MacroHistos/";
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

    TH1F* hjetCSV3_ggH115_Zbb		= (TH1F*) ggH115_Zbb->Get("jetCSV3");
    TH1F* hjetPt3_ggH115_Zbb		= (TH1F*) ggH115_Zbb->Get("jetPt3");
    TH1F* hnJets_ggH115_Zbb			= (TH1F*) ggH115_Zbb->Get("hnJets");
    TH1F* hnMuons_ggH115_Zbb		= (TH1F*) ggH115_Zbb->Get("hnMuons");
    TH1F* hjetPt1_ggH115_Zbb		= (TH1F*) ggH115_Zbb->Get("jetPt1");
    TH1F* hjetPt2_ggH115_Zbb		= (TH1F*) ggH115_Zbb->Get("jetPt2");
    TH1F* hPtb1_ggH115_Zbb			= (TH1F*) ggH115_Zbb->Get("hPtb1");
    TH1F* hPtb2_ggH115_Zbb			= (TH1F*) ggH115_Zbb->Get("hPtb2");
    TH1F* hPtjj_ggH115_Zbb			= (TH1F*) ggH115_Zbb->Get("hPtjj");
    TH1F* hPtmumu_ggH115_Zbb		= (TH1F*) ggH115_Zbb->Get("hPtmumu");
    TH1F* hPtmu1_ggH115_Zbb			= (TH1F*) ggH115_Zbb->Get("hPtmu1");
    TH1F* hPtmu2_ggH115_Zbb			= (TH1F*) ggH115_Zbb->Get("hPtmu2");
    TH1F* hCSV1_ggH115_Zbb			= (TH1F*) ggH115_Zbb->Get("hCSV1");
    TH1F* hCSV2_ggH115_Zbb			= (TH1F*) ggH115_Zbb->Get("hCSV2");
    TH1F* hdphiVH_ggH115_Zbb		= (TH1F*) ggH115_Zbb->Get("hdphiVH");
    TH1F* hdetaJJ_ggH115_Zbb		= (TH1F*) ggH115_Zbb->Get("hdetaJJ");
    TH1F* hNaj_ggH115_Zbb			= (TH1F*) ggH115_Zbb->Get("hNaj");
    TH1F* hMmumu_ggH115_Zbb			= (TH1F*) ggH115_Zbb->Get("hMmumu");
    TH1F* hMjj_ggH115_Zbb			= (TH1F*) ggH115_Zbb->Get("hMjj");
    TH1F* hCutFlow_ggH115_Zbb		= (TH1F*) ggH115_Zbb->Get("hCutFlow");

	// Histogram
	hjetCSV3_ggH115_Zbb->SetTitle("Third Jet ordered in CSV");
	hjetCSV3_ggH115_Zbb->GetXaxis()->SetRange(0,1);
	hjetCSV3_ggH115_Zbb->GetXaxis()->SetTitle("CSV");
	hjetCSV3_ggH115_Zbb->GetYaxis()->SetTitle(" ");
	hjetCSV3_ggH115_Zbb->SetLineColor(kBlue);
	hjetCSV3_ggH115_Zbb->Draw("hist");
	//gStyle->SetOptStat(kFALSE);
	plot = directory+"CSV3"+suffixps;
	c1->Print(plot);
	c1->Clear();
	
	hjetPt3_ggH115_Zbb->SetTitle("Jet with third highest Pt");
	hjetPt3_ggH115_Zbb->GetXaxis()->SetTitle("p_{T} (GeV)");
	hjetPt3_ggH115_Zbb->GetYaxis()->SetTitle(" ");
	hjetPt3_ggH115_Zbb->SetLineColor(kBlue);
	hjetPt3_ggH115_Zbb->Draw("hist");
	//gStyle->SetOptStat(kFALSE);
	plot = directory+"JetPt3"+suffixps;
	c1->Print(plot);
	c1->Clear();

	hnJets_ggH115_Zbb->SetTitle("Number of Good Jets");
	hnJets_ggH115_Zbb->GetXaxis()->SetTitle("Njets");
	hnJets_ggH115_Zbb->GetYaxis()->SetTitle(" ");
	hnJets_ggH115_Zbb->SetLineColor(kBlue);
	hnJets_ggH115_Zbb->Draw("hist");
	//gStyle->SetOptStat(kFALSE);
	plot = directory+"nJets"+suffixps;
	c1->Print(plot);
	c1->Clear();

	hnMuons_ggH115_Zbb->SetTitle("Number of Good Muonns");
	hnMuons_ggH115_Zbb->GetXaxis()->SetTitle("Nmuons");
	hnMuons_ggH115_Zbb->GetYaxis()->SetTitle(" ");
	hnMuons_ggH115_Zbb->SetLineColor(kBlue);
	hnMuons_ggH115_Zbb->Draw("hist");
	//gStyle->SetOptStat(kFALSE);
	plot = directory+"Nmuons"+suffixps;
	c1->Print(plot);
	c1->Clear();

	hjetPt2_ggH115_Zbb->SetTitle("Jet with second highest Pt");
	hjetPt2_ggH115_Zbb->GetXaxis()->SetTitle("p_{T} (GeV)");
	hjetPt2_ggH115_Zbb->GetYaxis()->SetTitle(" ");
	hjetPt2_ggH115_Zbb->SetLineColor(kBlue);
	hjetPt2_ggH115_Zbb->Draw("hist");
	//gStyle->SetOptStat(kFALSE);
	plot = directory+"JetPt2"+suffixps;
	c1->Print(plot);
	c1->Clear();

	hjetPt1_ggH115_Zbb->SetTitle("Jet with highest Pt");
	hjetPt1_ggH115_Zbb->GetXaxis()->SetTitle("p_{T} (GeV)");
	hjetPt1_ggH115_Zbb->GetYaxis()->SetTitle(" ");
	hjetPt1_ggH115_Zbb->SetLineColor(kBlue);
	hjetPt1_ggH115_Zbb->Draw("hist");
	//gStyle->SetOptStat(kFALSE);
	plot = directory+"JetPt1"+suffixps;
	c1->Print(plot);
	c1->Clear();
	
	hPtb1_ggH115_Zbb->SetTitle("Pt of Jet with highest CSV");
	hPtb1_ggH115_Zbb->GetXaxis()->SetTitle("p_{T} (GeV)");
	hPtb1_ggH115_Zbb->GetYaxis()->SetTitle(" ");
	hPtb1_ggH115_Zbb->SetLineColor(kBlue);
	hPtb1_ggH115_Zbb->Draw("hist");
	//gStyle->SetOptStat(kFALSE);
	plot = directory+"Ptb1"+suffixps;
	c1->Print(plot);
	c1->Clear();
	
	hPtb2_ggH115_Zbb->SetTitle("Pt of Jet with second highest CSV");
	hPtb2_ggH115_Zbb->GetXaxis()->SetTitle("p_{T} (GeV)");
	hPtb2_ggH115_Zbb->GetYaxis()->SetTitle(" ");
	hPtb2_ggH115_Zbb->SetLineColor(kBlue);
	hPtb2_ggH115_Zbb->Draw("hist");
	//gStyle->SetOptStat(kFALSE);
	plot = directory+"Ptb2"+suffixps;
	c1->Print(plot);
	c1->Clear();	

	hPtjj_ggH115_Zbb->SetTitle("Pt of dijet pair");
	hPtjj_ggH115_Zbb->GetXaxis()->SetTitle("p_{T} (GeV)");
	hPtjj_ggH115_Zbb->GetYaxis()->SetTitle(" ");
	hPtjj_ggH115_Zbb->SetLineColor(kBlue);
	hPtjj_ggH115_Zbb->Draw("hist");
	//gStyle->SetOptStat(kFALSE);
	plot = directory+"Ptjj"+suffixps;
	c1->Print(plot);
	c1->Clear();

	hPtmumu_ggH115_Zbb->SetTitle("Pt of dimuon pair");
	hPtmumu_ggH115_Zbb->GetXaxis()->SetTitle("p_{T} (GeV)");
	hPtmumu_ggH115_Zbb->GetYaxis()->SetTitle(" ");
	hPtmumu_ggH115_Zbb->SetLineColor(kBlue);
	hPtmumu_ggH115_Zbb->Draw("hist");
	//gStyle->SetOptStat(kFALSE);
	plot = directory+"Ptmumu"+suffixps;
	c1->Print(plot);
	c1->Clear();

	hPtmu1_ggH115_Zbb->SetTitle("Muon with highest Pt");
	hPtmu1_ggH115_Zbb->GetXaxis()->SetTitle("p_{T} (GeV)");
	hPtmu1_ggH115_Zbb->GetYaxis()->SetTitle(" ");
	hPtmu1_ggH115_Zbb->SetLineColor(kBlue);
	hPtmu1_ggH115_Zbb->Draw("hist");
	//gStyle->SetOptStat(kFALSE);
	plot = directory+"Ptmu1"+suffixps;
	c1->Print(plot);
	c1->Clear();

	hPtmu2_ggH115_Zbb->SetTitle("Muon with second highest P");
	hPtmu2_ggH115_Zbb->GetXaxis()->SetTitle("p_{T} (GeV)");
	hPtmu2_ggH115_Zbb->GetYaxis()->SetTitle(" ");
	hPtmu2_ggH115_Zbb->SetLineColor(kBlue);
	hPtmu2_ggH115_Zbb->Draw("hist");
	//gStyle->SetOptStat(kFALSE);
	plot = directory+"Ptmu2"+suffixps;
	c1->Print(plot);
	c1->Clear();

	hCSV1_ggH115_Zbb->SetTitle("First Jet ordered in CSV");
	hCSV1_ggH115_Zbb->GetXaxis()->SetRange(0,1);
	hCSV1_ggH115_Zbb->GetXaxis()->SetTitle("CSV");
	hCSV1_ggH115_Zbb->GetYaxis()->SetTitle(" ");
	hCSV1_ggH115_Zbb->SetLineColor(kBlue);
	hCSV1_ggH115_Zbb->Draw("hist");
	//gStyle->SetOptStat(kFALSE);
	plot = directory+"CSV1"+suffixps;
	c1->Print(plot);
	c1->Clear();

	hCSV2_ggH115_Zbb->SetTitle("Second Jet ordered in CSV");
	hCSV2_ggH115_Zbb->GetXaxis()->SetRange(0,1);
	hCSV2_ggH115_Zbb->GetXaxis()->SetTitle("CSV");
	hCSV2_ggH115_Zbb->GetYaxis()->SetTitle(" ");
	hCSV2_ggH115_Zbb->SetLineColor(kBlue);
	hCSV2_ggH115_Zbb->Draw("hist");
	//gStyle->SetOptStat(kFALSE);
	plot = directory+"CSV2"+suffixps;
	c1->Print(plot);
	c1->Clear();
	

	hdphiVH_ggH115_Zbb->SetTitle("Delta Phi of Vector Boson and Higgs");
	hdphiVH_ggH115_Zbb->GetXaxis()->SetTitle("#Delta#phi");
	hdphiVH_ggH115_Zbb->GetYaxis()->SetTitle(" ");
	hdphiVH_ggH115_Zbb->SetLineColor(kBlue);
	hdphiVH_ggH115_Zbb->Draw("hist");
	//gStyle->SetOptStat(kFALSE);
	plot = directory+"dphiVH"+suffixps;
	c1->Print(plot);
	c1->Clear();

	hdetaJJ_ggH115_Zbb->SetTitle("Delta Eta of dijet pair");
	hdetaJJ_ggH115_Zbb->GetXaxis()->SetTitle("#Delta#eta");
	hdetaJJ_ggH115_Zbb->GetYaxis()->SetTitle(" ");
	hdetaJJ_ggH115_Zbb->SetLineColor(kBlue);
	hdetaJJ_ggH115_Zbb->Draw("hist");
	//gStyle->SetOptStat(kFALSE);
	plot = directory+"DeltaEtaJJ"+suffixps;
	c1->Print(plot);
	c1->Clear();

	hNaj_ggH115_Zbb->SetTitle("Number of Jets in Addition to two b jets");
	hNaj_ggH115_Zbb->GetXaxis()->SetTitle("Njets");
	hNaj_ggH115_Zbb->GetYaxis()->SetTitle(" ");
	hNaj_ggH115_Zbb->SetLineColor(kBlue);
	hNaj_ggH115_Zbb->Draw("hist");
	//gStyle->SetOptStat(kFALSE);
	plot = directory+"Naj"+suffixps;
	c1->Print(plot);
	c1->Clear();

	hMmumu_ggH115_Zbb->SetTitle("Mass of dimuont pair");
	hMmumu_ggH115_Zbb->GetXaxis()->SetRange(30,130);
	hMmumu_ggH115_Zbb->GetXaxis()->SetTitle("Mass (GeV)");
	hMmumu_ggH115_Zbb->GetYaxis()->SetTitle(" ");
	hMmumu_ggH115_Zbb->SetLineColor(kBlue);
	hMmumu_ggH115_Zbb->Draw("hist");
	//gStyle->SetOptStat(kFALSE);
	plot = directory+"Mmumu"+suffixps;
	c1->Print(plot);
	c1->Clear();

	hMjj_ggH115_Zbb->SetTitle("Mass of dijet pair");
	hMjj_ggH115_Zbb->GetXaxis()->SetTitle("Mass (GeV)");
	hMjj_ggH115_Zbb->GetYaxis()->SetTitle(" ");
	hMjj_ggH115_Zbb->SetLineColor(kBlue);
	hMjj_ggH115_Zbb->Draw("hist");
	//gStyle->SetOptStat(kFALSE);
	plot = directory+"Mjj"+suffixps;
	c1->Print(plot);
	c1->Clear();
	
	 hCutFlow_ggH115_Zbb->GetXaxis()->SetBinLabel(1,"HLT");
	 hCutFlow_ggH115_Zbb->GetXaxis()->SetBinLabel(2,"PreSelection");
	 hCutFlow_ggH115_Zbb->GetXaxis()->SetBinLabel(3,"p_{T}(jj)");
	 hCutFlow_ggH115_Zbb->GetXaxis()->SetBinLabel(4,"p_{T}(Z)");
	 hCutFlow_ggH115_Zbb->GetXaxis()->SetBinLabel(5,"CSV1");
	 hCutFlow_ggH115_Zbb->GetXaxis()->SetBinLabel(6,"CSV2");
	 hCutFlow_ggH115_Zbb->GetXaxis()->SetBinLabel(7,"#Delta#phi");
	 hCutFlow_ggH115_Zbb->GetXaxis()->SetBinLabel(8,"N_{aj}");
	 hCutFlow_ggH115_Zbb->GetXaxis()->SetBinLabel(9,"M_{jj}");	
	hCutFlow_ggH115_Zbb->SetTitle("Cut Flow HZ");
	hCutFlow_ggH115_Zbb->GetYaxis()->SetTitle("Number of Events");
	hCutFlow_ggH115_Zbb->SetLineColor(kBlack);
	hCutFlow_ggH115_Zbb->Draw("hist");
	gStyle->SetOptStat(kFALSE);
	plot = directory+"CutFlow"+suffixps;
	c1->Print(plot);
	c1->Clear();
	
	
	
}
