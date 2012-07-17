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

#include "/home/hep/wilken/taus/CMSSW_4_4_2_patch8/src/UserCode/wilken/interface/xsecV21.h"

void CutFlow(){

  /* Macro to take first look at Zbb minitupes
   * hbb variable produced by Zbb_minituples.cc
   * Author:  Rachel Wilken - UNL
   */
   
//    TFile *ZH_M115_Zbb = TFile::Open("/home/hep/wilken/taus/CMSSW_4_4_2_patch8/src/UserCode/wilken/V21/ZH115_Fall.root");
    TFile *ZH_M115_Zbb = TFile::Open("/home/hep/wilken/taus/CMSSW_4_4_2_patch8/src/UserCode/wilken/V21/Filteredemu.root");
    TFile *DYJetsToLL_M50_Zbb = TFile::Open("/home/hep/wilken/taus/CMSSW_4_4_2_patch8/src/UserCode/wilken/V21/DY_M50.root");
    TFile *TTJets_Zbb = TFile::Open("/home/hep/wilken/taus/CMSSW_4_4_2_patch8/src/UserCode/wilken/V21/TTJets.root");
	TFile *WJets_Zbb = TFile::Open("/home/hep/wilken/taus/CMSSW_4_4_2_patch8/src/UserCode/wilken/V21/WJets.root");
	
	TFile *DYJetsToLL_PtZ_Zbb = TFile::Open("/home/hep/wilken/taus/CMSSW_4_4_2_patch8/src/UserCode/wilken/V21/DY_PtZ.root");
	

	TString suffixps = ".gif";
	TString directory = "StackedPlots_allEvts/";
	TString plot = directory+"somereallylongnameherehopethisarrayisbiggenough"+suffixps;
	// Files for histogram output --> set suffixps to desired file type:  e.g. .eps, .jpg, ...
	
	
	TCanvas *c1 = new TCanvas("c1","");
	//gStyle->SetOptStat("kTRUE");
	c1->SetFillColor(10);
	c1->SetFillColor(10);
	double normalization;
	
// ********************************************************************
// Read in Histograms
// ********************************************************************

     TH1F* hCutFlow_H115_Zbb		= (TH1F*) ZH_M115_Zbb->Get("hCutFlow");
    TH1F* hCutFlow_DYJetsToLL_M50_Zbb		= (TH1F*) DYJetsToLL_M50_Zbb->Get("hCutFlow");
    TH1F* hCutFlow_TTJets_Zbb		= (TH1F*) TTJets_Zbb->Get("hCutFlow");
    TH1F* hCutFlow_WJets_Zbb		= (TH1F*) WJets_Zbb->Get("hCutFlow");
	TH1F* hCutFlow_DYJetsToLL_PtZ_Zbb		= (TH1F*)DYJetsToLL_PtZ_Zbb->Get("hCutFlow");
	
	double lumi = 4.457;
 
	double H115_weight = lumi/lumiZH115;
//	double H115_weight = lumi/(2.28079234375000000e+05/xsecbfZH115);

    Double_t WJets_weight = lumi/(lumiWJ);
	Double_t TTJets_weight = lumi/(lumiTT);
	Double_t DYJetsToLL_M50_weight = lumi/(lumiZJL);
	Double_t DYJetsToLL_PtZ_weight = lumi/(lumiZJH);
	
	hCutFlow_DYJetsToLL_M50_Zbb		->Scale(DYJetsToLL_M50_weight);
	hCutFlow_DYJetsToLL_PtZ_Zbb		->Scale(DYJetsToLL_PtZ_weight);
	hCutFlow_TTJets_Zbb		->Scale(TTJets_weight);
	hCutFlow_WJets_Zbb		->Scale(WJets_weight);
	hCutFlow_H115_Zbb->Scale(H115_weight);
	
 
//Histograms
/*	
	hCutFlow_H115_Zbb->GetXaxis()->SetBinLabel(1,"Candidate");	
	 hCutFlow_H115_Zbb->GetXaxis()->SetBinLabel(2,"HLT");
	 hCutFlow_H115_Zbb->GetXaxis()->SetBinLabel(3,"PreSelection");
	 hCutFlow_H115_Zbb->GetXaxis()->SetBinLabel(4,"M(e#mu)");
	 hCutFlow_H115_Zbb->GetXaxis()->SetBinLabel(5,"#Delta#phi(Z,MET)");
	 hCutFlow_H115_Zbb->GetXaxis()->SetBinLabel(6,"#Delta#phi(H,Z)");
	 hCutFlow_H115_Zbb->GetXaxis()->SetBinLabel(7,"CSV0");
	hCutFlow_H115_Zbb->GetXaxis()->SetBinLabel(8,"CHF(b0)");
	hCutFlow_H115_Zbb->GetXaxis()->SetBinLabel(9,"#DeltaR(e,#mu)");
	hCutFlow_H115_Zbb->GetXaxis()->SetBinLabel(10,"ScalarSumPt");
	 hCutFlow_H115_Zbb->GetXaxis()->SetBinLabel(11,"M_{jj}");	
	hCutFlow_H115_Zbb->SetTitle("Cut Flow emu");
	hCutFlow_H115_Zbb->GetYaxis()->SetTitle("Number of Events");
	hCutFlow_H115_Zbb->SetLineColor(kBlack);
	hCutFlow_H115_Zbb->Draw("hist");
	gStyle->SetOptStat(kFALSE);
	plot = directory+"CutFlowSig"+suffixps;
	c1->Print(plot);
	c1->Clear();
	*/
	//gStyle->SetOptStat("kFALSE");
	hCutFlow_TTJets_Zbb->GetYaxis()->SetTitle("Number of Events");
	hCutFlow_TTJets_Zbb->SetTitle("Cut Flow ZH->emubb");
	hCutFlow_TTJets_Zbb->GetXaxis()->SetBinLabel(1,"Candidate");	
	hCutFlow_TTJets_Zbb->GetXaxis()->SetBinLabel(2,"HLT");
	hCutFlow_TTJets_Zbb->GetXaxis()->SetBinLabel(3,"PreSelection");
	hCutFlow_TTJets_Zbb->GetXaxis()->SetBinLabel(4,"EleFake");
	hCutFlow_TTJets_Zbb->GetXaxis()->SetBinLabel(5,"p_{T}(H)");
	hCutFlow_TTJets_Zbb->GetXaxis()->SetBinLabel(6,"M(e#mu)");
	hCutFlow_TTJets_Zbb->GetXaxis()->SetBinLabel(7,"CSV0");
	hCutFlow_TTJets_Zbb->GetXaxis()->SetBinLabel(8,"#Delta#phi(Z,MET)");
	hCutFlow_TTJets_Zbb->GetXaxis()->SetBinLabel(9,"M_{jj}");
	hCutFlow_TTJets_Zbb->GetXaxis()->SetBinLabel(10,"CHFb0");		
	hCutFlow_TTJets_Zbb->GetXaxis()->SetBinLabel(11,"Naj");		
	gStyle->SetOptStat(kFALSE);
	hCutFlow_TTJets_Zbb->SetMarkerColor(kBlue);
	hCutFlow_TTJets_Zbb->SetMarkerStyle(kOpenCross);
	hCutFlow_TTJets_Zbb->SetMarkerSize(2);
	hCutFlow_TTJets_Zbb->Draw("p");	
	hCutFlow_WJets_Zbb->SetMarkerColor(kGray+2);
	hCutFlow_WJets_Zbb->SetMarkerStyle(kOpenSquare);
	hCutFlow_WJets_Zbb->SetMarkerSize(2);
	hCutFlow_WJets_Zbb->Draw("same p");
	hCutFlow_DYJetsToLL_PtZ_Zbb->SetMarkerColor(kMagenta+2);
	hCutFlow_DYJetsToLL_PtZ_Zbb->SetMarkerStyle(kOpenCircle);
	hCutFlow_DYJetsToLL_PtZ_Zbb->SetMarkerSize(2);
	hCutFlow_DYJetsToLL_PtZ_Zbb->Draw("same p");
	hCutFlow_DYJetsToLL_M50_Zbb->SetMarkerColor(kRed);
	hCutFlow_DYJetsToLL_M50_Zbb->SetMarkerStyle(23); //kFullTriangleDown
	hCutFlow_DYJetsToLL_M50_Zbb->SetMarkerSize(2);
	hCutFlow_DYJetsToLL_M50_Zbb->Draw("same p");
	TLegend myLegend(0.65, 0.55, 0.89, 0.89);
	myLegend.SetTextSize(0.03);
	myLegend.AddEntry(hCutFlow_DYJetsToLL_M50_Zbb, "DYJetsToLL_M50", "p");	
	myLegend.AddEntry(hCutFlow_TTJets_Zbb, "TTJets", "p");	
	myLegend.AddEntry(hCutFlow_DYJetsToLL_PtZ_Zbb, "DYJetsToLL_PtZ", "p");	
	myLegend.AddEntry(hCutFlow_WJets_Zbb, "WJets", "p");		
	myLegend.Draw();	
	plot = directory+"CutFlowBkg"+suffixps;
	gStyle->SetOptStat(kFALSE);
	c1->Print(plot);
	hCutFlow_H115_Zbb->Scale(10000);
	hCutFlow_H115_Zbb->SetLineColor(kBlack);
	hCutFlow_H115_Zbb->Draw("hist same");
	gStyle->SetOptStat(kFALSE);
	myLegend.AddEntry(hCutFlow_H115_Zbb, "emu M115 x 1000", "l");
	myLegend.Draw();			
	gStyle->SetOptStat(kFALSE);
	plot = directory+"CutFlow"+suffixps;
	c1->Print(plot);
	c1->SetLogy();
	gStyle->SetOptStat(kFALSE);
	plot = directory+"CutFlowlogy_Bkg"+suffixps;	
	c1->Print(plot);
	c1->Clear();
	
	
}
