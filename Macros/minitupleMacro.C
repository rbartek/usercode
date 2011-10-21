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

    TFile *ggH115_Zbb = TFile::Open("/Users/rachelwilken/Documents/CMS/Zbb/RecreateAnalysis/AllCuts.root");


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

    TH1F* hMjj_ggH115_Zbb			= (TH1F*) ggH115_Zbb->Get("hMjj");

	// Histogram
	hMjj_ggH115_Zbb->SetTitle("Mass of dijet pair");
	hMjj_ggH115_Zbb->GetXaxis()->SetTitle("Mass (GeV)");
	hMjj_ggH115_Zbb->GetYaxis()->SetTitle(" ");
	hMjj_ggH115_Zbb->SetLineColor(kRed);
	hMjj_ggH115_Zbb->Draw("hist");
	//gStyle->SetOptStat(kFALSE);
	plot = directory+"Mjj"+suffixps;
	c1->Print(plot);
	c1->Clear();
	
	
	
	
	
}
