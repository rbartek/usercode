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


void BackgroundMacro_commonTuples(){

  /* Macro to take first look at Zbb minitupes
   * hbb variable produced by Zbb_minituples.cc
   * Author:  Rachel Wilken - UNL
   */

    TFile *ZZ_Zbb = TFile::Open("/Users/rachelwilken/Documents/CMS/Zbb/CommonNtuples/Halloween/ZZ.root");
    TFile *DY_Zbb = TFile::Open("/Users/rachelwilken/Documents/CMS/Zbb/CommonNtuples/Halloween/DY.root");
    TFile *ZJets_Zbb = TFile::Open("/Users/rachelwilken/Documents/CMS/Zbb/CommonNtuples/Halloween/ZJets.root");
    TFile *WZ_Zbb = TFile::Open("/Users/rachelwilken/Documents/CMS/Zbb/CommonNtuples/Halloween/WZ.root");
    TFile *TTJets_Zbb = TFile::Open("/Users/rachelwilken/Documents/CMS/Zbb/CommonNtuples/Halloween/TTJets.root");
	//TFile *ZZ_Zbb = TFile::Open("/afs/cern.ch/user/w/wilken/scratch0/CMSSW_4_2_5/src/UserCode/wilken/output.root");


TString suffixps = "_Bkg.gif";
TString directory = "Halloween/";
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

    TH1F* hnJets_ZZ_Zbb			= (TH1F*) ZZ_Zbb->Get("hnJets");
    TH1F* hnMuons_ZZ_Zbb		= (TH1F*) ZZ_Zbb->Get("hnMuons");
    TH1F* hPtb1_ZZ_Zbb			= (TH1F*) ZZ_Zbb->Get("hPtb1");
    TH1F* hPtb2_ZZ_Zbb			= (TH1F*) ZZ_Zbb->Get("hPtb2");
    TH1F* hPtjj_ZZ_Zbb			= (TH1F*) ZZ_Zbb->Get("hPtjj");
    TH1F* hPtmumu_ZZ_Zbb		= (TH1F*) ZZ_Zbb->Get("hPtmumu");
    TH1F* hPtbalZH_ZZ_Zbb		= (TH1F*) ZZ_Zbb->Get("hPtbalZH");
    TH1F* hPtmu1_ZZ_Zbb			= (TH1F*) ZZ_Zbb->Get("hPtmu1");
    TH1F* hEtamu2_ZZ_Zbb		= (TH1F*) ZZ_Zbb->Get("hEtamu1");
    TH1F* hEtamu1_ZZ_Zbb		= (TH1F*) ZZ_Zbb->Get("hEtamu2");
    TH1F* hPtmu2_ZZ_Zbb			= (TH1F*) ZZ_Zbb->Get("hPtmu2");
    TH1F* hCSV1_ZZ_Zbb			= (TH1F*) ZZ_Zbb->Get("hCSV1");
    TH1F* hCSV2_ZZ_Zbb			= (TH1F*) ZZ_Zbb->Get("hCSV2");
    TH1F* hdphiVH_ZZ_Zbb		= (TH1F*) ZZ_Zbb->Get("hdphiVH");
    TH1F* hdetaJJ_ZZ_Zbb		= (TH1F*) ZZ_Zbb->Get("hdetaJJ");
    TH1F* hdphiJJ_ZZ_Zbb		= (TH1F*) ZZ_Zbb->Get("hdphiJJ");
    TH1F* hNaj_ZZ_Zbb			= (TH1F*) ZZ_Zbb->Get("hNaj");
    TH1F* hMmumu_ZZ_Zbb			= (TH1F*) ZZ_Zbb->Get("hMmumu");
    TH1F* hMjj_ZZ_Zbb			= (TH1F*) ZZ_Zbb->Get("hMjj");
	TH1F* hMjj_aftercuts_ZZ_Zbb		= (TH1F*) ZZ_Zbb->Get("hMjj_aftercuts");
    TH1F* hCutFlow_ZZ_Zbb		= (TH1F*) ZZ_Zbb->Get("hCutFlow");
    TH2F* hDphiDetajj_ZZ_Zbb	= (TH2F*) ZZ_Zbb->Get("hDphiDetajj");
	
    TH1F* hnJets_DY_Zbb			= (TH1F*) DY_Zbb->Get("hnJets");
    TH1F* hnMuons_DY_Zbb		= (TH1F*) DY_Zbb->Get("hnMuons");
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
	
    TH1F* hnJets_ZJets_Zbb			= (TH1F*) ZJets_Zbb->Get("hnJets");
    TH1F* hnMuons_ZJets_Zbb		= (TH1F*) ZJets_Zbb->Get("hnMuons");
	TH1F* hPtb1_ZJets_Zbb		= (TH1F*) ZJets_Zbb->Get("hPtb1");
    TH1F* hPtb2_ZJets_Zbb		= (TH1F*) ZJets_Zbb->Get("hPtb2");
    TH1F* hPtjj_ZJets_Zbb		= (TH1F*) ZJets_Zbb->Get("hPtjj");
    TH1F* hPtmumu_ZJets_Zbb		= (TH1F*) ZJets_Zbb->Get("hPtmumu");
    TH1F* hPtbalZH_ZJets_Zbb		= (TH1F*) ZJets_Zbb->Get("hPtbalZH");
    TH1F* hPtmu1_ZJets_Zbb			= (TH1F*) ZJets_Zbb->Get("hPtmu1");
    TH1F* hPtmu2_ZJets_Zbb			= (TH1F*) ZJets_Zbb->Get("hPtmu2");
    TH1F* hEtamu1_ZJets_Zbb		= (TH1F*) ZJets_Zbb->Get("hEtamu1");
    TH1F* hEtamu2_ZJets_Zbb		= (TH1F*) ZJets_Zbb->Get("hEtamu2");
    TH1F* hCSV1_ZJets_Zbb			= (TH1F*) ZJets_Zbb->Get("hCSV1");
    TH1F* hCSV2_ZJets_Zbb			= (TH1F*) ZJets_Zbb->Get("hCSV2");
    TH1F* hdphiVH_ZJets_Zbb		= (TH1F*) ZJets_Zbb->Get("hdphiVH");
    TH1F* hdetaJJ_ZJets_Zbb		= (TH1F*) ZJets_Zbb->Get("hdetaJJ");
    TH1F* hdphiJJ_ZJets_Zbb		= (TH1F*) ZJets_Zbb->Get("hdphiJJ");
    TH1F* hNaj_ZJets_Zbb			= (TH1F*) ZJets_Zbb->Get("hNaj");
    TH1F* hMmumu_ZJets_Zbb			= (TH1F*) ZJets_Zbb->Get("hMmumu");
    TH1F* hMjj_ZJets_Zbb			= (TH1F*) ZJets_Zbb->Get("hMjj");
	TH1F* hMjj_aftercuts_ZJets_Zbb			= (TH1F*) ZJets_Zbb->Get("hMjj_aftercuts");
    TH1F* hCutFlow_ZJets_Zbb		= (TH1F*) ZJets_Zbb->Get("hCutFlow");
	TH2F* hDphiDetajj_ZJets_Zbb	= (TH2F*) ZJets_Zbb->Get("hDphiDetajj");
	
    TH1F* hnJets_WZ_Zbb			= (TH1F*) WZ_Zbb->Get("hnJets");
    TH1F* hnMuons_WZ_Zbb		= (TH1F*) WZ_Zbb->Get("hnMuons");
    TH1F* hPtb1_WZ_Zbb			= (TH1F*) WZ_Zbb->Get("hPtb1");
    TH1F* hPtb2_WZ_Zbb			= (TH1F*) WZ_Zbb->Get("hPtb2");
    TH1F* hPtjj_WZ_Zbb			= (TH1F*) WZ_Zbb->Get("hPtjj");
    TH1F* hPtmumu_WZ_Zbb		= (TH1F*) WZ_Zbb->Get("hPtmumu");
    TH1F* hPtbalZH_WZ_Zbb		= (TH1F*) WZ_Zbb->Get("hPtbalZH");
    TH1F* hPtmu1_WZ_Zbb			= (TH1F*) WZ_Zbb->Get("hPtmu1");
    TH1F* hEtamu2_WZ_Zbb			= (TH1F*) WZ_Zbb->Get("hEtamu1");
    TH1F* hEtamu1_WZ_Zbb			= (TH1F*) WZ_Zbb->Get("hEtamu2");
    TH1F* hPtmu2_WZ_Zbb			= (TH1F*) WZ_Zbb->Get("hPtmu2");
    TH1F* hCSV1_WZ_Zbb			= (TH1F*) WZ_Zbb->Get("hCSV1");
    TH1F* hCSV2_WZ_Zbb			= (TH1F*) WZ_Zbb->Get("hCSV2");
    TH1F* hdphiVH_WZ_Zbb		= (TH1F*) WZ_Zbb->Get("hdphiVH");
    TH1F* hdetaJJ_WZ_Zbb		= (TH1F*) WZ_Zbb->Get("hdetaJJ");
    TH1F* hdphiJJ_WZ_Zbb		= (TH1F*) WZ_Zbb->Get("hdphiJJ");
    TH1F* hNaj_WZ_Zbb			= (TH1F*) WZ_Zbb->Get("hNaj");
    TH1F* hMmumu_WZ_Zbb			= (TH1F*) WZ_Zbb->Get("hMmumu");
    TH1F* hMjj_WZ_Zbb			= (TH1F*) WZ_Zbb->Get("hMjj");
	TH1F* hMjj_aftercuts_WZ_Zbb			= (TH1F*) WZ_Zbb->Get("hMjj_aftercuts");
    TH1F* hCutFlow_WZ_Zbb		= (TH1F*) WZ_Zbb->Get("hCutFlow");
    TH2F* hDphiDetajj_WZ_Zbb		= (TH2F*) WZ_Zbb->Get("hDphiDetajj");
	
    TH1F* hnJets_TTJets_Zbb			= (TH1F*) TTJets_Zbb->Get("hnJets");
    TH1F* hnMuons_TTJets_Zbb		= (TH1F*) TTJets_Zbb->Get("hnMuons");
    TH1F* hPtb1_TTJets_Zbb			= (TH1F*) TTJets_Zbb->Get("hPtb1");
    TH1F* hPtb2_TTJets_Zbb			= (TH1F*) TTJets_Zbb->Get("hPtb2");
    TH1F* hPtjj_TTJets_Zbb			= (TH1F*) TTJets_Zbb->Get("hPtjj");
    TH1F* hPtmumu_TTJets_Zbb		= (TH1F*) TTJets_Zbb->Get("hPtmumu");
    TH1F* hPtbalZH_TTJets_Zbb		= (TH1F*) TTJets_Zbb->Get("hPtbalZH");
    TH1F* hPtmu1_TTJets_Zbb			= (TH1F*) TTJets_Zbb->Get("hPtmu1");
    TH1F* hEtamu2_TTJets_Zbb			= (TH1F*) TTJets_Zbb->Get("hEtamu1");
    TH1F* hEtamu1_TTJets_Zbb			= (TH1F*) TTJets_Zbb->Get("hEtamu2");
    TH1F* hPtmu2_TTJets_Zbb			= (TH1F*) TTJets_Zbb->Get("hPtmu2");
    TH1F* hCSV1_TTJets_Zbb			= (TH1F*) TTJets_Zbb->Get("hCSV1");
    TH1F* hCSV2_TTJets_Zbb			= (TH1F*) TTJets_Zbb->Get("hCSV2");
    TH1F* hdphiVH_TTJets_Zbb		= (TH1F*) TTJets_Zbb->Get("hdphiVH");
    TH1F* hdetaJJ_TTJets_Zbb		= (TH1F*) TTJets_Zbb->Get("hdetaJJ");
    TH1F* hdphiJJ_TTJets_Zbb		= (TH1F*) TTJets_Zbb->Get("hdphiJJ");
    TH1F* hNaj_TTJets_Zbb			= (TH1F*) TTJets_Zbb->Get("hNaj");
    TH1F* hMmumu_TTJets_Zbb			= (TH1F*) TTJets_Zbb->Get("hMmumu");
    TH1F* hMjj_TTJets_Zbb			= (TH1F*) TTJets_Zbb->Get("hMjj");
	TH1F* hMjj_aftercuts_TTJets_Zbb			= (TH1F*) TTJets_Zbb->Get("hMjj_aftercuts");
    TH1F* hCutFlow_TTJets_Zbb		= (TH1F*) TTJets_Zbb->Get("hCutFlow");
    TH2F* hDphiDetajj_TTJets_Zbb		= (TH2F*) TTJets_Zbb->Get("hDphiDetajj");
	
	
	

//	double TTTo2L2Nu2B_weight				=	(TTTo2L2Nu2B_Sigma*1000*TTTo2L2Nu2B_SkimEff)/TTTo2L2Nu2B_EventsCutbyCut;
	double H115_weight = 0.4107*0.704*0.101*10000/220000.0;
	
	//DYJetsToLL_PtZ-100_TuneZ2_7TeV-madgraph-tauola
	double DYJetsToLL_TuneZ2_M50_weight = 3151.864553*10000/ 36179400.0; //DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola
	//DYJetsToNuNu_PtZ-100_TuneZ2_7TeV-madgraph
	//QCD_Pt-80to120_TuneZ2_7TeV_pythia6
	//QCD_Pt-120to170_TuneZ2_7TeV_pythia6
	//QCD_Pt-170to300_TuneZ2_7TeV_pythia6
	//QCD_Pt-300to470_TuneZ2_7TeV_pythia6
	//QCD_Pt-470to600_TuneZ2_7TeV_pythia6
	//QCD_Pt-20_MuEnrichedPt-15_TuneZ2_7TeV-pythia6
	double TTJets_weight = 157.5*10000/3611944.0;  //TTJets_TuneZ2_7TeV-madgraph-tauola
	//T_TuneZ2_s-channel_7TeV-powheg-tauola
	//T_TuneZ2_t-channel_7TeV-powheg-tauola
	//T_TuneZ2_tW-channel-DR_7TeV-powheg-tauola
	//T_TuneZ2_tW-channel-DS_7TeV-powheg-tauola
	//Tbar_TuneZ2_s-channel_7TeV-powheg-tauola
	//Tbar_TuneZ2_t-channel_7TeV-powheg-tauola
	//Tbar_TuneZ2_tW-channel-DR_7TeV-powheg-tauola
	//Tbar_TuneZ2_tW-channel-DS_7TeV-powheg-tauola
	//WJetsToLNu_Pt-100_7TeV-herwigpp
	//WJetsToLNu_PtW-100_TuneZ2_7TeV-madgraph
	//WJetsToLNu_TuneZ2_7TeV-madgraph-tauola
	//WW_TuneZ2_7TeV_pythia6_tauola
	double WZ_weight = 18.2*10000/4145240.0;   //WZ_TuneZ2_7TeV_pythia6_tauola
	double ZJetsToLL_weight = 3048*10000/2659998.0;  //ZJetsToLL_Pt-100_7TeV-herwigpp
	//ZJetsToNuNu_Pt-100_7TeV-herwigpp
	double ZZ_weight = 5.9*10000/4157882.0; //ZZ_TuneZ2_7TeV_pythia6_tauola
	
	//scale histograms
     hnJets_ZZ_Zbb		->Scale(ZZ_weight);
     hnMuons_ZZ_Zbb		->Scale(ZZ_weight);
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
	 hDphiDetajj_ZZ_Zbb->Scale(ZZ_weight);
	
     hnJets_DY_Zbb		->Scale(DYJetsToLL_TuneZ2_M50_weight);
     hnMuons_DY_Zbb		->Scale(DYJetsToLL_TuneZ2_M50_weight);
     hPtb1_DY_Zbb		->Scale(DYJetsToLL_TuneZ2_M50_weight);
     hPtb2_DY_Zbb		->Scale(DYJetsToLL_TuneZ2_M50_weight);
     hPtjj_DY_Zbb		->Scale(DYJetsToLL_TuneZ2_M50_weight);
     hPtmumu_DY_Zbb		->Scale(DYJetsToLL_TuneZ2_M50_weight);
	 hPtbalZH_DY_Zbb->Scale(DYJetsToLL_TuneZ2_M50_weight);
     hPtmu1_DY_Zbb		->Scale(DYJetsToLL_TuneZ2_M50_weight);
     hPtmu2_DY_Zbb		->Scale(DYJetsToLL_TuneZ2_M50_weight);
	hEtamu1_DY_Zbb->Scale(ZZ_weight);
	 hEtamu2_DY_Zbb->Scale(ZZ_weight);
     hCSV1_DY_Zbb		->Scale(DYJetsToLL_TuneZ2_M50_weight);
     hCSV2_DY_Zbb		->Scale(DYJetsToLL_TuneZ2_M50_weight);
     hdphiVH_DY_Zbb		->Scale(DYJetsToLL_TuneZ2_M50_weight);
	hdetaJJ_DY_Zbb		->Scale(DYJetsToLL_TuneZ2_M50_weight);
	hdphiJJ_DY_Zbb		->Scale(DYJetsToLL_TuneZ2_M50_weight);
     hNaj_DY_Zbb		->Scale(DYJetsToLL_TuneZ2_M50_weight);
     hMmumu_DY_Zbb		->Scale(DYJetsToLL_TuneZ2_M50_weight);
     hMjj_DY_Zbb		->Scale(DYJetsToLL_TuneZ2_M50_weight);
	 hMjj_aftercuts_DY_Zbb->Scale(DYJetsToLL_TuneZ2_M50_weight);
     hCutFlow_DY_Zbb		->Scale(DYJetsToLL_TuneZ2_M50_weight);
	 hDphiDetajj_DY_Zbb->Scale(DYJetsToLL_TuneZ2_M50_weight);
	 
	hnJets_ZJets_Zbb		->Scale(ZJetsToLL_weight);
	hnMuons_ZJets_Zbb		->Scale(ZJetsToLL_weight);
	hPtb1_ZJets_Zbb		->Scale(ZJetsToLL_weight);
	hPtb2_ZJets_Zbb		->Scale(ZJetsToLL_weight);
	hPtjj_ZJets_Zbb		->Scale(ZJetsToLL_weight);
	hPtmumu_ZJets_Zbb		->Scale(ZJetsToLL_weight);
	hPtbalZH_ZJets_Zbb ->Scale(ZJetsToLL_weight);
	hPtmu1_ZJets_Zbb		->Scale(ZJetsToLL_weight);
	hPtmu2_ZJets_Zbb		->Scale(ZJetsToLL_weight);
	hEtamu1_ZJets_Zbb->Scale(ZJetsToLL_weight);
	hEtamu2_ZJets_Zbb->Scale(ZJetsToLL_weight);
	hCSV1_ZJets_Zbb		->Scale(ZJetsToLL_weight);
	hCSV2_ZJets_Zbb		->Scale(ZJetsToLL_weight);
	hdphiVH_ZJets_Zbb		->Scale(ZJetsToLL_weight);
	hdetaJJ_ZJets_Zbb		->Scale(ZJetsToLL_weight);
	hdphiJJ_ZJets_Zbb		->Scale(ZJetsToLL_weight);
	hNaj_ZJets_Zbb			->Scale(ZJetsToLL_weight);
	hMmumu_ZJets_Zbb		->Scale(ZJetsToLL_weight);
	hMjj_ZJets_Zbb			->Scale(ZJetsToLL_weight);
	hCutFlow_ZJets_Zbb		->Scale(ZJetsToLL_weight);
	hDphiDetajj_ZJets_Zbb->Scale(ZJetsToLL_weight);
	
	hnJets_WZ_Zbb		->Scale(WZ_weight);
	hnMuons_WZ_Zbb		->Scale(WZ_weight);
	hPtb1_WZ_Zbb		->Scale(WZ_weight);
	hPtb2_WZ_Zbb		->Scale(WZ_weight);
	hPtjj_WZ_Zbb		->Scale(WZ_weight);
	hPtmumu_WZ_Zbb		->Scale(WZ_weight);
	hPtbalZH_WZ_Zbb ->Scale(WZ_weight);
	hPtmu1_WZ_Zbb		->Scale(WZ_weight);
	hPtmu2_WZ_Zbb		->Scale(WZ_weight);
	hEtamu1_WZ_Zbb->Scale(WZ_weight);
	hEtamu2_WZ_Zbb->Scale(WZ_weight);
	hCSV1_WZ_Zbb		->Scale(WZ_weight);
	hCSV2_WZ_Zbb		->Scale(WZ_weight);
	hdphiVH_WZ_Zbb		->Scale(WZ_weight);
	hdetaJJ_WZ_Zbb		->Scale(WZ_weight);
	hdphiJJ_WZ_Zbb		->Scale(WZ_weight);
	hNaj_WZ_Zbb			->Scale(WZ_weight);
	hMmumu_WZ_Zbb		->Scale(WZ_weight);
	hMjj_WZ_Zbb			->Scale(WZ_weight);
	hCutFlow_WZ_Zbb		->Scale(WZ_weight);
	hDphiDetajj_WZ_Zbb->Scale(WZ_weight);
	
	hnJets_TTJets_Zbb		->Scale(TTJets_weight);
	hnMuons_TTJets_Zbb		->Scale(TTJets_weight);
	hPtb1_TTJets_Zbb		->Scale(TTJets_weight);
	hPtb2_TTJets_Zbb		->Scale(TTJets_weight);
	hPtjj_TTJets_Zbb		->Scale(TTJets_weight);
	hPtmumu_TTJets_Zbb		->Scale(TTJets_weight);
	hPtbalZH_TTJets_Zbb ->Scale(TTJets_weight);
	hPtmu1_TTJets_Zbb		->Scale(TTJets_weight);
	hPtmu2_TTJets_Zbb		->Scale(TTJets_weight);
	hEtamu1_TTJets_Zbb->Scale(TTJets_weight);
	hEtamu2_TTJets_Zbb->Scale(TTJets_weight);
	hCSV1_TTJets_Zbb		->Scale(TTJets_weight);
	hCSV2_TTJets_Zbb		->Scale(TTJets_weight);
	hdphiVH_TTJets_Zbb		->Scale(TTJets_weight);
	hdetaJJ_TTJets_Zbb		->Scale(TTJets_weight);
	hdphiJJ_TTJets_Zbb		->Scale(TTJets_weight);
	hNaj_TTJets_Zbb			->Scale(TTJets_weight);
	hMmumu_TTJets_Zbb		->Scale(TTJets_weight);
	hMjj_TTJets_Zbb			->Scale(TTJets_weight);
	hCutFlow_TTJets_Zbb		->Scale(TTJets_weight);
	hDphiDetajj_TTJets_Zbb->Scale(TTJets_weight);
	
	
	
	// ********************************************************************
	// Have to create output file here
	// ********************************************************************	
  	// Create the root file for the histograms
	TFile *theHistogramFile;
    theHistogramFile = new TFile("Backgrounds.root", "RECREATE");
    theHistogramFile->cd();
		
	//Stack plot
	THStack *hnJets_BkgStack = new THStack("hnJets_BkgStack","Stacked Background Njets");
	hnJets_ZZ_Zbb->GetXaxis()->SetTitle("Njets");
	hnJets_ZZ_Zbb->SetFillColor(kGray);
    hnJets_BkgStack->Add(hnJets_ZZ_Zbb);
	hnJets_DY_Zbb->SetFillColor(kRed);
    hnJets_BkgStack->Add(hnJets_DY_Zbb);
	hnJets_TTJets_Zbb->SetFillColor(kBlue);
    hnJets_BkgStack->Add(hnJets_TTJets_Zbb);
	hnJets_WZ_Zbb->SetFillColor(kGray+1);
    hnJets_BkgStack->Add(hnJets_WZ_Zbb);
	hnJets_ZJets_Zbb->SetFillColor(kMagenta+2);
    hnJets_BkgStack->Add(hnJets_ZJets_Zbb);
	//hjetCSV3_QCD_Zbb->SetFillColor(kMagenta); 
    //hjetCSV3_BkgStack->Add(hjetCSV3_QCD_Zbb);
	//hjetCSV3_SingleTop_Zbb->SetFillColor(kCyan); 
	//hjetCSV3_BkgStack->Add(hjetCSV3_SingleTop_Zbb);
	//hjetCSV3_Wbb_Zbb->SetFillColor(kGreen);	
	//hjetCSV3_Wcc_Zbb->SetFillColor(kGreen+3);	
	//hjetCSV3_Zbb_Zbb->SetFillColor(kYellow);	
	//hjetCSV3_Zcc_Zbb->SetFillColor(kOragne);	
	hnJets_BkgStack->SetTitle("Number of Good Jets");
	hnJets_BkgStack->Draw("hist");
	TLegend myLegend(0.7, 0.5, 0.89, 0.8);
	myLegend.SetTextSize(0.03);
	myLegend.AddEntry(hnJets_DY_Zbb, "DY", "f");	
	myLegend.AddEntry(hnJets_ZZ_Zbb, "ZZ", "f");	
	myLegend.AddEntry(hnJets_TTJets_Zbb, "TTJets", "f");	
	myLegend.AddEntry(hnJets_WZ_Zbb, "WZ", "f");	
	myLegend.AddEntry(hnJets_ZJets_Zbb, "ZJets", "f");	
	myLegend.Draw();		
	plot = directory+"nJets"+suffixps;
	c1->Print(plot);
	c1->Clear();	
	
	THStack *hnMuons_BkgStack = new THStack("hnMuons_BkgStack","Stacked Background Nmuons");
	hnMuons_ZZ_Zbb->GetXaxis()->SetTitle("Nmuons");
	hnMuons_ZZ_Zbb->SetFillColor(kGray);
    hnMuons_BkgStack->Add(hnMuons_ZZ_Zbb);
	hnMuons_DY_Zbb->SetFillColor(kRed);
    hnMuons_BkgStack->Add(hnMuons_DY_Zbb);
	hnMuons_TTJets_Zbb->SetFillColor(kBlue);
    hnMuons_BkgStack->Add(hnMuons_TTJets_Zbb);
	hnMuons_WZ_Zbb->SetFillColor(kGray+1);
    hnMuons_BkgStack->Add(hnMuons_WZ_Zbb);
	hnMuons_ZJets_Zbb->SetFillColor(kMagenta+2);
    hnMuons_BkgStack->Add(hnMuons_ZJets_Zbb);
	hnMuons_BkgStack->SetTitle("Number of Good Muons");
	hnMuons_BkgStack->Draw("hist");
	TLegend myLegend(0.7, 0.5, 0.89, 0.8);
	myLegend.SetTextSize(0.03);
	myLegend.AddEntry(hnMuons_DY_Zbb, "DY", "f");	
	myLegend.AddEntry(hnMuons_ZZ_Zbb, "ZZ", "f");	
	myLegend.AddEntry(hnMuons_TTJets_Zbb, "TTJets", "f");	
	myLegend.AddEntry(hnMuons_WZ_Zbb, "WZ", "f");	
	myLegend.AddEntry(hnMuons_ZJets_Zbb, "ZJets", "f");	
	myLegend.Draw();		
	plot = directory+"nMuons"+suffixps;
	c1->Print(plot);
	c1->Clear();	
	
	THStack *hPtb1_BkgStack = new THStack("hPtb1_BkgStack","Stacked Background Pt of Jet with highest CSV");
	hPtb1_ZZ_Zbb->GetXaxis()->SetTitle("p_{T} (GeV)");
	hPtb1_ZZ_Zbb->SetFillColor(kGray);
    hPtb1_BkgStack->Add(hPtb1_ZZ_Zbb);
	hPtb1_DY_Zbb->SetFillColor(kRed);
    hPtb1_BkgStack->Add(hPtb1_DY_Zbb);
	hPtb1_TTJets_Zbb->SetFillColor(kBlue);
    hPtb1_BkgStack->Add(hPtb1_TTJets_Zbb);
	hPtb1_WZ_Zbb->SetFillColor(kGray+1);
    hPtb1_BkgStack->Add(hPtb1_WZ_Zbb);
	hPtb1_ZJets_Zbb->SetFillColor(kMagenta+2);
    hPtb1_BkgStack->Add(hPtb1_ZJets_Zbb);
	hPtb1_BkgStack->SetTitle("Pt of Jet with highest CSV");
	hPtb1_BkgStack->Draw("hist");
	TLegend myLegend(0.7, 0.5, 0.89, 0.8);
	myLegend.SetTextSize(0.03);
	myLegend.AddEntry(hPtb1_DY_Zbb, "DY", "f");	
	myLegend.AddEntry(hPtb1_ZZ_Zbb, "ZZ", "f");	
	myLegend.AddEntry(hPtb1_TTJets_Zbb, "TTJets", "f");	
	myLegend.AddEntry(hPtb1_WZ_Zbb, "WZ", "f");	
	myLegend.AddEntry(hPtb1_ZJets_Zbb, "ZJets", "f");	
	myLegend.Draw();		
	plot = directory+"Ptb1"+suffixps;
	c1->Print(plot);
	c1->Clear();	
	
	THStack *hPtb2_BkgStack = new THStack("hPtb2_BkgStack","Stacked Background Pt of Jet with second highest CSV");
	hPtb2_ZZ_Zbb->GetXaxis()->SetTitle("p_{T} (GeV)");
	hPtb2_ZZ_Zbb->SetFillColor(kGray);
    hPtb2_BkgStack->Add(hPtb2_ZZ_Zbb);
	hPtb2_DY_Zbb->SetFillColor(kRed);
    hPtb2_BkgStack->Add(hPtb2_DY_Zbb);
	hPtb2_TTJets_Zbb->SetFillColor(kBlue);
    hPtb2_BkgStack->Add(hPtb2_TTJets_Zbb);
	hPtb2_WZ_Zbb->SetFillColor(kGray+1);
    hPtb2_BkgStack->Add(hPtb2_WZ_Zbb);
	hPtb2_ZJets_Zbb->SetFillColor(kMagenta+2);
    hPtb2_BkgStack->Add(hPtb2_ZJets_Zbb);
	hPtb2_BkgStack->SetTitle("Pt of Jet with second highest CSV");
	hPtb2_BkgStack->Draw("hist");
	TLegend myLegend(0.7, 0.5, 0.89, 0.8);
	myLegend.SetTextSize(0.03);
	myLegend.AddEntry(hPtb2_DY_Zbb, "DY", "f");	
	myLegend.AddEntry(hPtb2_ZZ_Zbb, "ZZ", "f");	
	myLegend.AddEntry(hPtb2_TTJets_Zbb, "TTJets", "f");	
	myLegend.AddEntry(hPtb2_WZ_Zbb, "WZ", "f");	
	myLegend.AddEntry(hPtb2_ZJets_Zbb, "ZJets", "f");	
	myLegend.Draw();		
	plot = directory+"Ptb2"+suffixps;
	c1->Print(plot);
	c1->Clear();	
	
	THStack *hPtjj_BkgStack = new THStack("hPtjj_BkgStack","Stacked Background Pt of dijet pair");
	hPtjj_ZZ_Zbb->GetXaxis()->SetTitle("p_{T} (GeV)");
	hPtjj_ZZ_Zbb->SetFillColor(kGray);
    hPtjj_BkgStack->Add(hPtjj_ZZ_Zbb);
	hPtjj_DY_Zbb->SetFillColor(kRed);
    hPtjj_BkgStack->Add(hPtjj_DY_Zbb);
	hPtjj_TTJets_Zbb->SetFillColor(kBlue);
    hPtjj_BkgStack->Add(hPtjj_TTJets_Zbb);
	hPtjj_WZ_Zbb->SetFillColor(kGray+1);
    hPtjj_BkgStack->Add(hPtjj_WZ_Zbb);
	hPtjj_ZJets_Zbb->SetFillColor(kMagenta+2);
    hPtjj_BkgStack->Add(hPtjj_ZJets_Zbb);
	hPtjj_BkgStack->SetTitle("Pt of dijet pair");
	hPtjj_BkgStack->Draw("hist");
	TLegend myLegend(0.7, 0.5, 0.89, 0.8);
	myLegend.SetTextSize(0.03);
	myLegend.AddEntry(hPtjj_DY_Zbb, "DY", "f");	
	myLegend.AddEntry(hPtjj_ZZ_Zbb, "ZZ", "f");	
	myLegend.AddEntry(hPtjj_TTJets_Zbb, "TTJets", "f");	
	myLegend.AddEntry(hPtjj_WZ_Zbb, "WZ", "f");	
	myLegend.AddEntry(hPtjj_ZJets_Zbb, "ZJets", "f");	
	myLegend.Draw();		
	plot = directory+"Ptjj"+suffixps;
	c1->Print(plot);
	c1->Clear();	
	
	THStack *hPtmumu_BkgStack = new THStack("hPtmumu_BkgStack","Stacked Background Pt of dimuon pair");
	hPtmumu_ZZ_Zbb->GetXaxis()->SetTitle("p_{T} (GeV)");
	hPtmumu_ZZ_Zbb->SetFillColor(kGray);
    hPtmumu_BkgStack->Add(hPtmumu_ZZ_Zbb);
	hPtmumu_DY_Zbb->SetFillColor(kRed);
    hPtmumu_BkgStack->Add(hPtmumu_DY_Zbb);
	hPtmumu_TTJets_Zbb->SetFillColor(kBlue);
    hPtmumu_BkgStack->Add(hPtmumu_TTJets_Zbb);
	hPtmumu_WZ_Zbb->SetFillColor(kGray+1);
    hPtmumu_BkgStack->Add(hPtmumu_WZ_Zbb);
	hPtmumu_ZJets_Zbb->SetFillColor(kMagenta+2);
    hPtmumu_BkgStack->Add(hPtmumu_ZJets_Zbb);
	hPtmumu_BkgStack->SetTitle("Pt of dimuon pair");
	hPtmumu_BkgStack->Draw("hist");
	TLegend myLegend(0.7, 0.5, 0.89, 0.8);
	myLegend.SetTextSize(0.03);
	myLegend.AddEntry(hPtmumu_DY_Zbb, "DY", "f");	
	myLegend.AddEntry(hPtmumu_ZZ_Zbb, "ZZ", "f");	
	myLegend.AddEntry(hPtmumu_TTJets_Zbb, "TTJets", "f");	
	myLegend.AddEntry(hPtmumu_WZ_Zbb, "WZ", "f");	
	myLegend.AddEntry(hPtmumu_ZJets_Zbb, "ZJets", "f");	
	myLegend.Draw();		
	plot = directory+"Ptmumu"+suffixps;
	c1->Print(plot);
	c1->Clear();	
	
	THStack *hPtmu1_BkgStack = new THStack("hPtmu1_BkgStack","Stacked Background Pt of Muon with highest CSV");
	hPtmu1_ZZ_Zbb->GetXaxis()->SetTitle("p_{T} (GeV)");
	hPtmu1_ZZ_Zbb->SetFillColor(kGray);
    hPtmu1_BkgStack->Add(hPtmu1_ZZ_Zbb);
	hPtmu1_DY_Zbb->SetFillColor(kRed);
    hPtmu1_BkgStack->Add(hPtmu1_DY_Zbb);
	hPtmu1_TTJets_Zbb->SetFillColor(kBlue);
    hPtmu1_BkgStack->Add(hPtmu1_TTJets_Zbb);
	hPtmu1_WZ_Zbb->SetFillColor(kGray+1);
    hPtmu1_BkgStack->Add(hPtmu1_WZ_Zbb);
	hPtmu1_ZJets_Zbb->SetFillColor(kMagenta+2);
    hPtmu1_BkgStack->Add(hPtmu1_ZJets_Zbb);
	hPtmu1_BkgStack->SetTitle("Pt of Muon with highest CSV");
	hPtmu1_BkgStack->Draw("hist");
	TLegend myLegend(0.7, 0.5, 0.89, 0.8);
	myLegend.SetTextSize(0.03);
	myLegend.AddEntry(hPtmu1_DY_Zbb, "DY", "f");	
	myLegend.AddEntry(hPtmu1_ZZ_Zbb, "ZZ", "f");	
	myLegend.AddEntry(hPtmu1_TTJets_Zbb, "TTJets", "f");	
	myLegend.AddEntry(hPtmu1_WZ_Zbb, "WZ", "f");	
	myLegend.AddEntry(hPtmu1_ZJets_Zbb, "ZJets", "f");	
	myLegend.Draw();		
	plot = directory+"Ptmu1"+suffixps;
	c1->Print(plot);
	c1->Clear();	
	
	THStack *hPtmu2_BkgStack = new THStack("hPtmu2_BkgStack","Stacked Background Pt of Muon with seccond highest CSV");
	hPtmu2_ZZ_Zbb->GetXaxis()->SetTitle("p_{T} (GeV)");
	hPtmu2_ZZ_Zbb->SetFillColor(kGray);
    hPtmu2_BkgStack->Add(hPtmu2_ZZ_Zbb);
	hPtmu2_DY_Zbb->SetFillColor(kRed);
    hPtmu2_BkgStack->Add(hPtmu2_DY_Zbb);
	hPtmu2_TTJets_Zbb->SetFillColor(kBlue);
    hPtmu2_BkgStack->Add(hPtmu2_TTJets_Zbb);
	hPtmu2_WZ_Zbb->SetFillColor(kGray+1);
    hPtmu2_BkgStack->Add(hPtmu2_WZ_Zbb);
	hPtmu2_ZJets_Zbb->SetFillColor(kMagenta+2);
    hPtmu2_BkgStack->Add(hPtmu2_ZJets_Zbb);
	hPtmu2_BkgStack->SetTitle("Pt of Muon with second highest Pt");
	hPtmu2_BkgStack->Draw("hist");
	TLegend myLegend(0.7, 0.5, 0.89, 0.8);
	myLegend.SetTextSize(0.03);
	myLegend.AddEntry(hPtmu2_DY_Zbb, "DY", "f");	
	myLegend.AddEntry(hPtmu2_ZZ_Zbb, "ZZ", "f");
	myLegend.AddEntry(hPtmu2_TTJets_Zbb, "TTJets", "f");	
	myLegend.AddEntry(hPtmu2_WZ_Zbb, "WZ", "f");	
	myLegend.AddEntry(hPtmu2_ZJets_Zbb, "ZJets", "f");		
	myLegend.Draw();		
	plot = directory+"Ptmu2"+suffixps;
	c1->Print(plot);
	c1->Clear();	
	
	THStack *hEtamu1_BkgStack = new THStack("hEtamu1_BkgStack","Stacked Background Eta of Muon with highest CSV");
	hEtamu1_ZZ_Zbb->GetXaxis()->SetTitle("#eta");
	hEtamu1_ZZ_Zbb->SetFillColor(kGray);
    hEtamu1_BkgStack->Add(hEtamu1_ZZ_Zbb);
	hEtamu1_DY_Zbb->SetFillColor(kRed);
    hEtamu1_BkgStack->Add(hEtamu1_DY_Zbb);
	hEtamu1_TTJets_Zbb->SetFillColor(kBlue);
    hEtamu1_BkgStack->Add(hEtamu1_TTJets_Zbb);
	hEtamu1_WZ_Zbb->SetFillColor(kGray+1);
    hEtamu1_BkgStack->Add(hEtamu1_WZ_Zbb);
	hEtamu1_ZJets_Zbb->SetFillColor(kMagenta+2);
    hEtamu1_BkgStack->Add(hEtamu1_ZJets_Zbb);
	hEtamu1_BkgStack->SetTitle("Eta of Muon with highest Pt");
	hEtamu1_BkgStack->Draw("hist");
	TLegend myLegend(0.7, 0.5, 0.89, 0.8);
	myLegend.SetTextSize(0.03);
	myLegend.AddEntry(hEtamu1_DY_Zbb, "DY", "f");	
	myLegend.AddEntry(hEtamu1_ZZ_Zbb, "ZZ", "f");	
	myLegend.AddEntry(hEtamu1_TTJets_Zbb, "TTJets", "f");	
	myLegend.AddEntry(hEtamu1_WZ_Zbb, "WZ", "f");	
	myLegend.AddEntry(hEtamu1_ZJets_Zbb, "ZJets", "f");	
	myLegend.Draw();		
	plot = directory+"Etamu1"+suffixps;
	c1->Print(plot);
	c1->Clear();	
	
	THStack *hEtamu2_BkgStack = new THStack("hEtamu2_BkgStack","Stacked Background Eta of Muon with seccond highest CSV");
	hEtamu2_ZZ_Zbb->GetXaxis()->SetTitle("#eta");
	hEtamu2_ZZ_Zbb->SetFillColor(kGray);
    hEtamu2_BkgStack->Add(hEtamu2_ZZ_Zbb);
	hEtamu2_DY_Zbb->SetFillColor(kRed);
    hEtamu2_BkgStack->Add(hEtamu2_DY_Zbb);
	hEtamu2_TTJets_Zbb->SetFillColor(kBlue);
    hEtamu2_BkgStack->Add(hEtamu2_TTJets_Zbb);
	hEtamu2_WZ_Zbb->SetFillColor(kGray+1);
    hEtamu2_BkgStack->Add(hEtamu2_WZ_Zbb);
	hEtamu2_ZJets_Zbb->SetFillColor(kMagenta+2);
    hEtamu2_BkgStack->Add(hEtamu2_ZJets_Zbb);
	hEtamu2_BkgStack->SetTitle("Eta of Muon with second highest Pt");
	hEtamu2_BkgStack->Draw("hist");
	TLegend myLegend(0.7, 0.5, 0.89, 0.8);
	myLegend.SetTextSize(0.03);
	myLegend.AddEntry(hEtamu2_DY_Zbb, "DY", "f");	
	myLegend.AddEntry(hEtamu2_ZZ_Zbb, "ZZ", "f");	
	myLegend.AddEntry(hEtamu2_TTJets_Zbb, "TTJets", "f");	
	myLegend.AddEntry(hEtamu2_WZ_Zbb, "WZ", "f");	
	myLegend.AddEntry(hEtamu2_ZJets_Zbb, "ZJets", "f");	
	myLegend.Draw();		
	plot = directory+"Etamu2"+suffixps;
	c1->Print(plot);
	c1->Clear();	
	
	THStack *hCSV2_BkgStack = new THStack("hCSV2_BkgStack","Stacked Background Second Jet ordered in CSV");
	hCSV2_ZZ_Zbb->GetXaxis()->SetTitle("CSV");
	hCSV2_DY_Zbb->SetFillColor(kRed);
    hCSV2_BkgStack->Add(hCSV2_DY_Zbb);
	hCSV2_ZJets_Zbb->SetFillColor(kMagenta+2);
    hCSV2_BkgStack->Add(hCSV2_ZJets_Zbb);
	hCSV2_TTJets_Zbb->SetFillColor(kBlue);
    hCSV2_BkgStack->Add(hCSV2_TTJets_Zbb);
	hCSV2_WZ_Zbb->SetFillColor(kGray+1);
    hCSV2_BkgStack->Add(hCSV2_WZ_Zbb);
	hCSV2_ZZ_Zbb->SetFillColor(kGray);
    hCSV2_BkgStack->Add(hCSV2_ZZ_Zbb);
	hCSV2_BkgStack->SetTitle("Second Jet ordered in CSV");
	hCSV2_BkgStack->Draw("hist");
	TLegend myLegend(0.7, 0.5, 0.89, 0.8);
	myLegend.SetTextSize(0.03);
	myLegend.AddEntry(hCSV2_DY_Zbb, "DY", "f");	
	myLegend.AddEntry(hCSV2_ZZ_Zbb, "ZZ", "f");	
	myLegend.AddEntry(hCSV2_TTJets_Zbb, "TTJets", "f");	
	myLegend.AddEntry(hCSV2_WZ_Zbb, "WZ", "f");	
	myLegend.AddEntry(hCSV2_ZJets_Zbb, "ZJets", "f");	
	myLegend.Draw();		
	plot = directory+"CSV2"+suffixps;
	c1->Print(plot);
	c1->Clear();	
	
	THStack *hCSV1_BkgStack = new THStack("hCSV1_BkgStack","Stacked Background First Jet ordered in CSV");
	hCSV1_ZZ_Zbb->GetXaxis()->SetTitle("CSV");
	hCSV1_ZZ_Zbb->SetFillColor(kGray);
    hCSV1_BkgStack->Add(hCSV1_ZZ_Zbb);
	hCSV1_DY_Zbb->SetFillColor(kRed);
    hCSV1_BkgStack->Add(hCSV1_DY_Zbb);
	hCSV1_TTJets_Zbb->SetFillColor(kBlue);
    hCSV1_BkgStack->Add(hCSV1_TTJets_Zbb);
	hCSV1_WZ_Zbb->SetFillColor(kGray+1);
    hCSV1_BkgStack->Add(hCSV1_WZ_Zbb);
	hCSV1_ZJets_Zbb->SetFillColor(kMagenta+2);
    hCSV1_BkgStack->Add(hCSV1_ZJets_Zbb);
	hCSV1_BkgStack->SetTitle("Jet with largest CSV");
	hCSV1_BkgStack->Draw("hist");
	TLegend myLegend(0.7, 0.5, 0.89, 0.8);
	myLegend.SetTextSize(0.03);
	myLegend.AddEntry(hCSV1_DY_Zbb, "DY", "f");	
	myLegend.AddEntry(hCSV1_ZZ_Zbb, "ZZ", "f");	
	myLegend.AddEntry(hCSV1_TTJets_Zbb, "TTJets", "f");	
	myLegend.AddEntry(hCSV1_WZ_Zbb, "WZ", "f");	
	myLegend.AddEntry(hCSV1_ZJets_Zbb, "ZJets", "f");	
	myLegend.Draw();		
	plot = directory+"CSV1"+suffixps;
	c1->Print(plot);
	c1->Clear();	
	
	THStack *hdphiVH_BkgStack = new THStack("hdphiVH_BkgStack","Stacked Background Delta Phi of Vector Boson and Higgs");
	hdphiVH_ZZ_Zbb->GetXaxis()->SetTitle("#Delta#phi");
	hdphiVH_ZZ_Zbb->SetFillColor(kGray);
    hdphiVH_BkgStack->Add(hdphiVH_ZZ_Zbb);
	hdphiVH_DY_Zbb->SetFillColor(kRed);
    hdphiVH_BkgStack->Add(hdphiVH_DY_Zbb);
	hdphiVH_TTJets_Zbb->SetFillColor(kBlue);
    hdphiVH_BkgStack->Add(hdphiVH_TTJets_Zbb);
	hdphiVH_WZ_Zbb->SetFillColor(kGray+1);
    hdphiVH_BkgStack->Add(hdphiVH_WZ_Zbb);
	hdphiVH_ZJets_Zbb->SetFillColor(kMagenta+2);
    hdphiVH_BkgStack->Add(hdphiVH_ZJets_Zbb);
	hdphiVH_BkgStack->SetTitle("Delta Phi of Vector Boson and Higgs");
	hdphiVH_BkgStack->Draw("hist");
	TLegend myLegend(0.3, 0.6, 0.49, 0.9);
	myLegend.SetTextSize(0.03);
	myLegend.AddEntry(hdphiVH_DY_Zbb, "DY", "f");	
	myLegend.AddEntry(hdphiVH_ZZ_Zbb, "ZZ", "f");	
	myLegend.AddEntry(hdphiVH_TTJets_Zbb, "TTJets", "f");	
	myLegend.AddEntry(hdphiVH_WZ_Zbb, "WZ", "f");	
	myLegend.AddEntry(hdphiVH_ZJets_Zbb, "ZJets", "f");	
	myLegend.Draw();		
	plot = directory+"dphiVH"+suffixps;
	c1->Print(plot);
	c1->Clear();

	THStack *hdEtaJJ_BkgStack = new THStack("hdEtaJJ_BkgStack","Stacked Background #Delta#eta of the dijet pair");
	hdetaJJ_ZZ_Zbb->GetXaxis()->SetTitle("#Delta#eta");
	hdetaJJ_ZZ_Zbb->SetFillColor(kGray);
    hdEtaJJ_BkgStack->Add(hdetaJJ_ZZ_Zbb);
	hdetaJJ_DY_Zbb->SetFillColor(kRed);
    hdEtaJJ_BkgStack->Add(hdetaJJ_DY_Zbb);
	hdetaJJ_TTJets_Zbb->SetFillColor(kBlue);
    hdEtaJJ_BkgStack->Add(hdetaJJ_TTJets_Zbb);
	hdetaJJ_WZ_Zbb->SetFillColor(kGray+1);
    hdEtaJJ_BkgStack->Add(hdetaJJ_WZ_Zbb);
	hdetaJJ_ZJets_Zbb->SetFillColor(kMagenta+2);
    hdEtaJJ_BkgStack->Add(hdetaJJ_ZJets_Zbb);
	hdEtaJJ_BkgStack->SetTitle("#Delta#eta of the dijet pair");
	hdEtaJJ_BkgStack->Draw("hist");
	TLegend myLegend(0.15, 0.7, 0.28, 0.93);
	myLegend.SetTextSize(0.03);
	myLegend.AddEntry(hdetaJJ_DY_Zbb, "DY", "f");	
	myLegend.AddEntry(hdetaJJ_ZZ_Zbb, "ZZ", "f");	
	myLegend.AddEntry(hdetaJJ_TTJets_Zbb, "TTJets", "f");	
	myLegend.AddEntry(hdetaJJ_WZ_Zbb, "WZ", "f");	
	myLegend.AddEntry(hdetaJJ_ZJets_Zbb, "ZJets", "f");	
	myLegend.Draw();		
	plot = directory+"DeltaEtaJJ"+suffixps;
	c1->Print(plot);
	c1->Clear();	
	
	THStack *hNaj_BkgStack = new THStack("hNaj_BkgStack","Stacked Background Additional Jets");
	hNaj_ZZ_Zbb->GetXaxis()->SetTitle("Njets");
	hNaj_ZZ_Zbb->SetFillColor(kGray);
    hNaj_BkgStack->Add(hNaj_ZZ_Zbb);
	hNaj_DY_Zbb->SetFillColor(kRed);
    hNaj_BkgStack->Add(hNaj_DY_Zbb);
	hNaj_TTJets_Zbb->SetFillColor(kBlue);
    hNaj_BkgStack->Add(hNaj_TTJets_Zbb);
	hNaj_WZ_Zbb->SetFillColor(kGray+1);
    hNaj_BkgStack->Add(hNaj_WZ_Zbb);
	hNaj_ZJets_Zbb->SetFillColor(kMagenta+2);
    hNaj_BkgStack->Add(hNaj_ZJets_Zbb);
	hNaj_BkgStack->SetTitle("Number of Jets in Addition to two b jets");
	hNaj_BkgStack->Draw("hist");
	TLegend myLegend(0.7, 0.5, 0.89, 0.8);
	myLegend.SetTextSize(0.03);
	myLegend.AddEntry(hNaj_DY_Zbb, "DY", "f");	
	myLegend.AddEntry(hNaj_ZZ_Zbb, "ZZ", "f");	
	myLegend.AddEntry(hNaj_TTJets_Zbb, "TTJets", "f");	
	myLegend.AddEntry(hNaj_WZ_Zbb, "WZ", "f");	
	myLegend.AddEntry(hNaj_ZJets_Zbb, "ZJets", "f");	
	myLegend.Draw();		
	plot = directory+"Naj"+suffixps;
	c1->Print(plot);
	c1->Clear();	
	
	THStack *hMmumu_BkgStack = new THStack("hMmumu_BkgStack","Stacked Background Mass of dimuon pair");
	hMmumu_ZZ_Zbb->GetXaxis()->SetTitle("Mass (GeV)");
	hMmumu_ZZ_Zbb->SetFillColor(kGray);
    hMmumu_BkgStack->Add(hMmumu_ZZ_Zbb);
	hMmumu_DY_Zbb->SetFillColor(kRed);
    hMmumu_BkgStack->Add(hMmumu_DY_Zbb);
	hMmumu_TTJets_Zbb->SetFillColor(kBlue);
    hMmumu_BkgStack->Add(hMmumu_TTJets_Zbb);
	hMmumu_WZ_Zbb->SetFillColor(kGray+1);
    hMmumu_BkgStack->Add(hMmumu_WZ_Zbb);
	hMmumu_ZJets_Zbb->SetFillColor(kMagenta+2);
    hMmumu_BkgStack->Add(hMmumu_ZJets_Zbb);
	hMmumu_BkgStack->SetTitle("Mass of dimuon pair");
	hMmumu_BkgStack->Draw("hist");
	TLegend myLegend(0.7, 0.5, 0.89, 0.8);
	myLegend.SetTextSize(0.03);
	myLegend.AddEntry(hMmumu_DY_Zbb, "DY", "f");	
	myLegend.AddEntry(hMmumu_ZZ_Zbb, "ZZ", "f");	
	myLegend.AddEntry(hMmumu_TTJets_Zbb, "TTJets", "f");	
	myLegend.AddEntry(hMmumu_WZ_Zbb, "WZ", "f");	
	myLegend.AddEntry(hMmumu_ZJets_Zbb, "ZJets", "f");	
	myLegend.Draw();		
	plot = directory+"Mmumu"+suffixps;
	c1->Print(plot);
	c1->Clear();	
	
	THStack *hMjj_BkgStack = new THStack("hMjj_BkgStack","Stacked Background Mass of dijet pair");
	hMjj_ZZ_Zbb->GetXaxis()->SetTitle("Mass (GeV)");
	hMjj_ZZ_Zbb->SetFillColor(kGray);
    hMjj_BkgStack->Add(hMjj_ZZ_Zbb);
	hMjj_DY_Zbb->SetFillColor(kRed);
    hMjj_BkgStack->Add(hMjj_DY_Zbb);
	hMjj_TTJets_Zbb->SetFillColor(kBlue);
    hMjj_BkgStack->Add(hMjj_TTJets_Zbb);
	hMjj_WZ_Zbb->SetFillColor(kGray+1);
    hMjj_BkgStack->Add(hMjj_WZ_Zbb);
	hMjj_ZJets_Zbb->SetFillColor(kMagenta+2);
    hMjj_BkgStack->Add(hMjj_ZJets_Zbb);
	hMjj_BkgStack->SetTitle("Mass of dijet pair");
	hMjj_BkgStack->Draw("hist");
	TLegend myLegend(0.7, 0.5, 0.89, 0.8);
	myLegend.SetTextSize(0.03);
	myLegend.AddEntry(hMjj_DY_Zbb, "DY", "f");	
	myLegend.AddEntry(hMjj_ZZ_Zbb, "ZZ", "f");	
	myLegend.AddEntry(hMjj_TTJets_Zbb, "TTJets", "f");	
	myLegend.AddEntry(hMjj_WZ_Zbb, "WZ", "f");	
	myLegend.AddEntry(hMjj_ZJets_Zbb, "ZJets", "f");	
	myLegend.Draw();		
	plot = directory+"Mjj"+suffixps;
	c1->Print(plot);
	c1->Clear();	

	THStack *hPtBalZH_BkgStack = new THStack("hPtBalZH_BkgStack","Stacked Background Pt Balance ZH");
	hPtbalZH_ZZ_Zbb->GetXaxis()->SetTitle("#Delta p_{T} (GeV)");
	hPtbalZH_ZZ_Zbb->SetFillColor(kGray);
    hPtBalZH_BkgStack->Add(hPtbalZH_ZZ_Zbb);
	hPtbalZH_DY_Zbb->SetFillColor(kRed);
    hPtBalZH_BkgStack->Add(hPtbalZH_DY_Zbb);
	hPtbalZH_TTJets_Zbb->SetFillColor(kBlue);
    hPtBalZH_BkgStack->Add(hPtbalZH_TTJets_Zbb);
	hPtbalZH_WZ_Zbb->SetFillColor(kGray+1);
    hPtBalZH_BkgStack->Add(hPtbalZH_WZ_Zbb);
	hPtbalZH_ZJets_Zbb->SetFillColor(kMagenta+2);
    hPtBalZH_BkgStack->Add(hPtbalZH_ZJets_Zbb);
	hPtBalZH_BkgStack->SetTitle("Pt balance of Z and H");
	hPtBalZH_BkgStack->Draw("hist");
	TLegend myLegend(0.7, 0.5, 0.89, 0.8);
	myLegend.SetTextSize(0.03);
	myLegend.AddEntry(hPtbalZH_DY_Zbb, "DY", "f");	
	myLegend.AddEntry(hPtbalZH_ZZ_Zbb, "ZZ", "f");	
	myLegend.AddEntry(hPtbalZH_TTJets_Zbb, "TTJets", "f");	
	myLegend.AddEntry(hPtbalZH_WZ_Zbb, "WZ", "f");	
	myLegend.AddEntry(hPtbalZH_ZJets_Zbb, "ZJets", "f");	
	myLegend.Draw();		
	plot = directory+"PtBalZH"+suffixps;
	c1->Print(plot);
	c1->Clear();	
	
	THStack *hdphiJJ_BkgStack = new THStack("hdphiJJ_BkgStack","Stacked Background #Delta#phi JJ");
	hdphiJJ_ZZ_Zbb->GetXaxis()->SetTitle("#Delta#phi");
	hdphiJJ_ZZ_Zbb->SetFillColor(kGray);
    hdphiJJ_BkgStack->Add(hdphiJJ_ZZ_Zbb);
	hdphiJJ_DY_Zbb->SetFillColor(kRed);
    hdphiJJ_BkgStack->Add(hdphiJJ_DY_Zbb);
	hdphiJJ_TTJets_Zbb->SetFillColor(kBlue);
    hdphiJJ_BkgStack->Add(hdphiJJ_TTJets_Zbb);
	hdphiJJ_WZ_Zbb->SetFillColor(kGray+1);
    hdphiJJ_BkgStack->Add(hdphiJJ_WZ_Zbb);
	hdphiJJ_ZJets_Zbb->SetFillColor(kMagenta+2);
    hdphiJJ_BkgStack->Add(hdphiJJ_ZJets_Zbb);
	hdphiJJ_BkgStack->SetTitle("#Delta#phi JJ");
	hdphiJJ_BkgStack->Draw("hist");
	TLegend myLegend(0.7, 0.5, 0.89, 0.8);
	myLegend.SetTextSize(0.03);
	myLegend.AddEntry(hdphiJJ_DY_Zbb, "DY", "f");	
	myLegend.AddEntry(hdphiJJ_ZZ_Zbb, "ZZ", "f");	
	myLegend.AddEntry(hdphiJJ_TTJets_Zbb, "TTJets", "f");	
	myLegend.AddEntry(hdphiJJ_WZ_Zbb, "DY", "f");	
	myLegend.AddEntry(hdphiJJ_ZJets_Zbb, "ZJets", "f");	
	myLegend.Draw();		
	plot = directory+"dphiJJ"+suffixps;
	c1->Print(plot);
	c1->Clear();
	
	THStack *hMjj_aftercuts_BkgStack = new THStack("hMjj_aftercuts_BkgStack","Stacked Background JJ Mass after Cuts");
	hMjj_aftercuts_ZZ_Zbb->GetXaxis()->SetTitle("Mass (GeV)");
	hMjj_aftercuts_ZZ_Zbb->SetFillColor(kGray);
    hMjj_aftercuts_BkgStack->Add(hMjj_aftercuts_ZZ_Zbb);
	hMjj_aftercuts_DY_Zbb->SetFillColor(kRed);
    hMjj_aftercuts_BkgStack->Add(hMjj_aftercuts_DY_Zbb);
	hMjj_aftercuts_TTJets_Zbb->SetFillColor(kBlue);
    hMjj_aftercuts_BkgStack->Add(hMjj_aftercuts_TTJets_Zbb);
	hMjj_aftercuts_WZ_Zbb->SetFillColor(kGray+1);
    hMjj_aftercuts_BkgStack->Add(hMjj_aftercuts_WZ_Zbb);
	hMjj_aftercuts_ZJets_Zbb->SetFillColor(kMagenta+2);
    hMjj_aftercuts_BkgStack->Add(hMjj_aftercuts_ZJets_Zbb);
	hMjj_aftercuts_BkgStack->SetTitle("JJ Mass After Cuts");
	hMjj_aftercuts_BkgStack->Draw("hist");
	TLegend myLegend(0.7, 0.5, 0.89, 0.8);
	myLegend.SetTextSize(0.03);
	myLegend.AddEntry(hMjj_aftercuts_DY_Zbb, "DY", "f");	
	myLegend.AddEntry(hMjj_aftercuts_ZZ_Zbb, "ZZ", "f");	
	myLegend.AddEntry(hMjj_aftercuts_TTJets_Zbb, "TTJets", "f");	
	myLegend.AddEntry(hMjj_aftercuts_WZ_Zbb, "WZ", "f");	
	myLegend.AddEntry(hMjj_aftercuts_ZJets_Zbb, "ZJets", "f");	
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
	hDphiDetajj_DY_Zbb->SetMarkerSize(.1);
    hDphiDetajj_BkgStack->Add(hDphiDetajj_DY_Zbb);
	hDphiDetajj_TTJets_Zbb->SetMarkerColor(kBlue);
	hDphiDetajj_TTJets_Zbb->SetMarkerStyle(20);
	hDphiDetajj_TTJets_Zbb->SetMarkerSize(.1);
    hDphiDetajj_BkgStack->Add(hDphiDetajj_TTJets_Zbb);
	hDphiDetajj_WZ_Zbb->SetMarkerColor(kOrange);
	hDphiDetajj_WZ_Zbb->SetMarkerStyle(20);
	hDphiDetajj_WZ_Zbb->SetMarkerSize(.1);
    hDphiDetajj_BkgStack->Add(hDphiDetajj_WZ_Zbb);
	hDphiDetajj_ZJets_Zbb->SetMarkerColor(kCyan);
	hDphiDetajj_ZJets_Zbb->SetMarkerStyle(20);
	hDphiDetajj_ZJets_Zbb->SetMarkerSize(.1);
    hDphiDetajj_BkgStack->Add(hDphiDetajj_ZJets_Zbb);
	hDphiDetajj_BkgStack->SetTitle(" #Delta#eta vs #Delta#phi of jets");
	hDphiDetajj_BkgStack->Draw("scat");
	TLegend myLegend(0.7, 0.5, 0.89, 0.8);
	myLegend.SetTextSize(0.03);
	myLegend.AddEntry(hDphiDetajj_DY_Zbb, "DY", "p");	
	myLegend.AddEntry(hDphiDetajj_ZZ_Zbb, "ZZ", "p");	
	myLegend.AddEntry(hDphiDetajj_TTJets_Zbb, "TTJets", "p");	
	myLegend.AddEntry(hDphiDetajj_WZ_Zbb, "WZ", "p");	
	myLegend.AddEntry(hDphiDetajj_ZJets_Zbb, "ZJets", "p");	
	myLegend.Draw();		
	plot = directory+"DphiDetajj"+suffixps;
	c1->Print(plot);
	c1->Clear();	
	hDphiDetajj_ZZ_Zbb->Add(hDphiDetajj_DY_Zbb);
	hDphiDetajj_ZZ_Zbb->Add(hDphiDetajj_TTJets_Zbb);
	hDphiDetajj_ZZ_Zbb->Add(hDphiDetajj_WZ_Zbb);
	hDphiDetajj_ZZ_Zbb->Add(hDphiDetajj_ZJets_Zbb);
	gStyle->SetPalette(1);
	hDphiDetajj_ZZ_Zbb->Draw("colz");
	plot = directory+"DphiDetajj_Colz"+suffixps;
	c1->Print(plot);
	c1->Clear();	
	
	
	
	hCutFlow_ZZ_Zbb->GetXaxis()->SetBinLabel(1,"PreSelection");
	hCutFlow_ZZ_Zbb->GetXaxis()->SetBinLabel(2,"p_{T}(jj)");
	hCutFlow_ZZ_Zbb->GetXaxis()->SetBinLabel(3,"p_{T}(Z)");
	hCutFlow_ZZ_Zbb->GetXaxis()->SetBinLabel(4,"CSV1");
	hCutFlow_ZZ_Zbb->GetXaxis()->SetBinLabel(5,"CSV2");
	hCutFlow_ZZ_Zbb->GetXaxis()->SetBinLabel(6,"#Delta#phi");
	hCutFlow_ZZ_Zbb->GetXaxis()->SetBinLabel(7,"N_{aj}");
	hCutFlow_ZZ_Zbb->GetXaxis()->SetBinLabel(8,"M_{jj}");		
	//THStack *hMjj_BkgStack = new THStack("hMjj_BkgStack","Stacked Background Cut Flow");
	gStyle->SetOptStat("kFALSE");
	hCutFlow_ZZ_Zbb->GetYaxis()->SetTitle("Number of Events");
	hCutFlow_ZZ_Zbb->SetTitle("Cut Flow HZ");
	//hCutFlow_ZZ_Zbb->SetAxisRange(0,9,"X");
	hCutFlow_ZZ_Zbb->GetXaxis()->SetRangeUser(0,9);
	hCutFlow_DY_Zbb->GetXaxis()->SetRangeUser(0,9);
	hCutFlow_TTJets_Zbb->GetXaxis()->SetRangeUser(0,9);
	hCutFlow_WZ_Zbb->GetXaxis()->SetRangeUser(0,9);
	hCutFlow_ZJets_Zbb->GetXaxis()->SetRangeUser(0,9);
	//hCutFlow_ZZ_Zbb->SetFillColor(kGray);
    //hMjj_BkgStack->Add(hCutFlow_ZZ_Zbb);
	hCutFlow_ZZ_Zbb->SetMarkerColor(kGreen);
	hCutFlow_ZZ_Zbb->SetMarkerStyle(3); //kStar
	hCutFlow_ZZ_Zbb->SetMarkerSize(2);
	
	//hCutFlow_DY_Zbb->SetFillColor(kRed);
    //hMjj_BkgStack->Add(hCutFlow_DY_Zbb);
	hCutFlow_DY_Zbb->SetMarkerColor(kRed);
	hCutFlow_DY_Zbb->SetMarkerStyle(23); //kFullTriangleDown
	hCutFlow_DY_Zbb->SetMarkerSize(2);
	hCutFlow_DY_Zbb->Draw("p");
	hCutFlow_ZZ_Zbb->Draw("same p");
	//hCutFlow_TTJets_Zbb->SetFillColor(kBlue);
    //hMjj_BkgStack->Add(hCutFlow_TTJets_Zbb);
	hCutFlow_TTJets_Zbb->SetMarkerColor(kBlue);
	hCutFlow_TTJets_Zbb->SetMarkerStyle(20);
	hCutFlow_TTJets_Zbb->SetMarkerSize(2);
	hCutFlow_TTJets_Zbb->Draw("same p");
	//hCutFlow_WZ_Zbb->SetFillColor(kGray+1);
    //hMjj_BkgStack->Add(hCutFlow_WZ_Zbb);
	hCutFlow_WZ_Zbb->SetMarkerColor(kOrange);
	hCutFlow_WZ_Zbb->SetMarkerStyle(kFullSquare);
	hCutFlow_WZ_Zbb->SetMarkerSize(2);
	hCutFlow_WZ_Zbb->Draw("same p");
	//hCutFlow_ZJets_Zbb->SetFillColor(kMagenta+2);
    //hMjj_BkgStack->Add(hCutFlow_ZJets_Zbb);
	hCutFlow_ZJets_Zbb->SetMarkerColor(kCyan);
	hCutFlow_ZJets_Zbb->SetMarkerStyle(kFullTriangleUp);
	hCutFlow_ZJets_Zbb->SetMarkerSize(2);
	hCutFlow_ZJets_Zbb->Draw("same p");
	TLegend myLegend(0.7, 0.5, 0.89, 0.8);
	myLegend.SetTextSize(0.03);
	myLegend.AddEntry(hCutFlow_DY_Zbb, "DY", "p");	
	myLegend.AddEntry(hCutFlow_ZZ_Zbb, "ZZ", "p");	
	myLegend.AddEntry(hCutFlow_TTJets_Zbb, "TTJets", "p");	
	myLegend.AddEntry(hCutFlow_WZ_Zbb, "WZ", "p");	
	myLegend.AddEntry(hCutFlow_ZJets_Zbb, "ZJets", "p");	
	myLegend.Draw();	
	plot = directory+"CutFlow"+suffixps;	
	//hMjj_BkgStack->SetTitle("Cut Flow HZ");
	//hMjj_BkgStack->Draw("hist");
	//plot = directory+"CutFlowStack"+suffixps;
	c1->Print(plot);
	c1->Clear();	

	theHistogramFile->Write();
	
//	hHighestPt_TTTo2L2Nu2B->Delete();
		
}
