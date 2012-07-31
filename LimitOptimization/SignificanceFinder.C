#include <TH2.h>
#include <TStyle.h>
#include <TFile.h>
#include "TH2F.h"
#include "THStack.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TRandom.h"
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <strstream>
#include <cstring>
#include <string>
//#include <iomanip.h>
//#include "TStopwatch.h"
#include "TTree.h"
#include "TROOT.h"
#include "TF2.h"
#include "TMath.h"
#include "TVirtualPad.h"
#include "TStyle.h"
#include "Riostream.h"
#include "TColor.h"
#include "TClass.h"
#include <vector>
#include "/home/hep/wilken/taus/CMSSW_4_4_2_patch8/src/UserCode/wilken/interface/xsecV21.h" 

using namespace std;

#include "TMVAClassification_BDT.class.C"
//#include "IClassifierReader.class.C"

void SignificanceFinder(){
	std::vector<string> BDTnames;
	BDTnames.push_back("Hmass");
	BDTnames.push_back("CSV0");
	BDTnames.push_back("Emumass");
	BDTnames.push_back("DeltaPhiHV");
	BDTnames.push_back("Mt");
	BDTnames.push_back("dPhiHMET");
	BDTnames.push_back("abs(Dphiemu)");
	BDTnames.push_back("abs(DphiZMET)");
	BDTnames.push_back("PtbalMETH");
	BDTnames.push_back("EtaStandDev");
	BDTnames.push_back("jetCHF0");
	BDTnames.push_back("ProjVisT");

	ReadBDT CalcBDT(BDTnames);
//	IClassifierReader CalcBDT;
    
	TFile *inputBtt  =   TFile::Open("/home/hep/wilken/taus/CMSSW_4_4_2_patch8/src/UserCode/wilken/V21_noWeights_92TTJets/TTJets.root");
	//	TFile *inputBzz  =   TFile::Open("/home/hep/wilken/taus/CMSSW_4_4_2_patch8/src/UserCode/wilken/V21_noWeights_92TTJets/ZZ.root"            );
	//	TFile *inputTt   =   TFile::Open("/home/hep/wilken/taus/CMSSW_4_4_2_patch8/src/UserCode/wilken/V21_noWeights_92TTJets/T_tchannel.root"         );
	//	TFile *inputtw   =   TFile::Open("/home/hep/wilken/taus/CMSSW_4_4_2_patch8/src/UserCode/wilken/V21_noWeights_92TTJets/T_tWchannel_DR.root"     );
	//	TFile *inputtbw   =  TFile::Open("/home/hep/wilken/taus/CMSSW_4_4_2_patch8/src/UserCode/wilken/V21_noWeights_92TTJets/Tbar_tWchannel_DR.root"  );
	//TFile *inputww   =   TFile::Open("/home/hep/wilken/taus/CMSSW_4_4_2_patch8/src/UserCode/wilken/V21_noWeights_92TTJets/WW.root"		      );
	//	TFile *inputwz   =   TFile::Open("/home/hep/wilken/taus/CMSSW_4_4_2_patch8/src/UserCode/wilken/V21_noWeights_92TTJets/WZ.root"                 );
	//	TFile *inputwj   =   TFile::Open("/home/hep/wilken/taus/CMSSW_4_4_2_patch8/src/UserCode/wilken/V21_noWeights_92TTJets/WJets.root"        );
	TFile *inputBzj   =  TFile::Open("/home/hep/wilken/taus/CMSSW_4_4_2_patch8/src/UserCode/wilken/V21_noWeights_92TTJets/DY_M50.root"   );
	//TFile *inputBzj2   = TFile::Open("/home/hep/wilken/taus/CMSSW_4_4_2_patch8/src/UserCode/wilken/V21_noWeights_92TTJets/DY_PtZ.root");
	
	//	TFile *inputS = TFile::Open("/home/hep/wilken/taus/CMSSW_4_4_2_patch8/src/UserCode/wilken/V21_noWeights_92TTJets/ZH115_1M.root");
	TFile *inputS = TFile::Open("/home/hep/wilken/taus/CMSSW_4_4_2_patch8/src/UserCode/wilken/V21_noWeights_92TTJets/Filteredemu.root");
	
	TTree *signal        = (TTree*) inputS->Get("FOM_tree"); 
	TTree *backgroundtt  = (TTree*) inputBtt ->Get("FOM_tree");   
	//	TTree *backgroundzz  = (TTree*) inputBzz ->Get("FOM_tree"); 
	TTree *backgroundzj  = (TTree*) inputBzj ->Get("FOM_tree"); 
	//TTree *backgroundzj2  = (TTree*) inputBzj2 ->Get("FOM_tree"); 
	//	TTree *backgroundTt = (TTree*) inputTt->Get("FOM_tree");
	//	TTree *backgroundtw  = (TTree*) inputtw ->Get("FOM_tree");
	//	TTree *backgroundtbw  = (TTree*) inputtbw ->Get("FOM_tree");
	
	//TTree *backgroundww  = (TTree*) inputww ->Get("FOM_tree");
	//	TTree *backgroundwz  = (TTree*) inputwz ->Get("FOM_tree");
	//	TTree *backgroundwj  = (TTree*) inputwj ->Get("FOM_tree");
	
	
	
	
	float Hmass, Emumass;
	float Hpt, Zpt;
	float CSV0, CSV1;
	float jetCHF0, jetCHF1;
	float DphiJJ, DetaJJ, delRjj;
	float Hphi, Zphi, DeltaPhiHV, Ht, MET, EventPt, PtbalZH, Mt, Mte, Mtmu, DphiZMET;
	float lep0pt, lep1pt, lep_pfCombRelIso0, lep_pfCombRelIso1, AngleEMU, delRemu;
	float EtaStandDev, RMS_eta, UnweightedEta;
	float Centrality, EvntShpAplanarity, EvntShpSphericity, EvntShpIsotropy, EvntShpCircularity;
	int nJets, eventFlavor, naJets, nSV;
	float PtbalMETH, PtbalZMET, DphiSecondMET, Detaemu, ScalarSumPt, Dphiemu, dPhiHMET, METsig;
	float ProjVisT, ProjMissT, ZmassSVD;
	
	signal->SetBranchAddress( "nJets", &nJets );
	signal->SetBranchAddress( "naJets", &naJets );
	signal->SetBranchAddress( "nSV", &nSV );	
	signal->SetBranchAddress( "Hmass", &Hmass );
	signal->SetBranchAddress( "Emumass", &Emumass );
	signal->SetBranchAddress( "Hpt", &Hpt );
	signal->SetBranchAddress( "Zpt", &Zpt );
	signal->SetBranchAddress( "CSV0", &CSV0 );
	signal->SetBranchAddress( "CSV1", &CSV1 );
	signal->SetBranchAddress( "jetCHF0", &jetCHF0 );
	signal->SetBranchAddress( "jetCHF1", &jetCHF1 );
	signal->SetBranchAddress( "DetaJJ", &DetaJJ );
	signal->SetBranchAddress( "DphiJJ", &DphiJJ );
	signal->SetBranchAddress( "delRjj", &delRjj );
	signal->SetBranchAddress( "Hphi", &Hphi );
	signal->SetBranchAddress( "Zphi", &Zphi );
	signal->SetBranchAddress( "DeltaPhiHV", &DeltaPhiHV );
	signal->SetBranchAddress( "Ht", &Ht );
	signal->SetBranchAddress( "MET", &MET );
	signal->SetBranchAddress( "EventPt", &EventPt );
	signal->SetBranchAddress( "PtbalZH", &PtbalZH );
	signal->SetBranchAddress( "Mt", &Mt);
	signal->SetBranchAddress( "Mte", &Mte );
	signal->SetBranchAddress( "Mtmu", &Mtmu );
	signal->SetBranchAddress( "DphiZMET", &DphiZMET );
	signal->SetBranchAddress( "Zphi", &Zphi );
	signal->SetBranchAddress( "Hphi", &Hphi );
	signal->SetBranchAddress( "lep_pfCombRelIso0", &lep_pfCombRelIso0 );
	signal->SetBranchAddress( "lep_pfCombRelIso1", &lep_pfCombRelIso1 );
	signal->SetBranchAddress( "lep0pt", &lep0pt );
	signal->SetBranchAddress( "lep1pt", &lep1pt );
	signal->SetBranchAddress( "AngleEMU", &AngleEMU );
	signal->SetBranchAddress( "delRemu", &delRemu );
	signal->SetBranchAddress( "RMS_eta", &RMS_eta );
	signal->SetBranchAddress( "UnweightedEta", &UnweightedEta );
	signal->SetBranchAddress( "EtaStandDev", &EtaStandDev );
	signal->SetBranchAddress( "Centrality", &Centrality );
	signal->SetBranchAddress( "EvntShpAplanarity", &EvntShpAplanarity );
	signal->SetBranchAddress( "EvntShpSphericity", &EvntShpSphericity );
	signal->SetBranchAddress( "EvntShpIsotropy", &EvntShpIsotropy );
	signal->SetBranchAddress( "EvntShpCircularity", &EvntShpCircularity );
	signal->SetBranchAddress( "ZmassSVD", &ZmassSVD );
	signal->SetBranchAddress( "PtbalMETH", &PtbalMETH );
	signal->SetBranchAddress( "PtbalZMET", &PtbalZMET );
	signal->SetBranchAddress( "DphiSecondMET", &DphiSecondMET );
	signal->SetBranchAddress( "delRjj", &delRjj );
	signal->SetBranchAddress( "Detaemu", &Detaemu );
	signal->SetBranchAddress( "Dphiemu", &Dphiemu );
	signal->SetBranchAddress( "ScalarSumPt", &ScalarSumPt );
	signal->SetBranchAddress( "dPhiHMET", &dPhiHMET );
	signal->SetBranchAddress( "METsig", &METsig );
	signal->SetBranchAddress( "ZmassSVD", &ZmassSVD );
	signal->SetBranchAddress( "ProjVisT", &ProjVisT );
	signal->SetBranchAddress( "ProjMissT", &ProjMissT );
	
	backgroundtt->SetBranchAddress( "nJets", &nJets );
	backgroundtt->SetBranchAddress( "naJets", &naJets );
	backgroundtt->SetBranchAddress( "nSV", &nSV );	
	backgroundtt->SetBranchAddress( "Hmass", &Hmass );
	backgroundtt->SetBranchAddress( "Emumass", &Emumass );
	backgroundtt->SetBranchAddress( "Hpt", &Hpt );
	backgroundtt->SetBranchAddress( "Zpt", &Zpt );
	backgroundtt->SetBranchAddress( "CSV0", &CSV0 );
	backgroundtt->SetBranchAddress( "CSV1", &CSV1 );
	backgroundtt->SetBranchAddress( "jetCHF0", &jetCHF0 );
	backgroundtt->SetBranchAddress( "jetCHF1", &jetCHF1 );
	backgroundtt->SetBranchAddress( "DetaJJ", &DetaJJ );
	backgroundtt->SetBranchAddress( "DphiJJ", &DphiJJ );
	backgroundtt->SetBranchAddress( "delRjj", &delRjj );
	backgroundtt->SetBranchAddress( "Hphi", &Hphi );
	backgroundtt->SetBranchAddress( "Zphi", &Zphi );
	backgroundtt->SetBranchAddress( "DeltaPhiHV", &DeltaPhiHV );
	backgroundtt->SetBranchAddress( "Ht", &Ht );
	backgroundtt->SetBranchAddress( "MET", &MET );
	backgroundtt->SetBranchAddress( "EventPt", &EventPt );
	backgroundtt->SetBranchAddress( "PtbalZH", &PtbalZH );
	backgroundtt->SetBranchAddress( "Mt", &Mt );
	backgroundtt->SetBranchAddress( "Mte", &Mte );
	backgroundtt->SetBranchAddress( "Mtmu", &Mtmu );
	backgroundtt->SetBranchAddress( "DphiZMET", &DphiZMET );
	backgroundtt->SetBranchAddress( "Zphi", &Zphi );
	backgroundtt->SetBranchAddress( "Hphi", &Hphi );
	backgroundtt->SetBranchAddress( "lep_pfCombRelIso0", &lep_pfCombRelIso0 );
	backgroundtt->SetBranchAddress( "lep_pfCombRelIso1", &lep_pfCombRelIso1 );
	backgroundtt->SetBranchAddress( "lep0pt", &lep0pt );
	backgroundtt->SetBranchAddress( "lep1pt", &lep1pt );
	backgroundtt->SetBranchAddress( "AngleEMU", &AngleEMU );
	backgroundtt->SetBranchAddress( "delRemu", &delRemu );
	backgroundtt->SetBranchAddress( "RMS_eta", &RMS_eta );
	backgroundtt->SetBranchAddress( "UnweightedEta", &UnweightedEta );
	backgroundtt->SetBranchAddress( "EtaStandDev", &EtaStandDev );
	backgroundtt->SetBranchAddress( "Centrality", &Centrality );
	backgroundtt->SetBranchAddress( "EvntShpAplanarity", &EvntShpAplanarity );
	backgroundtt->SetBranchAddress( "EvntShpSphericity", &EvntShpSphericity );
	backgroundtt->SetBranchAddress( "EvntShpIsotropy", &EvntShpIsotropy );
	backgroundtt->SetBranchAddress( "EvntShpCircularity", &EvntShpCircularity );
	backgroundtt->SetBranchAddress( "ZmassSVD", &ZmassSVD );
	backgroundtt->SetBranchAddress( "PtbalMETH", &PtbalMETH );
	backgroundtt->SetBranchAddress( "PtbalZMET", &PtbalZMET );
	backgroundtt->SetBranchAddress( "DphiSecondMET", &DphiSecondMET );
	backgroundtt->SetBranchAddress( "delRjj", &delRjj );
	backgroundtt->SetBranchAddress( "Detaemu", &Detaemu );
	backgroundtt->SetBranchAddress( "ScalarSumPt", &ScalarSumPt );
	backgroundtt->SetBranchAddress( "Dphiemu", &Dphiemu );
	backgroundtt->SetBranchAddress( "dPhiHMET", &dPhiHMET );
	backgroundtt->SetBranchAddress( "METsig", &METsig );
	backgroundtt->SetBranchAddress( "ZmassSVD", &ZmassSVD );
	backgroundtt->SetBranchAddress( "ProjVisT", &ProjVisT );
	backgroundtt->SetBranchAddress( "ProjMissT", &ProjMissT );
	
	backgroundzj->SetBranchAddress( "nJets", &nJets );
	backgroundzj->SetBranchAddress( "naJets", &naJets );
	backgroundzj->SetBranchAddress( "nSV", &nSV );	
	backgroundzj->SetBranchAddress( "Hmass", &Hmass );
	backgroundzj->SetBranchAddress( "Emumass", &Emumass );
	backgroundzj->SetBranchAddress( "Hpt", &Hpt );
	backgroundzj->SetBranchAddress( "Zpt", &Zpt );
	backgroundzj->SetBranchAddress( "CSV0", &CSV0 );
	backgroundzj->SetBranchAddress( "CSV1", &CSV1 );
	backgroundzj->SetBranchAddress( "jetCHF0", &jetCHF0 );
	backgroundzj->SetBranchAddress( "jetCHF1", &jetCHF1 );
	backgroundzj->SetBranchAddress( "DetaJJ", &DetaJJ );
	backgroundzj->SetBranchAddress( "DphiJJ", &DphiJJ );
	backgroundzj->SetBranchAddress( "delRjj", &delRjj );
	backgroundzj->SetBranchAddress( "Hphi", &Hphi );
	backgroundzj->SetBranchAddress( "Zphi", &Zphi );
	backgroundzj->SetBranchAddress( "DeltaPhiHV", &DeltaPhiHV );
	backgroundzj->SetBranchAddress( "Ht", &Ht );
	backgroundzj->SetBranchAddress( "MET", &MET );
	backgroundzj->SetBranchAddress( "EventPt", &EventPt );
	backgroundzj->SetBranchAddress( "PtbalZH", &PtbalZH );
	backgroundzj->SetBranchAddress( "Mt", &Mt );
	backgroundzj->SetBranchAddress( "Mte", &Mte );
	backgroundzj->SetBranchAddress( "Mtmu", &Mtmu );
	backgroundzj->SetBranchAddress( "DphiZMET", &DphiZMET );
	backgroundzj->SetBranchAddress( "Zphi", &Zphi );
	backgroundzj->SetBranchAddress( "Hphi", &Hphi );
	backgroundzj->SetBranchAddress( "lep_pfCombRelIso0", &lep_pfCombRelIso0 );
	backgroundzj->SetBranchAddress( "lep_pfCombRelIso1", &lep_pfCombRelIso1 );
	backgroundzj->SetBranchAddress( "lep0pt", &lep0pt );
	backgroundzj->SetBranchAddress( "lep1pt", &lep1pt );
	backgroundzj->SetBranchAddress( "AngleEMU", &AngleEMU );
	backgroundzj->SetBranchAddress( "delRemu", &delRemu );
	backgroundzj->SetBranchAddress( "RMS_eta", &RMS_eta );
	backgroundzj->SetBranchAddress( "UnweightedEta", &UnweightedEta );
	backgroundzj->SetBranchAddress( "EtaStandDev", &EtaStandDev );
	backgroundzj->SetBranchAddress( "Centrality", &Centrality );
	backgroundzj->SetBranchAddress( "EvntShpAplanarity", &EvntShpAplanarity );
	backgroundzj->SetBranchAddress( "EvntShpSphericity", &EvntShpSphericity );
	backgroundzj->SetBranchAddress( "EvntShpIsotropy", &EvntShpIsotropy );
	backgroundzj->SetBranchAddress( "EvntShpCircularity", &EvntShpCircularity );
	backgroundzj->SetBranchAddress( "ZmassSVD", &ZmassSVD );
	backgroundzj->SetBranchAddress( "PtbalMETH", &PtbalMETH );
	backgroundzj->SetBranchAddress( "PtbalZMET", &PtbalZMET );
	backgroundzj->SetBranchAddress( "DphiSecondMET", &DphiSecondMET );
	backgroundzj->SetBranchAddress( "delRjj", &delRjj );
	backgroundzj->SetBranchAddress( "Detaemu", &Detaemu );
	backgroundzj->SetBranchAddress( "ScalarSumPt", &ScalarSumPt );
	backgroundzj->SetBranchAddress( "Dphiemu", &Dphiemu );
	backgroundzj->SetBranchAddress( "dPhiHMET", &dPhiHMET );
	backgroundzj->SetBranchAddress( "METsig", &METsig );
	backgroundzj->SetBranchAddress( "ZmassSVD", &ZmassSVD );
	backgroundzj->SetBranchAddress( "ProjVisT", &ProjVisT );
	backgroundzj->SetBranchAddress( "ProjMissT", &ProjMissT );
	
	double lumi = 4.457;
	Double_t  ZH_M115_weight = lumi/(lumiZH115);
	Double_t  TTJets_weight = lumi/(lumiTT) ;
	Double_t  DYJetsToLL_M50_weight = lumi/(lumiZJL);
	
	TString suffixps = "_Sig.gif";
	TString directory = "/home/hep/wilken/taus/CMSSW_4_4_2_patch8/src/UserCode/wilken/LimitOptimization/";
	TString plot = directory+"somereallylongnameherehopethisarrayisbiggenough"+suffixps;	
	
	double Nsig = 0, NTTJ = 0, NZJ = 0;
	float maxSig = 0.0, corrCut = 0.0, significance, SoB;
	
	const int steps = 10;
	float cutmin = 0;
	float cutmax = 3;
	float range = -99.99; 
	range = cutmax-cutmin;
	float stepsize = range/float(steps);
	cout << "range " <<range<< endl;
	float cut = -99.99;
	std::vector<double> BDTvariables;
	double BDTvalue = -99.99;

	TCanvas *c1 = new TCanvas("c1","");
	gStyle->SetOptStat("kTRUE");
	c1->SetFillColor(10);
	c1->SetFillColor(10);
	TFile *theHistogramFile;
    theHistogramFile = new TFile("BDTHistos.root", "RECREATE");
    theHistogramFile->cd();
	TH1F* hSig= new TH1F	("hSig", "BDT Signal", 20, -1.0, 0.5);
    TH1F* hTTJets= new TH1F	("hTTJets", "BDT TTJets", 20, -1.0, 0.5);
    TH1F* hDYM50= new TH1F	("hDYM50", "BDT DYM50", 20, -1.0, 0.5);

	
	for (float icut = 0; icut<steps; icut++){
	hSig->Clear();
	hTTJets->Clear();
	hDYM50->Clear();
		cut = icut*stepsize+cutmin;
		//clear root file
		
		Nsig = 0;
		for (Long64_t ievt=0; ievt<signal->GetEntries();ievt++) {
			signal->GetEntry(ievt);
			if((lep0pt > 20 || lep1pt >20) && (lep0pt > 10 && lep1pt >10)){
				/*				if (jetCHF0> 0.32 && CSV0 > 0.679 && Hmass> 65.5 && Hmass < 135 && delRemu > 0.4 && delRemu < 2 
				 && DeltaPhiHV > 2.3 && Emumass > 10 && Emumass < 65 && Hpt > 23 && AngleEMU > 0.2 && Centrality < 0.42 && PtbalZH > -33
				 && (ZmassSVD<0 || (ZmassSVD > 33 && ZmassSVD < 96)) && PtbalMETH > -20 && MET >10 && ScalarSumPt > 103
				 && fabs(DphiSecondMET)<1.276 && dPhiHMET > 2 && jetCHF1 > 0.03 && METsig > 4 && Ht > 16.5 && EventPt > 7){ */
				if (Hmass > 45 && Hmass < 150 && CSV0 > 0.244 && Emumass > 10 && Emumass < 70 
					&& delRemu > 0.4 && jetCHF0 > 0.2 && AngleEMU > 0.25 && DeltaPhiHV>2 && (ProjVisT-(0.25*ProjVisT)) > 10){
					if(fabs(DphiZMET) <cut) {
						BDTvariables.push_back(Hmass);
						BDTvariables.push_back(Hmass);
						BDTvariables.push_back(CSV0);
						BDTvariables.push_back(Emumass);
						BDTvariables.push_back(DeltaPhiHV);
						BDTvariables.push_back(Mt);
						BDTvariables.push_back(dPhiHMET);
						BDTvariables.push_back(abs(Dphiemu));
						BDTvariables.push_back(abs(DphiZMET));
						BDTvariables.push_back(PtbalMETH);
						BDTvariables.push_back(EtaStandDev);
						BDTvariables.push_back(jetCHF0);
						BDTvariables.push_back(ProjVisT);					
						BDTvalue = CalcBDT.GetMvaValue(BDTvariables);
						hSig->Fill(BDTvalue);
						Nsig++;
					}
				}}//signal trigger requirement
		}//end signal event loop
		
		NTTJ = 0;
		for (Long64_t ievt=0; ievt<backgroundtt->GetEntries();ievt++) {
			backgroundtt->GetEntry(ievt);
			if((lep0pt > 20 || lep1pt >20) && (lep0pt > 10 && lep1pt >10)){
				/*				if (jetCHF0> 0.32 && CSV0 > 0.679 && Hmass> 65.5 && Hmass < 135 && delRemu > 0.4 && delRemu < 2 
				 && DeltaPhiHV > 2.3 && Emumass > 10 && Emumass < 65 && Hpt > 23 && AngleEMU > 0.2 && Centrality < 0.42 && PtbalZH > -33
				 && (ZmassSVD<0 || (ZmassSVD > 33 && ZmassSVD < 96)) && PtbalMETH > -20 && MET >10 && ScalarSumPt > 103
				 && fabs(DphiSecondMET)<1.276 && dPhiHMET > 2 && jetCHF1 > 0.03 && METsig > 4 && Ht > 16.5 && EventPt > 7){ */
				if (Hmass > 45 && Hmass < 150 && CSV0 > 0.244 && Emumass > 10 && Emumass < 70 
					&& delRemu > 0.4 && jetCHF0 > 0.2 && AngleEMU > 0.25 && DeltaPhiHV>2 && (ProjVisT-(0.25*ProjVisT)) > 10){
					if(fabs(DphiZMET) <cut){
						BDTvariables.push_back(Hmass);
						BDTvariables.push_back(Hmass);
						BDTvariables.push_back(CSV0);
						BDTvariables.push_back(Emumass);
						BDTvariables.push_back(DeltaPhiHV);
						BDTvariables.push_back(Mt);
						BDTvariables.push_back(dPhiHMET);
						BDTvariables.push_back(abs(Dphiemu));
						BDTvariables.push_back(abs(DphiZMET));
						BDTvariables.push_back(PtbalMETH);
						BDTvariables.push_back(EtaStandDev);
						BDTvariables.push_back(jetCHF0);
						BDTvariables.push_back(ProjVisT);					
						BDTvalue = CalcBDT.GetMvaValue(BDTvariables);
						hTTJets->Fill(BDTvalue);
						NTTJ++;
					}
				}}//TTJets trigger requirement
		}//end TTJets event loop
		
		NZJ = 0;
		for (Long64_t ievt=0; ievt<backgroundzj->GetEntries();ievt++) {
			backgroundzj->GetEntry(ievt);
			if((lep0pt > 20 || lep1pt >20) && (lep0pt > 10 && lep1pt >10)){
				/*				if (jetCHF0> 0.32 && CSV0 > 0.679 && Hmass> 65.5 && Hmass < 135 && delRemu > 0.4 && delRemu < 2 
				 && DeltaPhiHV > 2.3 && Emumass > 10 && Emumass < 65 && Hpt > 23 && AngleEMU > 0.2 && Centrality < 0.42 && PtbalZH > -33
				 && (ZmassSVD<0 || (ZmassSVD > 33 && ZmassSVD < 96)) && PtbalMETH > -20 && MET >10 && ScalarSumPt > 103
				 && fabs(DphiSecondMET)<1.276 && dPhiHMET > 2 && jetCHF1 > 0.03 && METsig > 4 && Ht > 16.5 && EventPt > 7){ */
				if (Hmass > 45 && Hmass < 150 && CSV0 > 0.244 && Emumass > 10 && Emumass < 70 
					&& delRemu > 0.4 && jetCHF0 > 0.2 && AngleEMU > 0.25 && DeltaPhiHV>2 && (ProjVisT-(0.25*ProjVisT)) > 10){
					if(fabs(DphiZMET) <cut){
						BDTvariables.push_back(Hmass);
						BDTvariables.push_back(Hmass);
						BDTvariables.push_back(CSV0);
						BDTvariables.push_back(Emumass);
						BDTvariables.push_back(DeltaPhiHV);
						BDTvariables.push_back(Mt);
						BDTvariables.push_back(dPhiHMET);
						BDTvariables.push_back(abs(Dphiemu));
						BDTvariables.push_back(abs(DphiZMET));
						BDTvariables.push_back(PtbalMETH);
						BDTvariables.push_back(EtaStandDev);
						BDTvariables.push_back(jetCHF0);
						BDTvariables.push_back(ProjVisT);					
						BDTvalue = CalcBDT.GetMvaValue(BDTvariables);
						cout << "DYM50 BDT value " <<BDTvalue  << endl;
						hDYM50->Fill(BDTvalue);
						NZJ++;
					}
				}}//DY_M50 trigger requirement
		}//end TTJets event loop
		cout << "Cut " <<cut << " Number of Events Signal: " <<Nsig << " TTJets " << NTTJ<< " DY " << NZJ<< endl;
		
		float background = (NTTJ*TTJets_weight)+(NZJ*DYJetsToLL_M50_weight);
		//		float background = (NTTJ*TTJets_weight);
		significance = (Nsig*ZH_M115_weight)/(1.5+sqrt(background)+.2*background);
		SoB = (Nsig*ZH_M115_weight)/background;
		if (significance>maxSig) {
			maxSig = significance;
			corrCut = cut;
		}
		//Calculate limit here using combine function
		hSig->Write();
		hTTJets->Write();
		hDYM50->Write();
		
		THStack *histBdt_BkgStack = new THStack("histBdt_BkgStack","Stacked Background BDT");
		hDYM50->SetFillColor(kYellow);
		hTTJets->SetFillColor(kBlue);
		histBdt_BkgStack->Add(hDYM50);
		histBdt_BkgStack->Add(hTTJets);
		histBdt_BkgStack->Draw();
		hSig->Draw("same hist");
		TLegend myLegend(0.7, 0.5, 0.89, 0.8);
		myLegend.SetTextSize(0.03);
		myLegend.AddEntry(hSig, "ZH_M115 x 100", "l");	
		myLegend.AddEntry(hTTJets, "DYJetsToLL_M50", "f");	
		myLegend.AddEntry(hDYM50, "TTJets", "f");	
		myLegend.Draw();		
		//plot = directory+"DphiZMET"+cut.c_str()+suffixps;
		plot = TString::Format(directory+"%f.gif",cut).Data();
		c1->Print(plot);
		c1->Clear();		
		
		cout << "significance is " <<significance << " at " << cut  << " S/B " << SoB << endl;
	}//icut for loop
	
	cout << "The max significance is " <<maxSig << " at " << corrCut<< endl;
		
	c1->Close();
	hSig->Delete();
	hTTJets->Delete();
	hDYM50->Delete();
	
	
}
