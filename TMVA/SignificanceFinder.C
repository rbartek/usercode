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
#include "TTree.h"
#include "TROOT.h"
#include "TF2.h"
#include "TMath.h"
#include "TVirtualPad.h"
#include "TStyle.h"
#include "Riostream.h"
#include "TColor.h"
#include "TClass.h"
//#include <vector>
#include "/home/hep/wilken/taus/CMSSW_4_4_2_patch8/src/UserCode/wilken/interface/xsecV21.h" 

using namespace std;


void SignificanceFinder(){
    
	TFile *inputBtt  =   TFile::Open("/home/hep/wilken/taus/CMSSW_4_4_2_patch8/src/UserCode/wilken/V21/TTJets.root");
//	TFile *inputBzz  =   TFile::Open("/home/hep/wilken/taus/CMSSW_4_4_2_patch8/src/UserCode/wilken/V21/ZZ.root"            );
//	TFile *inputTt   =   TFile::Open("/home/hep/wilken/taus/CMSSW_4_4_2_patch8/src/UserCode/wilken/V21/T_tchannel.root"         );
//	TFile *inputtw   =   TFile::Open("/home/hep/wilken/taus/CMSSW_4_4_2_patch8/src/UserCode/wilken/V21/T_tWchannel_DR.root"     );
//	TFile *inputtbw   =  TFile::Open("/home/hep/wilken/taus/CMSSW_4_4_2_patch8/src/UserCode/wilken/V21/Tbar_tWchannel_DR.root"  );
	TFile *inputww   =   TFile::Open("/home/hep/wilken/taus/CMSSW_4_4_2_patch8/src/UserCode/wilken/V21/WW.root"		      );
//	TFile *inputwz   =   TFile::Open("/home/hep/wilken/taus/CMSSW_4_4_2_patch8/src/UserCode/wilken/V21/WZ.root"                 );
//	TFile *inputwj   =   TFile::Open("/home/hep/wilken/taus/CMSSW_4_4_2_patch8/src/UserCode/wilken/V21/WJets.root"        );
	TFile *inputBzj   =  TFile::Open("/home/hep/wilken/taus/CMSSW_4_4_2_patch8/src/UserCode/wilken/V21/DY_M50.root"   );
	TFile *inputBzj2   = TFile::Open("/home/hep/wilken/taus/CMSSW_4_4_2_patch8/src/UserCode/wilken/V21/DY_PtZ.root");
	
//	TFile *inputS = TFile::Open("/home/hep/wilken/taus/CMSSW_4_4_2_patch8/src/UserCode/wilken/V21/ZH115_1M.root");
	TFile *inputS = TFile::Open("/home/hep/wilken/taus/CMSSW_4_4_2_patch8/src/UserCode/wilken/V21/Filteredemu.root");
	
	TTree *signal        = (TTree*) inputS->Get("FOM_tree"); 
	TTree *backgroundtt  = (TTree*) inputBtt ->Get("FOM_tree");   
//	TTree *backgroundzz  = (TTree*) inputBzz ->Get("FOM_tree"); 
	TTree *backgroundzj  = (TTree*) inputBzj ->Get("FOM_tree"); 
	TTree *backgroundzj2  = (TTree*) inputBzj2 ->Get("FOM_tree"); 
//	TTree *backgroundTt = (TTree*) inputTt->Get("FOM_tree");
//	TTree *backgroundtw  = (TTree*) inputtw ->Get("FOM_tree");
	
//	TTree *backgroundtbw  = (TTree*) inputtbw ->Get("FOM_tree");
	
	
	TTree *backgroundww  = (TTree*) inputww ->Get("FOM_tree");
//	TTree *backgroundwz  = (TTree*) inputwz ->Get("FOM_tree");
//	TTree *backgroundwj  = (TTree*) inputwj ->Get("FOM_tree");
	
	
	TString suffixps = "_Sig.gif";
	TString directory = "/home/hep/wilken/taus/CMSSW_4_4_2_patch8/src/UserCode/wilken/V21/";
	TString plot = directory+"somereallylongnameherehopethisarrayisbiggenough"+suffixps;
	
	
	float Hmass, Emumass;
	float Hpt, Zpt;
	float CSV0, CSV1;
	float jetCHF0, jetCHF1;
	float DphiJJ, DetaJJ, delRjj;
	float Hphi, Zphi, DeltaPhiHV, Ht, MET, EventPt, PtbalZH, Mte, Mtmu, DphiZMET;
	float lep0pt, lep1pt, lep_pfCombRelIso0, lep_pfCombRelIso1, AngleEMU, delRemu;
	float EtaStandDev, RMS_eta, UnweightedEta;
	float Centrality, EvntShpAplanarity, EvntShpSphericity, EvntShpIsotropy, EvntShpCircularity;
	int nJets, eventFlavor, naJets, nSV;
	float ZmassSVD, PtbalMETH, PtbalZMET, DphiSecondMET, Detaemu, ScalarSumPt, Dphiemu, dPhiHMET, METsig;
	
	double lumi = 4.457;
	Double_t  ZH_M115_weight = lumi/(lumiZH115/2.0);
	Double_t  TTJets_weight = lumi/(lumiTT/2.0) ;
	Double_t  DYJetsToLL_M50_weight = lumi/(lumiZJL/2.0);
	
    TH2F* hOptiSig= new TH2F	("hOptiSig", "cut vs Significance", 1001, 0, 1.05, 1001, 0.00001, 0.001);
	double Nsig = 0, NTTJ = 0, NZJ = 0;
	float maxSig = 0.0, corrCut = 0.0, significance;
	
	const int steps = 100;
	float cutmin = 0.5;
	float cutmax = 1.5;
	float range = -99.99; 
	range = cutmax-cutmin;
	float stepsize = range/float(steps);
	cout << "range " <<range<< endl;
	float cut = -99.99;
	
	for (float icut = 0; icut<steps; icut++){
		cut = icut*stepsize+cutmin;
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

		
		Nsig = 0;
		for (Long64_t ievt=0; ievt<signal->GetEntries();ievt++) {
			signal->GetEntry(ievt);
			if((lep0pt > 20 || lep1pt >20) && (lep0pt > 10 && lep1pt >10)){
				/*				if (jetCHF0> 0.32 && CSV0 > 0.679 && Hmass> 65.5 && Hmass < 135 && delRemu > 0.4 && delRemu < 2 
				 && DeltaPhiHV > 2.3 && Emumass > 10 && Emumass < 65 && Hpt > 23 && AngleEMU > 0.2 && Centrality < 0.42 && PtbalZH > -33
				 && (ZmassSVD<0 || (ZmassSVD > 33 && ZmassSVD < 96)) && PtbalMETH > -20 && MET >10 && ScalarSumPt > 103
				 && fabs(DphiSecondMET)<1.276 && dPhiHMET > 2 && jetCHF1 > 0.03 && METsig > 4 && Ht > 16.5 && EventPt > 7){ */
				if (lep0pt > 10 && lep1pt> 10 && Hmass< 150 && Hmass > 45 && CSV0 > 0.244 && Emumass > 10 && Emumass < 70 && Hpt > 20 ){
					if(fabs(DphiZMET)<cut) Nsig++;
			}}//signal trigger requirement
		}//end signal event loop
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


		
		NTTJ = 0;
		for (Long64_t ievt=0; ievt<backgroundtt->GetEntries();ievt++) {
			backgroundtt->GetEntry(ievt);
			if((lep0pt > 20 || lep1pt >20) && (lep0pt > 10 && lep1pt >10)){
				/*				if (jetCHF0> 0.32 && CSV0 > 0.679 && Hmass> 65.5 && Hmass < 135 && delRemu > 0.4 && delRemu < 2 
				 && DeltaPhiHV > 2.3 && Emumass > 10 && Emumass < 65 && Hpt > 23 && AngleEMU > 0.2 && Centrality < 0.42 && PtbalZH > -33
				 && (ZmassSVD<0 || (ZmassSVD > 33 && ZmassSVD < 96)) && PtbalMETH > -20 && MET >10 && ScalarSumPt > 103
				 && fabs(DphiSecondMET)<1.276 && dPhiHMET > 2 && jetCHF1 > 0.03 && METsig > 4 && Ht > 16.5 && EventPt > 7){ */
				if (lep0pt > 10 && lep1pt> 10 && Hmass< 150 && Hmass > 80 && (ZmassSVD<0 ||(ZmassSVD >30 && ZmassSVD<100)) && CSV0 > 0.244 && Emumass > 10 && Emumass < 85 && Hpt > 20 ){
					if(fabs(DphiZMET)<cut) NTTJ++;
			}}//TTJets trigger requirement
		}//end TTJets event loop
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

		NZJ = 0;
		for (Long64_t ievt=0; ievt<backgroundzj->GetEntries();ievt++) {
			backgroundzj->GetEntry(ievt);
			if((lep0pt > 20 || lep1pt >20) && (lep0pt > 10 && lep1pt >10)){
/*				if (jetCHF0> 0.32 && CSV0 > 0.679 && Hmass> 65.5 && Hmass < 135 && delRemu > 0.4 && delRemu < 2 
					&& DeltaPhiHV > 2.3 && Emumass > 10 && Emumass < 65 && Hpt > 23 && AngleEMU > 0.2 && Centrality < 0.42 && PtbalZH > -33
					&& (ZmassSVD<0 || (ZmassSVD > 33 && ZmassSVD < 96)) && PtbalMETH > -20 && MET >10 && ScalarSumPt > 103
					&& fabs(DphiSecondMET)<1.276 && dPhiHMET > 2 && jetCHF1 > 0.03 && METsig > 4 && Ht > 16.5 && EventPt > 7){ */
				if (lep0pt > 10 && lep1pt> 10 && Hmass< 150 && Hmass > 80 && (ZmassSVD<0 ||(ZmassSVD >30 && ZmassSVD<100)) && CSV0 > 0.244 && Emumass > 10 && Emumass < 85 && Hpt > 20 ){
			if(fabs(DphiZMET)<cut) NZJ++;
			}}//DY_M50 trigger requirement
		}//end TTJets event loop
		cout << "Cut " <<cut << " Number of Events Signal: " <<Nsig << " TTJets " << NTTJ<< " DY " << NZJ<< endl;

		float background = (NTTJ*TTJets_weight)+(NZJ*DYJetsToLL_M50_weight);
//		float background = (NTTJ*TTJets_weight);
		significance = (Nsig*ZH_M115_weight)/(1.5+sqrt(background)+.2*background);
		SoB = (Nsig*ZH_M115_weight)/background;
		hOptiSig->Fill(cut,significance);
		if (significance>maxSig) {
		maxSig = significance;
		corrCut = cut;
		}
		cout << "significance is " <<significance << " at " << cut  << " S/B " << SoB << endl;
	}//icut for loop
	
	cout << "The max significance is " <<maxSig << " at " << corrCut<< endl;
		
	TCanvas *c1 = new TCanvas("c1","");
	gStyle->SetOptStat("kTRUE");
	c1->SetFillColor(10);
	c1->SetFillColor(10);
	
		
		hOptiSig->GetXaxis()->SetTitle("DphiZMETwithsomeCuts");
		hOptiSig->GetYaxis()->SetTitle("S/(1.5+#sqrtB+.2B)");
		hOptiSig->SetMarkerStyle(20);
		hOptiSig->SetMarkerSize(1);
		hOptiSig->Draw();
		plot = directory+"DphiZMETwithsomeCuts"+suffixps;
		c1->Print(plot);
		c1->Clear();	
	
	c1->Close();
	
}
