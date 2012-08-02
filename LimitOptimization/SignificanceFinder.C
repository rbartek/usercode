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

void makeDataCard(float VHrate, float TTrate, float ZjLFrate, float Thiscut, string dir, string shapefile);

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
	string directory = "/home/hep/wilken/taus/CMSSW_4_4_2_patch8/src/UserCode/wilken/LimitOptimization/";
	TString plot = directory+"somereallylongnameherehopethisarrayisbiggenough"+suffixps;	
	
	double Nsig = 0, NTTJ = 0, NZJ = 0;
	float maxSig = 0.0, corrCut = 0.0, significance, SoB;
	
	const int steps = 1;
	float cutmin = 1.5;
	float cutmax = 3;
	float range = -99.99; 
	range = cutmax-cutmin;
	float stepsize = range/float(steps);
	cout << "range " <<range<< endl;
	float cut = -99.99;
	//std::vector<double> BDTvariables;
	double BDTvalue = -99.99;
	
	TCanvas *c1 = new TCanvas("c1","");
	gStyle->SetOptStat("kTRUE");
	c1->SetFillColor(10);
	c1->SetFillColor(10);
/*	TFile *theHistogramFile;
	string filename = "HistoForShapeLimit.root";
    theHistogramFile = new TFile(filename, "RECREATE");
    theHistogramFile->cd();*/
	int nbins = 20;
	TH1F* VH= new TH1F	("VH", "BDT Signal", nbins, -1.0, 0);
    TH1F* TT= new TH1F	("TT", "BDT TTJets", nbins, -1.0, 0);
    TH1F* ZjLF= new TH1F	("ZjLF", "BDT DYM50", nbins, -1.0, 0);
	TH1F* VH_CMSstatUp= new TH1F	("VH_CMSstatUp", "BDT Signal plus sigma", nbins, -1.0, 0);
    TH1F* TT_CMSstatUp= new TH1F	("TT_CMSstatUp", "BDT TTJets plus sigma", nbins, -1.0, 0);
    TH1F* ZjLF_CMSstatUp= new TH1F	("ZjLF_CMSstatUp", "BDT DYM50 plus sigma", nbins, -1.0, 0);
	TH1F* VH_CMSstatDown= new TH1F	("VH_CMSstatDown", "BDT Signal minus sigma", nbins, -1.0, 0);
    TH1F* TT_CMSstatDown= new TH1F	("TT_CMSstatDown", "BDT TTJets minus sigma", nbins, -1.0, 0);
    TH1F* ZjLF_CMSstatDown= new TH1F	("ZjLF_CMSstatDown", "BDT DYM50 minus sigma", nbins, -1.0, 0);
    TH1F* data_obs= new TH1F	("data_obs", "null histo for data", nbins, -1.0, 0);
	
	for (float icut = 0; icut<steps; icut++){
		cut = icut*stepsize+cutmin;
		string filename;
		filename = TString::Format("HistoForShapeLimit%0.3f.root",cut).Data();
		TFile *theHistogramFile= new TFile(TString::Format("%s/%s",directory.c_str(),filename.c_str()).Data(), "RECREATE", "histogram file",0);
		 theHistogramFile->cd();
		
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
						std::vector<double> BDTvariables;
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
						VH->Fill(BDTvalue);
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
						std::vector<double> BDTvariables;
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
						TT->Fill(BDTvalue);
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
						std::vector<double> BDTvariables;
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
						//cout << "input values";
						//for (unsigned int i = 0 ; i < BDTvariables.size(); i++)  cout << " " << BDTvariables[i];
						//cout << endl;
						//cout << "BDTvalue " << BDTvalue << endl;
						ZjLF->Fill(BDTvalue);
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
		
		VH->Scale(ZH_M115_weight);
		TT->Scale(TTJets_weight);
		ZjLF->Scale(DYJetsToLL_M50_weight);
		VH->Write();
		TT->Write();
		ZjLF->Write();
		data_obs->Write();
		makeDataCard(VH->Integral(0,-1), TT->Integral(0,-1), ZjLF->Integral(0,-1), cut, directory, filename);

		for (int i = 0; i<nbins; i++){
			double binContent = -99.99;
			double sigma =-99.99;
			binContent = VH->GetBinContent(i);
			sigma = VH->GetBinError(i);
			VH_CMSstatUp->SetBinContent(i,binContent+sigma);
			if((binContent-sigma) < 0) sigma  = 0;
			VH_CMSstatDown->SetBinContent(i,binContent-sigma);
			binContent = TT->GetBinContent(i);
			sigma = TT->GetBinError(i);
			TT_CMSstatUp->SetBinContent(i,binContent+sigma);
			if((binContent-sigma) < 0) sigma  = 0;
			TT_CMSstatDown->SetBinContent(i,binContent-sigma);
			binContent = ZjLF->GetBinContent(i);
			sigma = ZjLF->GetBinError(i);
			ZjLF_CMSstatUp->SetBinContent(i,binContent+sigma);
			if((binContent-sigma) < 0) sigma  = 0;
			ZjLF_CMSstatDown->SetBinContent(i,binContent-sigma);
		}
		VH_CMSstatUp->Write();
		TT_CMSstatUp->Write();
		ZjLF_CMSstatUp->Write();
		VH_CMSstatDown->Write();
		TT_CMSstatDown->Write();
		ZjLF_CMSstatDown->Write();
		
		
		THStack *histBdt_BkgStack = new THStack("histBdt_BkgStack","Stacked Background BDT");
		ZjLF->SetFillColor(kYellow);
		TT->SetFillColor(kBlue);
		histBdt_BkgStack->Add(ZjLF);
		histBdt_BkgStack->Add(TT);
		histBdt_BkgStack->Draw();
		VH->SetLineWidth(3);
		VH->Scale(100);
		VH->Draw("same hist");
		TLegend myLegend(0.7, 0.5, 0.89, 0.8);
		myLegend.SetTextSize(0.03);
		myLegend.AddEntry(VH, "ZH_M115 x 100", "l");	
		myLegend.AddEntry(TT, "DYJetsToLL_M50", "f");	
		myLegend.AddEntry(ZjLF, "TTJets", "f");	
		myLegend.Draw();		
		//plot = directory+"DphiZMET"+cut.c_str()+suffixps;
		plot = TString::Format("%sBDT%0.02f.gif",directory.c_str(),cut).Data();
		c1->Print(plot);
		c1->Clear();		
		
		cout << "significance is " <<significance << " at " << cut  << " S/B " << SoB << endl;
		VH->Reset("ICES");
		TT->Reset("ICES");
		ZjLF->Reset("ICES");
		VH_CMSstatUp->Reset("ICES");
		TT_CMSstatUp->Reset("ICES");
		ZjLF_CMSstatUp->Reset("ICES");
		VH_CMSstatDown->Reset("ICES");
		TT_CMSstatDown->Reset("ICES");
		ZjLF_CMSstatDown->Reset("ICES");
		theHistogramFile->Delete();
		
	}//icut for loop
	
	cout << "The max significance is " <<maxSig << " at " << corrCut<< endl;
	
	c1->Close();
	VH->Delete();
	TT->Delete();
	ZjLF->Delete();
	VH_CMSstatUp->Delete();
	TT_CMSstatUp->Delete();
	ZjLF_CMSstatUp->Delete();
	VH_CMSstatDown->Delete();
	TT_CMSstatDown->Delete();
	ZjLF_CMSstatDown->Delete();
	data_obs->Delete();
	
	
	
}

void makeDataCard(float VHrate, float TTrate, float ZjLFrate, float Thiscut, string dir, string shapefile){
ofstream myfile (TString::Format("%s/ShapeDataCard%0.03fcut.txt",dir.c_str(),Thiscut).Data());
	myfile<<"imax 1  number of channels"<<endl;
	myfile<<"jmax 2  number of backgrounds"<<endl;;
	myfile<<"kmax *  number of nuisance parameters (sources of systematical uncertainties)"<<endl;;
	myfile<<TString::Format("shapes * * %s $PROCESS $PROCESS_$SYSTEMATIC", shapefile.c_str()).Data()<<endl;
	myfile<<"bin Zemu" <<endl;
	myfile<<"observation 0"<<endl;
	myfile<<"bin Zemu Zemu Zemu"<<endl;
	myfile<<"process \t VH \t ZjLF \t TT "<<endl;
	myfile<<"process  \t 0 \t 1 \t 2 "<<endl;
	myfile<<TString::Format("rate \t %0.3f \t %0.3f \t %0.3f ",VHrate,ZjLFrate,TTrate).Data()<<endl;
	myfile<<"\n";
	myfile<<"lumi \t lnN \t 1.022 \t - \t - "<<endl; 
	myfile<<"pdf_qqbar \t lnN \t 1.01 \t - \t - "<<endl; 
	myfile<<"QCDscale_VH \t lnN \t 1.04 \t - \t - "<<endl; 
	myfile<<"CMS_vhbb_ZjLF_SF \t lnN \t - \t 1.06 \t - "<<endl; 
	myfile<<"CMS_vhbb_ZjLF_ex \t lnN \t - \t 1.05 \t -" << endl; 
	myfile<<"CMS_vhbb_TT_SF \t lnN  \t - \t -  \t 1.21 "<<endl; 
	myfile<<"CMS_vhbb_TT_ex \tlnN \t - \t - \t 1.05"<<endl;
	myfile<<"CMS_trigger_m  \t lnN  \t 1.01  \t - \t - "<<endl; 
	myfile<<"CMS_trigger_3  \t lnN  \t 1.02  \t - \t - "<<endl; 
	myfile<<"CMS_eff_m  \t lnN  \t 1.04 \t - \t - "<<endl; 
	myfile<<"CMS_eff_e  \t lnN  \t 1.04 \t - \t - "<<endl; 
	myfile<<"CMS_eff_b  \t lnN  \t 1.11 \t 1.07 \t 1 "<<endl; 
	myfile<<"CMS_fake_b  \t lnN  \t 1.05 \t 1.12 \t 1 "<<endl; 
	myfile<<"CMS_scale_j  \t lnN  \t 1.03 \t -  \t -  "<<endl;  
	myfile<<"CMS_res_j  \t lnN  \t 1.05 \t 1.04 \t 1.04 "<<endl;  
	myfile<<"CMSstat \t shape \t 1.0 \t - \t -"<< endl;
	myfile<<"CMSstat \t shape \t - \t 1.0 \t -"<< endl;
	myfile<<"CMSstat \t shape \t - \t - \t 1.0"<< endl;
	myfile.close();
	


	


}