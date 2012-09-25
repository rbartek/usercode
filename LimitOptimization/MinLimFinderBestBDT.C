#include <TH2.h>
#include <TStyle.h>
#include <TFile.h>
#include "TH2F.h"
#include "THStack.h"
#include "TGraph.h"
#include "TMultiGraph.h";
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
#include "TSystem.h"
#include "/home/hep/wilken/taus/CMSSW_4_4_2_patch8/src/UserCode/wilken/interface/xsecV21.h" 

using namespace std;

#include "TMVAClassification_BDT.class.C"
//#include "IClassifierReader.class.C"

void makeDataCard(float VHrate, float ZjLFrate, float ZjHFrate, float TTrate, float sToprate, float VVrate, float Thiscut, string dir, string shapefile);

void MinLimFinderBestBDT(){
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
	
    
	TFile *inputBtt  =   TFile::Open("/home/hep/wilken/taus/CMSSW_4_4_2_patch8/src/UserCode/wilken/V21_bugSingleTop/TTJets.root");
	TFile *inputBzz  =   TFile::Open("/home/hep/wilken/taus/CMSSW_4_4_2_patch8/src/UserCode/wilken/V21_bugSingleTop/ZZ.root"            );
	//TFile *inputTt   =   TFile::Open("/home/hep/wilken/taus/CMSSW_4_4_2_patch8/src/UserCode/wilken/V21_bugSingleTop/T_tchannel.root"         );
	TFile *inputtw   =   TFile::Open("/home/hep/wilken/taus/CMSSW_4_4_2_patch8/src/UserCode/wilken/V21_bugSingleTop/T_tW.root"     );
	TFile *inputtbw   =  TFile::Open("/home/hep/wilken/taus/CMSSW_4_4_2_patch8/src/UserCode/wilken/V21_bugSingleTop/Tbar_tW.root"  );
	TFile *inputww   =   TFile::Open("/home/hep/wilken/taus/CMSSW_4_4_2_patch8/src/UserCode/wilken/V21_bugSingleTop/WW.root"		      );
	TFile *inputwz   =   TFile::Open("/home/hep/wilken/taus/CMSSW_4_4_2_patch8/src/UserCode/wilken/V21_bugSingleTop/WZ.root"                 );
	//	TFile *inputwj   =   TFile::Open("/home/hep/wilken/taus/CMSSW_4_4_2_patch8/src/UserCode/wilken/V21_bugSingleTop/WJets.root"        );
	TFile *inputBzjLF   =  TFile::Open("/home/hep/wilken/taus/CMSSW_4_4_2_patch8/src/UserCode/wilken/V21_bugSingleTop/DY_LF_M50.root"   );
	TFile *inputBzjHF   =  TFile::Open("/home/hep/wilken/taus/CMSSW_4_4_2_patch8/src/UserCode/wilken/V21_bugSingleTop/DY_HF_M50.root"   );
	TFile *inputBzj2LF   = TFile::Open("/home/hep/wilken/taus/CMSSW_4_4_2_patch8/src/UserCode/wilken/V21_bugSingleTop/DY_LF_PtZ.root");
	TFile *inputBzj2HF   = TFile::Open("/home/hep/wilken/taus/CMSSW_4_4_2_patch8/src/UserCode/wilken/V21_bugSingleTop/DY_HF_PtZ.root");
	//	TFile *inputS = TFile::Open("/home/hep/wilken/taus/CMSSW_4_4_2_patch8/src/UserCode/wilken/V21_bugSingleTop/ZH115_1M.root");
	TFile *inputS = TFile::Open("/home/hep/wilken/taus/CMSSW_4_4_2_patch8/src/UserCode/wilken/V21_bugSingleTop/Filteredemu.root");
	
	TTree *signal        = (TTree*) inputS->Get("FOM_tree"); 
	TTree *backgroundtt  = (TTree*) inputBtt ->Get("FOM_tree");   
	TTree *backgroundzz  = (TTree*) inputBzz ->Get("FOM_tree"); 
	TTree *backgroundzjLF  = (TTree*) inputBzjLF ->Get("FOM_tree"); 
	TTree *backgroundzjHF  = (TTree*) inputBzjHF ->Get("FOM_tree"); 
	TTree *backgroundzj2LF  = (TTree*) inputBzj2LF ->Get("FOM_tree"); 
	TTree *backgroundzj2HF  = (TTree*) inputBzj2HF ->Get("FOM_tree"); 
	//TTree *backgroundTt = (TTree*) inputTt->Get("FOM_tree");
	TTree *backgroundtw  = (TTree*) inputtw ->Get("FOM_tree");
	TTree *backgroundtbw  = (TTree*) inputtbw ->Get("FOM_tree");
	TTree *backgroundww  = (TTree*) inputww ->Get("FOM_tree");
	TTree *backgroundwz  = (TTree*) inputwz ->Get("FOM_tree");
	//	TTree *backgroundwj  = (TTree*) inputwj ->Get("FOM_tree");
	
	vector <TTree*> sample;
	sample.push_back(signal); //0
	sample.push_back(backgroundtt); //1
	sample.push_back(backgroundzjLF); //2
	sample.push_back(backgroundzj2LF); //3
	sample.push_back(backgroundzjHF); //4
	sample.push_back(backgroundzj2HF);//5 
	sample.push_back(backgroundtw);//6 
	sample.push_back(backgroundtbw);//7
	sample.push_back(backgroundww);//8 
	sample.push_back(backgroundwz);//9
	sample.push_back(backgroundzz);//10 
	
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
	float Trigweight, B2011PUweight, A2011PUweight;
	
	for (unsigned int i = 0; i<sample.size(); i++){
		sample[i]->SetBranchAddress( "nJets", &nJets );
		sample[i]->SetBranchAddress( "naJets", &naJets );
		sample[i]->SetBranchAddress( "nSV", &nSV );	
		sample[i]->SetBranchAddress( "Hmass", &Hmass );
		sample[i]->SetBranchAddress( "Emumass", &Emumass );
		sample[i]->SetBranchAddress( "Hpt", &Hpt );
		sample[i]->SetBranchAddress( "Zpt", &Zpt );
		sample[i]->SetBranchAddress( "CSV0", &CSV0 );
		sample[i]->SetBranchAddress( "CSV1", &CSV1 );
		sample[i]->SetBranchAddress( "jetCHF0", &jetCHF0 );
		sample[i]->SetBranchAddress( "jetCHF1", &jetCHF1 );
		sample[i]->SetBranchAddress( "DetaJJ", &DetaJJ );
		sample[i]->SetBranchAddress( "DphiJJ", &DphiJJ );
		sample[i]->SetBranchAddress( "delRjj", &delRjj );
		sample[i]->SetBranchAddress( "Hphi", &Hphi );
		sample[i]->SetBranchAddress( "Zphi", &Zphi );
		sample[i]->SetBranchAddress( "DeltaPhiHV", &DeltaPhiHV );
		sample[i]->SetBranchAddress( "Ht", &Ht );
		sample[i]->SetBranchAddress( "MET", &MET );
		sample[i]->SetBranchAddress( "EventPt", &EventPt );
		sample[i]->SetBranchAddress( "PtbalZH", &PtbalZH );
		sample[i]->SetBranchAddress( "Mt", &Mt);
		sample[i]->SetBranchAddress( "Mte", &Mte );
		sample[i]->SetBranchAddress( "Mtmu", &Mtmu );
		sample[i]->SetBranchAddress( "DphiZMET", &DphiZMET );
		sample[i]->SetBranchAddress( "Zphi", &Zphi );
		sample[i]->SetBranchAddress( "Hphi", &Hphi );
		sample[i]->SetBranchAddress( "lep_pfCombRelIso0", &lep_pfCombRelIso0 );
		sample[i]->SetBranchAddress( "lep_pfCombRelIso1", &lep_pfCombRelIso1 );
		sample[i]->SetBranchAddress( "lep0pt", &lep0pt );
		sample[i]->SetBranchAddress( "lep1pt", &lep1pt );
		sample[i]->SetBranchAddress( "AngleEMU", &AngleEMU );
		sample[i]->SetBranchAddress( "delRemu", &delRemu );
		sample[i]->SetBranchAddress( "RMS_eta", &RMS_eta );
		sample[i]->SetBranchAddress( "UnweightedEta", &UnweightedEta );
		sample[i]->SetBranchAddress( "EtaStandDev", &EtaStandDev );
		sample[i]->SetBranchAddress( "Centrality", &Centrality );
		sample[i]->SetBranchAddress( "EvntShpAplanarity", &EvntShpAplanarity );
		sample[i]->SetBranchAddress( "EvntShpSphericity", &EvntShpSphericity );
		sample[i]->SetBranchAddress( "EvntShpIsotropy", &EvntShpIsotropy );
		sample[i]->SetBranchAddress( "EvntShpCircularity", &EvntShpCircularity );
		sample[i]->SetBranchAddress( "ZmassSVD", &ZmassSVD );
		sample[i]->SetBranchAddress( "PtbalMETH", &PtbalMETH );
		sample[i]->SetBranchAddress( "PtbalZMET", &PtbalZMET );
		sample[i]->SetBranchAddress( "DphiSecondMET", &DphiSecondMET );
		sample[i]->SetBranchAddress( "delRjj", &delRjj );
		sample[i]->SetBranchAddress( "Detaemu", &Detaemu );
		sample[i]->SetBranchAddress( "Dphiemu", &Dphiemu );
		sample[i]->SetBranchAddress( "ScalarSumPt", &ScalarSumPt );
		sample[i]->SetBranchAddress( "dPhiHMET", &dPhiHMET );
		sample[i]->SetBranchAddress( "METsig", &METsig );
		sample[i]->SetBranchAddress( "ZmassSVD", &ZmassSVD );
		sample[i]->SetBranchAddress( "ProjVisT", &ProjVisT );
		sample[i]->SetBranchAddress( "ProjMissT", &ProjMissT );
		sample[i]->SetBranchAddress( "Trigweight", &Trigweight );
		sample[i]->SetBranchAddress( "B2011PUweight", &B2011PUweight );
		sample[i]->SetBranchAddress( "A2011PUweight", &A2011PUweight );
		sample[i]->SetBranchAddress( "eventFlavor", &eventFlavor );
	}
	
	
	double lumi = 4.457;
	vector <double> lumiScale;
	lumiScale.push_back(lumi/(lumiZH115)); //0 signal
	//lumiScale.push_back(lumi/(2.28079234375000000e+05/xsecbfZH115)); //0 signal Fall11 ZH_ZToLL
	lumiScale.push_back(lumi/(lumiTT)); //1 TTJets
	lumiScale.push_back(lumi/( lumiZJL)); //2 DYM50 LF
	lumiScale.push_back(lumi/( lumiZJH)); //3 DYPtZ LF
	lumiScale.push_back(lumi/( lumiZJL)); //4 DYM50 HF
	lumiScale.push_back(lumi/( lumiZJH)); //5 DYPtZ HF
	lumiScale.push_back(lumi/( lumiTtW)); //6 T_tw
	lumiScale.push_back(lumi/( lumiTtWb)); //7 Tbar_tW
	lumiScale.push_back(lumi/( lumiWW)); //8 WW
	lumiScale.push_back(lumi/( lumiWZ)); //9 WZ
	lumiScale.push_back(lumi/( lumiZZ)); //10 ZZ 
	
	TString suffixps = "_Sig.gif";
	string directory = "/home/hep/wilken/limit/CMSSW_4_2_8/src/BestBDTforward/";
	TString plot = directory+"somereallylongnameherehopethisarrayisbiggenough"+suffixps;	
	
	float minLimit = 9999.99, corrCut = -99.99, significance, SoB, GlobalMin = 9999.99;
	
	vector <string> SampleName;
	SampleName.push_back("VH"); //0 signal
	SampleName.push_back("TT"); //1 TTJets
	SampleName.push_back("ZjLF"); //2 DYM50 LF
	SampleName.push_back("DYPtZLF"); //3 DYPtZ LF
	SampleName.push_back("ZjHF"); //4 DYM50 HF
	SampleName.push_back("DYPtZHF"); //5 DYPtZ HF
	SampleName.push_back("s_Top"); //6 T_tw
	SampleName.push_back("Tbar_tW"); //7 Tbar_tW
	SampleName.push_back("WW"); //8 WW
	SampleName.push_back("WZ"); //9 WZ
	SampleName.push_back("VV"); //10 ZZ 
	
	TCanvas *c1 = new TCanvas("c1","");
	gStyle->SetOptStat("kTRUE");
	c1->SetFillColor(10);
	c1->SetFillColor(10);	
		
	float BestCutValue[] = {80, //0	Higgs mass min
		280, //1	Higgs mass max
		0.4, //2	CSV0
		6, //3		Z mass min
		155, //4		Z mass max
	0.45, //5 Delta Phi (Z,MET)
	0.2, //6 Charged Hadron Fraction
	 2.612, //7 Delta Phi (H,V)
	 10, //8 Pzeta Cut
	 95, //9 Mt electron
105};//10 Mt muon
	
	string VarNames[] = {"HmassMin","HmassMax","CSV0","ZmassMin","ZmassMax","DeltaPhiZMET","CHFb0","DeltaPhiHV","Pzeta","Mte","Mtmu"};
	float cutmin[] = {60, 140, 0, 0, 90, 0.35, 0.1, 0.6, -20, 55, 55};
	float cutmax[] = {110, 290, 0.5, 20, 190, 1.35, 0.3, 3.1, 30, 155, 155};
//	float csvsteps[] = {0.1, 0.244, 0.5, 0.679, 0.898};
	float csvsteps[] = {0.0, 0.122, 0.244, 0.372, 0.5};
	
	float range = -99.99; 
	float cut = -99.99;
	double BDTvalue = -99.99;
	
	//	for (int Var_iter = 10; Var_iter>-1;Var_iter--){
	for (int Var_iter = 0; Var_iter<11;Var_iter++){
		int steps = 10;
		range = cutmax[Var_iter]-cutmin[Var_iter];
		float stepsize = range/float(steps);
//		if(Var_iter==2) steps = 5;
	
		vector <int> count;
		count.push_back(0); //0 signal
		count.push_back(0); //1 TTJets
		count.push_back(0); //2 DYM50 LF
		count.push_back(0); //3 DYPtZ LF
		count.push_back(0); //4 DYM50 HF
		count.push_back(0); //5 DYPtZ HF
		count.push_back(0); //6 T_tw
		count.push_back(0); //7 Tbar_tW
		count.push_back(0); //8 WW
		count.push_back(0); //9 WZ
		count.push_back(0); //10 ZZ 
		
	
	int nbins = 18;
	float bl=-1.0, bh=-0.2;
	TH1F* EventsBDT[sample.size()];
	TH1F* h_stat_up[sample.size()];
	TH1F* h_stat_down[sample.size()];
	for(unsigned int j=0;j<sample.size();j++){
		EventsBDT[j]=new TH1F(TString::Format("%s",SampleName[j].c_str()).Data(),
							  TString::Format("BDT Value %s",SampleName[j].c_str()).Data(),
							  nbins,bl,bh);
		h_stat_up[j]=new TH1F(TString::Format("%s_CMSstatUp",SampleName[j].c_str()).Data(),
							  TString::Format("BDT %s plus sigma",SampleName[j].c_str()).Data(),
							  nbins,bl,bh);
		h_stat_down[j]=new TH1F(TString::Format("%s_CMSstatDown",SampleName[j].c_str()).Data(),
								TString::Format("BDT %s minus sigma",SampleName[j].c_str()).Data(),
								nbins,bl,bh);
	}
	
	TH1F* data_obs= new TH1F	("data_obs", "null histo for data", nbins,bl,bh);
	
	float PUweight2011=1.0;
	TGraph *Limits= new TGraph();
	TGraph *Punzi= new TGraph();
	TGraph *SigOverSqrtB= new TGraph();
	
 for (int icut = 0; icut<steps; icut++){
			 cut = icut*stepsize+cutmin[Var_iter];
			//if(Var_iter==2) cut = csvsteps[icut];
			string filename;
			
			filename = TString::Format("HistoForShapeLimit%0.3f.root",cut).Data();
			TFile *theHistogramFile= new TFile(TString::Format("%s%s",directory.c_str(),filename.c_str()).Data(), "RECREATE", "histogram file",0);
			theHistogramFile->cd();
			
			//clear root file
			for(unsigned int j =0; j<sample.size(); j++){
				EventsBDT[j]->Reset("ICES");
				h_stat_up[j]->Reset("ICES");
				h_stat_down[j]->Reset("ICES");
				
				EventsBDT[j]->Sumw2();
				h_stat_up[j]->Sumw2();
				h_stat_down[j]->Sumw2();
			}
			for(unsigned int isample = 0; isample<count.size(); isample++){
				count[isample]=0;
				for (Long64_t ievt=0; ievt<sample[isample]->GetEntries();ievt++) {
					sample[isample]->GetEntry(ievt);
					float Pzeta = ProjMissT-(0.25*ProjVisT) ;
					if((lep0pt > 20 || lep1pt >20) && (lep0pt > 10 && lep1pt >10)&&delRemu > 0.4){
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
						PUweight2011 = (2.219*A2011PUweight + 2.238*B2011PUweight)/4.457;
						switch (Var_iter){
							case 0:
								if (Hmass < BestCutValue[1] && CSV0 > BestCutValue[2] && Emumass > BestCutValue[3] && Emumass < BestCutValue[4] 
								&& fabs(DphiZMET)  < BestCutValue[5] && jetCHF0 > BestCutValue[6] && DeltaPhiHV > BestCutValue[7] 
								&& Pzeta > BestCutValue[8] && Mte < BestCutValue[9] && Mtmu < BestCutValue[10] ){
									if(Hmass>cut) {
										EventsBDT[isample]->Fill(BDTvalue,lumiScale[isample]*PUweight2011);
										count[isample]++;
									}}
								break;																								
							case 1:
								if (Hmass > BestCutValue[0] && CSV0 > BestCutValue[2] && Emumass > BestCutValue[3] && Emumass < BestCutValue[4] 
									&& fabs(DphiZMET)  < BestCutValue[5] && jetCHF0 > BestCutValue[6] && DeltaPhiHV > BestCutValue[7] 
									&& Pzeta > BestCutValue[8] && Mte < BestCutValue[9] && Mtmu < BestCutValue[10] ){
									if(Hmass<cut) {
										EventsBDT[isample]->Fill(BDTvalue,lumiScale[isample]*PUweight2011);
										count[isample]++;
									}}
								break;																
							case 2:
								if (Hmass > BestCutValue[0] && Hmass < BestCutValue[1] && Emumass > BestCutValue[3] && Emumass < BestCutValue[4] 
									&& fabs(DphiZMET)  < BestCutValue[5] && jetCHF0 > BestCutValue[6] && DeltaPhiHV > BestCutValue[7] 
									&& Pzeta > BestCutValue[8] && Mte < BestCutValue[9] && Mtmu < BestCutValue[10] ){
									if(CSV0 >cut) {
										EventsBDT[isample]->Fill(BDTvalue,lumiScale[isample]*PUweight2011);
										count[isample]++;
									}}
								break;								
							case 3:
								if (Hmass > BestCutValue[0] && Hmass < BestCutValue[1] && CSV0 > BestCutValue[2] && Emumass < BestCutValue[4] 
									&& fabs(DphiZMET)  < BestCutValue[5] && jetCHF0 > BestCutValue[6] && DeltaPhiHV > BestCutValue[7] 
									&& Pzeta > BestCutValue[8] && Mte < BestCutValue[9] && Mtmu < BestCutValue[10] ){
									if(Emumass >cut) {
										EventsBDT[isample]->Fill(BDTvalue,lumiScale[isample]*PUweight2011);
										count[isample]++;
									}}
								break;								
							case 4:
								if (Hmass > BestCutValue[0] && Hmass < BestCutValue[1] && CSV0 > BestCutValue[2] && Emumass > BestCutValue[3] 
									&& fabs(DphiZMET)  < BestCutValue[5] && jetCHF0 > BestCutValue[6] && DeltaPhiHV > BestCutValue[7] 
									&& Pzeta > BestCutValue[8] && Mte < BestCutValue[9] && Mtmu < BestCutValue[10] ){
									if(Emumass <cut) {
										EventsBDT[isample]->Fill(BDTvalue,lumiScale[isample]*PUweight2011);
										count[isample]++;
									}}
								break;
							case 5:
								if (Hmass > BestCutValue[0] && Hmass < BestCutValue[1] && CSV0 > BestCutValue[2] && Emumass > BestCutValue[3] 
									&& Emumass < BestCutValue[4] && jetCHF0 > BestCutValue[6] && DeltaPhiHV > BestCutValue[7] 
									&& Pzeta > BestCutValue[8] && Mte < BestCutValue[9] && Mtmu < BestCutValue[10] ){
									if(fabs(DphiZMET) <cut) {
										EventsBDT[isample]->Fill(BDTvalue,lumiScale[isample]*PUweight2011);
										count[isample]++;
									}}
								break;
							case 6:
								if (Hmass > BestCutValue[0] && Hmass < BestCutValue[1] && CSV0 > BestCutValue[2] && Emumass > BestCutValue[3] 
								&& Emumass < BestCutValue[4] && fabs(DphiZMET)  < BestCutValue[5] && DeltaPhiHV > BestCutValue[7] 
									&& Pzeta > BestCutValue[8] && Mte < BestCutValue[9] && Mtmu < BestCutValue[10] ){
									if(jetCHF0 > cut) {
										EventsBDT[isample]->Fill(BDTvalue,lumiScale[isample]*PUweight2011);
										count[isample]++;
									}}
								break;
							case 7:
								if (Hmass > BestCutValue[0] && Hmass < BestCutValue[1] && CSV0 > BestCutValue[2] && Emumass > BestCutValue[3] 
								&& Emumass < BestCutValue[4] && fabs(DphiZMET)  < BestCutValue[5] && jetCHF0 > BestCutValue[6]
									&& Pzeta > BestCutValue[8] && Mte < BestCutValue[9] && Mtmu < BestCutValue[10] ){
									if(DeltaPhiHV>cut) {
										EventsBDT[isample]->Fill(BDTvalue,lumiScale[isample]*PUweight2011);
										count[isample]++;
									}}
								break;
							case 8:
								if (Hmass > BestCutValue[0] && Hmass < BestCutValue[1] && CSV0 > BestCutValue[2] && Emumass > BestCutValue[3] 
								&& Emumass < BestCutValue[4] && fabs(DphiZMET)  < BestCutValue[5] && jetCHF0 > BestCutValue[6] 
								&& DeltaPhiHV > BestCutValue[7] && Mte < BestCutValue[9] && Mtmu < BestCutValue[10] ){
									if(Pzeta>cut) {
										EventsBDT[isample]->Fill(BDTvalue,lumiScale[isample]*PUweight2011);
										count[isample]++;
									}}
								break;
							case 9:
								if (Hmass > BestCutValue[0] && Hmass < BestCutValue[1] && CSV0 > BestCutValue[2] && Emumass > BestCutValue[3] 
								&& Emumass < BestCutValue[4] && fabs(DphiZMET)  < BestCutValue[5] && jetCHF0 > BestCutValue[6] 
								&& DeltaPhiHV > BestCutValue[7] && Pzeta > BestCutValue[8] && Mtmu < BestCutValue[10] ){
									if(Mte < cut) {
										EventsBDT[isample]->Fill(BDTvalue,lumiScale[isample]*PUweight2011);
										count[isample]++;
									}}
								break;
							case 10:
								if (Hmass > BestCutValue[0] && Hmass < BestCutValue[1] && CSV0 > BestCutValue[2] && Emumass > BestCutValue[3] 
									&& Emumass < BestCutValue[4] && fabs(DphiZMET)  < BestCutValue[5] && jetCHF0 > BestCutValue[6] 
									&& DeltaPhiHV > BestCutValue[7] && Pzeta > BestCutValue[8] && Mte < BestCutValue[9] ){
									if(Mtmu<cut) {
										EventsBDT[isample]->Fill(BDTvalue,lumiScale[isample]*PUweight2011);
										count[isample]++;
									}}
								break;
							default: 
								cout << "BUG IN SWITCH" << endl;
						}//switch
					}//signal trigger requirement
				}//end signal event loop
			}//end loop on all samples
			
			
			cout << "Cut " <<cut << " Number of Events Signal: " <<count[0] << " TTJets " << count[1]<< " DY M50 LF " << count[2]<< endl;
			
			float background = 0.00;
			//		background = (NTTJ*lumiScale[1])+(NZJ*lumiScale[2])+(NZJ2*lumiScale[3]);
			for (unsigned int k = 1; k<lumiScale.size(); k++) background=background+(lumiScale[k]*count[k]);
			significance = (count[0]*lumiScale[0])/(1.5+sqrt(background)+.2*background);
			SoB = (count[0]*lumiScale[0])/background;
			//Calculate limit here using combine function
			
			//		VH->Scale(ZH_M115_weight);
			//		TT->Scale(TTJets_weight);
			//		ZjLF->Scale(DYJetsToLL_M50_weight);
			
			for(unsigned int j =0; j<sample.size(); j++){
				for (int i = 1; i<(nbins+1); i++){
					double binContent = -99.99;
					double sigma =-99.99;
					binContent = EventsBDT[j]->GetBinContent(i);
					sigma = EventsBDT[j]->GetBinError(i);
					h_stat_up[j]->SetBinContent(i,binContent+sigma);
					//if((binContent-sigma) < 0) sigma  = 0;
					h_stat_down[j]->SetBinContent(i,binContent-sigma);
				}//end loop over bins
			}//end loop over samples
			
			EventsBDT[2]->Add(EventsBDT[3]);
			h_stat_up[2]->Add(h_stat_up[3]);
			h_stat_down[2]->Add(h_stat_down[3]);
			//Heavy Flavor	
			EventsBDT[4]->Add(EventsBDT[5]);
			h_stat_up[4]->Add(h_stat_up[5]);
			h_stat_down[4]->Add(h_stat_down[5]);
			//Single Top
			EventsBDT[6]->Add(EventsBDT[7]);
			h_stat_up[6]->Add(h_stat_up[7]);
			h_stat_down[6]->Add(h_stat_down[7]);
			//DiBoson
			EventsBDT[10]->Add(EventsBDT[8]);
			h_stat_up[10]->Add(h_stat_up[8]);
			h_stat_down[10]->Add(h_stat_down[8]);
			EventsBDT[10]->Add(EventsBDT[9]);
			h_stat_up[10]->Add(h_stat_up[9]);
			h_stat_down[10]->Add(h_stat_down[9]);
			
			
			makeDataCard(EventsBDT[0]->Integral(), EventsBDT[2]->Integral(), EventsBDT[4]->Integral(), EventsBDT[1]->Integral(), EventsBDT[6]->Integral(), EventsBDT[10]->Integral(), cut, directory, filename);

			
			for(unsigned int k =0; k<sample.size(); k++){
				if(k==3||k==5||k==7||k==8||k==9) continue;
				EventsBDT[k]->Write();
				h_stat_up[k]->Write();
				h_stat_down[k]->Write();
			}
			data_obs->Write();

			
			THStack *histBdt_BkgStack = new THStack("histBdt_BkgStack","Stacked Background BDT");
			EventsBDT[2]->SetFillColor(kYellow);
			EventsBDT[1]->SetFillColor(kBlue);
			histBdt_BkgStack->Add(EventsBDT[2]);
			histBdt_BkgStack->Add(EventsBDT[1]);
			histBdt_BkgStack->Draw("hist");
			EventsBDT[0]->SetLineWidth(3);
			EventsBDT[0]->Scale(100);
			EventsBDT[0]->Draw("same hist");
			TLegend myLegend(0.7, 0.5, 0.89, 0.8);
			myLegend.SetTextSize(0.03);
			myLegend.AddEntry(EventsBDT[0], "ZH_M115 x 100", "l");	
			myLegend.AddEntry(EventsBDT[1], "DYJetsToLL LF", "f");	
			myLegend.AddEntry(EventsBDT[2], "TTJets", "f");	
			myLegend.Draw();		
			//plot = directory+"DphiZMET"+cut.c_str()+suffixps;
			plot = TString::Format("%sBDT%0.02f.gif",directory.c_str(),cut).Data();
			c1->Print(plot);
			c1->Clear();	
			
			theHistogramFile->Close();		
			theHistogramFile->Delete();		
			
			
			string combiner = TString::Format("combine -M ProfileLikelihood %sShapeDataCard%0.03fcut.txt -m 115 -t 1000 | tee %soutput%0.03fcut.log",directory.c_str(),cut,directory.c_str(),cut).Data();
			gSystem->Exec(combiner.c_str());
			string GetMedianLimit  = TString::Format("grep 'median expected' %soutput%0.03fcut.log  | gawk '{print $6}' > tmpMedianLimit.out",directory.c_str(),cut).Data();
			gSystem->Exec(GetMedianLimit.c_str());
			ifstream FileStreamForLimit("tmpMedianLimit.out");
			double medianExpLim = 99999.99;
			FileStreamForLimit>>medianExpLim;
			cout << "limit is " <<medianExpLim << " at " << cut  << " S/B " << SoB << " for " << VarNames[Var_iter]<< endl;
			FileStreamForLimit.close();
			
			cout << "Seting Point in Limits " << icut<<","<<cut<<","<<medianExpLim<<endl;
			Limits->SetPoint(icut,cut,medianExpLim);
			SigOverSqrtB->SetPoint(icut,cut,35000*(count[0]*lumiScale[0])/sqrt(background));
			Punzi->SetPoint(icut,cut,significance*200000);
			
			if (medianExpLim<minLimit) {
				minLimit = medianExpLim;
				corrCut = cut;
			}			
			
			
		}//icut for loop
		
		
		plot=directory+"LimitvsCut"+VarNames[Var_iter]+".root";
		TFile *TGraphFile= new TFile(plot, "RECREATE");
		TGraphFile->cd();
		double x, y;
		Limits->GetPoint(0,x,y);
		cout << "Geting point in Limits " <<  x <<","<< y << endl;
		Limits->Draw("AC*");
			plot=directory+"LimitvsCut"+VarNames[Var_iter]+suffixps;
			c1->Print(plot);
			c1->Clear();
			
//		SigOverSqrtB->Draw("AC*");
//		Punzi->Draw("AC*");
		
		Limits->Write();
		SigOverSqrtB->Write();
		Punzi->Write();
		
		
		TMultiGraph *mg = new TMultiGraph();
		mg->SetTitle("Punzi, Signal/SQRT(B), Limit");
			Limits->SetMarkerStyle(kStar);
		SigOverSqrtB->SetLineColor(kBlue);
		SigOverSqrtB->SetMarkerColor(kBlue);
		SigOverSqrtB->SetMarkerStyle(kOpenTriangleUp);
		Punzi->SetLineColor(kRed);
		Punzi->SetMarkerColor(kRed);
		Punzi->SetMarkerStyle(kOpenSquare);
		mg->Add(SigOverSqrtB);
		mg->Add(Punzi);	
		mg->Add(Limits);	
		mg->Draw("ACP");	
		//		c1->BuildLegend();
		plot=directory+"PunziSosqrtBLimCompare"+VarNames[Var_iter]+suffixps;
		c1->Print(plot);
		
		
		TGraphFile->Delete();
		mg->Delete();
		for(unsigned int k =0; k<sample.size(); k++){
			EventsBDT[k]->Delete();
			h_stat_up[k]->Delete();
			h_stat_down[k]->Delete();
		}
		data_obs->Delete();
		
		
		cout << "The best limit is " <<minLimit << " at " << corrCut<< " for " << VarNames[Var_iter]<<endl;
			if (minLimit<GlobalMin) BestCutValue[Var_iter] = corrCut;
		corrCut = -99.99;
		minLimit = 9999.99;
	}//loop on variables
	
	
	c1->Close();
	
	cout << "Ending Cut Values: " << endl;
	cout <<  VarNames[0] << " " << BestCutValue[0] << endl;
	cout <<  VarNames[1] << " " << BestCutValue[1] << endl;
	cout <<  VarNames[2] << " " << BestCutValue[2] << endl;
	cout <<  VarNames[3] << " " << BestCutValue[3] << endl;
	cout <<  VarNames[4] << " " << BestCutValue[4] << endl;
	cout <<  VarNames[5] << " " << BestCutValue[5] << endl;
	cout <<  VarNames[6] << " " << BestCutValue[6] << endl;
	cout <<  VarNames[7] << " " << BestCutValue[7] << endl;
	cout <<  VarNames[8] << " " << BestCutValue[8] << endl;
	cout <<  VarNames[9] << " " << BestCutValue[9] << endl;
	cout <<  VarNames[10] << " " << BestCutValue[10] << endl;
	
	
	
}

void makeDataCard(float VHrate, float ZjLFrate, float ZjHFrate, float TTrate, float sToprate, float VVrate, float Thiscut, string dir, string shapefile){
	ofstream myfile (TString::Format("%s/ShapeDataCard%0.03fcut.txt",dir.c_str(),Thiscut).Data());
	myfile<<"imax 1  number of channels"<<endl;
	myfile<<"jmax 5  number of backgrounds"<<endl;;
	myfile<<"kmax *  number of nuisance parameters (sources of systematical uncertainties)"<<endl;
	myfile<<TString::Format("shapes * * %s%s $PROCESS $PROCESS_$SYSTEMATIC",dir.c_str(),shapefile.c_str()).Data()<<endl;
	myfile<<"bin Zemu"<<endl;;
	myfile<<"observation 0"<<endl;
	myfile<<"bin \t Zemu \t Zemu \t Zemu  \t Zemu \t Zemu \t Zemu"<<endl;
	myfile<<"process \t VH \t ZjLF \t ZjHF \t TT \t s_Top \t VV"<<endl;
	myfile<<"process  \t 0 \t 1 \t 2  \t 3 \t 4  \t 5 "<<endl;
	myfile<<TString::Format("rate \t %0.3f \t %0.3f \t %0.3f \t %0.3f \t %0.3f \t %0.3f",VHrate,ZjLFrate,ZjHFrate,TTrate,sToprate,VVrate).Data()<<endl;
	myfile<<"\n";
	myfile<<"lumi \t lnN \t 1.022 \t - \t - \t - \t 1.022 \t 1.022"<<endl; 
	myfile<<"pdf_qqbar \t lnN \t 1.01 \t - \t - \t - \t - \t 1.01"<<endl; 
	myfile<<"pdf_gg \t lnN \t - \t - \t - \t - \t 1.01 \t -"<<endl; 
	myfile<<"QCDscale_VH \t lnN \t 1.04 \t - \t - \t - \t - \t -"<<endl; 
	myfile<<"QCDscale_VV \t lnN \t - \t - \t - \t - \t - \t 1.04"<<endl; 
	myfile<<"QCDscale_ttbar \t lnN \t - \t - \t - \t - \t 1.06 \t -"<<endl; 
	myfile<<"CMS_vhbb_ST \t lnN  \t - \t -  \t - \t - \t 1.29 \t -"<<endl; 
	myfile<<"CMS_vhbb_VV \t lnN  \t - \t -  \t - \t - \t - \t 1.3"<<endl; 
	myfile<<"CMS_vhbb_TT_SF \t lnN  \t - \t 1.002  \t 0.974 \t 1.123 \t - \t -"<<endl; 
	myfile<<"CMS_vhbb_TT_ex \tlnN \t - \t - \t 1.05 \t - \t - \t -"<<endl;
	myfile<<"CMS_vhbb_ZjHF_SF \t lnN \t - \t 0.873 \t 1.202 \t 0.957 \t - \t -"<<endl; 
	myfile<<"CMS_vhbb_ZjHF_ex \t lnN \t - \t - \t 1.05 \t - \t - \t -" << endl; 
	myfile<<"CMS_vhbb_ZjLF_SF \t lnN \t - \t 1.198 \t 0.875 \t 1.003 \t - \t -"<<endl; 
	myfile<<"CMS_vhbb_ZjLF_ex \t lnN \t - \t 1.05 \t - \t - \t - \t -" << endl; 
	myfile<<"CMS_trigger_m  \t lnN  \t 1.01  \t - \t - \t - \t - \t -"<<endl; 
	myfile<<"CMS_trigger_3  \t lnN  \t 1.02  \t - \t - \t - \t - \t -"<<endl; 
	myfile<<"CMS_eff_m  \t lnN  \t 1.04 \t - \t - \t - \t 1.04 \t 1.04"<<endl; 
	myfile<<"CMS_eff_e  \t lnN  \t 1.04 \t - \t - \t - \t 1.04 \t 1.04"<<endl; 
	myfile<<"CMS_eff_b  \t lnN  \t 1.11 \t 1.07 \t 1 \t 1 \t 1.05 \t 1.11"<<endl; 
	myfile<<"CMS_fake_b  \t lnN  \t 1.05 \t 1.12 \t 1 \t 1 \t 1.15 \t 1"<<endl; 
	myfile<<"CMS_scale_j  \t lnN  \t 1.03 \t -  \t -  \t - \t 1.03 \t 1.03"<<endl;  
	myfile<<"CMS_res_j  \t lnN  \t 1.03 \t 1.04 \t 1.04 \t 1.04 \t 1.04 \t 1.05"<<endl;  
	myfile<<"CMSstat \t shape \t 1.0 \t - \t - \t - \t - \t -"<< endl;
	myfile<<"CMSstat \t shape \t - \t 1.0 \t - \t - \t - \t -"<< endl;
	myfile<<"CMSstat \t shape \t - \t - \t 1.0 \t - \t - \t -"<< endl;
	myfile<<"CMSstat \t shape \t - \t - \t - \t 1.0 \t - \t -"<< endl;
	myfile<<"CMSstat \t shape \t - \t - \t - \t - \t 1.0 \t -"<< endl;
	myfile<<"CMSstat \t shape \t - \t - \t - \t - \t - \t 1.0"<< endl;
	
	myfile.close();
	
}