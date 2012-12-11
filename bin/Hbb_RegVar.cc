/**
 @brief Example of an analysis code to read tree ntuple
 and create histogram
 Usage:
 
 \verbatim
 mkdir [directoryName]
 Hbb_RegVar <inputfile> [outputfilename] [directoryName]
 \endverbatim
 
 @param inputfile Either a ROOT file or an ASCII file containing list of
 ROOT files.
 
 @param outputfile Name of the ROOT file which contains the histogram.
 Defaulted to 'output.root'
 
 @author Rachel Wilken <rachel.wilken@cern.ch>
 
 @date Fri Nov 9 2011
 
 */

#include "Hbb_TauemuTuples.h"
#include <iostream>
#include <fstream>

using namespace std;

double deltaPhi(double phi1, double phi2) 
{ 
    double PI = 3.14159265;
    double result = phi1 - phi2;
    while (result > PI) result -= 2*PI;
    while (result <= -PI) result += 2*PI;
    return abs(result);
}

double evalEt( double pt, double eta, double phi, double e){
	TLorentzVector j;
	j.SetPtEtaPhiE(pt,eta,phi, e );
	return j.Et(); 
	
}


double evalMt( double pt, double eta, double phi, double e){
	TLorentzVector j;
	j.SetPtEtaPhiE(pt,eta,phi, e );
	return j.Mt(); 
	
}

float min(float csv1, float  csv2)
{
	if(csv1 > csv2) return csv2;
	else return csv1;
}

float resolutionBias(float eta)
{
	// return 0;//Nominal!
	if(eta< 1.1) return 0.05;
	if(eta< 2.5) return 0.10;
	if(eta< 5) return 0.30;
	return 0;
}

double evalJERBias( double ptreco, double ptgen, double eta){
	double cor =1;   
	if ((fabs(ptreco - ptgen)/ ptreco)<0.5) { //Limit the effect to the core 
		cor = (ptreco +resolutionBias(eta) *(ptreco-ptgen))/ptreco;   
	}
	return ptreco*cor;
}



int main(int argc, char** argv) {
	
	bool debug = false;	
	
    // There must be at least one argument, the input file.
    if (argc < 2) {
        std::cout << "Input filename is not specified! Exiting!" << std::endl;
        exit(1);
    }
	
    // First argument: a ROOT file or an ASCII file containing list of
    // ROOT files.
    std::string ifilename(argv[1]);
    std::string ofilename("output.root");
    std::string directory("plots");
	
    if (argc == 3) {
        // If there is a second argument, use it for the output filename.
        ofilename = std::string(argv[2]);
    }
	if (argc == 4) {
        // If there is a third argument, use it for the plots drectory
		ofilename = std::string(argv[2]);
		directory = std::string(argv[3]);
	}
	
    // Create the treeReader object.
   treeReader sample(ifilename, std::string("tree"));
	
    // If the input file(s) doesn't contain any event, exit.
    if (!(sample.readEvent(0))) {
        return 0;
    }
	
    // Set the sumw2 option for all histogram
    TH1::SetDefaultSumw2(true);
	
	BTagShapeNew* btagNew;
	btagNew = new BTagShapeNew("csvdiscrNew.root");
	btagNew->computeFunctions();
	if (debug){
		cout << "Done computing btagSF " << endl;
		cout << "b quark SF for CVSL: " << btagNew->ib->Eval(0.244) << "   CVSM: " << btagNew->ib->Eval(0.5) << "  CSVT: " << btagNew->ib->Eval(0.898) << endl;
		cout << "c quark SF for CVSL: " << btagNew->ic->Eval(0.244) << "   CVSM: " << btagNew->ic->Eval(0.5) << "  CSVT: " << btagNew->ic->Eval(0.898) << endl;
		cout << "light quark SF for CVSL: " << btagNew->il->Eval(0.244) << "   CVSM: " << btagNew->il->Eval(0.5) << "  CSVT: " << btagNew->il->Eval(0.898) << endl << endl;
	}	
	
    // Create the output ROOT file
    TFile ofile(ofilename.c_str(),"RECREATE");
    ofile.cd();
	LumiWeight = SetWeight(ifilename);
	bool isZjets = false;
	bool isM50sample = false;
	isdata = false;
	isDATA = false;
	bool LFfile = false;
	bool HFfile = false;
	if (findString(ifilename, "DY")) isZjets = true;
	if (findString(ifilename, "M50")) isM50sample = true;
	if (findString(ifilename, "M-50")) isM50sample = true;
	if (findString(ifilename, "M_50")) isM50sample = true;
	if (findString(ifilename, "ZJets")) isZjets = true;
	if (findString(ifilename, "data")) {
		isdata = true;
		isDATA = true;
	}
	if (findString(ifilename, "Data")) {
		isdata = true;
		isDATA = true;
	}
	if (findString(ifilename, "Run")) {
		isdata = true;
		isDATA = true;
	}
	if (findString(ofilename, "LF")) LFfile = true;
	if (findString(ofilename, "HF")) HFfile = true;
	
	
    // Check how many events and file(s) are analyzed
    long int nEvent     = sample.nEvent();
    size_t   nFile      = sample.nFile();
	
    std::cout << "There are " << nEvent << " events in the chain."
	<< std::endl ;
    std::cout << "There are " << nFile  << " files in the chain."
	<< std::endl ;
    std::cout << "The tree name  is: " << sample.treeName()  << std::endl;
    std::cout << "The tree title is: " << sample.treeTitle() << std::endl;
	
	TTree *FOM_tree = new TTree("FOM_tree","Tree for Significance Optimization");
	FOM_tree->Branch("nJets",&nJets, "nJets/I");
	FOM_tree->Branch("Naj",&Naj, "Naj/I");
	FOM_tree->Branch("Nab",&Nab, "Nab/I");
	FOM_tree->Branch("naJets",&naJets, "naJets/I");
	FOM_tree->Branch("eventFlavor",&eventFlavor, "eventFlavor/I");
	FOM_tree->Branch("CSV0",&CSV0, "CSV0/F");
	FOM_tree->Branch("CSV1",&CSV1, "CSV1/F");
	FOM_tree->Branch("Emumass",&Emumass, "Emumass/F");
	FOM_tree->Branch("Hmass",&Hmass, "Hmass/F");
	FOM_tree->Branch("oldHmass",&oldHmass, "oldHmass/F");
	FOM_tree->Branch("RegHmass",&RegHmass, "RegHmass/F");
	FOM_tree->Branch("DeltaPhiHV",&DeltaPhiHV, "DeltaPhiHV/F");
	FOM_tree->Branch("Hpt",&Hpt, "Hpt/F");
	FOM_tree->Branch("Zpt",&Zpt, "Zpt/F");
	FOM_tree->Branch("lep0pt",&lep0pt, "lep0pt/F");
	FOM_tree->Branch("ScalarSumPt",&ScalarSumPt, "ScalarSumPt/F");
	FOM_tree->Branch("ScalarSumJetPt",&ScalarSumJetPt, "ScalarSumJetPt/F");
	FOM_tree->Branch("ScalarSumHiggsJetPt",&ScalarSumHiggsJetPt, "ScalarSumHiggsJetPt/F");
	FOM_tree->Branch("Ht",&Ht, "Ht/F");
	FOM_tree->Branch("EtaStandDev",&EtaStandDev, "EtaStandDev/F");
	FOM_tree->Branch("UnweightedEta",&UnweightedEta, "UnweightedEta/F");
	FOM_tree->Branch("EvntShpCircularity",&EvntShpCircularity, "EvntShpCircularity/F");
	FOM_tree->Branch("alpha_j",&alpha_j, "alpha_j/F");
	FOM_tree->Branch("qtb1",&qtb1, "qtb1/F");
	FOM_tree->Branch("nSV",&nSV, "nSV/I");
	FOM_tree->Branch("Trigweight",&Trigweight, "Trigweight/F");
	FOM_tree->Branch("B2011PUweight",&B2011PUweight, "B2011PUweight/F");
	FOM_tree->Branch("A2011PUweight",&A2011PUweight, "A2011PUweight/F");
	FOM_tree->Branch("btag2CSF",&btag2CSF, "btag2CSF/F");
	FOM_tree->Branch("DetaJJ",&DetaJJ, "DetaJJ/F");
	FOM_tree->Branch("jetCHF0",&jetCHF[0], "jetCHF0/F");
	FOM_tree->Branch("jetCHF1",&jetCHF[1], "jetCHF1/F");
	FOM_tree->Branch("jetPt0",&jetPt[0], "jetPt0/F");
	FOM_tree->Branch("jetPt1",&jetPt[1], "jetPt1/F");
	FOM_tree->Branch("jetEta0",&jetEta[0], "jetEta0/F");
	FOM_tree->Branch("jetEta1",&jetEta[1], "jetEta1/F");
	FOM_tree->Branch("CSVNewShape0",&CSVNewShape[0], "CSVNewShape0/F");
	FOM_tree->Branch("CSVNewShape1",&CSVNewShape[1], "CSVNewShape1/F");
	FOM_tree->Branch("lep1pt",&leptonPt[1], "lep1pt/F");
	FOM_tree->Branch("lep_pfCombRelIso0",&lep_pfCombRelIso[0], "lep_pfCombRelIso0/F");
	FOM_tree->Branch("lep_pfCombRelIso1",&lep_pfCombRelIso[1], "lep_pfCombRelIso1/F");
	FOM_tree->Branch("DphiJJ",&DphiJJ, "DphiJJ/F");
	FOM_tree->Branch("RMS_eta",&RMS_eta, "RMS_eta/F");
	FOM_tree->Branch("PtbalZH",&PtbalZH, "PtbalZH/F");
	FOM_tree->Branch("EventPt",&EventPt, "EventPt/F");
	FOM_tree->Branch("EventMass",&EventMass, "EventMass/F");
	FOM_tree->Branch("AngleHemu",&AngleHemu, "AngleHemu/F");
	FOM_tree->Branch("Centrality",&Centrality, "Centrality/F");
	FOM_tree->Branch("MET",&MET, "MET/F");
	FOM_tree->Branch("EvntShpAplanarity",&EvntShpAplanarity, "EvntShpAplanarity/F");
	FOM_tree->Branch("EvntShpSphericity",&EvntShpSphericity, "EvntShpSphericity/F");
	FOM_tree->Branch("EvntShpIsotropy",&EvntShpIsotropy, "EvntShpIsotropy/F");
	FOM_tree->Branch("Zphi",&Zphi, "Zphi/F");
	FOM_tree->Branch("Hphi",&Hphi, "Hphi/F");
	FOM_tree->Branch("SV_mass",&SV_mass, "SV_mass/F");
	FOM_tree->Branch("Mte",&Mte, "Mte/F");
	FOM_tree->Branch("Mtmu",&Mtmu, "Mtmu/F");
	FOM_tree->Branch("delPullAngle",&delPullAngle, "delPullAngle/F");
	FOM_tree->Branch("delPullAngle2",&delPullAngle2, "delPullAngle2/F");
	FOM_tree->Branch("Mt",&Mt, "Mt/F");
	FOM_tree->Branch("dPhiHMET",&dPhiHMET, "dPhiHMET/F");
	FOM_tree->Branch("DeltaPhijetMETmin",&DeltaPhijetMETmin, "DeltaPhijetMETmin/F");
	FOM_tree->Branch("DeltaPhijetMETZtaumin",&DeltaPhijetMETZtaumin, "DeltaPhijetMETZtaumin/F");
	FOM_tree->Branch("AngleEMU",&AngleEMU, "AngleEMU/F");
	FOM_tree->Branch("AaronEleMissE",&AaronEleMissE, "AaronEleMissE/F");
	FOM_tree->Branch("AaronMuMissE",&AaronMuMissE, "AaronMuMissE/F");
	FOM_tree->Branch("Dphiemu",&Dphiemu, "Dphiemu/F");
	FOM_tree->Branch("delRjj",&delRjj, "delRjj/F");
	FOM_tree->Branch("Detaemu",&Detaemu, "Detaemu/F");
	FOM_tree->Branch("EleMissE",&EleMissE, "EleMissE/F");
	FOM_tree->Branch("DphiEleMET",&DphiEleMET, "DphiEleMET/F");
	FOM_tree->Branch("dphiMuMET",&dphiMuMET, "dphiMuMET/F");
	FOM_tree->Branch("PtbalMETH",&PtbalMETH, "PtbalMETH/F");
	FOM_tree->Branch("topPt",&topPt, "topPt/F");
	FOM_tree->Branch("MassEleb0",&MassEleb0, "MassEleb0/F");
	FOM_tree->Branch("MassMub0",&MassMub0, "MassMub0/F");
	FOM_tree->Branch("MassEleb1",&MassEleb1, "MassEleb1/F");
	FOM_tree->Branch("MassMub1",&MassMub1, "MassMub1/F");
	FOM_tree->Branch("METsig",&METsig, "METsig/F");
	FOM_tree->Branch("delRemu",&delRemu, "delRemu/F");
	FOM_tree->Branch("PtbalZMET",&PtbalZMET, "PtbalZMET/F");
	FOM_tree->Branch("DphiZMET",&DphiZMET, "DphiZMET/F");
	FOM_tree->Branch("Zmass",&Zmass, "Zmass/F");
	FOM_tree->Branch("ZmassSVD",&ZmassSVD, "ZmassSVD/F");
	FOM_tree->Branch("ZmassSVDnegSol",&ZmassSVDnegSol, "ZmassSVDnegSol/F");
	FOM_tree->Branch("ZmassNegInclu",&ZmassNegInclu, "ZmassNegInclu/F");
	FOM_tree->Branch("DphiSecondMET",&DphiSecondMET, "DphiSecondMET/F");
	FOM_tree->Branch("DphiLeadMET",&DphiLeadMET, "DphiLeadMET/F");
	FOM_tree->Branch("topMass",&topMass, "topMass/F");
	FOM_tree->Branch("ProjVisT",&ProjVisT, "ProjVisT/F");
	FOM_tree->Branch("ProjMissT",&ProjMissT, "ProjMissT/F");

	
	
	
    // Here one can declare histograms
    // In compiled C++, it possible not to use pointer type (TH1F*) and
    // use the non-pointer type (TH1F) of ROOT object.
	
	//Declare Histograms
	
    TH1F hallhJet_pt  ("hallhJet_pt","Pt of all jets in event",		150, 0.0, 150);
    TH1F hEleIDReco  ("hEleIDReco","Electron ID/Reco weight after EleFakeCuts ",		101, -0.5, 2.0);
	TH1F hCutFlow	("hCutFlow",  "Selection",					11, 0.5, 11.5);
	
	//NotThisCut
	TH1F hMemu_NotThisCut		("hMemu_NotThisCut",  "Invariant Mass of two Leptons", 20, 0, 100);
	TH1F hMjj_NotThisCut		("hMjj_NotThisCut",  "Invariant Mass of two Jets", 25, 0, 250);
	TH1F hDphiZMET_NotThisCut		("hDphiZMET_NotThisCut", "Delta phi between Z and MET",  24, 0, 4.71238898);
	TH1F hCSV0_NotThisCut	("hCSV0_NotThisCut", "CSV BTag Shape", 30, 0, 1.5);
	TH1F hPzeta_NotThisCut		("hPzeta_NotThisCut", "P_#zeta Cut", 25, -100, 150);
	TH1F hdelRemu_NotThisCut		("hdelRemu_NotThisCut", "Delta R emu", 20, 0, 5);
	TH1F hDeltaPhiHV_NotThisCut ("hDeltaPhiHV_NotThisCut", "Delta phi between Z and Higgs", 24, 0, 4.71238898);
	
 	
	if (debug) std::cout << "all histograms declared " << std::endl;
	
    // Loop over all events.
    // For other methods to access event/navigate through the sample,
    // see the documentation of RootTreeReader class.
	float  N_Vtype =0.0, Ntrigger = 0.0, Npreselect =0.0, N_Mjj =0.0, N_DphiZMET =0.0, NMemu =0.0, N_Naj =0.0, N_CSV0 =0.0, N_EfakeCuts =0.0, N_jetCHF0 = 0.0;
	float N_Mt = 0.0, N_Pzeta =0.0, N_DeltaPhiHV = 0.0;
	int event =0, Ntree =0, N_isomu = 0, N_tightDoubleele = 0, N_loosedoubleEle = 0, N_sigleEleWP80 = 0, N_singleEleDiJet = 0;
	float N_TopCR = 0.0, N_SingleTopCR = 0.0, N_LFCR = 0.0, N_HFCR =0.0;
	float FailedJetID=0.0, NTMVAtree =0.0, NBDTtree =0.0, NBDTbtree = 0.0;
	int  NNegEleMissE = 0, NPosEleMissE = 0, NPosMuMissE = 0, NNegMuMissE = 0;
	int NBothMissENeg = 0, NMixedEleMissENeg = 0, NMixedMuonMissENeg =0;
	float SFUnc_TOPCR = 0.0, SFUnc_SingleTOP = 0.0, SFUnc_LFCR = 0.0, SFUnc_HFCR = 0.0;
	bool firstevent = true;
	
	//	double Emumass = 91.1976;
	int Ntrainingsample = 0, Nmatch = 0, NSameCandidate=0, NgoodGenInfo = 0, Njet1Bbar =0, Njet2Bbar =0, NORMatch = 0, Nmatch4 = 0, Nmatch15 = 0;
	
    do {
		//event = event + 1*PUweight2011;
		event++;
		
		//	if (isM50sample && sample.genZpt < 100) cout << "event cut out to remove overlap with DY_PtZ" << endl;
		if (((isM50sample && sample.genZpt < 100) || !isM50sample) && ((HFfile && sample.eventFlav == 5)||!HFfile) && ((LFfile && sample.eventFlav!=5)||!LFfile)){
			//cout << "Vtype is " << sample.Vtype << endl;
			B2011PUweight = 1.0, A2011PUweight = 1.0;
			B2011PUweight = sample.PUweight2011B;
			A2011PUweight = sample.PUweight;
			float PUweight2011 = (2.219*A2011PUweight + (lumi-2.219)*B2011PUweight)/lumi;
			//cout << "pile up weight: " << PUweight2011 << endl;
			isdata = false;
			Trigweight = 1.0;
			//if (!isDATA) Trigweight = sample.weightTrig;
			weight = 1.0;
			
			
			if(sample.H_mass>0.0){
			isdata = false;
				for (int i=0; i < 5; ++i) {
					CSVNewShape[i] = -99.99;
					jetPt[i] = -99.99;
					jetEta[i] = -99.99;
					jetPhi[i] = -99.99;
					jetCSV[i] = -99.99;
					jetCHF[i] = -99.99;
					RegjetPt[i] = -99.99;
					jetPtRaw[i] = -99.99;
					jetE[i] = -99.99;
					jetVtx3dL[i] = -99.99;
					jetVtx3deL[i] = -99.99;
					jetVtxPt[i] = -99.99;
					jetVtxMass[i] = -99.99;
					jetPtLeadTrack[i] = -99.99;
					jetNconstintuents[i] = -99.99;
					jetCEF[i] = -99.99;
					jetNCH[i] = -99.99;
					jetJECUnc[i] = -99.99;
					jetEt[i] = -99.99;
					jetMt[i] = -99.99;
					jetPtRawJER[i] = -99.99;
					jetGenPt[i] = -99.99;
					SVnumTracks[i] = -99.99;
					jetNHF[i] = -99.99;
					SVmass[i] = -99.99;
				}
				double CSVshapeNew = -99.99;
				for (int k=0;k<sample.nhJets;k++){
					if (debug) cout << "for "<< k << "th pass" << endl;
					if (debug)	cout << "CSV discriminator: " << sample.hJet_csv[k] << endl;
					if(sample.hJet_csv[k]<=0 || sample.hJet_csv[k]>=1) CSVshapeNew=sample.hJet_csv[k];
					else if(sample.hJet_flavour[k]==0) CSVshapeNew=sample.hJet_csv[k];
					else if(fabs(sample.hJet_flavour[k])==5) CSVshapeNew=btagNew->ib->Eval(sample.hJet_csv[k]);
					else if(fabs(sample.hJet_flavour[k])==4) CSVshapeNew=btagNew->ic->Eval(sample.hJet_csv[k]);
					else if(fabs(sample.hJet_flavour[k])!=5 && fabs(sample.hJet_flavour[k])!=4)  CSVshapeNew=btagNew->il->Eval(sample.hJet_csv[k]);
					jetPt[k] = sample.hJet_pt[k];
					jetEta[k] = sample.hJet_eta[k];
					jetPhi[k] = sample.hJet_phi[k];
					jetCSV[k] = sample.hJet_csv[k];
					jetCHF[k] = sample.hJet_chf[k];
					jetPtRaw[k] = sample.hJet_ptRaw[k];
					jetE[k] = sample.hJet_e[k];
					jetVtx3dL[k] = sample.hJet_vtx3dL[k];
					jetVtx3deL[k] = sample.hJet_vtx3deL[k];
					jetVtxPt[k] = sample.hJet_vtxPt[k];
					jetVtxMass[k] = sample.hJet_vtxMass[k];
					jetPtLeadTrack[k] = sample.hJet_ptLeadTrack[k];
					jetNconstintuents[k] = sample.hJet_nconstituents[k];
					jetCEF[k] = sample.hJet_cef[k];
					jetNCH[k] = sample.hJet_nch[k];
					jetJECUnc[k] = sample.hJet_JECUnc[k];
					jetEt[k] = evalEt(jetPt[k],jetEta[k],jetPhi[k],jetE[k]);
					jetMt[k] = evalMt(jetPt[k],jetEta[k],jetPhi[k],jetE[k]);
					jetPtRawJER[k] = evalJERBias(jetPtRaw[k], sample.hJet_genPt[k], jetEta[k]);
					jetGenPt[k] = sample.hJet_genPt[k];
					SVnumTracks[k] = sample.hJet_numTracksSV[k];
					jetNHF[k] = sample.hJet_nhf[k];
					CSVNewShape[k] = CSVshapeNew;
					CSVshapeNew = -99.99;
				}// end k for loop
				for (int z=0;z<sample.nSvs;z++) SVmass[z]=sample.Sv_massSv[z];
				RegjetPt[0] = -99.99;
				RegjetPt[1] = -99.99;	
				RegHmass = -99.99;
				oldHmass = sample.H_mass;
				Hmass = sample.H_mass;
				ScalarSumJetPt = jetPt[0]+jetPt[1];
				float MinDphiaJet = 99.99;
				for (int a=0;a<sample.naJets;a++){
					ScalarSumJetPt=ScalarSumJetPt+sample.aJet_pt[a];
					if(sample.aJet_pt[a]>30&&abs(sample.aJet_eta[a])<2.5){
						if(MinDphiaJet>deltaPhi(sample.MET_phi,sample.aJet_phi[a]))
							MinDphiaJet =deltaPhi(sample.MET_phi,sample.aJet_phi[a]);
					}}
				ScalarSumHiggsJetPt = jetPt[0]+jetPt[1];
				DetaJJ = sample.hJet_eta[0]-sample.hJet_eta[1];
				Hpt = sample.H_pt;
				Rho25 = sample.rho25;
				DeltaPhijetMETmin = sample.minDeltaPhijetMET;
				MinMET = min(sample.METtype1corr_et,sample.METnoPUCh_et);
				MindPhiMEThJetOR30aJet = min(sample.minDeltaPhijetMET,MinDphiaJet);
				
				bool bmatch = false;
				bool goodGenInfo = true;
				float DelR0b = 99.99, DelR0Bbar = 99.99, DelR1b = 99.99, DelR1Bbar = 99.99;
				if ((sample.genB_pt-sample.genB_eta)<0.-1 || (sample.genBbar_pt-sample.genBbar_eta)<0.01) goodGenInfo = false;
				if (goodGenInfo){
					DelR0b = sqrt((sample.hJet_genEta[0]-sample.genB_eta)*(sample.hJet_genEta[0]-sample.genB_eta)+(sample.hJet_genPhi[0]-sample.genB_phi)*(sample.hJet_genPhi[0]-sample.genB_phi));
					DelR0Bbar = sqrt((sample.hJet_genEta[0]-sample.genBbar_eta)*(sample.hJet_genEta[0]-sample.genBbar_eta)+(sample.hJet_genPhi[0]-sample.genBbar_phi)*(sample.hJet_genPhi[0]-sample.genBbar_phi));
					DelR1b = sqrt((sample.hJet_genEta[1]-sample.genB_eta)*(sample.hJet_genEta[1]-sample.genB_eta)+(sample.hJet_genPhi[1]-sample.genB_phi)*(sample.hJet_genPhi[1]-sample.genB_phi));
					DelR1Bbar = sqrt((sample.hJet_genEta[1]-sample.genBbar_eta)*(sample.hJet_genEta[1]-sample.genBbar_eta)+(sample.hJet_genPhi[1]-sample.genBbar_phi)*(sample.hJet_genPhi[1]-sample.genBbar_phi));
				
				
				//both match
				if ( (DelR0b<0.3||DelR0Bbar<0.3)&&(DelR1b<0.3||DelR1Bbar<0.3) ){				//reasonable match for each jet
				Nmatch++;
					if (DelR0b<DelR0Bbar&&DelR1Bbar<DelR1b) {
					bmatch = true; //hJet[0] matched to b, hJet[1] matched to b-bar
					Njet2Bbar++;}
					if (DelR0b>DelR0Bbar&&DelR1Bbar>DelR1b){
					 bmatch = true; //hJet[0] matched to b-bar, hJet[[0] matched to b
					Njet1Bbar++;}
				if ((DelR0b>0.3&&DelR1b>0.3)||(DelR0Bbar>0.3&&DelR1Bbar>0.3)) {
					bmatch = false; //both jets matched to same b candidate
					NSameCandidate++;
					}
				}
				if ((DelR0b<0.3||DelR0Bbar<0.3)||(DelR1b<0.3||DelR1Bbar<0.3)) NORMatch++;
				if ( (DelR0b<0.4||DelR0Bbar<0.4)&&(DelR1b<0.4||DelR1Bbar<0.4) ) Nmatch4++;
				if ( (DelR0b<0.15||DelR0Bbar<0.15)&&(DelR1b<0.15||DelR1Bbar<0.15) ) Nmatch15++;
					NgoodGenInfo++;
				}
				
			//	JetDistributions("Zbb", weight);
			//	if (sample.hJet_csv[0]>0.244&&sample.hJet_csv[1]>0.244)		JetDistributions("ZbbCSVL", weight);
			//	if (sample.hJet_csv[0]>0.898&&sample.hJet_csv[1]>0.898)		JetDistributions("ZbbCSVT", weight);
			//	if (sample.V_mass > 75 && sample.V_mass < 105 && (sample.Vtype==1||sample.Vtype==0)) JetDistributions("ZbbZmass", weight);
				if (sample.hJet_pt[0]>15&&sample.hJet_pt[1]>15 && (sample.Vtype!=5) && bmatch && sample.H_mass<300) {
				JetDistributions("TrainingSample", weight);
				TH2FDistributions("TrainingSample",weight);
				Ntrainingsample++;
				}
				if(sample.hJet_pt[0]>15&&sample.hJet_pt[1]>15 &&sample.H_mass<300){
					if(jetGenPt[0]>50)JetDistributions("HighEnergy0jet",weight);
					if(jetGenPt[0]<50)JetDistributions("LowEnergy0jet",weight);
					if(jetGenPt[1]>50)JetDistributions("HighEnergy1jet",weight);
					if(jetGenPt[1]<50)JetDistributions("LowEnergy1jet",weight);
				}
			}//if H_mass
			
			weight = LumiWeight*PUweight2011;
			isdata = false;
			Trigweight = 1.0;
			//if (!isDATA) Trigweight = sample.weightTrig;
			weight = Trigweight*weight;
			
			if (sample.Vtype == 5 && sample.EVENT_json == 1 ){
				if (!(event%500))  std::cout << "entered event loop " << event << std::endl;
				isdata= false;
				
				std::vector< std::pair<size_t,double> > indexedJetPt;
				std::vector<size_t> PtSortedJetIndex;
				std::vector< std::pair<size_t,double> > indexedJetCSV;
				std::vector<size_t> CSVSortedJetIndex;
				std::vector< std::pair<size_t,double> > indexedPt;
				std::vector<size_t> PtSortedIndex;
				std::vector< std::pair<size_t,double> > indexedMuPt;
				std::vector<size_t> MuonPtSortedIndex;
				
				// Analysis loop.
				// One can access the ntuple leaves directly from sample object
				
				for (int i=0; i < 5; ++i) {
					CSVNewShape[i] = -99.99;
					jetPt[i] = -99.99;
					jetEta[i] = -99.99;
					jetPhi[i] = -99.99;
					jetCSV[i] = -99.99;
					jetCHF[i] = -99.99;
					RegjetPt[i] = -99.99;
					jetPtRaw[i] = -99.99;
					jetE[i] = -99.99;
					jetVtx3dL[i] = -99.99;
					jetVtx3deL[i] = -99.99;
					jetVtxPt[i] = -99.99;
					jetVtxMass[i] = -99.99;
					jetPtLeadTrack[i] = -99.99;
					jetNconstintuents[i] = -99.99;
					jetCEF[i] = -99.99;
					jetNCH[i] = -99.99;
					jetJECUnc[i] = -99.99;
					jetEt[i] = -99.99;
					jetMt[i] = -99.99;
					jetPtRawJER[i] = -99.99;
					jetGenPt[i] = -99.99;
					leptonPt[i] = -99.99;
					leptonEta[i] = -99.99;
					leptonPhi[i] = -99.99;
					lep_pfCombRelIso[i] = -99.99;
					lep_id95[i] = -99.99;
					lep_flavor[i] = -99;
					SortedMuonPt[i] = -99.99;
				}
				
				nJets =0, nSV =-99, nMuons = 0,  nElectrons = 0,  nLeptons = 0, Na_lep= 0, nPV= -99, MET= -99.99, Naj= 0, Nab =0, eventFlavor = -99;
				CSV0 = -1.0, CSV1 = -1.0, Emumass = -99.99, Hmass = -99.99, oldHmass = -99.99, RegHmass = -99.99, DeltaPhiHV = -99.99, Hpt = -99.99, Zpt = -99.99;
				lep0pt = -99.99, ScalarSumPt = -99.99, EtaStandDev = -99.99, UnweightedEta = -99.99, EvntShpCircularity = -99.99;
				alpha_j = -99.99, qtb1 = 0.0, DphiJJ = -99.99,  btag2CSF = 1.0;
				RMS_eta = -99.99, PtbalZH= -99.99, EventPt= -99.99, AngleHemu= -99.99, Centrality = -99.99;
				qtlep1 = 0.0, alpha_lep = -99.99;
				EvntShpAplanarity = -99.99, EvntShpSphericity = -99.99, EvntShpIsotropy = -99.99;
				DetaJJ = -99.99;
				Zphi = -99.99, Hphi =-99.99;
				dPhiHMET = -99.99, Mt = -99.99, DeltaPhijetMETmin = -99.99, DeltaPhijetMETZtaumin = -99.99;
				SV_mass = -99.99, Mte = -99.99, Mtmu = -99.99, Ht = -99.99, delPullAngle = -99.99, delPullAngle2 = -99.99;
				delRjj = -99.99, Detaemu = -99.99, DphiEleMET = -99.99, dphiMuMET = -99.99, PtbalMETH = -99.99, naJets = 0;
				topMass = -99.99, topPt = -99.99, topWmass = -99.99;
				MassEleb0 = -99.99, MassMub0 = -99.99, MassEleb1 = -99.99;
				MassMub1 = -99.99, METsig = -99.99;
				delRemu = -99.99, PtbalZMET = -99.99, DphiZMET = -99.99;
				DphiLeadMET = -99.99, DphiSecondMET = -99.99;
				ProjVisT= -99.99, ProjMissT= -99.99;
				ScalarSumJetPt = 0.0, ScalarSumHiggsJetPt = -99.99, EventMass = -99.99;
				
				
				double StandDevEta[4];
				double AverageEta = -99.99;
				TLorentzVector FirstJet;
				TLorentzVector SecondJet;
				TLorentzVector Muon;
				TLorentzVector Electron;
				//TLorentzVector METvector;
				TLorentzVector Higgs;
				TLorentzVector ZBoson;
				TLorentzVector HVsystem;
				TLorentzVector DummyVector;
				double pllep1 = 0.0, pllep2 = 0.0 , qtlep2 = 0.0;
				double plb1 = 0.0, plb2 = 0.0 , qtb2 = 0.0;
				
				AngleEMU = -99.99, CosThetaEle = -99.99, CosThetaMu = -99.99;
				EleMissE = -99.99, MuonMissE =-99.99, Dphiemu = -99.99;
				Zmass = -99.99, ZmassSVD = -99.99, ZmassSVDnegSol = -99.99, ZmassNegInclu = -99.99;
				AaronEleMissE = -599.99, AaronMuMissE = -599.99;
				
				double CSVshapeNew = -99.99;
				for (int k=0;k<sample.nhJets;k++){
					if (debug) cout << "for "<< k << "th pass" << endl;
					indexedJetPt.push_back(std::pair<size_t,double>(k,(double) sample.hJet_pt[k]));
					indexedPt.push_back(std::pair<size_t,double>(k,(double) sample.hJet_pt[k]));
					hallhJet_pt.Fill(sample.hJet_pt[k]); 
					nJets++;
					//indexedJetCSV.push_back(std::pair<size_t,double>(k,(double) sample.hJet_csv[k]));
					if (debug)	cout << "CSV discriminator: " << sample.hJet_csv[k] << endl;
					if(sample.hJet_csv[k]<=0 || sample.hJet_csv[k]>=1) CSVshapeNew=sample.hJet_csv[k];
					else if(sample.hJet_flavour[k]==0) CSVshapeNew=sample.hJet_csv[k];
					else if(fabs(sample.hJet_flavour[k])==5) CSVshapeNew=btagNew->ib->Eval(sample.hJet_csv[k]);
					else if(fabs(sample.hJet_flavour[k])==4) CSVshapeNew=btagNew->ic->Eval(sample.hJet_csv[k]);
					else if(fabs(sample.hJet_flavour[k])!=5 && fabs(sample.hJet_flavour[k])!=4)  CSVshapeNew=btagNew->il->Eval(sample.hJet_csv[k]);
					if (!(event%500)) cout << "The orignial CSV value was " << sample.hJet_csv[k] << "  the corrected value is " << CSVshapeNew << endl;
					indexedJetCSV.push_back(std::pair<size_t,double>(k,(double) CSVshapeNew));
					CSVshapeNew = -99.99;
				}// end k for loop
				
				float MinDphiaJet = 99.99;
				for (int a=0;a<sample.naJets;a++){
					hallhJet_pt.Fill(sample.aJet_pt[a]);
					nJets++; 
					if(sample.aJet_csv[a]<=0 || sample.aJet_csv[a]>=1) CSVshapeNew=sample.aJet_csv[a];
					else if(sample.aJet_flavour[a]==0) CSVshapeNew=sample.aJet_csv[a];
					else if(fabs(sample.aJet_flavour[a])==5) CSVshapeNew=btagNew->ib->Eval(sample.aJet_csv[a]);
					else if(fabs(sample.aJet_flavour[a])==4) CSVshapeNew=btagNew->ic->Eval(sample.aJet_csv[a]);
					else if(fabs(sample.aJet_flavour[a])!=5 && fabs(sample.aJet_flavour[a])!=4)  CSVshapeNew=btagNew->il->Eval(sample.aJet_csv[a]);				
					if (a<3){
						jetPt[sample.nhJets+a] = sample.aJet_pt[a];
						jetEta[sample.nhJets+a] = sample.aJet_eta[a];
						jetPhi[sample.nhJets+a] = sample.aJet_phi[a];
						jetCSV[sample.nhJets+a] = sample.aJet_csv[a];
						jetCHF[sample.nhJets+a] = sample.aJet_chf[a];
						CSVNewShape[sample.nhJets+a] = CSVshapeNew;
						//jetPtRaw[sample.nhJets+a] = sample.aJet_ptRaw[a];
						jetE[sample.nhJets+a] = sample.aJet_e[a];
						jetVtx3dL[sample.nhJets+a] = sample.aJet_vtx3dL[a];
						jetVtx3deL[sample.nhJets+a] = sample.aJet_vtx3deL[a];
						//jetVtxPt[sample.nhJets+a] = sample.aJet_vtxPt[a];
						jetVtxMass[sample.nhJets+a] = sample.aJet_vtxMass[a];
						//jetPtLeadTrack[sample.nhJets+a] = sample.aJet_ptLeadTrack[a];
						jetNconstintuents[sample.nhJets+a] = sample.aJet_nconstituents[a];
						jetCEF[sample.nhJets+a] = sample.aJet_cef[a];
						jetNCH[sample.nhJets+a] = sample.aJet_nch[a];
						jetJECUnc[sample.nhJets+a] = sample.aJet_JECUnc[a];
						jetEt[sample.nhJets+a] = evalEt(sample.aJet_pt[a],sample.aJet_eta[a],sample.aJet_phi[a],sample.aJet_e[a]);
						jetMt[sample.nhJets+a] = evalMt(sample.aJet_pt[a],sample.aJet_eta[a],sample.aJet_phi[a],sample.aJet_e[a]);
						//jetPtRawJER[sample.nhJets+a] = evalJERBias(sample.aJet_ptRaw[a], sample.aJet_genPt[a], sample.aJet_eta[a]);
						jetGenPt[sample.nhJets+a]= sample.aJet_genPt[a];
					}
					if(sample.aJet_pt[a]>30&&abs(sample.aJet_eta[a])<2.5) {
					if(MinDphiaJet>deltaPhi(sample.MET_phi,sample.aJet_phi[a])) MinDphiaJet =deltaPhi(sample.MET_phi,sample.aJet_phi[a]);
					}
					if ( fabs(sample.aJet_eta[a]) < 2.4 && (sample.aJet_pt[a] > 20)) Naj++;
					if (CSVshapeNew > 0.4 ) Nab++;
					indexedJetPt.push_back(std::pair<size_t,double>(sample.nhJets+a,(double) sample.aJet_pt[a]));	
					CSVshapeNew = -99.99;
					ScalarSumJetPt = ScalarSumJetPt + sample.aJet_pt[a];
				}
				indexedPt.push_back(std::pair<size_t,double>(2,(double) sample.vLepton_pt[0]));
				indexedPt.push_back(std::pair<size_t,double>(3,(double) sample.vLepton_pt[1]));
				
				
				if (firstevent) for (size_t i = 0 ; (i != indexedPt.size()) ; ++i) {cout << "Unsorted pt of objects: " << indexedPt[i].second << endl; }
				
				std::sort(indexedJetPt.begin(),indexedJetPt.end(),::IndexedQuantityGreaterThan<double>);
				std::sort(indexedJetCSV.begin(),indexedJetCSV.end(),::IndexedQuantityGreaterThan<double>);
				std::sort(indexedPt.begin(),indexedPt.end(),::IndexedQuantityGreaterThan<double>);
				
				for (size_t i = 0 ; (i != indexedJetPt.size()) ; ++i) {   PtSortedJetIndex.push_back(indexedJetPt[i].first);        }
				for (size_t i = 0 ; (i != indexedJetCSV.size()) ; ++i) {   CSVSortedJetIndex.push_back(indexedJetCSV[i].first);        }
				for (size_t i = 0 ; (i != indexedPt.size()) ; ++i) {   PtSortedIndex.push_back(indexedPt[i].first);        }
				
				if (firstevent) for (size_t i = 0 ; (i != indexedPt.size()) ; ++i) {cout << "Sorted pt of objects: " << indexedPt[i].second << endl; }
				firstevent = false;
				CSV0 = indexedJetCSV[0].second; 
				jetPt[0] = sample.hJet_pt[0];
				jetEta[0] = sample.hJet_eta[0];
				jetPhi[0] = sample.hJet_phi[0];
				jetCSV[0] = sample.hJet_csv[0];
				jetCHF[0] = sample.hJet_chf[0];
				jetPtRaw[0] = sample.hJet_ptRaw[0];
				jetE[0] = sample.hJet_e[0];
				jetVtx3dL[0] = sample.hJet_vtx3dL[0];
				jetVtx3deL[0] = sample.hJet_vtx3deL[0];
				jetVtxPt[0] = sample.hJet_vtxPt[0];
				jetVtxMass[0] = sample.hJet_vtxMass[0];
				jetPtLeadTrack[0] = sample.hJet_ptLeadTrack[0];
				jetNconstintuents[0] = sample.hJet_nconstituents[0];
				jetCEF[0] = sample.hJet_cef[0];
				jetNCH[0] = sample.hJet_nch[0];
				jetJECUnc[0] = sample.hJet_JECUnc[0];
				jetEt[0] = evalEt(jetPt[0],jetEta[0],jetPhi[0],jetE[0]);
				jetMt[0] = evalMt(jetPt[0],jetEta[0],jetPhi[0],jetE[0]);
				jetPtRawJER[0] = evalJERBias(jetPtRaw[0], sample.hJet_genPt[0], jetEta[0]);
				jetGenPt[0]=sample.hJet_genPt[0];
				
				//if (indexedJetCSV[0].first == 0) RegjetPt[0] = sample.hJet_genPtReg0;
				//if (indexedJetCSV[0].first == 1) RegjetPt[0] = sample.hJet_genPtReg1;
				if(jetPt[0] == sample.hJet_pt[0]){
				RegjetPt[0] = -99.99;
				RegjetPt[1] = -99.99;
			}
				CSVNewShape[0] = indexedJetCSV[0].second;
				if (sample.nhJets > 1) { 
					CSV1 = indexedJetCSV[1].second;	
					jetPt[1] = sample.hJet_pt[1];
					jetEta[1] = sample.hJet_eta[1];
					jetPhi[1] = sample.hJet_phi[1];
					jetCSV[1] = sample.hJet_csv[1];
					jetCHF[1] = sample.hJet_chf[1];
					jetPtRaw[1] = sample.hJet_ptRaw[1];
					jetE[1] = sample.hJet_e[1];
					jetVtx3dL[1] = sample.hJet_vtx3dL[1];
					jetVtx3deL[1] = sample.hJet_vtx3deL[1];
					jetVtxPt[1] = sample.hJet_vtxPt[1];
					jetVtxMass[1] = sample.hJet_vtxMass[1];
					jetPtLeadTrack[1] = sample.hJet_ptLeadTrack[1];
					jetNconstintuents[1] = sample.hJet_nconstituents[1];
					jetCEF[1] = sample.hJet_cef[1];
					jetNCH[1] = sample.hJet_nch[1];
					jetJECUnc[1] = sample.hJet_JECUnc[1];
					jetEt[1] = evalEt(jetPt[1],jetEta[1],jetPhi[1],jetE[1]);
					jetMt[1] = evalMt(jetPt[1],jetEta[1],jetPhi[1],jetE[1]);
					jetPtRawJER[1] = evalJERBias(jetPtRaw[1], sample.hJet_genPt[1], jetEta[1]);
					jetGenPt[1]=sample.hJet_genPt[1];

					//if (indexedJetCSV[1].first == 0) RegjetPt[1] = sample.hJet_genPtReg0;
					//if (indexedJetCSV[1].first == 1) RegjetPt[1] = sample.hJet_genPtReg1;					
					
					DetaJJ = sample.hJet_eta[0]-sample.hJet_eta[1];				
					Hmass = sample.H_mass;
					oldHmass = sample.H_mass;
					RegHmass = -99.99;
					GenHiggsMass = sample.genH_mass;
					Hpt = sample.H_pt;
					if(sample.svmass<999)DitauMass = sample.svmass;
					ScalarSumHiggsJetPt = sample.hJet_pt[1] + sample.hJet_pt[0];
					ScalarSumJetPt = ScalarSumJetPt+ sample.hJet_pt[1] + sample.hJet_pt[0];
					Rho25 = sample.rho25;
					MindPhiMEThJetOR30aJet = min(sample.minDeltaPhijetMET,MinDphiaJet);
					MinMET = min(sample.METtype1corr_et,sample.METnoPUCh_et);
					//cout << "weight of the Jet Histograms: " << weight << endl;
					JetDistributions("allEvts", weight);
				}// end njet requirement			
				
				nSV = sample.nSvs;
				nPV = sample.nPVs;
				MET = sample.METtype1corr_et;
				METsig = sample.METtype1corr_sig;
				eventFlavor = sample.eventFlav;
				SV_mass = sample.Sv_massSv[0];
				naJets = sample.naJets;
				
				btag2CSF = sample.btag2CSF;
				if (debug) cout << "halfway through filling histos" << endl;
				
				for (int z = 0; z< sample.nvlep;z++){
					leptonPt[z] = sample.vLepton_pt[z];
					leptonEta[z] = sample.vLepton_eta[z];
					leptonPhi[z] = sample.vLepton_phi[z];
					lep_pfCombRelIso[z] = sample.vLepton_pfCombRelIso[z];
					lep_id95[z] = sample.vLepton_id95[z];
					if (abs(sample.aLepton_type[z] == 13)) {
						nMuons++;
						indexedMuPt.push_back(std::pair<size_t,double>(z,(double) sample.vLepton_pt[z]));
					}
					if (abs(sample.aLepton_type[z] == 11)) nElectrons++;
					nLeptons++;
				}
				for (int zz = 0; zz< sample.nalep;zz++){
					if (abs(sample.aLepton_type[zz] == 13)) {
						nMuons++;
						indexedMuPt.push_back(std::pair<size_t,double>(zz+1,(double) sample.vLepton_pt[zz+1]));
					}
					if (abs(sample.aLepton_type[zz] == 11)) nElectrons++;
					if (zz<3){
						if ((sample.aLepton_pt[zz] == sample.vLepton_pt[0]) || (sample.aLepton_pt[zz] == sample.vLepton_pt[1])){
							cout << "Aditional Letpon " << zz+sample.nvlep << " pt is the same as vector lepton" << endl;
						} else {
							leptonPt[zz+sample.nvlep] = sample.aLepton_pt[zz];
							nLeptons++;
						}
						leptonPhi[zz+sample.nvlep] = sample.aLepton_phi[zz];
						leptonEta[zz+sample.nvlep] = sample.aLepton_eta[zz];
						lep_pfCombRelIso[zz+sample.nvlep] = sample.aLepton_pfCombRelIso[zz];
						Na_lep++;
					}
				}
				if (sample.nvlep > 1) {
					Emumass = sample.V_mass;
					Zpt = sample.V_pt;
					DphiZMET = deltaPhi(sample.METtype1corr_phi,sample.V_phi);
				}
				
				if(firstevent) for (size_t i = 0 ; (i != indexedMuPt.size()) ; ++i) {cout << "Unsorted pt of muons: " << indexedMuPt[i].second << endl; }
				
				std::sort(indexedMuPt.begin(),indexedMuPt.end(),::IndexedQuantityGreaterThan<double>);
				
				if (firstevent) for (size_t i = 0 ; (i != indexedMuPt.size()) ; ++i) {   MuonPtSortedIndex.push_back(indexedMuPt[i].first);        }
				
				if(firstevent) for (size_t i = 0 ; (i != indexedMuPt.size()) ; ++i) {cout << "Sorted pt of muons: " << indexedMuPt[i].second << endl; }
				for (size_t k = 0; k < indexedMuPt.size(); k++){
					SortedMuonPt[k] = sample.hJet_pt[indexedMuPt[k].first];
				}
				
				LeptonDistributions("allEvts", weight);
				dPhiHMET = sample.HMETdPhi;
				DeltaPhiHV = sample.HVdPhi;			
				Mt = sample.VMt;
				DeltaPhijetMETmin = sample.minDeltaPhijetMET; 
				DeltaPhijetMETZtaumin = sample.minDeltaPhijetMETZtau;
				//Mte = sample.VMte;
				//Mtmu = sample.VMtmu;
				Ht = sample.MHT_mht;
				delPullAngle = sample.deltaPullAngle;
				delPullAngle2 = sample.deltaPullAngle2;
				topMass = sample.top_mass;
				//cout << "top mass " << topMass << endl;
				topPt = sample.top_pt;
				topWmass = sample.top_wMass;
				
				if( (sample.nvlep > 1) && (sample.nhJets > 1)) {
					RMS_eta = sample.vLepton_eta[1]*sample.vLepton_eta[1]+sample.vLepton_eta[0]*sample.vLepton_eta[0]+sample.hJet_eta[0]*sample.hJet_eta[0]+sample.hJet_eta[1]*sample.hJet_eta[1];
					RMS_eta = sqrt(RMS_eta/4);
					StandDevEta[0] =sample.vLepton_eta[1];
					StandDevEta[1] =sample.hJet_eta[0];
					StandDevEta[2] =sample.vLepton_eta[0];
					StandDevEta[3] =sample.hJet_eta[1];
					EtaStandDev = TMath::RMS(4,StandDevEta);
					AverageEta = sample.vLepton_eta[0]+sample.vLepton_eta[1]+sample.hJet_eta[0];
					AverageEta = (AverageEta+sample.hJet_eta[1])/4;
					UnweightedEta = (sample.vLepton_eta[0]-AverageEta)*(sample.vLepton_eta[0]-AverageEta);
					UnweightedEta = UnweightedEta + (sample.vLepton_eta[1]-AverageEta)*(sample.vLepton_eta[1]-AverageEta);
					UnweightedEta = UnweightedEta + (sample.hJet_eta[0]-AverageEta)*(sample.hJet_eta[0]-AverageEta);
					UnweightedEta = UnweightedEta + (sample.hJet_eta[1]-AverageEta)*(sample.hJet_eta[1]-AverageEta);
					
					PtbalZH = (Hpt-Zpt);
					PtbalMETH = Hpt-MET;
					PtbalZMET = Zpt - MET;
					lep0pt = sample.vLepton_pt[0];
					
					FirstJet.SetPtEtaPhiM(sample.hJet_pt[0],sample.hJet_eta[0],sample.hJet_phi[0],4.2/1000.0);
					SecondJet.SetPtEtaPhiM(sample.hJet_pt[1],sample.hJet_eta[1],sample.hJet_phi[1],4.2/1000.0);
					Muon.SetPtEtaPhiM(sample.vLepton_pt[1],sample.vLepton_eta[1],sample.vLepton_phi[1],105.65836668/1000.0);
					Electron.SetPtEtaPhiM(sample.vLepton_pt[0],sample.vLepton_eta[0],sample.vLepton_phi[0],0.510998928/1000.0);
					
					DummyVector = Electron+FirstJet;
					MassEleb0 = DummyVector.M();
					DummyVector = Muon+FirstJet;
					MassMub0 = DummyVector.M();
					DummyVector = Electron+SecondJet;
					MassEleb1 = DummyVector.M();
					DummyVector = Muon+SecondJet;
					MassMub1 = DummyVector.M();
					
					Higgs = FirstJet+SecondJet;
					ZBoson = Muon+Electron;
					
					Zphi = sample.H_phi;
					Hphi = sample.V_phi;
					
					Dphiemu = Electron.DeltaPhi(Muon);
					Detaemu = sample.vLepton_eta[1]-sample.vLepton_eta[0];
					delRemu = Electron.DeltaR(Muon);
					DphiEleMET = deltaPhi(sample.vLepton_phi[0],sample.METtype1corr_phi);
					dphiMuMET = deltaPhi(sample.vLepton_phi[1],sample.METtype1corr_phi);
					if (sample.vLepton_pt[0] > sample.vLepton_pt[1]) {
						DphiLeadMET = DphiEleMET;
						DphiSecondMET = dphiMuMET;
					}
					if (sample.vLepton_pt[1] > sample.vLepton_pt[0]) {
						DphiLeadMET = dphiMuMET;
						DphiSecondMET = DphiEleMET;
					}
					
					//Bisector
					CosThetaEle = Electron.CosTheta();
					CosThetaMu = Muon.CosTheta();
					float phielectron = sample.vLepton_phi[0];
					float phimuon = sample.vLepton_phi[1];
					TVector3 m,e,metv;
					m.SetPtEtaPhi(sample.vLepton_pt[1],0,sample.vLepton_phi[1]);
					e.SetPtEtaPhi(sample.vLepton_pt[0],0,sample.vLepton_phi[0]);
					metv.SetPtEtaPhi(MET,0,sample.METtype1corr_phi);
					TVector3 bisector(m.Unit() + e.Unit());
					bisector = bisector.Unit();
					ProjVisT  = (m+e).Dot(bisector);
					ProjMissT  =  metv.Dot(bisector);
					//Double_t pzeta = projMet - 0.85*projVis;
					
					//Transverse Mass because Ntupler bug
					//m=Muon.Vect();
					//e=Electron.Vect();
					
					//m.SetZ(0);
					//e.SetZ(0);
					//metv.SetZ(0);
					float dotprod = e.Dot(metv);
					float cosphi = dotprod/(e.Mag()*metv.Mag());
					float et = 2*sample.vLepton_pt[0]*MET;
					Mte = sqrt(et*(1-cosphi));
					
					dotprod = m.Dot(metv);
					cosphi = dotprod/(m.Mag()*metv.Mag());
					et = 2*sample.vLepton_pt[1]*MET;
					Mtmu = sqrt(et*(1-cosphi));
					
					
					float METx = MET*cos(sample.MET_phi);
					float METy = MET*sin(sample.MET_phi);
					if (debug) cout << "MET is " << MET << " METx " << METx << " METy " << METy << endl;
					AngleEMU = Electron.Angle(Muon.Vect());
					
					//Try Calculating with matricies
					TMatrixDSym A(2);
					A(0,0) = CosThetaEle*cos(phielectron);
					A(1,0) = CosThetaMu*cos(phimuon);
					A(0,1) = CosThetaEle*sin(phielectron);
					A(1,1) = CosThetaMu*sin(phimuon);
					TMatrixD InverseA = A.Invert();
					
					//Single value decomposition method
					TVectorD b;
					bool ok;
					double ab[] = {METx,METy};
					b.Use(2,ab);
					TDecompSVD svd(A);
					const TVectorD c_svd = svd.Solve(b,ok);
					//cout << "Missing Energy via SVD: Electron " << c_svd(0)<< " Muon  " << c_svd(1)<< endl;
					if(c_svd(0) < 0 || c_svd(1) < 0){
						//cout << "Negative missing energy from SVD method. ";
						//cout << "Angle between electron and muon " << AngleEMU << endl;
						//cout << "Phi between electron and muon " << Dphiemu << endl;
						ZmassSVDnegSol = sqrt((2*1.77682*1.77682)+(2*(c_svd(0)+Electron.E())*(c_svd(1)+Muon.E())*(1-cos(AngleEMU))));
					} else {
						ZmassSVD = sqrt((2*1.77682*1.77682)+(2*(c_svd(0)+Electron.E())*(c_svd(1)+Muon.E())*(1-cos(AngleEMU))));
					}				
					
					//calculating from Aaron's piece of paper
					double DetA = (CosThetaEle*cos(phielectron)*CosThetaMu*sin(phimuon)) - (CosThetaMu*cos(phimuon)*CosThetaEle*sin(phielectron));
					AaronEleMissE = (1/DetA)*((METx*CosThetaMu*sin(phimuon))-(METy*CosThetaMu*cos(phimuon)));
					AaronMuMissE = (1/DetA)*(((-1)*METx*CosThetaEle*sin(phielectron))+(METy*CosThetaEle*cos(phielectron)));
					if (AaronEleMissE> 0 && AaronMuMissE>0)Zmass = sqrt((2*1.77682*1.77682)+(2*(AaronEleMissE+Electron.E())*(AaronMuMissE+Muon.E())*(1-cos(AngleEMU))));
					//cout << "Missing Energy via Aaron linear algebra: Electron " << AaronEleMissE << " Muon " << AaronMuMissE << endl;
					
					if(AaronMuMissE > 0 && AaronEleMissE > 0){
						Zmass = sqrt((2*1.77682*1.77682)+(2*(AaronEleMissE+Electron.E())*(AaronMuMissE+Muon.E())*(1-cos(AngleEMU))));
					}
					ZmassNegInclu = sqrt((2*1.77682*1.77682)+(2*(AaronEleMissE+Electron.E())*(AaronMuMissE+Muon.E())*(1-cos(AngleEMU))));
					if(AaronMuMissE < 0 || AaronEleMissE < 0){
						if (AaronMuMissE < 0 && AaronEleMissE < 0) NBothMissENeg++;
						if (AaronMuMissE > 0 && AaronEleMissE < 0) NMixedEleMissENeg++;
						if (AaronMuMissE < 0 && AaronEleMissE > 0) NMixedMuonMissENeg++;
					}				
					if (AaronEleMissE < 0) NNegEleMissE++;
					if (AaronEleMissE > 0) NPosEleMissE++;
					if (AaronMuMissE > 0) NPosMuMissE++;
					if (AaronMuMissE < 0) NNegMuMissE++;
					
					ScalarSumPt = FirstJet.Pt()+SecondJet.Pt()+Muon.Pt()+Electron.Pt();
					DphiJJ = FirstJet.DeltaPhi(SecondJet);
					delRjj = FirstJet.DeltaR(SecondJet);
					
					Centrality = (indexedPt[2].second+indexedPt[3].second)/ScalarSumPt;
					AngleHemu = Higgs.Angle(ZBoson.Vect());
					HVsystem = 	Higgs+ZBoson;
					EventPt = HVsystem.Pt();
					EventMass = HVsystem.M();
					
					//Armenteros-Podolansky Plots
					TVector3 VectZ = ZBoson.Vect();
					TVector3 Unit_Z = VectZ.Unit();
					pllep1 = Unit_Z.Dot(Muon.Vect());
					pllep2 = Unit_Z.Dot(Electron.Vect());
					alpha_lep = (pllep1-pllep2)/(pllep1+pllep2);
					qtlep1= Muon.Perp(ZBoson.Vect()); 
					qtlep2 = Electron.Perp(ZBoson.Vect());
					
					TVector3 VectH = Higgs.Vect();
					TVector3 Unit_H = VectH.Unit();
					plb1 = Unit_H.Dot(FirstJet.Vect());
					plb2 = Unit_H.Dot(SecondJet.Vect());
					if (FirstJet.Pt()>SecondJet.Pt()){
						alpha_j = (plb1-plb2)/(plb1+plb2);
					}
					if (FirstJet.Pt()<SecondJet.Pt()){
						alpha_j = (plb2-plb1)/(plb1+plb2);
					}
					qtb1= FirstJet.Perp(Higgs.Vect()); 
					qtb2 = SecondJet.Perp(Higgs.Vect());
					TH2FDistributions("allEvts", weight);
					
					//EventShapeVariables
					std::vector<math::XYZVector> particles;
					TVector3 TVect_particle = Muon.Vect();
					math::XYZVector Vparticle1(TVect_particle.X(), TVect_particle.Y(), TVect_particle.Z());
					particles.push_back(Vparticle1);
					TVect_particle = Electron.Vect();
					math::XYZVector Vparticle2(TVect_particle.X(), TVect_particle.Y(), TVect_particle.Z());
					particles.push_back(Vparticle2);
					TVect_particle = FirstJet.Vect();
					math::XYZVector Vparticle3(TVect_particle.X(), TVect_particle.Y(), TVect_particle.Z());
					particles.push_back(Vparticle3);
					TVect_particle = SecondJet.Vect();
					math::XYZVector Vparticle4(TVect_particle.X(), TVect_particle.Y(), TVect_particle.Z());
					particles.push_back(Vparticle4);
					
					EventShapeVariables eventshape(particles);
					EvntShpAplanarity = eventshape.aplanarity();
					EvntShpSphericity = eventshape.sphericity();
					EvntShpCircularity = eventshape.circularity();
					EvntShpIsotropy = eventshape.isotropy();
					
					EventShapeDistributions("allEvts", weight);
				}// two leptons and two jets requirement
				
				
				if (debug) cout << "done filling histograms for event: " << event << endl;
				
				N_Vtype++;
				if (sample.triggerFlags[0]){//IsoMuHLT
					N_isomu++;
				}
				if (sample.triggerFlags[6]){
					N_tightDoubleele++;
				}
				if (sample.triggerFlags[5]){
					N_loosedoubleEle++;
				}
				if (sample.triggerFlags[4]){
					N_sigleEleWP80++;
				}
				if (sample.triggerFlags[17]){
					N_singleEleDiJet++;
				}
				
				
				
				
				float HmassMaxCut = 170, HmassMinCut = 85, CSV0cut = 0.4, DphiZMETcut = 0.85, dphiHVcut = 1.6;
				float ZmassMinCut = 35, ZmassMaxCut = 170, CHFb0cut = 0.2, PzetaCut = 10;
				float MtmuCut =115, MteCut = 95;
				
				
				float Pzeta = ProjMissT-(0.25*ProjVisT);
				if ((sample.vLepton_pt[1] > 20 || sample.vLepton_pt[0] > 20) && (sample.vLepton_pt[1]>10 && sample.vLepton_pt[0]> 10)){
					//	if (sample.triggerFlags[39]||sample.triggerFlags[40]||sample.triggerFlags[41]){
					Ntrigger = Ntrigger + 1*PUweight2011;
					//cout << "PUweight2011 inside trigger requirement" << PUweight2011 << endl;
					JetDistributions("HLT", weight);
					LeptonDistributions("HLT", weight);
					TH2FDistributions("HLT", weight);
					EventShapeDistributions("HLT", weight);
					EventDistributions("HLT", weight);				
					if ( jetPt[0] > 20 && jetPt[1] > 20 && oldHmass<300 &&
						fabs(sample.vLepton_eta[0]) < 2.5 && fabs(sample.vLepton_eta[1]) < 2.4 && fabs(sample.hJet_eta[0]) < 2.5 &&
						fabs(sample.hJet_eta[1]) < 2.5 && sample.hJet_id[0]==1 && sample.hJet_id[1]==1 && sample.hbhe){
						Npreselect = Npreselect + 1*PUweight2011;
						isdata = false;
						JetDistributions("PreSelect", weight);
						LeptonDistributions("PreSelect", weight);
						TH2FDistributions("PreSelect", weight);
						EventShapeDistributions("PreSelect", weight);
						EventDistributions("PreSelect", weight);
						if (isDATA) isdata = true;
						if (fabs(DphiZMET) < DphiZMETcut && (oldHmass>=HmassMinCut)&&(oldHmass<=HmassMaxCut) && CSV0>CSV0cut && DeltaPhiHV > dphiHVcut && Pzeta > PzetaCut && Emumass< 68 && DitauMass > ZmassMinCut && DitauMass < ZmassMaxCut && jetCHF[0]>CHFb0cut) hdelRemu_NotThisCut.Fill(delRemu, weight);
						if(delRemu>0.4){
							N_EfakeCuts = N_EfakeCuts + 1*PUweight2011;
							isdata = false;
							JetDistributions("delRemu", weight);
							LeptonDistributions("delRemu", weight);
							TH2FDistributions("delRemu", weight);
							EventShapeDistributions("delRemu", weight);
							EventDistributions("delRemu", weight);
							FOM_tree->Fill();
							if (fabs(DphiZMET) < DphiZMETcut && (oldHmass>=HmassMinCut)&&(oldHmass<=HmassMaxCut) && CSV0>CSV0cut && DeltaPhiHV > dphiHVcut && Pzeta > PzetaCut && jetCHF[0]>CHFb0cut) hMemu_NotThisCut.Fill(DitauMass, weight);
							if (Emumass< 68 && DitauMass > ZmassMinCut && DitauMass < ZmassMaxCut) { 
								NMemu = NMemu + 1*PUweight2011;
								if(isDATA)isdata=true;
								JetDistributions("MemuCut", weight);
								LeptonDistributions("MemuCut", weight);
								TH2FDistributions("MemuCut", weight);
								EventShapeDistributions("MemuCut", weight);
								EventDistributions("MemuCut", weight); 
								if (fabs(DphiZMET) < DphiZMETcut && fabs(DphiSecondMET) < 1.5 && (oldHmass>HmassMinCut)&&(oldHmass<=HmassMaxCut) && DeltaPhiHV > dphiHVcut && CSV0<CSV0cut && Naj < 2 && nSV ==0){	
									N_LFCR = N_LFCR + 1*PUweight2011;
									isdata = false;
									SFUnc_LFCR = SFUnc_LFCR + weight*weight;
								}// LF Control Region	
								if (fabs(DphiZMET) < DphiZMETcut && (oldHmass>=HmassMinCut)&&(oldHmass<=HmassMaxCut) && DeltaPhiHV > dphiHVcut && jetCHF[0]>CHFb0cut && Pzeta > PzetaCut) hCSV0_NotThisCut.Fill(CSV0, weight);							
								if(CSV0>CSV0cut){
									N_CSV0 = N_CSV0 + 1*PUweight2011;
									if(isDATA)isdata=true;
									JetDistributions("CSV0", weight);
									LeptonDistributions("CSV0", weight);
									TH2FDistributions("CSV0", weight);
									EventShapeDistributions("CSV0", weight);
									EventDistributions("CSV0", weight);							
									if (fabs(DphiZMET) > DphiZMETcut && (oldHmass>HmassMinCut)&&(oldHmass<=HmassMaxCut) && DeltaPhiHV > dphiHVcut && jetCHF[0]>CHFb0cut ){
										if (CSV1 > 0.5 ){
											N_TopCR = N_TopCR + 1*PUweight2011;
											isdata = false;
											JetDistributions("TopCR", weight);
											LeptonDistributions("TopCR", weight);
											TH2FDistributions("TopCR", weight);
											EventShapeDistributions("TopCR", weight);
											EventDistributions("TopCR", weight);
											SFUnc_TOPCR = SFUnc_TOPCR + weight*weight;
										}//Top CR
										if (CSV0 > 0.679 && EvntShpAplanarity < 0.1 && CSV1 < 0.5  && Nab < 1){
											N_SingleTopCR = N_SingleTopCR + 1*PUweight2011;
											isdata = false;
											JetDistributions("SingleTopCR", weight);
											LeptonDistributions("SingleTopCR", weight);
											TH2FDistributions("SingleTopCR", weight);
											EventShapeDistributions("SingleTopCR", weight);
											EventDistributions("SingleTopCR", weight);
											SFUnc_SingleTOP = SFUnc_SingleTOP + weight*weight;
										}//Single Top CR
									}// Top and Single Top Orthogonal to Signal, reverse DphiZMET cut																		
									if ((oldHmass>=HmassMinCut)&&(oldHmass<=HmassMaxCut) && DeltaPhiHV > dphiHVcut && Pzeta > PzetaCut && jetCHF[0]>CHFb0cut) hDphiZMET_NotThisCut.Fill(DphiZMET, weight);
									if (fabs(DphiZMET) < DphiZMETcut){
										N_DphiZMET = N_DphiZMET + 1*PUweight2011;
										if(isDATA)isdata=true;
										JetDistributions("DphiZMET", weight);
										LeptonDistributions("DphiZMET", weight);
										TH2FDistributions("DphiZMET", weight);
										EventShapeDistributions("DphiZMET", weight);
										EventDistributions("DphiZMET", weight);
										if ((oldHmass>=HmassMinCut)&&(oldHmass<=HmassMaxCut) && Pzeta > PzetaCut && jetCHF[0]>CHFb0cut) hDeltaPhiHV_NotThisCut.Fill(DeltaPhiHV, weight);
										if (DeltaPhiHV > dphiHVcut){
											N_DeltaPhiHV = N_DeltaPhiHV + 1*PUweight2011;
											if(isDATA)isdata=true;
											JetDistributions("DeltaPhiHV", weight);
											LeptonDistributions("DeltaPhiHV", weight);
											TH2FDistributions("DeltaPhiHV", weight);
											EventShapeDistributions("DeltaPhiHV", weight);
											EventDistributions("DeltaPhiHV", weight);
											if ((oldHmass>=HmassMinCut)&&(oldHmass<=HmassMaxCut) && jetCHF[0]>CHFb0cut) hPzeta_NotThisCut.Fill(Pzeta, weight);
											if (Pzeta > PzetaCut){
												N_Pzeta = N_Pzeta + 1*PUweight2011;
												if(isDATA)isdata=true;
												JetDistributions("Pzeta", weight);
												LeptonDistributions("Pzeta", weight);
												TH2FDistributions("Pzeta", weight);
												EventShapeDistributions("Pzeta", weight);
												EventDistributions("Pzeta", weight);												
												if (Naj<2) hMjj_NotThisCut.Fill(oldHmass, weight);
												if((oldHmass>=HmassMinCut)&&(oldHmass<=HmassMaxCut)){
													N_Mjj = N_Mjj + 1*PUweight2011; 
													if(isDATA)isdata=true;
													JetDistributions("Mjj", weight);
													LeptonDistributions("Mjj", weight);
													TH2FDistributions("Mjj", weight);
													EventShapeDistributions("Mjj", weight);
													EventDistributions("Mjj", weight);
													if (Naj<2){
														N_Naj = N_Naj + 1*PUweight2011;
														if(isDATA)isdata=true;
														JetDistributions("Naj", weight);
														LeptonDistributions("Naj", weight);
														TH2FDistributions("Naj", weight);
														EventShapeDistributions("Naj", weight);
														EventDistributions("Naj", weight);
													}	//Naj
													if (jetCHF[0]>CHFb0cut){
														N_jetCHF0 = N_jetCHF0 + 1*PUweight2011;
														if(isDATA)isdata=true;
														JetDistributions("CHF0", weight);
														LeptonDistributions("CHF0", weight);
														TH2FDistributions("CHF0", weight);
														EventShapeDistributions("CHF0", weight);
														EventDistributions("CHF0", weight);
														if (Mte < MteCut && Mtmu < MtmuCut){
															if(isDATA)isdata=true;
															N_Mt = N_Mt + 1*PUweight2011;
															JetDistributions("Mt", weight);
															LeptonDistributions("Mt", weight);
															TH2FDistributions("Mt", weight);
															EventShapeDistributions("Mt", weight);
															EventDistributions("Mt", weight);
														}	//Mte and Mtmu
													}	//CHF0
												}//Mjj							
											} //PzetaCut
										}//DeltaPhiHV
									}// Delta Phi Z, MET (Z is emu vectors only)
								}//CSV
							}//emu mass requirement
						}//EleFakeCuts
					} else {
						FailedJetID = FailedJetID + 1*PUweight2011;
					}//Jet ID and eta requirement
				}//end trigger emulation
				isdata = false;
				EventDistributions("allEvts", weight);
			}//end requirement Zemu event
		}//end isM50sample genZpt cut
	} while (sample.nextEvent());
	
	
	
	std::cout << endl << endl;
	std::cout << "Number of events " << event << endl;
	std::cout << "Vtype 5 " << N_Vtype << endl;	
	std::cout << "emu trigger " << Ntrigger << endl;
	std::cout << "PreSelection: " << Npreselect << endl;
	std::cout << "EleFakeCuts: " << N_EfakeCuts << endl;
	std::cout << "Memu: " << NMemu << endl;	
	std::cout << "CSV0: " << N_CSV0 << endl;
	std::cout << "DphiZMET: " << N_DphiZMET << endl;
	std::cout << "DeltaPhiHV: " << N_DeltaPhiHV << endl;
	std::cout << "Pzeta: " << N_Pzeta << endl;
	std::cout << "Mjj: " << N_Mjj << endl;
	std::cout << "Naj: " << N_Naj << endl;
	std::cout << "CHF0: " << N_jetCHF0 << endl;
	std::cout << "Mte&Mtmu: " << N_Mt << endl;
	
	
	
	SFUnc_TOPCR = sqrt(SFUnc_TOPCR);
	SFUnc_SingleTOP = sqrt(SFUnc_SingleTOP);
	SFUnc_LFCR = sqrt(SFUnc_LFCR);
	SFUnc_HFCR = sqrt(SFUnc_HFCR);
	
	std::cout << endl << endl;
	std::cout << "Top CR: " << N_TopCR << " " << SFUnc_TOPCR << endl;
	std::cout << "Single Top CR: " << N_SingleTopCR<< " "<< SFUnc_SingleTOP << endl;
	std::cout << "LF CR: " << N_LFCR<< " " << SFUnc_LFCR << endl;
	std::cout << "HF CR: " << N_HFCR<< " "<< SFUnc_HFCR << endl;
	
	
	std::cout << endl << endl;
	std::cout << "Number of Events that failed JetID " << FailedJetID << endl;
	std::cout << "BDT tree: " << NBDTtree << endl;
	std::cout << "BDT  bjet tree: " << NBDTbtree << endl;
	std::cout << "TMVA tree: " << NTMVAtree << endl;
	
	std::cout << endl << endl;
	std::cout << "Number of Events with good Gen Info and Hmass > 0 " << NgoodGenInfo << endl;
	std::cout << "Number of Events with nJet[1] matched to b-bar " << Njet2Bbar << endl;
	std::cout << "Number of Events with nJet[0] matched to b-bar " << Njet1Bbar << endl;
	std::cout << "Number of Events with one OR both jets matched " << NORMatch << endl;
	std::cout << "Number of Events with both jets matched to same candidate " << NSameCandidate << endl;
	std::cout << "Number of Events with both jets matched to b candidate delR<0.3 " << Nmatch << endl;
	std::cout << "Number of Events with both jets matched to b candidate delR<0.4 " << Nmatch4 << endl;
	std::cout << "Number of Events with both jets matched to b candidate delR<0.15 " << Nmatch15 << endl;
	std::cout << "Number of Events with one jet matched  " << NORMatch-Nmatch << endl;
	std::cout << "Number of Events with good matches  " << Njet2Bbar+Njet1Bbar << endl;
	std::cout << "Number of Events in training sample  " << Ntrainingsample << endl;
	
	
	ofstream myfile;
	myfile.open(TString::Format("%sDataEntry.txt",directory.c_str()).Data());	
	myfile<<TString::Format("\t %s",directory.c_str()).Data()<<endl;
	//	myfile<<TString::Format("Number of events \t %f \t %0.5f",event,event*LumiWeight).Data()<<endl;
	//	myfile<<TString::Format("Vtype 5 \t %f \t %0.5f",N_Vtype,N_Vtype*LumiWeight).Data()<<endl;
	myfile<<TString::Format("emu trigger \t %0.5f \t %0.5f",Ntrigger,Ntrigger*LumiWeight).Data()<<endl;
	myfile<<TString::Format("PreSelection \t %0.5f \t %0.5f",Npreselect,Npreselect*LumiWeight).Data()<<endl;
	myfile<<TString::Format("EleFakeCuts \t %0.5f \t %0.5f",N_EfakeCuts,N_EfakeCuts*LumiWeight).Data()<<endl;
	myfile<<TString::Format("Memu \t %0.5f \t %0.5f",NMemu,NMemu*LumiWeight).Data()<<endl;
	myfile<<TString::Format("CSV0 \t %0.5f \t %0.5f",N_CSV0,N_CSV0*LumiWeight).Data()<<endl;
	myfile<<TString::Format("DphiZMET \t %0.5f \t %0.5f",N_DphiZMET,N_DphiZMET*LumiWeight).Data()<<endl;
	myfile<<TString::Format("DphiZH \t %0.5f \t %0.5f",N_DeltaPhiHV,N_DeltaPhiHV*LumiWeight).Data()<<endl;
	myfile<<TString::Format("Pzeta \t %0.5f \t %0.5f",N_Pzeta,N_Pzeta*LumiWeight).Data()<<endl;
	myfile<<TString::Format("Mjj \t %0.5f \t %0.5f",N_Mjj,N_Mjj*LumiWeight).Data()<<endl;
	myfile<<TString::Format("Naj \t %0.5f \t %0.5f",N_Naj,N_Naj*LumiWeight).Data()<<endl;
	myfile<<TString::Format("Mte&Mtmu \t %0.5f \t %0.5f",N_Mt,N_Mt*LumiWeight).Data()<<endl;
	myfile<<TString::Format("CHF0 \t %0.5f \t %0.5f",N_jetCHF0,N_jetCHF0*LumiWeight).Data()<<endl;
	myfile<<endl;
	myfile<<TString::Format("Top CR	\t %0.5f \t %0.5f \t %0.5f \t %0.5f",N_TopCR,N_TopCR*LumiWeight,SFUnc_TOPCR*SFUnc_TOPCR,SFUnc_TOPCR).Data()<<endl;
	myfile<<TString::Format("Single Top CR \t %0.5f \t %0.5f \t %0.5f \t %0.5f",N_SingleTopCR,N_SingleTopCR*LumiWeight,SFUnc_SingleTOP*SFUnc_SingleTOP,SFUnc_SingleTOP).Data()<<endl;
	myfile<<TString::Format("LF CR \t %0.5f \t %0.5f \t %0.5f \t %0.5f",N_LFCR,N_LFCR*LumiWeight,SFUnc_LFCR*SFUnc_LFCR,SFUnc_LFCR).Data()<<endl;
	myfile<<TString::Format("HF CR \t %0.5f \t %0.5f \t %0.5f \t %0.5f",N_HFCR,N_HFCR*LumiWeight,SFUnc_HFCR*SFUnc_HFCR,SFUnc_HFCR).Data()<<endl;
	
	myfile.close();
	
	
	
	// Here one can create canvas, draw and print the histogram.
	TCanvas c1("c1","c1");
	c1.cd();
	c1.SetFillColor(kWhite);
	string suffixps = ".gif";
	
	c1.Clear(); // don't create a new canvas
	hallhJet_pt.Draw();
	c1.Print((directory+"/PtAllJets"+suffixps).c_str());
	
	c1.Clear(); // don't create a new canvas
	hMemu_NotThisCut.Draw();
	c1.Print((directory+"/Memu_NotThisCut"+suffixps).c_str());
	
	c1.Clear(); // don't create a new canvas
	hMjj_NotThisCut.Draw();
	c1.Print((directory+"/Mjj_NotThisCut"+suffixps).c_str());
	
	c1.Clear(); // don't create a new canvas
	hDphiZMET_NotThisCut.Draw();
	c1.Print((directory+"/DphiZMET_NotThisCut"+suffixps).c_str());
	
	c1.Clear(); // don't create a new canvas
	hCSV0_NotThisCut.Draw();
	c1.Print((directory+"/CSV0_NotThisCut"+suffixps).c_str());
	
	c1.Clear(); // don't create a new canvas
	hdelRemu_NotThisCut.Draw();
	c1.Print((directory+"/delRemu_NotThisCut"+suffixps).c_str());
	
	c1.Clear(); // don't create a new canvas
	hDeltaPhiHV_NotThisCut.Draw();
	c1.Print((directory+"/DeltaPhiHV_NotThisCut"+suffixps).c_str());
	
	c1.Clear(); // don't create a new canvas
	hPzeta_NotThisCut.Draw();
	c1.Print((directory+"/Pzeta_NotThisCut"+suffixps).c_str());		
	
	
	hCutFlow.SetBinContent(1,N_Vtype );
	hCutFlow.SetBinContent(2,Ntrigger );
	hCutFlow.SetBinContent(3,Npreselect );
	hCutFlow.SetBinContent(4,N_EfakeCuts );
	hCutFlow.SetBinContent(5,NMemu );
	hCutFlow.SetBinContent(6,N_CSV0 );
	hCutFlow.SetBinContent(7,N_DphiZMET );
	hCutFlow.SetBinContent(8,N_DeltaPhiHV );
	hCutFlow.SetBinContent(9,N_Pzeta );
	hCutFlow.SetBinContent(10,N_Mjj );
	hCutFlow.SetBinContent(11,N_Naj );
	hCutFlow.Draw();
	c1.Print((directory+"/CutFlow"+suffixps).c_str());
		
	////////
	// write histos to file:
	////////
	for (map<std::string,TH1*>::iterator it=histmap.begin(); it!=histmap.end();it++) {
		(*it).second->Write();
		c1.Clear();
		(*it).second->Draw();
		c1.Print((directory+"/"+(*it).first+suffixps).c_str());
		if (debug)		cout << "This is the histmap string: " << (*it).first << endl;
	}
	for (map<string,TH2*>::iterator it2=bidimhistmap.begin(); it2!=bidimhistmap.end();it2++) {
		(*it2).second->Write();
		c1.Clear();
		(*it2).second->Draw();
		c1.Print((directory+"/"+(*it2).first+suffixps).c_str());
		delete (*it2).second;
	}
	
	
	// Write and Close the output file.
	ofile.Write();
	ofile.Close();
	
	return 0;
}


double SetWeight( std::string filename){
	double SampleWeight = 1.0;
	if (findString(filename, "ZH115_1M")){ 
		SampleWeight = lumi/lumiZH115;
	cout << "found ZH_ZToLL_HToBB_M-115 string" << endl;}
	if (findString(filename, "ZH115")){ SampleWeight = lumi/(2.19999000000000000e+05/xsecbfZH115);}
	if (findString(filename, "ZH125Fall11")){ SampleWeight = lumi/(2.16759000000000000e+05/xsecbfZH125);}
	if (findString(filename, "Filtered115emu")){ SampleWeight = lumi/lumiZH115;}
	if (findString(filename, "Filtered125emu")){ SampleWeight = lumi/lumiZH125;}
	if (findString(filename, "DYPtZ")){ SampleWeight = lumi/lumiZJH;}
	if (findString(filename, "DYM50")){ SampleWeight = lumi/lumiZJL;}
	if (findString(filename, "120to170")){ SampleWeight = lumi/lumiQCD120;}
	if (findString(filename, "170to300")){ SampleWeight = lumi/lumiQCD170;}
	if (findString(filename, "300to470")){ SampleWeight = lumi/lumiQCD300;}
	if (findString(filename, "470to600")){ SampleWeight = 70.21630282*100000/ 3598283.0;}
	if (findString(filename, "80to120")){ SampleWeight = lumi/(6397439.5000/(823744.5*1000));}
	if (findString(filename, "TTJets")){ SampleWeight = lumi/(lumiTT);}
	if (findString(filename, "T_TuneZ2_s")){ SampleWeight = lumi/lumiTs;}
	if (findString(filename, "T_TuneZ2_t-channel")){ SampleWeight = lumi/lumiTt;}
	if (findString(filename, "T_tW")){ SampleWeight = lumi/lumiTtW;}
	if (findString(filename, "T_TuneZ2_tW-channel-DS")){ SampleWeight = lumi/(828002.00/xsecbfTtW) ;}
	if (findString(filename, "Tbar_TuneZ2_s")){ SampleWeight = lumi/lumiTsb;}
	if (findString(filename, "Tbar_TuneZ2_t-channel")){ SampleWeight = lumi/lumiTtb;}
	if (findString(filename, "Tbar_tW")){ SampleWeight = lumi/lumiTtWb;}
	if (findString(filename, "Tbar_TuneZ2_tW-channel-DS")){ SampleWeight = lumi/(727800.6250000/xsecbfTtWb);}
	if (findString(filename, "WJetsToLNu_Pt-100")){ SampleWeight = lumi/(18859930.0/(685.69*1000) );}
	if (findString(filename, "WJets.root")){ SampleWeight = lumi/lumiWJ;}
	if (findString(filename, "WW")){ SampleWeight = lumi/lumiWW;}
	if (findString(filename, "WZ")){ SampleWeight = lumi/lumiWZ;}
	if (findString(filename, "ZJetsToLL")){ SampleWeight = lumi/(2349387.0/(3048*1000));}
	if (findString(filename, "ZJetsToNuNu")){ SampleWeight = lumi/(3643062.7500/(1.3*25.8*1000));}
	if (findString(filename, "ZZ")){ 
		SampleWeight = lumi/lumiZZ;
		cout << "found ZZ string" << endl;
	}
	
	return SampleWeight;
}

bool findString(std::string strToSearch, std::string strPattern){
	size_t found;
	
	bool foundStr = false;
	
	found=strToSearch.find(strPattern);
    if (found!=string::npos)
		foundStr = true;
	
	return foundStr;		
}

void fillhisto(string histname,double val, double hweight, string title, int nbins, double min, double max) {
	if (!histmap[histname]) {
		if (nbins==0) {
			cout << " ERROR: must define histogram " << histname << " properties on first call of fillhisto(...)" << endl;
			exit(1);
		}
		
		histmap[histname]=new TH1F(histname.c_str(),title.c_str(),nbins,min,max);
		
	}
	histmap[histname]->Fill(val, hweight);
}
void fill2Dhisto(string histname,double valx, double valy, double weight, string title, int nxbins, double xmin, double xmax, int nybins, double ymin, double ymax) {
	if (!bidimhistmap[histname]) {
		if (nxbins==0 || nybins==0) {
			//cout << " ERROR: must define bi dim histogram " << histname << " properties on first call of fillbidimhisto(...)" << endl;
			exit(1);
		}
		bidimhistmap[histname]=new TH2F(histname.c_str(),title.c_str(),nxbins,xmin,xmax,nybins,ymin,ymax);
		
	}
	bidimhistmap[histname]->Fill(valx,valy,weight);
}// fill2Dhisto


void JetDistributions(string cut, double ph_weight){
	
	for (int i=0; i != nJets && i < 5; ++i) {
		string jetpT = Form("hPtb%i_", i) + cut;
		string jetNewCSV = Form("hCSVNewShape%i_", i) + cut;
		string jetNewCSVAC = Form("hCSVNewShape%iAC_", i) + cut;
		string jeteta = Form("hEtab%i_", i) + cut;
		string jetphi = Form("hPhib%i_", i) + cut;
		string jetcsv = Form("hCSV%i_", i) + cut;
		string jetchf = Form("hCHFb%i_", i) + cut;
		string jetptRaw = Form("hPtRawb%i_",i)+cut;
		string jetEnergy = Form("hEnergyb%i_",i)+cut;
		string jetvtx3dL = Form("hVtx3dLb%i_",i)+cut;
		string jetvtx3deL = Form("hVtx3dLErrorb%i_",i)+cut;
		string jetvtxPt  = Form("hVtxPtb%i_",i)+cut;
		string jetvtxMass  = Form("hVtxMassb%i_",i)+cut;
		string jetptLeadTrack = Form("hptLeadTrackb%i_",i)+cut;
		string jetNconstints = Form("hNconstintuentsb%i_",i)+cut;
		string jetcef = Form("hCEFb%i_",i)+cut;
		string jetnch = Form("hNCHb%i_",i)+cut;
		string jetjecUnc = Form("hJECUncb%i_",i)+cut;
		string jetet = Form("hEtb%i_",i)+cut;
		string jetmt = Form("hMtb%i_",i)+cut;
		string jetptRawJER = Form("hptRawJERb%i_",i)+cut;
		string jetgenPt = Form("hGenPtb%i_",i)+cut;
		string svnumTracks = Form("hSVnumTracksb%i_",i)+cut;
		string jetnhf = Form("hNHFb%i_",i)+cut;
		string RecoRes = Form("hRecoResb%i_",i)+cut;
		string RecoDiff = Form("hRecoDiffb%i_",i)+cut;
		string RecoRatio = Form("hRecoRatiob%i_",i)+cut;
		
		
		fillhisto(jetpT, jetPt[i], ph_weight, "jet pT", 10, 0.0, 200);
		fillhisto(jetNewCSV, CSVNewShape[i], ph_weight, "CSV BTag Shape", 30, 0, 1.5);
		fillhisto(jetNewCSVAC, CSVNewShape[i], ph_weight, "CSV BTag Shape", 40, 0.244, 1.244);
		fillhisto(jeteta, jetEta[i], ph_weight, "jet #eta", 13, -3, 3.5);
		fillhisto(jetphi, jetPhi[i], ph_weight, "jet #phi", 20, -3.14159265, 4.7123889);
		fillhisto(jetcsv, jetCSV[i], ph_weight, "jet CSV", 30, 0, 1.5);
		if(jetCHF[i]> 0.0000001)fillhisto(jetchf, jetCHF[i], ph_weight, "charged Hadron Energy Fraction", 20, 0.0, 1.2);
		fillhisto(jetptRaw, jetPtRaw[i], ph_weight, "jet pT Raw", 10, 0.0, 200);
		fillhisto(jetEnergy, jetE[i], ph_weight, "jet energy", 10, 0.0, 200);
		fillhisto(jetvtx3dL, jetVtx3dL[i], ph_weight, "3D vertex length", 30, 0.0, 3);
		fillhisto(jetvtx3deL, jetVtx3deL[i], ph_weight, "error on vertex 3D length", 20, 0.0, .2);
		fillhisto(jetvtxPt, jetVtxPt[i], ph_weight, "vertex pT", 10, 0.0, 100);
		fillhisto(jetvtxMass, jetVtxMass[i], ph_weight, "vertex Mass", 8, 0.0, 8);
		fillhisto(jetptLeadTrack, jetPtLeadTrack[i], ph_weight, "pt of leading track", 10, 0.0, 100);
		fillhisto(jetNconstints, jetNconstintuents[i], ph_weight, "Number of Constintuents", 20,0,60);
		if(jetCEF[i]> 0.0000001)fillhisto(jetcef, jetCEF[i], ph_weight, "charged EM Energy Fraction", 20, 0.0, 1.2);
		fillhisto(jetnch, jetNCH[i], ph_weight, "Number of Charged Hadrons",25,0,50);
		fillhisto(jetjecUnc, jetJECUnc[i], ph_weight, "Jet Energy Correction Uncertainty", 20, 0, 0.1);
		fillhisto(jetet, jetEt[i], ph_weight, "Et", 10, 0.0, 200);
		fillhisto(jetmt, jetMt[i], ph_weight, "Mt", 10, 0.0, 200);
		fillhisto(jetptRawJER, jetPtRawJER[i], ph_weight, "Pt Raw with Jet Energy Resolution", 10, 0.0, 200);
		fillhisto(jetgenPt, jetGenPt[i], ph_weight, "Generated Pt", 10, 0.0, 200);
		fillhisto(svnumTracks, SVnumTracks[i], ph_weight, "number of tracks in secondary vertex", 20,0,60);
		if(jetNHF[i]> 0.0000001)fillhisto(jetnhf, jetNHF[i], ph_weight, "Neutral Hadron energy fraction", 20, 0.0, 1.2);
		fillhisto(RecoRes, (jetGenPt[i]-jetPt[i])/jetGenPt[i], ph_weight, "Reco Ptgen-Pt/Genpt", 10, -1,1);
		fillhisto(RecoDiff, jetGenPt[i]-jetPt[i], ph_weight, "Reco Ptgen-Pt", 40, -60,60);
		fillhisto(RecoRatio, jetPt[i]/jetGenPt[i], ph_weight, "Reco pt/gen ", 20, 0, 2);

	}
	
	string mjj = "hMjj_" + cut;
	string mjjReg = "hMjjReg_" + cut;
	string mjjOld = "hMjjOld_" + cut;
	string mjjAC = "hMjjAC_" + cut;
	string ptjj = "hPtjj_" +cut;
	string detajj = "hdetaJJ_" +cut;
	string scalarSumHiggsJetPt = "hScalarSumHiggsJetPt_" +cut;
	string scalarSumJetPt = "hScalarSumJetPt_" +cut;
	string Regressedjetpt0 = "hRegJetPt0_"+cut;
	string Regressedjetpt1 = "hRegJetPt1_"+cut;
	string rho25 = "hRho25_"+cut;
	string hminDeltaPhijetMET = "hminDeltaPhijetMET_" +cut;
	string hminMET = "hMinMET_"+cut;
	string hminDeltaPhijetMETOR30aJET = "hMindPhiMEThJetOR30aJet"+cut;
	string SVmass0 = "hSVmass0"+cut;
	string SVmass1 = "hSVmass1"+cut;
	string SVmass2 = "hSVmass2"+cut;
	string RegRes0 = "hRegResb0_"+cut;
	string RegRes1 = "hRegResb1_"+cut;
	string RegDiff0 = "hRegDiffb0_"+cut;
	string RegDiff1 = "hRegDiffb1_"+cut;
	string RegRatio0 = "hRegRatiob0_"+cut;
	string RegRatio1 = "hRegRatiob1_"+cut;
	string RegRes = "hRegressionResolutionMjj_"+cut;
	string RecoRes = "hRecoResolutionMjj_"+cut;
	string RegDiff = "hRegressionDiffMjj_"+cut;
	string RecoDiff = "hRecoDiffMjj_"+cut;
	string RegRatio = "hRegressionRatioMjj_"+cut;
	string RecoRatio = "hRecoRatioMjj_"+cut;
	string Recob0Res = "hRecoResolutionb0_"+cut;
	string Recob1Res = "hRecoResolutionb1_"+cut;
	string Regb0Res = "hRegResolutionb0_"+cut;
	string Regb1Res = "hRegResolutionb1_"+cut;
	

	
	if (!isdata || (Hmass < 90 || Hmass > 150)) fillhisto(mjjReg, RegHmass, ph_weight, "Regressed Mass", 25, 0, 250);
	if (!isdata || (Hmass < 90 || Hmass > 150)) fillhisto(mjjOld, oldHmass, ph_weight, "Old Higgs Mass", 25, 0, 250);
	if (!isdata || (Hmass < 90 || Hmass > 150)) fillhisto(mjj, Hmass, ph_weight, "Invariant Mass of two Jets", 25, 0, 250);
	if (!isdata || (Hmass < 90 || Hmass > 150)) fillhisto(mjjAC, Hmass, ph_weight, "Invariant Mass of two Jets", 20, 75, 175);
	fillhisto(ptjj, Hpt, ph_weight, "Pt of two b jets with highest CSV", 25, 0, 250);
	fillhisto(detajj, fabs(DetaJJ), ph_weight, "Delta eta between two jets", 10, 0, 5);
	fillhisto(scalarSumHiggsJetPt, ScalarSumHiggsJetPt, ph_weight, "scalar sum higgs jet pt", 25, 25, 375);
	fillhisto(scalarSumJetPt, ScalarSumJetPt, ph_weight, "scalar sum alljet pt", 25, 30, 300);
	fillhisto(Regressedjetpt0, RegjetPt[0], ph_weight, "Regressed jet0 pT", 10, 0.0, 200);
	fillhisto(Regressedjetpt1, RegjetPt[1], ph_weight, "Regressed jet1 pT", 10, 0.0, 200);
	fillhisto(rho25, Rho25, ph_weight, "Energy density eta<2.5", 10, 0.0, 20);	
	fillhisto(hminDeltaPhijetMET,DeltaPhijetMETmin, ph_weight, "Delta phi between MET and nearest jet", 16, 0, 3.14159265);
	fillhisto(hminMET, MinMET, ph_weight, "min(METtype1corr,METnoPUCh)",		13, 0.0, 260);
	fillhisto(hminDeltaPhijetMETOR30aJET,MindPhiMEThJetOR30aJet, ph_weight, "min(DeltaPhijetMET,Dphiajetpt30MET)", 16, 0, 3.14159265);
	fillhisto(SVmass0,SVmass[0], ph_weight, "secondary vertex mass 0",8, 0.0, 8);
	fillhisto(SVmass1,SVmass[1], ph_weight, "secondary vertex mass 1",8, 0.0, 8);
	fillhisto(SVmass2,SVmass[2], ph_weight, "secondary vertex mass 2",8, 0.0, 8);
	fillhisto(RegRes0, (jetGenPt[0]-RegjetPt[0])/jetGenPt[0], ph_weight, "Regressed Ptgen-Pt/Genpt leading jet", 10, -1,1);
	fillhisto(RegRes1, (jetGenPt[1]-RegjetPt[1])/jetGenPt[1], ph_weight, "Regressed Ptgen-Pt/Genpt subleading jet", 10, -1,1);
	fillhisto(RegRatio0, jetGenPt[0]/RegjetPt[0], ph_weight, "Regressed pt/gen leading jet", 20, 0, 2);
	fillhisto(RegRatio1, jetGenPt[1]/RegjetPt[1], ph_weight, "Regressed pt/gen subleading jet", 20, 0, 2);
	fillhisto(RegDiff0, jetGenPt[0]-RegjetPt[0], ph_weight, "Regressed Ptgen-Pt leading jet", 40, -60,60);
	fillhisto(RegDiff1, jetGenPt[1]-RegjetPt[1], ph_weight, "Regressed Ptgen-Pt subleading jet", 40, -60,60);
			  if (!isdata) fillhisto(RegRes, (GenHiggsMass-RegHmass)/GenHiggsMass, ph_weight, "Resoltuion Higgs Regressed", 10, -1, 1);
			  if (!isdata) fillhisto(RecoRes, (GenHiggsMass-oldHmass)/GenHiggsMass, ph_weight, "Resolution Higgs Reco", 10, -1, 1);
			  if (!isdata) fillhisto(RegDiff, GenHiggsMass-RegHmass, ph_weight, "Diff Higgs Regressed", 80, -90, 90);
			  if (!isdata) fillhisto(RecoDiff, GenHiggsMass-oldHmass, ph_weight, "Diff Higgs Reco", 80, -90, 90);
			  if (!isdata) fillhisto(RegRatio, RegHmass/GenHiggsMass, ph_weight, "Ratio Higgs Regressed", 20, 0, 2);
			  if (!isdata) fillhisto(RecoRatio, oldHmass/GenHiggsMass, ph_weight, "Ratio Higgs Reco", 20, 0, 2);

}// JetDistributions


void LeptonDistributions(string cut, double ph_weight){
	
	for (int i=0; i != nLeptons && i < 10; ++i) {
		string leptonpT = Form("hPtlep%i_", i) + cut;
		string leptoneta = Form("hEtalep%i_", i) + cut;
		string leptonphi = Form("hPhilep%i_", i) + cut;
		string leptonmupt = Form("hPtMu%i_", i) + cut;
		string leppfCombRelIso = Form("hPFRelIso%i_", i) + cut;
		
		fillhisto(leptonpT, leptonPt[i], ph_weight, "Lepton pT", 10, 0.0, 100);
		fillhisto(leptoneta, leptonEta[i], ph_weight, "Lepton #eta", 13, -3, 3.5);
		fillhisto(leptonphi, leptonPhi[i], ph_weight, "Lepton #phi", 20, -3.14159265, 4.7123889);
		fillhisto(leptonmupt, SortedMuonPt[i], ph_weight, "Lepton pT", 30, 0.0, 150);
		if (lep_pfCombRelIso[i] > 0 ) fillhisto(leppfCombRelIso, lep_pfCombRelIso[i], ph_weight, "PF Rel Iso of Lepton", 25, 0, 0.5);
		
	}
	
	string memu = "hMemu_" + cut;
	string memuAC = "hMemuAC_" + cut;
	string ptemu = "hPtemu_" +cut;
	string LeadLepPt = "hLeadLepPt_" +cut;
	string SecLepPt = "hSecLepPt_" +cut;
	
	fillhisto(memu, Emumass, ph_weight, "Invariant Mass of two Leptons", 20, 0, 100);
	fillhisto(memuAC, Emumass, ph_weight, "Invariant Mass of two Leptons ", 40, 10, 80);
	fillhisto(ptemu, Zpt, ph_weight, "Pt of two Leptons", 10, 0, 150);
	if(leptonPt[0]> leptonPt[1]){
		fillhisto(LeadLepPt, leptonPt[0], ph_weight, "Leading Lepton Pt", 10, 0.0, 100);
		fillhisto(SecLepPt, leptonPt[1], ph_weight, "Second Lepton Pt", 10, 0.0, 100);
	}else{
		fillhisto(LeadLepPt, leptonPt[1], ph_weight, "Leading Lepton Pt", 10, 0.0, 100);
		fillhisto(SecLepPt, leptonPt[0], ph_weight, "Second Lepton Pt", 10, 0.0, 100);
	}
	
}// LeptonDistributions	

void TH2FDistributions(string cut, double ph_weight){
	
	string hRecoPt0vsGenPt = "hRecoPt0vsGenPt_"+cut;
	string hRecoPt1vsGenPt = "hRecoPt1vsGenPt_"+cut;
	string hRegPt0vsGenPt = "hRegPt0vsGenPt_"+cut;
	string hRegPt1vsGenPt = "hRegPt1vsGenPt_"+cut;
	string hPt0vsgenReco = "hPt0vsgenReco_"+cut;
	string hPt1vsgenReco = "hPt1vsgenReco_"+cut;
	string hRawPt0vsgenReco = "hRawPt0vsgenReco_"+cut;
	string hRawPt1vsgenReco = "hRawPt1vsgenReco_"+cut;
	string hEta0vsgenReco = "hEta0vsgenReco_"+cut;
	string hEta1vsgenReco = "hEta1vsgenReco_"+cut;
	string hEt0vsgenReco = "hEt0vsgenReco_"+cut;
	string hEt1vsgenReco = "hEt1vsgenReco_"+cut;
	string hMt0vsgenReco = "hMt0vsgenReco_"+cut;
	string hMt1vsgenReco = "hMt1vsgenReco_"+cut;
	string hPtleadTrack0vsgenReco = "hPtleadTrack0vsgenReco_"+cut;
	string hPtleadTrack1vsgenReco = "hPtleadTrack1vsgenReco_"+cut;
	string hCHF0vsgenReco = "hCHF0vsgenReco_"+cut;
	string hCHF1vsgenReco = "hCHF1vsgenReco_"+cut;
	string hCEF0vsgenReco = "hCEF0vsgenReco_"+cut;
	string hCEF1vsgenReco = "hCEF1vsgenReco_"+cut;
	string hNCH0vsgenReco = "hNCH0vsgenReco_"+cut;
	string hNCH1vsgenReco = "hNCH1vsgenReco_"+cut;
	string hJECUnc0vsgenReco = "hJECUnc0vsgenReco_"+cut;
	string hJECUnc1vsgenReco = "hJECUnc1vsgenReco_"+cut;
	string hvtxPt0vsgenReco = "hvtxPt0vsgenReco_"+cut;
	string hvtxPt1vsgenReco = "hvtxPt1vsgenReco_"+cut;
	string hvtxMass0vsgenReco = "hvtxMass0vsgenReco_"+cut;
	string hvtxMass1vsgenReco = "hvtxMass1vsgenReco_"+cut;
	string hvtx3dL0vsgenReco = "hvtx3dL0vsgenReco_"+cut;
	string hvtx3dL1vsgenReco = "hvtx3dL1vsgenReco_"+cut;
	string hvtx3deL0vsgenReco = "hvtx3deL0vsgenReco_"+cut;
	string hvtx3deL1vsgenReco = "hvtx3deL1vsgenReco_"+cut;
	string hEnergy0vsgenReco = "hEnergy0vsgenReco_"+cut;
	string hEnergy1vsgenReco = "hEnergy1vsgenReco_"+cut;
	string hNconsint0vsgenReco = "hNconsint0vsgenReco_"+cut;
	string hNconsint1vsgenReco = "hNconsint1vsgenReco_"+cut;
	string hRawJER0vsgenReco = "hRawJER0vsgenReco_"+cut;
	string hRawJER1vsgenReco = "hRawJER1vsgenReco_"+cut;
	string hRho25vs0genReco = "hRho25vs0genReco_"+cut;
	string hRho25vs1genReco = "hRho25vs1genReco_"+cut;
	string hCSV0vsgenReco = "hCSV0vsgenReco_"+cut;
	string hCSV1vsgenReco = "hCSV1vsgenReco_"+cut;

	
	fill2Dhisto(hRecoPt0vsGenPt, jetGenPt[0], jetPt[0], ph_weight, "RecoPt vs GenPt jet0", 100, 0, 200, 100, 0.0, 200);
	fill2Dhisto(hRecoPt1vsGenPt, jetGenPt[1], jetPt[1], ph_weight, "RecoPt vs GenPt jet1", 100, 0, 200, 100, 0.0, 200);
	fill2Dhisto(hRegPt0vsGenPt, jetGenPt[0], RegjetPt[0], ph_weight, "RegPt vs GenPt jet0", 100, 0, 200, 100, 0.0, 200);
	fill2Dhisto(hRegPt1vsGenPt, jetGenPt[1], RegjetPt[1], ph_weight, "RegPt vs GenPt jet1", 100, 0, 200, 100, 0.0, 200);
	fill2Dhisto(hPt0vsgenReco, jetGenPt[0]-jetPt[0], jetPt[0] , ph_weight, "Reco Pt vs gen-Reco jet0 ", 120, -60, 60, 75, 0.0, 150);
	fill2Dhisto(hPt1vsgenReco, jetGenPt[1]-jetPt[1], jetPt[1] , ph_weight, "Reco Pt vs gen-Reco jet1 ", 120, -60, 60, 75, 0.0, 150);
	fill2Dhisto(hRawPt0vsgenReco, jetGenPt[0]-jetPt[0], jetPtRaw[0] , ph_weight, "RawPt vs gen-Reco jet0 ", 120, -60, 60, 75, 0.0, 150);
	fill2Dhisto(hRawPt1vsgenReco, jetGenPt[1]-jetPt[1], jetPtRaw[1] , ph_weight, "RawPt vs gen-Reco jet1 ", 120, -60, 60, 75, 0.0, 150);
	fill2Dhisto(hEta0vsgenReco, jetGenPt[0]-jetPt[0], jetEta[0], ph_weight, "#eta vs gen-Reco jet0 ", 120, -60, 60,	36, -3, 3);
	fill2Dhisto(hEta1vsgenReco, jetGenPt[1]-jetPt[1], jetEta[1] , ph_weight, "#eta vs GenPt jet1", 120, -60, 60, 36, -3, 3);
	fill2Dhisto(hEt0vsgenReco, jetGenPt[0]-jetPt[0], jetEt[0], ph_weight, "E_{T} vs gen-Reco jet0 ", 120, -60, 60, 75, 0.0, 150);
	fill2Dhisto(hEt1vsgenReco, jetGenPt[1]-jetPt[1], jetEt[1] , ph_weight, "E_{T} vs GenPt jet1", 120, -60, 60, 75, 0.0, 150);
	fill2Dhisto(hMt0vsgenReco, jetGenPt[0]-jetPt[0], jetMt[0], ph_weight, "M_{T} vs gen-Reco jet0 ", 120, -60, 60, 75, 0.0, 150);
	fill2Dhisto(hMt1vsgenReco, jetGenPt[1]-jetPt[1], jetMt[1] , ph_weight, "M_{T} vs GenPt jet1", 120, -60, 60, 75, 0.0, 150);
	fill2Dhisto(hPtleadTrack0vsgenReco, jetGenPt[0]-jetPt[0], jetPtLeadTrack[0], ph_weight, "pt of leading track vs gen-Reco jet0 ", 120, -60, 60, 101, 0.0, 100	);
	fill2Dhisto(hPtleadTrack1vsgenReco, jetGenPt[1]-jetPt[1], jetPtLeadTrack[1] , ph_weight, "pt of leading track vs GenPt jet1", 120, -60, 60,	101, 0.0, 100);
	fill2Dhisto(hCHF0vsgenReco, jetGenPt[0]-jetPt[0], jetCHF[0], ph_weight, "Charged Hadron Fraction vs gen-Reco jet0 ", 120, -60, 60, 101, 0.0, 1.0);	
	fill2Dhisto(hCHF1vsgenReco, jetGenPt[1]-jetPt[1], jetCHF[1], ph_weight, "Charged Hadron Fraction vs GenPt jet1", 120, -60, 60, 101, 0.0, 1.0);
	fill2Dhisto(hCEF0vsgenReco, jetGenPt[0]-jetPt[0], jetCEF[0], ph_weight, "charged EM Energy Fraction vs gen-Reco jet0 ", 120, -60, 60, 101, 0.0, 1.0);
	fill2Dhisto(hCEF1vsgenReco, jetGenPt[1]-jetPt[1], jetCEF[1] , ph_weight, "charged EM Energy Fraction vs gen-Reco jet1 ", 120, -60, 60, 101, 0.0, 1.0);
	fill2Dhisto(hNCH0vsgenReco, jetGenPt[0]-jetPt[0], jetNCH[0], ph_weight, "Number of Charged Hadrons vs gen-Reco jet0 ", 120, -60, 60,50,0,50);
	fill2Dhisto(hNCH1vsgenReco, jetGenPt[1]-jetPt[1], jetNCH[1] , ph_weight, "Number of Charged Hadrons vs gen-Reco jet1 ", 120, -60, 60,50,0,50);
	fill2Dhisto(hJECUnc0vsgenReco, jetGenPt[0]-jetPt[0], jetJECUnc[0], ph_weight, "Jet Energy Correction Uncertainty vs gen-Reco jet0 ", 120, -60, 60, 51, 0, 0.05);
	fill2Dhisto(hJECUnc1vsgenReco, jetGenPt[1]-jetPt[1], jetJECUnc[1] , ph_weight, "Jet Energy Correction Uncertainty vs gen-Reco jet1 ", 120, -60, 60, 51, 0, 0.05);
	fill2Dhisto(hvtxPt0vsgenReco, jetGenPt[0]-jetPt[0], jetVtxPt[0], ph_weight, "vertex pT vs gen-Reco jet0 ", 120, -60, 60, 50, 0.0, 50);
	fill2Dhisto(hvtxPt1vsgenReco, jetGenPt[1]-jetPt[1], jetVtxPt[1] , ph_weight, "vertex pT vs gen-Reco jet1 ", 120, -60, 60, 50, 0.0, 50);
	fill2Dhisto(hvtxMass0vsgenReco, jetGenPt[0]-jetPt[0], jetVtxMass[0], ph_weight, "vertex Mass vs gen-Reco jet0 ", 120, -60, 60, 50, 0.0, 5);
	fill2Dhisto(hvtxMass1vsgenReco, jetGenPt[1]-jetPt[1], jetVtxMass[1] , ph_weight, "vertex Mass vs gen-Reco jet1 ", 120, -60, 60, 50, 0.0, 5);
	fill2Dhisto(hvtx3dL0vsgenReco, jetGenPt[0]-jetPt[0], jetVtx3dL[0], ph_weight, "3D vertex length vs gen-Reco jet0 ", 120, -60, 60, 51, 0.0, 3);
	fill2Dhisto(hvtx3dL1vsgenReco, jetGenPt[1]-jetPt[1], jetVtx3dL[1] , ph_weight, "3D vertex length vs gen-Reco jet1 ", 120, -60, 60, 51, 0.0, 3);
	fill2Dhisto(hvtx3deL0vsgenReco, jetGenPt[0]-jetPt[0], jetVtx3deL[0], ph_weight,"error on vertex 3D length vs gen-Reco jet0 ", 120, -60, 60, 51, 0.0, .2);
	fill2Dhisto(hvtx3deL1vsgenReco, jetGenPt[1]-jetPt[1], jetVtx3deL[1] , ph_weight, "error on vertex 3D length vs gen-Reco jet1 ", 120, -60, 60, 51, 0.0, .2);
	fill2Dhisto(hEnergy0vsgenReco, jetGenPt[0]-jetPt[0], jetE[0], ph_weight, "jet energy vs gen-Reco jet0 ", 120, -60, 60, 75, 0.0, 150);
	fill2Dhisto(hEnergy1vsgenReco, jetGenPt[1]-jetPt[1], jetE[1] , ph_weight, "jet energy vs gen-Reco jet1 ", 120, -60, 60, 75, 0.0, 150);
	fill2Dhisto(hNconsint0vsgenReco, jetGenPt[0]-jetPt[0], jetNconstintuents[0], ph_weight, "Number of Constintuents vs gen-Reco jet0 ", 120, -60, 60, 60,0,60);
	fill2Dhisto(hNconsint1vsgenReco, jetGenPt[1]-jetPt[1], jetNconstintuents[1] , ph_weight,"Number of Constintuents vs gen-Reco jet1 ", 120, -60, 60, 60,0,60);
	fill2Dhisto(hRawJER0vsgenReco, jetGenPt[0]-jetPt[0], jetPtRawJER[0], ph_weight, "Pt Raw with Jet Energy Resolution vs gen-Reco jet0 ", 120, -60, 60, 101, 0.0, 100);
	fill2Dhisto(hRawJER1vsgenReco, jetGenPt[1]-jetPt[1], jetPtRawJER[1] , ph_weight, "Pt Raw with Jet Energy Resolution vs gen-Reco jet1 ", 120, -60, 60, 101, 0.0, 100);
	fill2Dhisto(hRho25vs0genReco, jetGenPt[0]-jetPt[0], Rho25, ph_weight, "#rho |#eta| < 2.5 vs gen-Reco jet0 ", 120, -60, 60,	51, 0.0, 10);	
	fill2Dhisto(hRho25vs1genReco, jetGenPt[1]-jetPt[1], Rho25, ph_weight, "#rho |#eta| < 2.5 vs gen-Reco jet1 ", 120, -60, 60,	51, 0.0, 10);
	fill2Dhisto(hCSV0vsgenReco, jetGenPt[0]-jetPt[0], CSVNewShape[0], ph_weight, "csv vs gen-Reco jet0 ", 120, -60, 60, 101, 0.0, 1.0);	
	fill2Dhisto(hCSV1vsgenReco, jetGenPt[1]-jetPt[1], CSVNewShape[1], ph_weight, "csv vs GenPt jet1", 120, -60, 60, 101, 0.0, 1.0);

}// TH2FDistributions		

void EventShapeDistributions(string cut, double ph_weight){
	
	string hPtbalZH = "hPtbalZH_" + cut;
	string hPtbalMETH = "hPtbalMETH_" + cut;
	string hPtbalZMET = "hPtbalZMET_" + cut;
	string hdphiVH = "hdphiVH_" +cut;
	string hdphiVHAC = "hdphiVHAC_" +cut;
	string hRMSeta = "hRMSeta_" +cut;
	string hStaDeveta = "hStaDeveta_" +cut;
	string hUnweightedEta = "hUnweightedEta_" +cut;
	string hScalarSumPt = "hScalarSumPt_" +cut;
	string hCentrality = "hCentrality_" +cut;
	string hEventPt = "hEventPt_" +cut;
	string hEventMass = "hEventMass_" +cut;
	string hAngleHemu = "hAngleHemu_" +cut;
	string hSphericity = "hSphericity_" +cut;
	string hAplanarity = "hAplanarity_" +cut;
	string hCircularity = "hCircularity_" +cut;
	string hIsotropy = "hIsotropy_" +cut;
	string dphijj = "hdphiJJ_vect_" +cut;
	string hHt = "hHt_" +cut;
	string hZmass = "hZmass_" +cut;
	string hZmassSVD = "hZmassSVD_" +cut;
	string hZmassSVDnegSol = "hZmassSVDnegSol_" +cut;
	string hZmassSVDAC = "hZmassSVDAC_" +cut;
	string hZmassNegInclu = "hZmassNegInclu_" +cut;
	string hAngleEMU = "hAngleEMU_" +cut;
	string hCosThetaEle = "hCosThetaEle_" +cut;
	string hCosThetaMu = "hCosThetaMu_" +cut;
	string hEleMissE = "hEleMissE_" +cut;
	string hMuonMissE = "hMuonMissE_" +cut;
	string hDphiemu = "hDphiemu_" +cut;
	string hDetaemu = "hDetaemu_" +cut;
	string hMassEleb0 = "hMassEleb0_" +cut;
	string hMassMub0 = "hMassMub0_" +cut;
	string hMassEleb1 = "hMassEleb1_" +cut;
	string hMassMub1 = "hMassMub1_" +cut;
	string hDphiEleMET = "hDphiEleMET_" +cut;
	string hdphiMuMET = "hdphiMuMET_" +cut;
	string hDphiLeadMET = "hDphiLeadMET_" +cut;
	string hDphiSecondMET = "hDphiSecondMET_" +cut;
	string hDphiZMET = "hDphiZMET_" +cut;
	string hDphiZMETAC = "hDphiZMETAC_" +cut;
	string hdelRjj = "hdelRjj_" +cut;
	string hdelRemu = "hdelRemu_" +cut;
	string hProjVisT = "hProjVisT_" +cut;
	string hProjMissT = "hProjMissT_" +cut;
	string hPzetaCut = "hPzetaCut_" +cut;
	
	
	fillhisto(hPtbalZH, PtbalZH, ph_weight, "Pt balance of Z and H", 20, -75, 175);
	fillhisto(hPtbalMETH, PtbalMETH, ph_weight, "Pt balance of MET and H", 20, -75, 175);
	fillhisto(hPtbalZMET, PtbalZMET, ph_weight, "Pt balance of Z and MET", 25, -100, 125);
	fillhisto(hdphiVH, DeltaPhiHV, ph_weight, "Delta phi between Z and Higgs", 24, 0, 4.71238898);
	fillhisto(hdphiVHAC, DeltaPhiHV, ph_weight, "Delta phi between Z and Higgs", 44, 0.785398, 4.71238898);
	fillhisto(hRMSeta, RMS_eta, ph_weight, "RMS Eta", 22, 0, 2.2);
	fillhisto(hStaDeveta, EtaStandDev, ph_weight, "Standard Deviation Eta",		40, 0, 2);
	fillhisto(hUnweightedEta, UnweightedEta, ph_weight, "Unweighted Eta",		45, 0, 7);
	fillhisto(hScalarSumPt, ScalarSumPt, ph_weight, "scalar sum of pt of four particles", 25, 75, 375);
	fillhisto(hCentrality, Centrality, ph_weight, "Centrality", 30, 0.0, 0.70);
	fillhisto(hEventPt, EventPt, ph_weight, "Pt of HV system", 25, 0.0, 250);
	fillhisto(hEventMass, EventMass, ph_weight, "Mass of HV system", 30, 100, 400);
	fillhisto(hAngleHemu, AngleHemu, ph_weight, "Angle between H and Z", 15, 0, 3.5);
	fillhisto(hSphericity, EvntShpSphericity, ph_weight,"EventShapeVariables sphericity", 25, 0, 1);
	fillhisto(hAplanarity, EvntShpAplanarity, ph_weight, "EventShapeVariables Aplanarity", 30, 0.0, .3);
	fillhisto(hCircularity, EvntShpCircularity, ph_weight,  "EventShapeVariables circularity", 20, 0.0, 1.2);
	fillhisto(hIsotropy, EvntShpIsotropy, ph_weight,  "EventShapeVariables isotropy", 40, 0.05, 1.0);
	fillhisto(dphijj, fabs(DphiJJ), ph_weight, "Delta phi between two jets",  24, 0, 4.71238898);
	if (Ht>1)fillhisto(hHt, Ht, ph_weight, "from MHT class MHT_ht", 20, 0, 200);
	if(Zmass>-98) fillhisto(hZmass, Zmass, ph_weight, "Invariant Mass of two Leptons corrected", 20, 0, 100);
	if (ZmassSVD>-98) fillhisto(hZmassSVD, ZmassSVD, ph_weight, "Invariant Mass of two Leptons corrected SVD", 20, 0, 100);
	if (ZmassSVD>-98) fillhisto(hZmassSVDAC, ZmassSVD, ph_weight, "Invariant Mass of two Leptons corrected SVD", 20, 0, 100);
	if (ZmassSVDnegSol>-98) fillhisto(hZmassSVDnegSol, ZmassSVDnegSol, ph_weight, "Invariant Mass of two Leptons corrected SVD", 20, 0, 100);
	if(ZmassNegInclu>-98)fillhisto(hZmassNegInclu, ZmassNegInclu, ph_weight, "Invariant Mass of two Leptons corrected Matrix", 20, 0, 100);
	fillhisto(hAngleEMU, AngleEMU, ph_weight, "Angle between electron and muon", 15, 0, 3.5);
	fillhisto(hCosThetaMu,CosThetaMu, ph_weight, "Cos Theta Muon", 30, -1, 2.0);
	fillhisto(hCosThetaEle,CosThetaEle, ph_weight, "Cos Theta Electron", 30, -1, 2.0);
	fillhisto(hEleMissE, AaronEleMissE, ph_weight, "Missing Energy electron", 40, -400, 400);
	fillhisto(hMuonMissE, AaronMuMissE, ph_weight, "Missing Energy muon", 40, -400, 400);
	fillhisto(hDphiemu, fabs(Dphiemu), ph_weight, "Delta phi between e and muon", 24, 0, 4.71238898);
	fillhisto(hDetaemu, fabs(Detaemu), ph_weight, "Delta eta between e and muon", 10, 0, 5);
	if(MassEleb0>1)fillhisto(hMassEleb0, MassEleb0, ph_weight, "Invariant Mass of Electron and b0", 20, 0, 200);
	if(MassMub0>1)fillhisto(hMassMub0, MassMub0, ph_weight, "Invariant Mass of Muon and b0", 20, 0, 200);
	if(MassEleb1>1)fillhisto(hMassEleb1, MassEleb1, ph_weight, "Invariant Mass of Electron and b1", 20, 0, 200);
	if(MassMub1>1)fillhisto(hMassMub1, MassMub1, ph_weight, "Invariant Mass of Muon and b1", 20, 0, 200);
	fillhisto(hDphiEleMET, fabs(DphiEleMET), ph_weight, "Delta phi between Electron and MET",  24, 0, 4.71238898);
	fillhisto(hdphiMuMET,fabs( dphiMuMET), ph_weight, "Delta phi between Muon and MET",  24, 0, 4.71238898);
	fillhisto(hDphiLeadMET, fabs(DphiLeadMET), ph_weight, "Delta phi between leading lepton and MET",  24, 0, 4.71238898);
	fillhisto(hDphiSecondMET, fabs(DphiSecondMET), ph_weight, "Delta phi between second lepton and MET",  24, 0, 4.71238898);
	fillhisto(hDphiZMET, fabs(DphiZMET), ph_weight, "Delta phi between Z and MET",  24, 0, 4.71238898);
	fillhisto(hDphiZMETAC, fabs(DphiZMET), ph_weight, "Delta phi between Z and MET",  24, 0, 3.14159265);
	fillhisto(hdelRjj, delRjj, ph_weight, "Delta R jj", 20, 0, 5);
	fillhisto(hdelRemu, delRemu, ph_weight, "Delta R emu", 20, 0, 5);
	fillhisto(hProjVisT, ProjVisT, ph_weight, "Transverse componenet of Projection of Z onto bisector", 20, 0, 200);
	fillhisto(hProjMissT, ProjMissT, ph_weight, "Transeverse compenent of Projection of MET onto e,mu bisector", 30, -100, 200);
	fillhisto(hPzetaCut, ProjMissT-(0.25*ProjVisT), ph_weight, "P_#zeta Cut", 25, -50, 150);
	
	
	
	
	
}// EventShapeDistributions


void EventDistributions(string cut, double ph_weight){
	
	string hnJets = "hnJets_" + cut;
	string hnMuons = "hnMuons_" +cut;
	string hnElectron = "hnElectrons_" +cut;
	string hnLeptons = "hnLeptons_" +cut;
	string hnSV = "hnSV_" +cut;
	string hnPV = "hnPV_" +cut;
	string hnaJets = "hnaJets_" +cut;
	string hMET = "hMET_" +cut;
	string hMETsig = "hMETsig_" +cut;
	string hNaj = "hNaj_" +cut;
	string hNab = "hNab_" +cut;
	string Nalep = "hNalep_" +cut;
	string hHMETdPhi = "hHMETdPhi_" +cut;
	string hVMt = "hVMt_" +cut;
	string hminDeltaPhijetMETZtau = "hminDeltaPhijetMETZtau_" +cut;
	string hSVmass = "hSVmass_" +cut;
	string hVMte = "hVMte_" +cut;
	string hVMtmu = "hVMtmu_" +cut;
	string hdeltaPullAngle = "hdeltaPullAngle_" +cut;
	string hdeltaPullAngle2 = "hdeltaPullAngle2_" +cut;
	string htopMass = "htopMass_" +cut;
	string htopPt = "htopPt_" +cut;
	string htopWmass = "htopWmass_" +cut;
	
	
	fillhisto(hnJets, nJets, ph_weight, "Number of Good Jets", 13, -0.5, 12.5);
	fillhisto(hnMuons, nMuons, ph_weight, "Number of Good Muons", 6, -0.5, 5.5);
	fillhisto(hnElectron, nElectrons, ph_weight, "Number of Good Electrons", 6, -0.5, 5.5);
	fillhisto(hnLeptons, nLeptons, ph_weight, "Number of Good Leptons", 6, -0.5, 5.5);
	fillhisto(hnSV, nSV, ph_weight, "Number of Secondary Verticies", 4, -0.5, 3.5);
	fillhisto(hnPV, nPV, ph_weight, "Number of Primary Verticies", 26, -0.5, 25.5);
	fillhisto(hnaJets, naJets, ph_weight, "Number of Aditional Jets", 11, -0.5, 10.5);
	fillhisto(hMET, MET, ph_weight, "Missing Et",		13, 0.0, 260);
	fillhisto(hMETsig, METsig, ph_weight, "Missing Et significance",		15, 0, 75);
	fillhisto(hNaj, Naj, ph_weight, "Number of Additional Jets",		7, -0.5, 6.5);
	fillhisto(hNab, Nab, ph_weight, "Number of Additional b jets",		5, -0.5, 4.5);
	fillhisto(Nalep, Na_lep, ph_weight, "Number of Additional Leptons",		5, -1.5, 4.5);
	fillhisto(hHMETdPhi, dPhiHMET, ph_weight, "Delta phi between MET and Higgs", 24, 0, 4.71238898);
	fillhisto(hVMt, Mt, ph_weight, "VMt", 16, 0, 160);
	fillhisto(hminDeltaPhijetMETZtau,DeltaPhijetMETZtaumin, ph_weight, "Delta phi between MET+Zpt and nearest jet", 24, 0, 4.71238898);
	fillhisto(hSVmass, SV_mass, ph_weight, "mass of Secondary Vertex",		20, 0.0, 5);
	fillhisto(hVMte, Mte, ph_weight, "VMte",  16, 0, 160);
	fillhisto(hVMtmu, Mtmu, ph_weight, "VMtmu",  16, 0, 160);
	fillhisto(hdeltaPullAngle,delPullAngle, ph_weight, "Delta pull Angle", 42, -3.25, 5.25);
	fillhisto(hdeltaPullAngle2,delPullAngle2, ph_weight, "Delta Pull Angle 2", 42, -3.25, 5.25);
	fillhisto(htopMass,topMass, ph_weight, "Top Mass", 30, 75, 375);
	fillhisto(htopPt,topPt, ph_weight, "Pt of Top", 40, 0, 200);
	if (topWmass > 82 || topWmass < 80) fillhisto(htopWmass,topWmass, ph_weight, "Mass of W coming from Top", 25, 75, 125);
	
	
	
}// EventDistributions	
