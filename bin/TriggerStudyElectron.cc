/**
 @brief Example of an analysis code to read tree ntuple
 and create histogram
 Usage:
 
 \verbatim
 mkdir [directoryName]
 TriggerStudyElectron <inputfile> [outputfilename] [directoryName]
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
    return result;
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
    std::string ofilename("trigEle.root");
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
	bool isDATA = false;
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
	
	float ptMin = -99.99, ptMax = -99.99, etaMin = -99.99, etaMax = -99.99, scale = -99.99, error = -99.99;
	
	TTree *tree = new TTree("tree","Tree for Trigger Weight");
	tree->Branch("ptMin",&ptMin, "ptMin/F");
	tree->Branch("ptMax",&ptMax, "ptMax/F");
	tree->Branch("etaMin",&etaMin, "etaMin/F");
	tree->Branch("etaMax",&etaMax, "etaMax/F");
	tree->Branch("scale",&scale, "scale/F");
	tree->Branch("error",&error, "error/F");
	
	
    // Here one can declare histograms
    // In compiled C++, it possible not to use pointer type (TH1F*) and
    // use the non-pointer type (TH1F) of ROOT object.
	
	//Declare Histograms
	
	double etas[5] = {0.,0.8,1.44,1.57,2.5};
	double pts[] = {10.,15.,20.,25.,30.,35.,40.,45.,50.,200.};
	
	TH1F hallhJet_pt   ("hallhJet_pt","Pt of all jets in event",		150, 0.0, 150);
    TH1F hEleIDReco  ("hEleIDReco","Muon ID/Reco weight after EleFakeCuts ",		101, -0.5, 2.0);
	TH1F *hPtTurnOnCurve;
	TH1F *hEtaTurnOnCurve;
	TH1F *hPtPass = new TH1F	("hPtPass",  "hPt Pass Electron",					9, pts);
	TH1F *hEtaPass = new TH1F	("hEtaPass",  "hEta Pass Electron",	4, etas);
	TH1F *hPtAll = new TH1F	("hPtAll",  "hPt All Electron",					9, pts);
	TH1F *hEtaAll = new TH1F	("hEtaAll",  "hEta All Electron",	4, etas);
	hPtPass->Sumw2();
	hEtaPass->Sumw2();
	hPtAll->Sumw2();
	hEtaAll->Sumw2();
	
	TH1F *hPtTurnOnCurve0_08;
	TH1F *hPtTurnOnCurve08_144;
	TH1F *hPtTurnOnCurve157_25;
	TH1F *hPtPass0_08 = new TH1F	("hPtPass0_08",  "hPt Pass Electron 0 < eta < 0.8",					9, pts);
	TH1F *hPtPass08_144 = new TH1F	("hPtPass08_144",  "hPt Pass Electron 0.8 < ela < 1.44",					9, pts);
	TH1F *hPtPass157_25 = new TH1F	("hPtPass157_25",  "hPt Pass Electron 1.57 < eta < 2.5",					9, pts);
	TH1F *hPtAll0_08 = new TH1F	("hPtAll0_08",  "hPt All Electron  0 < eta < 0.8",					9, pts);
	TH1F *hPtAll08_144 = new TH1F	("hPtAll08_144",  "hPt All Electron 0.8 < ela < 1.44",					9, pts);
	TH1F *hPtAll157_25 = new TH1F	("hPtAll157_25",  "hPt All Electron 1.57 < eta < 2.5",					9, pts);
	
	
	if (debug) std::cout << "all histograms declared " << std::endl;
	
    // Loop over all events.
    // For other methods to access event/navigate through the sample,
    // see the documentation of RootTreeReader class.
	float  N_Vtype =0.0, Npreselect =0.0, N_EfakeCuts =0.0;
	int event =0;
	float FailedJetID=0.0;
	bool firstevent = true;
	
	//	double Emumass = 91.1976;
	
    do {
		//event = event + 1*PUweight2011;
		event++;
		isdata = false;
		//	if (isM50sample && sample.genZpt < 100) cout << "event cut out to remove overlap with DY_PtZ" << endl;
		if (((isM50sample && sample.genZpt < 100) || !isM50sample) && ((HFfile && sample.eventFlav == 5)||!HFfile) && ((LFfile && sample.eventFlav!=5)||!LFfile)){
			//cout << "Vtype is " << sample.Vtype << endl;
			if (sample.Vtype == 5 && sample.EVENT_json == 1 ){
				if (!(event%500))  std::cout << "entered event loop " << event << std::endl;
				B2011PUweight = 1.0, A2011PUweight = 1.0;
				B2011PUweight = sample.PUweight2011B;
				A2011PUweight = sample.PUweight;
				float PUweight2011 = (2.219*A2011PUweight + 2.238*B2011PUweight)/4.457;
				//cout << "pile up weight: " << PUweight2011 << endl;
				weight = LumiWeight*PUweight2011;
				Trigweight = 1.0;
				//if (!isDATA) Trigweight = sample.weightTrig;
				weight = Trigweight*weight;
				
				
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
					leptonPt[i] = -99.99;
					leptonEta[i] = -99.99;
					leptonPhi[i] = -99.99;
					lep_pfCombRelIso[i] = -99.99;
					lep_id95[i] = -99.99;
					lep_flavor[i] = -99;
					SortedMuonPt[i] = -99.99;
				}
				
				nJets =0, nSV =-99, nMuons = 0,  nElectrons = 0,  nLeptons = 0, Na_lep= 0, nPV= -99, MET= -99.99, Naj= 0, Nab =0, eventFlavor = -99;
				CSV0 = -1.0, CSV1 = -1.0, Emumass = -99.99, Hmass = -99.99, DeltaPhiHV = -99.99, Hpt = -99.99, Zpt = -99.99;
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
					}
					if ( fabs(sample.aJet_eta[a]) < 2.4 && (sample.aJet_pt[a] > 20)) Naj++;
					if (CSVshapeNew > 0.5 ) Nab++;
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
				jetPt[0] = sample.hJet_pt[indexedJetCSV[0].first];
				jetEta[0] = sample.hJet_eta[indexedJetCSV[0].first];
				jetPhi[0] = sample.hJet_phi[indexedJetCSV[0].first];
				jetCSV[0] = sample.hJet_csv[indexedJetCSV[0].first];
				jetCHF[0] = sample.hJet_chf[indexedJetCSV[0].first];	
				CSVNewShape[0] = indexedJetCSV[0].second;
				if (sample.nhJets > 1) { 
					CSV1 = indexedJetCSV[1].second;	
					jetPt[1] = sample.hJet_pt[indexedJetCSV[1].first];
					jetEta[1] = sample.hJet_eta[indexedJetCSV[1].first];
					jetPhi[1] = sample.hJet_phi[indexedJetCSV[1].first];
					jetCSV[1] = sample.hJet_csv[indexedJetCSV[1].first];
					jetCHF[1] = sample.hJet_chf[indexedJetCSV[1].first];	
					CSVNewShape[1] = indexedJetCSV[1].second;
					
					DetaJJ = sample.hJet_eta[indexedJetCSV[0].first]-sample.hJet_eta[indexedJetCSV[1].first];				
					Hmass = sample.H_mass;
					Hpt = sample.H_pt;
					ScalarSumHiggsJetPt = sample.hJet_pt[indexedJetCSV[1].first] + sample.hJet_pt[indexedJetCSV[0].first];
					ScalarSumJetPt = ScalarSumJetPt+ sample.hJet_pt[indexedJetCSV[1].first] + sample.hJet_pt[indexedJetCSV[0].first];
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
					RMS_eta = sample.vLepton_eta[1]*sample.vLepton_eta[1]+sample.vLepton_eta[0]*sample.vLepton_eta[0]+sample.hJet_eta[indexedJetCSV[0].first]*sample.hJet_eta[indexedJetCSV[0].first]+sample.hJet_eta[indexedJetCSV[1].first]*sample.hJet_eta[indexedJetCSV[1].first];
					RMS_eta = sqrt(RMS_eta/4);
					StandDevEta[0] =sample.vLepton_eta[1];
					StandDevEta[1] =sample.hJet_eta[indexedJetCSV[0].first];
					StandDevEta[2] =sample.vLepton_eta[0];
					StandDevEta[3] =sample.hJet_eta[indexedJetCSV[1].first];
					EtaStandDev = TMath::RMS(4,StandDevEta);
					AverageEta = sample.vLepton_eta[0]+sample.vLepton_eta[1]+sample.hJet_eta[indexedJetCSV[0].first];
					AverageEta = (AverageEta+sample.hJet_eta[indexedJetCSV[1].first])/4;
					UnweightedEta = (sample.vLepton_eta[0]-AverageEta)*(sample.vLepton_eta[0]-AverageEta);
					UnweightedEta = UnweightedEta + (sample.vLepton_eta[1]-AverageEta)*(sample.vLepton_eta[1]-AverageEta);
					UnweightedEta = UnweightedEta + (sample.hJet_eta[indexedJetCSV[0].first]-AverageEta)*(sample.hJet_eta[indexedJetCSV[0].first]-AverageEta);
					UnweightedEta = UnweightedEta + (sample.hJet_eta[indexedJetCSV[1].first]-AverageEta)*(sample.hJet_eta[indexedJetCSV[1].first]-AverageEta);
					
					PtbalZH = (Hpt-Zpt);
					PtbalMETH = Hpt-MET;
					PtbalZMET = Zpt - MET;
					lep0pt = sample.vLepton_pt[0];
					
					FirstJet.SetPtEtaPhiM(sample.hJet_pt[indexedJetCSV[0].first],sample.hJet_eta[indexedJetCSV[0].first],sample.hJet_phi[indexedJetCSV[0].first],4.2/1000.0);
					SecondJet.SetPtEtaPhiM(sample.hJet_pt[indexedJetCSV[1].first],sample.hJet_eta[indexedJetCSV[1].first],sample.hJet_phi[indexedJetCSV[1].first],4.2/1000.0);
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
				
				
				
				//	if (sample.triggerFlags[39]||sample.triggerFlags[40]||sample.triggerFlags[41]){
				FOM_tree->Fill();
				if ( jetPt[0] > 20 && jetPt[1] > 20 && 
					fabs(sample.vLepton_eta[0]) < 2.5 && fabs(sample.vLepton_eta[1]) < 2.4 && fabs(sample.hJet_eta[0]) < 2.5 &&
					fabs(sample.hJet_eta[1]) < 2.5 && sample.hJet_id[0]==1 && sample.hJet_id[1]==1 && sample.hbhe){
					Npreselect = Npreselect + 1*PUweight2011;
					JetDistributions("PreSelect", weight);
					LeptonDistributions("PreSelect", weight);
					TH2FDistributions("PreSelect", weight);
					EventShapeDistributions("PreSelect", weight);
					EventDistributions("PreSelect", weight);
					if (isDATA) isdata = true;
					if(delRemu>0.4){
						if(sample.vLepton_pt[1]>20){
							//	"HLT_IsoMu17_v.*" , #0
							//							"HLT_IsoMu17_eta2p1_DiCentralJet30_v.*", #24
							if (!(sample.triggerFlags[0]||sample.triggerFlags[24])) cout << "NO Muon TRIGGER"<< endl;
							N_EfakeCuts = N_EfakeCuts + 1*PUweight2011;
							hPtAll->Fill(sample.vLepton_pt[0]);
							hEtaAll->Fill(fabs(sample.vLepton_eta[0]));
							if (sample.triggerFlags[39]){
								hPtPass->Fill(sample.vLepton_pt[0]);
								hEtaPass->Fill(fabs(sample.vLepton_eta[0]));
							}
							for(int i = 0; i<4; i++){
								if(sample.vLepton_eta[0] > etas[i] && sample.vLepton_eta[0] < etas[i+1]){
									switch (i){
										case 0:
											hPtAll0_08->Fill(sample.vLepton_pt[0]);
											if (sample.triggerFlags[39]) hPtPass0_08->Fill(sample.vLepton_pt[0]);
											break;
										case 1:
											hPtAll08_144->Fill(sample.vLepton_pt[0]);
											if (sample.triggerFlags[39]) hPtPass08_144->Fill(sample.vLepton_pt[0]);
											break;
										case 2:
											break;
										case 3:
											hPtAll157_25->Fill(sample.vLepton_pt[0]);
											if (sample.triggerFlags[39]) hPtPass157_25->Fill(sample.vLepton_pt[0]);
											break;
										default:
											cout << "ERROR ERROR ERROR BINS ARE WRONG" << endl;
									}//end switch
								}//end if statement
							}//end for loop
						}
						JetDistributions("delRemu", weight);
						LeptonDistributions("delRemu", weight);
						TH2FDistributions("delRemu", weight);
						EventShapeDistributions("delRemu", weight);
						EventDistributions("delRemu", weight);
					}//EleFakeCuts
				} else {
					FailedJetID = FailedJetID + 1*PUweight2011;
				}//Jet ID and eta requirement
				isdata = false;
				EventDistributions("allEvts",weight);
			}//end requirement Zemu event
		}//end isM50sample genZpt cut
	} while (sample.nextEvent());
	
	
	
	std::cout << endl << endl;
	std::cout << "Number of events " << event << endl;
	std::cout << "Vtype 5 " << N_Vtype << endl;	
	std::cout << "PreSelection: " << Npreselect << endl;
	std::cout << "EleFakeCuts: " << N_EfakeCuts << endl;
	
	std::cout << endl << endl;
	std::cout << "Number of Events that failed JetID " << FailedJetID << endl;
	
	ofstream myfile;
	myfile.open(TString::Format("%sElectronTrigEff.txt",directory.c_str()).Data());	
	myfile<<TString::Format("\t %s",directory.c_str()).Data()<<endl;
	//	myfile<<TString::Format("Number of events \t %f \t %0.5f",event,event*LumiWeight).Data()<<endl;
	//	myfile<<TString::Format("Vtype 5 \t %f \t %0.5f",N_Vtype,N_Vtype*LumiWeight).Data()<<endl;
	myfile<<TString::Format("pass pt bins \t %0.0f \t %0.0f  \t %0.0f \t %0.0f  \t %0.0f \t %0.0f  \t %0.0f ",hPtPass->GetBinContent(1),hPtPass->GetBinContent(2),hPtPass->GetBinContent(3)
							,hPtPass->GetBinContent(4),hPtPass->GetBinContent(5),hPtPass->GetBinContent(6),hPtPass->GetBinContent(7)).Data()<<endl;
	myfile<<TString::Format("all pt bins \t %0.0f \t %0.0f  \t %0.0f \t %0.0f  \t %0.0f \t %0.0f  \t %0.0f ",hPtAll->GetBinContent(1),hPtAll->GetBinContent(2),hPtAll->GetBinContent(3)
							,hPtAll->GetBinContent(4),hPtAll->GetBinContent(5),hPtAll->GetBinContent(6),hPtAll->GetBinContent(7)).Data()<<endl;
	myfile<<TString::Format("eff pt bin \t %0.5f \t %0.5f  \t %0.5f \t %0.5f  \t %0.5f \t %0.5f  \t %0.5f ",hPtPass->GetBinContent(1)/hPtAll->GetBinContent(1),
							hPtPass->GetBinContent(2)/hPtAll->GetBinContent(2),hPtPass->GetBinContent(3)/hPtAll->GetBinContent(3),hPtPass->GetBinContent(4)/hPtAll->GetBinContent(4),
							hPtPass->GetBinContent(5)/hPtAll->GetBinContent(5),hPtPass->GetBinContent(6)/hPtAll->GetBinContent(6),hPtPass->GetBinContent(7)/hPtAll->GetBinContent(7)).Data()<<endl;
	myfile<<TString::Format("pass eta bins \t %0.0f \t %0.0f  \t %0.0f \t %0.0f ",hEtaPass->GetBinContent(1),hEtaPass->GetBinContent(2),hEtaPass->GetBinContent(3),hEtaPass->GetBinContent(4)).Data()<<endl;
	myfile<<TString::Format("all eta bins \t %0.0f \t %0.0f  \t %0.0f  \t %0.0f",hEtaAll->GetBinContent(1),hEtaAll->GetBinContent(2),hEtaAll->GetBinContent(3),hEtaAll->GetBinContent(4)).Data()<<endl;
	myfile<<TString::Format("eff eta bins \t %0.5f \t %0.5f \t %0.5f \t %0.5f",hEtaPass->GetBinContent(1)/hEtaAll->GetBinContent(1),
							hEtaPass->GetBinContent(2)/hEtaAll->GetBinContent(2),hEtaPass->GetBinContent(3)/hEtaAll->GetBinContent(3),hEtaPass->GetBinContent(4)/hEtaAll->GetBinContent(4)).Data()<<endl;
	
	hPtTurnOnCurve = new TH1F	("hPtTurnOnCurve",  "hPtTurnOnCurve Electron",		9, pts);
	hEtaTurnOnCurve = new TH1F	("hEtaTurnOnCurve",  "hEtaTurnOnCurve Electron",	4, etas);
	
	hPtTurnOnCurve->Sumw2();
	hEtaTurnOnCurve->Sumw2();
	hPtTurnOnCurve->Divide(hPtPass,hPtAll,1,1,"B");
	hEtaTurnOnCurve->Divide(hEtaPass,hEtaAll,1,1,"B");
	
	
	myfile<<TString::Format("error pt bins \t %0.5f \t %0.5f  \t %0.5f \t %0.5f  \t %0.5f \t %0.5f  \t %0.5f ",hPtTurnOnCurve->GetBinError(1),hPtTurnOnCurve->GetBinError(2),hPtTurnOnCurve->GetBinError(3)
							,hPtTurnOnCurve->GetBinError(4),hPtTurnOnCurve->GetBinError(5),hPtTurnOnCurve->GetBinError(6),hPtTurnOnCurve->GetBinError(7)).Data()<<endl;
	myfile<<TString::Format("error eta bins \t %0.5f \t %0.5f  \t %0.5f \t %0.5f ",hEtaTurnOnCurve->GetBinError(1),hEtaTurnOnCurve->GetBinError(2),hEtaTurnOnCurve->GetBinError(3),hEtaTurnOnCurve->GetBinError(4)).Data()<<endl;
	
	
	myfile<<endl;
	myfile.close();
	
	hPtTurnOnCurve0_08 = new TH1F	("hPtTurnOnCurve0_08",  "hPtTurnOnCurve Electron 0 < eta < 0.8",		9, pts);
	hPtTurnOnCurve0_08->Sumw2();
	hPtTurnOnCurve0_08->Divide(hPtPass0_08,hPtAll0_08,1,1,"B");
	hPtTurnOnCurve08_144 = new TH1F	("hPtTurnOnCurve08_144",  "hPtTurnOnCurve Electron 0.8 < ela < 1.44",		9, pts);
	hPtTurnOnCurve08_144->Sumw2();
	hPtTurnOnCurve08_144->Divide(hPtPass08_144,hPtAll08_144,1,1,"B");
	hPtTurnOnCurve157_25 = new TH1F	("hPtTurnOnCurve157_25",  "hPtTurnOnCurve Electron 1.57 < eta < 2.5",		9, pts);
	hPtTurnOnCurve157_25->Sumw2();
	hPtTurnOnCurve157_25->Divide(hPtPass157_25,hPtAll157_25,1,1,"B");
	
			for(int i = 0; i<4; i++){
					for(int j = 0; j<9; j++){
							ptMin = pts[j];
							ptMax = pts[j+1];
							etaMin = etas[i];
							etaMax  = etas[i+1];
							switch (i){
								case 0:
									scale = hPtTurnOnCurve0_08->GetBinContent(j);
									error = hPtTurnOnCurve0_08->GetBinError(j);
									break;
								case 1:
									scale = hPtTurnOnCurve08_144->GetBinContent(j);
									error = hPtTurnOnCurve08_144->GetBinError(j);
									break;
								case 2:
									scale = -99.99;
									break;
								case 3:
									scale = hPtTurnOnCurve157_25->GetBinContent(j);
									error = hPtTurnOnCurve157_25->GetBinError(j);
									break;
								default:
									cout << "ERROR ERROR ERROR BINS ARE WRONG" << endl;
							}//end switch
							if ((scale +99.99)>1) tree->Fill();
					}// end pt for loop
			}//end eta for loop
	
	// Here one can create canvas, draw and print the histogram.
	TCanvas c1("c1","c1");
	c1.cd();
	c1.SetFillColor(kWhite);
	string suffixps = ".gif";
	
	c1.Clear(); // don't create a new canvas
	hallhJet_pt.Draw();
	c1.Print((directory+"/PtAllJets"+suffixps).c_str());
	
	c1.Clear(); // don't create a new canvas
	hPtTurnOnCurve->Draw();
	c1.Print((directory+"/PtTurnOnCurve"+suffixps).c_str());
	
	c1.Clear(); // don't create a new canvas
	hEtaTurnOnCurve->Draw();
	c1.Print((directory+"/EtaTurnOnCurve"+suffixps).c_str());
	
	c1.Clear(); // don't create a new canvas
	hPtPass->Draw();
	c1.Print((directory+"/PtPass"+suffixps).c_str());
	
	c1.Clear(); // don't create a new canvas
	hEtaPass->Draw();
	c1.Print((directory+"/EtaPass"+suffixps).c_str());
	
	c1.Clear(); // don't create a new canvas
	hPtAll->Draw();
	c1.Print((directory+"/PtAll"+suffixps).c_str());
	
	FOM_tree->Write();
	tree->Write();
	
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
	double lumi = 4.457;
	if (findString(filename, "ZH115_1M")){ 
		SampleWeight = lumi/lumiZH115;
	cout << "found ZH_ZToLL_HToBB_M-115 string" << endl;}
	if (findString(filename, "ZH115_Fall11")){ SampleWeight = lumi/(2.28079234375000000e+05/xsecbfZH115);}
	if (findString(filename, "Filteredemu")){ SampleWeight = lumi/lumiZH115;}
	if (findString(filename, "DY_PtZ")){ SampleWeight = lumi/lumiZJH;}
	if (findString(filename, "DY_M50")){ SampleWeight = lumi/lumiZJL;}
	if (findString(filename, "120to170")){ SampleWeight = lumi/lumiQCD120;}
	if (findString(filename, "170to300")){ SampleWeight = lumi/lumiQCD170;}
	if (findString(filename, "300to470")){ SampleWeight = lumi/lumiQCD300;}
	if (findString(filename, "470to600")){ SampleWeight = 70.21630282*100000/ 3598283.0;}
	if (findString(filename, "80to120")){ SampleWeight = lumi/(6397439.5000/(823744.5*1000));}
	if (findString(filename, "TTJets")){ SampleWeight = lumi/(lumiTT);}
	if (findString(filename, "T_TuneZ2_s")){ SampleWeight = lumi/lumiTs;}
	if (findString(filename, "T_TuneZ2_t-channel")){ SampleWeight = lumi/lumiTt;}
	if (findString(filename, "T_TuneZ2_tW-channel-DR")){ SampleWeight = lumi/lumiTtW;}
	if (findString(filename, "T_TuneZ2_tW-channel-DS")){ SampleWeight = lumi/(828002.00/xsecbfTtW) ;}
	if (findString(filename, "Tbar_TuneZ2_s")){ SampleWeight = lumi/lumiTsb;}
	if (findString(filename, "Tbar_TuneZ2_t-channel")){ SampleWeight = lumi/lumiTtb;}
	if (findString(filename, "Tbar_TuneZ2_tW-channel-DR")){ SampleWeight = lumi/lumiTtWb;}
	if (findString(filename, "Tbar_TuneZ2_tW-channel-DS")){ SampleWeight = lumi/(727800.6250000/xsecbfTtWb);}
	if (findString(filename, "WJetsToLNu_Pt-100")){ SampleWeight = lumi/(18859930.0/(685.69*1000) );}
	if (findString(filename, "WJetsToLNu_TuneZ2")){ SampleWeight = lumi/lumiWJ;}
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
	
	for (int i=0; i != nJets && i < 10; ++i) {
		string jetpT = Form("hPtb%i_", i) + cut;
		string jetNewCSV = Form("hCSVNewShape%i_", i) + cut;
		string jetNewCSVAC = Form("hCSVNewShape%iAC_", i) + cut;
		string jeteta = Form("hEtab%i_", i) + cut;
		string jetphi = Form("hPhib%i_", i) + cut;
		string jetcsv = Form("hCSV%i_", i) + cut;
		string jetchf = Form("hCHFb%i_", i) + cut;
		
		fillhisto(jetpT, jetPt[i], ph_weight, "jet pT", 10, 0.0, 200);
		fillhisto(jetNewCSV, CSVNewShape[i], ph_weight, "CSV BTag Shape", 30, 0, 1.5);
		fillhisto(jetNewCSVAC, CSVNewShape[i], ph_weight, "CSV BTag Shape", 40, 0.244, 1.244);
		fillhisto(jeteta, jetEta[i], ph_weight, "jet #eta", 13, -3, 3.5);
		fillhisto(jetphi, jetPhi[i], ph_weight, "jet #phi", 20, -3.14159265, 4.7123889);
		fillhisto(jetcsv, jetCSV[i], ph_weight, "jet CSV", 30, 0, 1.5);
		if(jetCHF[i]> 0.0000001)fillhisto(jetchf, jetCHF[i], ph_weight, "charged Hadron Energy Fraction", 30, 0.0, 1.2);
	}
	
	string mjj = "hMjj_" + cut;
	string mjjAC = "hMjjAC_" + cut;
	string ptjj = "hPtjj_" +cut;
	string ptjjAC = "hPtjjAC_" +cut;
	string detajj = "hdetaJJ_" +cut;
	string scalarSumHiggsJetPt = "hScalarSumHiggsJetPt_" +cut;
	string scalarSumJetPt = "hScalarSumJetPt_" +cut;
	
	if (!isdata || (Hmass < 80 || Hmass > 150)) fillhisto(mjj, Hmass, ph_weight, "Invariant Mass of two Jets", 25, 0, 250);
	if (!isdata || (Hmass < 80 || Hmass > 150)) fillhisto(mjjAC, Hmass, ph_weight, "Invariant Mass of two Jets", 20, 70, 170);
	fillhisto(ptjj, Hpt, ph_weight, "Pt of two b jets with highest CSV", 14, 0, 210);
	fillhisto(ptjjAC, Hpt, ph_weight, "Pt of two b jets with highest CSV", 10, 15, 115);
	fillhisto(detajj, fabs(DetaJJ), ph_weight, "Delta eta between two jets", 10, 0, 5);
	fillhisto(scalarSumHiggsJetPt, ScalarSumHiggsJetPt, ph_weight, "scalar sum higgs jet pt", 25, 30, 300);
	fillhisto(scalarSumJetPt, ScalarSumJetPt, ph_weight, "scalar sum alljet pt", 25, 30, 300);
	
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
	
	string hDphiDetajj = "hDphiDetajj_" + cut;
	string hqtvsalphaJJ = "hqtvsalphaJJ_" +cut;
	string hqtvsalphaZ = "hqtvsalphaZ_" +cut;
	string hDphiHZDphiHMET = "hDphiHZDphiHMET_" + cut;
	string hDphiHZDphiZMET = "hDphiHZDphiZMET_" + cut;
	string hDphiZMETDphiHMET = "hDphiZMETDphiHMET_" + cut;
	string hZmassvsDphiemu = "hZmassvsDphiemu_" + cut;
	string hMtevsMtmu = "hMtevsMtmu_" + cut;
	string hMETvsPtZ = "hMETvsPtZ_" + cut;
	string hProjBis = "hProjBis_" + cut;
	string hZptvsHpt = "hZptvsHpt_" + cut;
	string hABCD = "hABCD_" + cut;
	
	fill2Dhisto(hDphiDetajj, DphiJJ, DetaJJ, ph_weight, "#Delta#phi vs #Delta#eta JJ", 75, -3.5, 6, 75, -5, 5);
	fill2Dhisto(hqtvsalphaJJ, alpha_j, qtb1, ph_weight, "Armenteros-Podolansky Plot Jets", 70, -2, 2, 75, 0, 75);
	fill2Dhisto(hqtvsalphaZ, alpha_lep, qtlep1, ph_weight,  "Armenteros-Podolansky Plot Z", 70, -2, 2, 75, 0, 75);
	fill2Dhisto(hDphiHZDphiHMET, DeltaPhiHV, dPhiHMET, ph_weight, "#Delta#phi HZ vs #Delta#phi HMET", 75, 0, 4.5, 65, 0, 3.25);
	fill2Dhisto(hDphiHZDphiZMET, DeltaPhiHV, fabs(DphiZMET), ph_weight, "#Delta#phi HZ vs #Delta#phi ZMET", 75, 0, 4.5, 65, 0, 3.25);
	fill2Dhisto(hDphiZMETDphiHMET, fabs(DphiZMET), dPhiHMET, ph_weight, "#Delta#phi ZMET vs #Delta#phi HMET", 75, 0, 4.5, 65, 0, 3.25);
	if (Zmass > 0 )fill2Dhisto(hZmassvsDphiemu, Dphiemu, Zmass, ph_weight, "#Delta#phi emu vs Zmass", 75, -3.25, 5.5, 75, 0, 700);
	fill2Dhisto(hMtevsMtmu, Mte, Mtmu, ph_weight, "Mtmu vs Mte", 101, -0.1, 200, 101, -0.1, 200);
	fill2Dhisto(hMETvsPtZ, MET, Zpt, ph_weight, "MET vs Pt(e,mu)", 50, 0.0, 200, 150, 0, 300);
	fill2Dhisto(hProjBis, ProjVisT, ProjMissT, ph_weight, "ProjVisT vs ProjMissT", 100, 0, 300, 100, -200, 200);
	fill2Dhisto(hZptvsHpt, Zpt, Hpt, ph_weight, "Zpt vs Hpt", 150, 0, 300, 150, 0, 300);
	fill2Dhisto(hABCD, CSVNewShape[1]+CSVNewShape[0], DeltaPhijetMETmin, ph_weight, "ABCD method for controlling QCD", 50, 0, 2, 50, 0, 3);
	
	
	
	
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
	
	
	fillhisto(hPtbalZH, PtbalZH, ph_weight, "Pt balance of Z and H", 83, -99.5, 149.5);
	fillhisto(hPtbalMETH, PtbalMETH, ph_weight, "Pt balance of MET and H", 83, -99.5, 149.5);
	fillhisto(hPtbalZMET, PtbalZMET, ph_weight, "Pt balance of Z and MET", 83, -99.5, 149.5);
	fillhisto(hdphiVH, DeltaPhiHV, ph_weight, "Delta phi between Z and Higgs", 24, 0, 4.71238898);
	fillhisto(hdphiVHAC, DeltaPhiHV, ph_weight, "Delta phi between Z and Higgs", 44, 0.785398, 4.71238898);
	fillhisto(hRMSeta, RMS_eta, ph_weight, "RMS Eta", 22, 0, 2.2);
	fillhisto(hStaDeveta, EtaStandDev, ph_weight, "Standard Deviation Eta",		50, 0, 2);
	fillhisto(hUnweightedEta, UnweightedEta, ph_weight, "Unweighted Eta",		60, 0, 7);
	fillhisto(hScalarSumPt, ScalarSumPt, ph_weight, "scalar sum of pt of four particles", 25, 30, 300);
	fillhisto(hCentrality, Centrality, ph_weight, "Centrality", 30, 0.0, 0.70);
	fillhisto(hEventPt, EventPt, ph_weight, "Pt of HV system", 50, 0.0, 250);
	fillhisto(hEventMass, EventMass, ph_weight, "Mass of HV system", 100, 100, 700);
	fillhisto(hAngleHemu, AngleHemu, ph_weight, "Angle between H and Z", 15, 0, 3.5);
	fillhisto(hSphericity, EvntShpSphericity, ph_weight,"EventShapeVariables sphericity", 45, 0.0, .9);
	fillhisto(hAplanarity, EvntShpAplanarity, ph_weight, "EventShapeVariables Aplanarity", 30, 0.0, .3);
	fillhisto(hCircularity, EvntShpCircularity, ph_weight,  "EventShapeVariables circularity", 20, 0.0, 1.2);
	fillhisto(hIsotropy, EvntShpIsotropy, ph_weight,  "EventShapeVariables isotropy", 50, 0.1, 1.4);
	fillhisto(dphijj, fabs(DphiJJ), ph_weight, "Delta phi between two jets",  24, 0, 4.71238898);
	if (Ht>1)fillhisto(hHt, Ht, ph_weight, "from MHT class MHT_ht", 50, 0, 250);
	if(Zmass>-98) fillhisto(hZmass, Zmass, ph_weight, "Invariant Mass of two Leptons corrected", 20, 0, 100);
	if (ZmassSVD>-98) fillhisto(hZmassSVD, ZmassSVD, ph_weight, "Invariant Mass of two Leptons corrected SVD", 20, 0, 100);
	if (ZmassSVD>-98) fillhisto(hZmassSVDAC, ZmassSVD, ph_weight, "Invariant Mass of two Leptons corrected SVD", 20, 0, 100);
	if (ZmassSVDnegSol>-98) fillhisto(hZmassSVDnegSol, ZmassSVDnegSol, ph_weight, "Invariant Mass of two Leptons corrected SVD", 20, 0, 100);
	if(ZmassNegInclu>-98)fillhisto(hZmassNegInclu, ZmassNegInclu, ph_weight, "Invariant Mass of two Leptons corrected Matrix", 20, 0, 100);
	fillhisto(hAngleEMU, AngleEMU, ph_weight, "Angle between electron and muon", 15, 0, 3.5);
	fillhisto(hCosThetaMu,CosThetaMu, ph_weight, "Cos Theta Muon", 30, -1, 2.0);
	fillhisto(hCosThetaEle,CosThetaEle, ph_weight, "Cos Theta Electron", 30, -1, 2.0);
	fillhisto(hEleMissE, AaronEleMissE, ph_weight, "Missing Energy electron", 50, -400, 400);
	fillhisto(hMuonMissE, AaronMuMissE, ph_weight, "Missing Energy muon", 50, -400, 400);
	fillhisto(hDphiemu, fabs(Dphiemu), ph_weight, "Delta phi between e and muon", 24, 0, 4.71238898);
	fillhisto(hDetaemu, fabs(Detaemu), ph_weight, "Delta eta between e and muon", 10, 0, 5);
	if(MassEleb0>1)fillhisto(hMassEleb0, MassEleb0, ph_weight, "Invariant Mass of Electron and b0", 50, 0, 250);
	if(MassMub0>1)fillhisto(hMassMub0, MassMub0, ph_weight, "Invariant Mass of Muon and b0", 50, 0, 250);
	if(MassEleb1>1)fillhisto(hMassEleb1, MassEleb1, ph_weight, "Invariant Mass of Electron and b1", 50, 0, 250);
	if(MassMub1>1)fillhisto(hMassMub1, MassMub1, ph_weight, "Invariant Mass of Muon and b1", 50, 0, 250);
	fillhisto(hDphiEleMET, fabs(DphiEleMET), ph_weight, "Delta phi between Electron and MET",  24, 0, 4.71238898);
	fillhisto(hdphiMuMET,fabs( dphiMuMET), ph_weight, "Delta phi between Muon and MET",  24, 0, 4.71238898);
	fillhisto(hDphiLeadMET, fabs(DphiLeadMET), ph_weight, "Delta phi between leading lepton and MET",  24, 0, 4.71238898);
	fillhisto(hDphiSecondMET, fabs(DphiSecondMET), ph_weight, "Delta phi between second lepton and MET",  24, 0, 4.71238898);
	fillhisto(hDphiZMET, fabs(DphiZMET), ph_weight, "Delta phi between Z and MET",  24, 0, 4.71238898);
	fillhisto(hDphiZMETAC, fabs(DphiZMET), ph_weight, "Delta phi between Z and MET",  24, 0, 3.14159265);
	fillhisto(hdelRjj, delRjj, ph_weight, "Delta R jj", 20, 0, 5);
	fillhisto(hdelRemu, delRemu, ph_weight, "Delta R emu", 20, 0, 5);
	fillhisto(hProjVisT, ProjVisT, ph_weight, "Transverse componenet of Projection of Z onto bisector", 100, -200, 200);
	fillhisto(hProjMissT, ProjMissT, ph_weight, "Transeverse compenent of Projection of MET onto e,mu bisector", 100, -200, 200);
	fillhisto(hPzetaCut, ProjMissT-(0.25*ProjVisT), ph_weight, "P_#zeta Cut", 75, -200, 100);
	
	
	
	
	
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
	string hminDeltaPhijetMET = "hminDeltaPhijetMET_" +cut;
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
	fillhisto(hminDeltaPhijetMET,DeltaPhijetMETmin, ph_weight, "Delta phi between MET and nearest jet", 16, 0, 3.14159265);
	fillhisto(hminDeltaPhijetMETZtau,DeltaPhijetMETZtaumin, ph_weight, "Delta phi between MET+Zpt and nearest jet", 24, 0, 4.71238898);
	fillhisto(hSVmass, SV_mass, ph_weight, "mass of Secondary Vertex",		20, 0.0, 5);
	fillhisto(hVMte, Mte, ph_weight, "VMte",  16, 0, 160);
	fillhisto(hVMtmu, Mtmu, ph_weight, "VMtmu",  16, 0, 160);
	fillhisto(hdeltaPullAngle,delPullAngle, ph_weight, "Delta pull Angle", 42, -3.25, 5.25);
	fillhisto(hdeltaPullAngle2,delPullAngle2, ph_weight, "Delta Pull Angle 2", 42, -3.25, 5.25);
	fillhisto(htopMass,topMass, ph_weight, "Top Mass", 50, 75, 375);
	fillhisto(htopPt,topPt, ph_weight, "Pt of Top", 50, 0, 200);
	if (topWmass > 82 || topWmass < 80) fillhisto(htopWmass,topWmass, ph_weight, "Mass of W coming from Top", 25, 75, 125);
	
	
	
}// EventDistributions	
