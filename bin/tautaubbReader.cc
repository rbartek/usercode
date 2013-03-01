/**
 @brief Example of an analysis code to read tree ntuple
 and create histogram
 Usage:
 
 \verbatim
 mkdir [directoryName]
 Hbb_TauemuTuples <inputfile> [outputfilename] [directoryName]
 \endverbatim
 
 @param inputfile Either a ROOT file or an ASCII file containing list of
 ROOT files.
 
 @param outputfile Name of the ROOT file which contains the histogram.
 Defaulted to 'output.root'
 
 @author Rachel Wilken <rachel.wilken@cern.ch>
 
 @date Fri Dec 21 2012
 
 */

#include "tautaubbReader.h"
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

float max(float csv1, float  csv2)
{
	if(csv1 < csv2) return csv2;
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

float corrCSV(BTagShapeNew* btag, float csv,int flav){
	if(csv < 0.) return csv;
	if(csv > 1.) return csv;
	if(flav == 0) return csv;
	if(fabs(flav) == 5) return  btag->ib->Eval(csv);
	if(fabs(flav) == 4) return  btag->ic->Eval(csv);
	if(fabs(flav) != 4  and fabs(flav) != 5) return  btag->il->Eval(csv);
	return -10000;
}

float JERSys(bool isUp, float eta, float pt, float genpt){
	float inner, outer;
	float rPt = -99.99;
	if (isUp) inner = 0.06, outer = 0.1;
	if (!isUp) inner = -0.06, outer = -0.1;
	if (fabs(eta)<1.1){
		rPt = pt + (pt-genpt)*inner;
	} else {
		rPt = pt + (pt-genpt)*outer;
	}
	return rPt;
} 

float AssignWeight(TTree* mytree, float pt1, float ABSeta){
	float ptMin,ptMax,etaMin,etaMax,scale,error;
	float s1 = 1.0;
	if(!mytree) return s1;
	int count = 0;
	mytree->SetBranchAddress("ptMin",&ptMin);
	mytree->SetBranchAddress("ptMax",&ptMax);
	mytree->SetBranchAddress("etaMin",&etaMin);
	mytree->SetBranchAddress("etaMax",&etaMax);
	mytree->SetBranchAddress("scale",&scale);
	mytree->SetBranchAddress("error",&error);
    float lastPtBin = 200;
    for(int jentry = 0; jentry < mytree->GetEntries(); jentry++)
		{
			mytree->GetEntry(jentry);
			if(ptMax==lastPtBin) ptMax=1e99;
			if((pt1 > ptMin) && (pt1 < ptMax) && (ABSeta > etaMin) && (ABSeta < etaMax))
				{
					s1 = scale;
					count++;
				}
		}
	return s1;
}

void HiggsCandBuilder(int *X, int *Y, double CSV0, double CSV1, int posCSV0, int posCSV1, int pospt0, int pospt1){
	//choose Higgs Jets with new Algo
	int CSV44mid1 = 99, CSV44mid2 = 99;
	if(CSV0>0.4&&CSV1>0.4){
		CSV44mid1 = posCSV0;
		CSV44mid2 = posCSV1;							
	}else if (CSV0>0.4&&CSV1<0.4){
		CSV44mid1 = posCSV0;
		if(posCSV0!=pospt0)CSV44mid2 = pospt0;
		if(posCSV0==pospt0)CSV44mid2 = pospt1;	
	}else if(CSV0<0.4&&CSV1>0.4){
		CSV44mid1 = posCSV0;
		CSV44mid2 = posCSV1;							
	}else{
		CSV44mid1 = pospt0;
		CSV44mid2 = pospt1;
	}
	(*X) = CSV44mid1;
	(*Y) = CSV44mid2;
}	

int main(int argc, char** argv) {
	
	bool debug = true;	
	
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
	
	TFile *MuonTrigWeightFilelow = TFile::Open("/home/hep/wilken/CMSSW_5_3_3_patch2/src/VHbbAnalysis/VHbbDataFormats/bin/triggerRootFiles/DoubleMu8.TrigEff.2012A.root");
	TFile *EleTrigWeightFilelow = TFile::Open("/home/hep/wilken/CMSSW_5_3_3_patch2/src/VHbbAnalysis/VHbbDataFormats/bin/triggerRootFiles/DoubleEle8.TrigEff.wp95.2012A.root");
	TFile *MuonTrigWeightFilehigh = TFile::Open("/home/hep/wilken/CMSSW_5_3_3_patch2/src/VHbbAnalysis/VHbbDataFormats/bin/triggerRootFiles/DoubleMu17.TrigEff.2012A.root");
	TFile *EleTrigWeightFilehigh = TFile::Open("/home/hep/wilken/CMSSW_5_3_3_patch2/src/VHbbAnalysis/VHbbDataFormats/bin/triggerRootFiles/DoubleEle17.TrigEff.wp95.2012A.root");
	TFile *MuonIDWeightFile = TFile::Open("/home/hep/wilken/CMSSW_5_3_3_patch2/src/VHbbAnalysis/VHbbDataFormats/bin/triggerRootFiles/MuRecoId.ScaleFactor.2012A.root");
	TFile *EleIDWeightFile = TFile::Open("/home/hep/wilken/CMSSW_5_3_3_patch2/src/VHbbAnalysis/VHbbDataFormats/bin/triggerRootFiles/EleRecoId.ScaleFactor.wp95.2012A.root");
	TFile *lowptMuonIDWeightFile = TFile::Open("/home/hep/wilken/CandReaderForTaus/src/UserCode/wilken/ScaleFactorsMuon_lowptBoth.root");
	TFile *lowptEleIDWeightFile = TFile::Open("/home/hep/wilken/CandReaderForTaus/src/UserCode/wilken/ScaleFactorsElectronWP95_lowpt.root");
	
	TTree *MuonTrigWeightTreelow = (TTree*) MuonTrigWeightFilelow->Get("tree"); 
	TTree *EleTrigWeightTreelow = (TTree*) EleTrigWeightFilelow->Get("tree"); 
	TTree *MuonTrigWeightTreehigh = (TTree*) MuonTrigWeightFilehigh->Get("tree"); 
	TTree *EleTrigWeightTreehigh = (TTree*) EleTrigWeightFilehigh->Get("tree"); 
	TTree *MuonIDWeightTree = (TTree*) MuonIDWeightFile->Get("tree"); 
	TTree *EleIDWeightTree = (TTree*) EleIDWeightFile->Get("tree"); 
	TTree *lowptMuonIDWeightTree = (TTree*) lowptMuonIDWeightFile->Get("tree"); 
	TTree *lowptEleIDWeightTree = (TTree*) lowptEleIDWeightFile->Get("tree"); 
	
	
	BTagShapeNew* btagNew;
	btagNew = new BTagShapeNew("csvdiscrNew.root");
	btagNew->computeFunctions();
	if (debug){
		cout << "Done computing btagSF " << endl;
		cout << "b quark SF for CVSL: " << btagNew->ib->Eval(0.244) << "   CVSM: " << btagNew->ib->Eval(0.5) << "  CSVT: " << btagNew->ib->Eval(0.898) << endl;
		cout << "c quark SF for CVSL: " << btagNew->ic->Eval(0.244) << "   CVSM: " << btagNew->ic->Eval(0.5) << "  CSVT: " << btagNew->ic->Eval(0.898) << endl;
		cout << "light quark SF for CVSL: " << btagNew->il->Eval(0.244) << "   CVSM: " << btagNew->il->Eval(0.5) << "  CSVT: " << btagNew->il->Eval(0.898) << endl << endl;
	}	
	BTagShapeNew* btagUp;
	btagUp = new BTagShapeNew("csvdiscrNew.root");
	btagUp->computeFunctions(+1.,0.);
	BTagShapeNew* btagDown;
	btagDown = new BTagShapeNew("csvdiscrNew.root");
	btagDown->computeFunctions(-1.,0.);
	BTagShapeNew* btagFUp;
	btagFUp = new BTagShapeNew("csvdiscrNew.root");
	btagFUp->computeFunctions(0.,+1.);
	BTagShapeNew* btagFDown;
	btagFDown = new BTagShapeNew("csvdiscrNew.root");
	btagFDown->computeFunctions(0.,-1.);
	
	
    // Create the output ROOT file
    TFile ofile(ofilename.c_str(),"RECREATE");
    ofile.cd();
	LumiWeight = SetWeight(ifilename);
	bool isZjets = false;
	bool isM50sample = false;
	isdata = false;
	isDATA = false;
	bool LFfile = false;
	bool TTfile = false;
	bool HFfile = false;
	if (findString(ifilename, "DY")) isZjets = true;
	if (findString(ifilename, "M50")) isM50sample = true;
	if (findString(ifilename, "M-50")) isM50sample = true;
	if (findString(ifilename, "M_50")) isM50sample = true;
	if (findString(ifilename, "ZJets")) isZjets = true;
	if (findString(ifilename, "Run")||findString(ifilename, "Prompt")||findString(ifilename, "Aug")||
		findString(ifilename, "Jul")||findString(ifilename, "recovered")) {
		isdata = true;
		isDATA = true;
	}
	if (findString(ofilename, "TT")) TTfile = true;
	if (findString(ofilename, "LF")) LFfile = true;
	if (findString(ofilename, "HF")) HFfile = true;
	
	if (debug) cout << "HFfile is " << HFfile << " LumiWeight is " << LumiWeight << " isDATA " << isDATA << endl;
	
	
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
	//FOM_tree->Branch("RegHmass",&RegHmass, "RegHmass/F");
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
	FOM_tree->Branch("PUweight2012",&PUweight2012, "PUweight2012/F");
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
	FOM_tree->Branch("DitauMass",&DitauMass, "DitauMass/F");
	FOM_tree->Branch("pZeta25",&pZeta25, "pZeta25/F");
	FOM_tree->Branch("pZeta45",&pZeta45, "pZeta45/F");
	FOM_tree->Branch("pZeta65",&pZeta65, "pZeta65/F");
	FOM_tree->Branch("pZeta85",&pZeta85, "pZeta85/F");	
	FOM_tree->Branch("CSV0up",&CSVup[0], "CSV0up/F");
	FOM_tree->Branch("CSV1up",&CSVup[1], "CSV1up/F");
	FOM_tree->Branch("CSV0down",&CSVdown[0], "CSV0down/F");
	FOM_tree->Branch("CSV1down",&CSVdown[1], "CSV1down/F");
	FOM_tree->Branch("csv0Fup",&csvFup[0], "csv0Fup/F");
	FOM_tree->Branch("csv1Fup",&csvFup[1], "csv1Fup/F");
	FOM_tree->Branch("csv0Fdown",&csvFdown[0], "csv0Fdown/F");
	FOM_tree->Branch("csv1Fdown",&csvFdown[1], "csv1Fdown/F");
	FOM_tree->Branch("JER0Eup",&JER_e_up[0], "JER0Eup/F");
	FOM_tree->Branch("JER1Eup",&JER_e_up[1], "JER1Eup/F");
	FOM_tree->Branch("JER0Edown",&JER_e_down[0], "JER0Edown/F");
	FOM_tree->Branch("JER1Edown",&JER_e_down[1], "JER1Edown/F");
	FOM_tree->Branch("JER0PtUp",&JER_pt_up[0], "JER0PtUp/F");
	FOM_tree->Branch("JER1PtUp",&JER_pt_up[1], "JER1PtUp/F");
	FOM_tree->Branch("JER0PtDown",&JER_pt_down[0], "JER0PtDown/F");
	FOM_tree->Branch("JER1PtDown",&JER_pt_down[1], "JER1PtDown/F");
	FOM_tree->Branch("JES0Eup",&JES_e_up[0], "JES0Eup/F");
	FOM_tree->Branch("JES1Eup",&JES_e_up[1], "JES1Eup/F");
	FOM_tree->Branch("JES0Edown",&JES_e_down[0], "JES0Edown/F");
	FOM_tree->Branch("JES1Edown",&JES_e_down[1], "JES1Edown/F");
	FOM_tree->Branch("JES0PtUp",&JES_pt_up[0], "JES0PtUp/F");
	FOM_tree->Branch("JES1PtUp",&JES_pt_up[1], "JES1PtUp/F");
	FOM_tree->Branch("JES0PtDown",&JES_pt_down[0], "JES0PtDown/F");
	FOM_tree->Branch("JES1PtDown",&JES_pt_down[1], "JES1PtDown/F");
	FOM_tree->Branch("JERHmassUp",&JERHmassUp, "JERHmassUp/F");
	FOM_tree->Branch("JERHmassDown",&JERHmassDown, "JERHmassDown/F");
	FOM_tree->Branch("JESHmassUp",&JESHmassUp, "JESHmassUp/F");
	FOM_tree->Branch("JESHmassDown",&JESHmassDown, "JESHmassDown/F");
	FOM_tree->Branch("CSVHmassUp",&CSVHmassUp, "CSVHmassUp/F");
	FOM_tree->Branch("CSVHmassFUp",&CSVHmassFUp, "CSVHmassFUp/F");
	FOM_tree->Branch("CSVHmassDown",&CSVHmassDown, "CSVHmassDown/F");
	FOM_tree->Branch("CSVHmassFDown",&CSVHmassFDown, "CSVHmassFDown/F");
	
	
	TTree *TMVA_tree = new TTree("TMVA_tree","Tree for TMVA input");
	TMVA_tree->Branch("nJets",&nJets, "nJets/I");
	TMVA_tree->Branch("Naj",&Naj, "Naj/I");
	TMVA_tree->Branch("Nab",&Nab, "Nab/I");
	TMVA_tree->Branch("naJets",&naJets, "naJets/I");
	TMVA_tree->Branch("eventFlavor",&eventFlavor, "eventFlavor/I");
	TMVA_tree->Branch("CSV0",&CSV0, "CSV0/F");
	TMVA_tree->Branch("CSV1",&CSV1, "CSV1/F");
	TMVA_tree->Branch("Emumass",&Emumass, "Emumass/F");
	TMVA_tree->Branch("Hmass",&Hmass, "Hmass/F");
	TMVA_tree->Branch("oldHmass",&oldHmass, "oldHmass/F");
	//TMVA_tree->Branch("RegHmass",&RegHmass, "RegHmass/F");
	TMVA_tree->Branch("DeltaPhiHV",&DeltaPhiHV, "DeltaPhiHV/F");
	TMVA_tree->Branch("Hpt",&Hpt, "Hpt/F");
	TMVA_tree->Branch("Zpt",&Zpt, "Zpt/F");
	TMVA_tree->Branch("lep0pt",&lep0pt, "lep0pt/F");
	TMVA_tree->Branch("ScalarSumPt",&ScalarSumPt, "ScalarSumPt/F");
	TMVA_tree->Branch("ScalarSumJetPt",&ScalarSumJetPt, "ScalarSumJetPt/F");
	TMVA_tree->Branch("ScalarSumHiggsJetPt",&ScalarSumHiggsJetPt, "ScalarSumHiggsJetPt/F");
	TMVA_tree->Branch("Ht",&Ht, "Ht/F");
	TMVA_tree->Branch("EtaStandDev",&EtaStandDev, "EtaStandDev/F");
	TMVA_tree->Branch("UnweightedEta",&UnweightedEta, "UnweightedEta/F");
	TMVA_tree->Branch("EvntShpCircularity",&EvntShpCircularity, "EvntShpCircularity/F");
	TMVA_tree->Branch("alpha_j",&alpha_j, "alpha_j/F");
	TMVA_tree->Branch("qtb1",&qtb1, "qtb1/F");
	TMVA_tree->Branch("nSV",&nSV, "nSV/I");
	TMVA_tree->Branch("Trigweight",&Trigweight, "Trigweight/F");
	TMVA_tree->Branch("PUweight2012",&PUweight2012, "PUweight2012/F");
	TMVA_tree->Branch("btag2CSF",&btag2CSF, "btag2CSF/F");
	TMVA_tree->Branch("DetaJJ",&DetaJJ, "DetaJJ/F");
	TMVA_tree->Branch("jetCHF0",&jetCHF[0], "jetCHF0/F");
	TMVA_tree->Branch("jetCHF1",&jetCHF[1], "jetCHF1/F");
	TMVA_tree->Branch("jetPt0",&jetPt[0], "jetPt0/F");
	TMVA_tree->Branch("jetPt1",&jetPt[1], "jetPt1/F");
	TMVA_tree->Branch("jetEta0",&jetEta[0], "jetEta0/F");
	TMVA_tree->Branch("jetEta1",&jetEta[1], "jetEta1/F");
	TMVA_tree->Branch("CSVNewShape0",&CSVNewShape[0], "CSVNewShape0/F");
	TMVA_tree->Branch("CSVNewShape1",&CSVNewShape[1], "CSVNewShape1/F");
	TMVA_tree->Branch("lep1pt",&leptonPt[1], "lep1pt/F");
	TMVA_tree->Branch("lep_pfCombRelIso0",&lep_pfCombRelIso[0], "lep_pfCombRelIso0/F");
	TMVA_tree->Branch("lep_pfCombRelIso1",&lep_pfCombRelIso[1], "lep_pfCombRelIso1/F");
	TMVA_tree->Branch("DphiJJ",&DphiJJ, "DphiJJ/F");
	TMVA_tree->Branch("RMS_eta",&RMS_eta, "RMS_eta/F");
	TMVA_tree->Branch("PtbalZH",&PtbalZH, "PtbalZH/F");
	TMVA_tree->Branch("EventPt",&EventPt, "EventPt/F");
	TMVA_tree->Branch("EventMass",&EventMass, "EventMass/F");
	TMVA_tree->Branch("AngleHemu",&AngleHemu, "AngleHemu/F");
	TMVA_tree->Branch("Centrality",&Centrality, "Centrality/F");
	TMVA_tree->Branch("MET",&MET, "MET/F");
	TMVA_tree->Branch("EvntShpAplanarity",&EvntShpAplanarity, "EvntShpAplanarity/F");
	TMVA_tree->Branch("EvntShpSphericity",&EvntShpSphericity, "EvntShpSphericity/F");
	TMVA_tree->Branch("EvntShpIsotropy",&EvntShpIsotropy, "EvntShpIsotropy/F");
	TMVA_tree->Branch("Zphi",&Zphi, "Zphi/F");
	TMVA_tree->Branch("Hphi",&Hphi, "Hphi/F");
	TMVA_tree->Branch("SV_mass",&SV_mass, "SV_mass/F");
	TMVA_tree->Branch("Mte",&Mte, "Mte/F");
	TMVA_tree->Branch("Mtmu",&Mtmu, "Mtmu/F");
	TMVA_tree->Branch("delPullAngle",&delPullAngle, "delPullAngle/F");
	TMVA_tree->Branch("delPullAngle2",&delPullAngle2, "delPullAngle2/F");
	TMVA_tree->Branch("Mt",&Mt, "Mt/F");
	TMVA_tree->Branch("dPhiHMET",&dPhiHMET, "dPhiHMET/F");
	TMVA_tree->Branch("DeltaPhijetMETmin",&DeltaPhijetMETmin, "DeltaPhijetMETmin/F");
	TMVA_tree->Branch("DeltaPhijetMETZtaumin",&DeltaPhijetMETZtaumin, "DeltaPhijetMETZtaumin/F");
	TMVA_tree->Branch("AngleEMU",&AngleEMU, "AngleEMU/F");
	TMVA_tree->Branch("AaronEleMissE",&AaronEleMissE, "AaronEleMissE/F");
	TMVA_tree->Branch("AaronMuMissE",&AaronMuMissE, "AaronMuMissE/F");
	TMVA_tree->Branch("Dphiemu",&Dphiemu, "Dphiemu/F");
	TMVA_tree->Branch("delRjj",&delRjj, "delRjj/F");
	TMVA_tree->Branch("Detaemu",&Detaemu, "Detaemu/F");
	TMVA_tree->Branch("DphiEleMET",&DphiEleMET, "DphiEleMET/F");
	TMVA_tree->Branch("dphiMuMET",&dphiMuMET, "dphiMuMET/F");
	TMVA_tree->Branch("PtbalMETH",&PtbalMETH, "PtbalMETH/F");
	TMVA_tree->Branch("topPt",&topPt, "topPt/F");
	TMVA_tree->Branch("MassEleb0",&MassEleb0, "MassEleb0/F");
	TMVA_tree->Branch("MassMub0",&MassMub0, "MassMub0/F");
	TMVA_tree->Branch("MassEleb1",&MassEleb1, "MassEleb1/F");
	TMVA_tree->Branch("MassMub1",&MassMub1, "MassMub1/F");
	TMVA_tree->Branch("METsig",&METsig, "METsig/F");
	TMVA_tree->Branch("delRemu",&delRemu, "delRemu/F");
	TMVA_tree->Branch("PtbalZMET",&PtbalZMET, "PtbalZMET/F");
	TMVA_tree->Branch("DphiZMET",&DphiZMET, "DphiZMET/F");
	TMVA_tree->Branch("Zmass",&Zmass, "Zmass/F");
	TMVA_tree->Branch("ZmassSVD",&ZmassSVD, "ZmassSVD/F");
	TMVA_tree->Branch("ZmassSVDnegSol",&ZmassSVDnegSol, "ZmassSVDnegSol/F");
	TMVA_tree->Branch("ZmassNegInclu",&ZmassNegInclu, "ZmassNegInclu/F");
	TMVA_tree->Branch("DphiSecondMET",&DphiSecondMET, "DphiSecondMET/F");
	TMVA_tree->Branch("DphiLeadMET",&DphiLeadMET, "DphiLeadMET/F");
	TMVA_tree->Branch("topMass",&topMass, "topMass/F");
	TMVA_tree->Branch("ProjVisT",&ProjVisT, "ProjVisT/F");
	TMVA_tree->Branch("ProjMissT",&ProjMissT, "ProjMissT/F");
	TMVA_tree->Branch("DitauMass",&DitauMass, "DitauMass/F");
	TMVA_tree->Branch("pZeta25",&pZeta25, "pZeta25/F");
	TMVA_tree->Branch("pZeta45",&pZeta45, "pZeta45/F");
	TMVA_tree->Branch("pZeta65",&pZeta65, "pZeta65/F");
	TMVA_tree->Branch("pZeta85",&pZeta85, "pZeta85/F");	
	TMVA_tree->Branch("CSV0up",&CSVup[0], "CSV0up/F");
	TMVA_tree->Branch("CSV1up",&CSVup[1], "CSV1up/F");
	TMVA_tree->Branch("CSV0down",&CSVdown[0], "CSV0down/F");
	TMVA_tree->Branch("CSV1down",&CSVdown[1], "CSV1down/F");
	TMVA_tree->Branch("csv0Fup",&csvFup[0], "csv0Fup/F");
	TMVA_tree->Branch("csv1Fup",&csvFup[1], "csv1Fup/F");
	TMVA_tree->Branch("csv0Fdown",&csvFdown[0], "csv0Fdown/F");
	TMVA_tree->Branch("csv1Fdown",&csvFdown[1], "csv1Fdown/F");
	TMVA_tree->Branch("JER0Eup",&JER_e_up[0], "JER0Eup/F");
	TMVA_tree->Branch("JER1Eup",&JER_e_up[1], "JER1Eup/F");
	TMVA_tree->Branch("JER0Edown",&JER_e_down[0], "JER0Edown/F");
	TMVA_tree->Branch("JER1Edown",&JER_e_down[1], "JER1Edown/F");
	TMVA_tree->Branch("JER0PtUp",&JER_pt_up[0], "JER0PtUp/F");
	TMVA_tree->Branch("JER1PtUp",&JER_pt_up[1], "JER1PtUp/F");
	TMVA_tree->Branch("JER0PtDown",&JER_pt_down[0], "JER0PtDown/F");
	TMVA_tree->Branch("JER1PtDown",&JER_pt_down[1], "JER1PtDown/F");
	TMVA_tree->Branch("JES0Eup",&JES_e_up[0], "JES0Eup/F");
	TMVA_tree->Branch("JES1Eup",&JES_e_up[1], "JES1Eup/F");
	TMVA_tree->Branch("JES0Edown",&JES_e_down[0], "JES0Edown/F");
	TMVA_tree->Branch("JES1Edown",&JES_e_down[1], "JES1Edown/F");
	TMVA_tree->Branch("JES0PtUp",&JES_pt_up[0], "JES0PtUp/F");
	TMVA_tree->Branch("JES1PtUp",&JES_pt_up[1], "JES1PtUp/F");
	TMVA_tree->Branch("JES0PtDown",&JES_pt_down[0], "JES0PtDown/F");
	TMVA_tree->Branch("JES1PtDown",&JES_pt_down[1], "JES1PtDown/F");
	TMVA_tree->Branch("JERHmassUp",&JERHmassUp, "JERHmassUp/F");
	TMVA_tree->Branch("JERHmassDown",&JERHmassDown, "JERHmassDown/F");
	TMVA_tree->Branch("JESHmassUp",&JESHmassUp, "JESHmassUp/F");
	TMVA_tree->Branch("JESHmassDown",&JESHmassDown, "JESHmassDown/F");
	TMVA_tree->Branch("CSVHmassUp",&CSVHmassUp, "CSVHmassUp/F");
	TMVA_tree->Branch("CSVHmassFUp",&CSVHmassFUp, "CSVHmassFUp/F");
	TMVA_tree->Branch("CSVHmassDown",&CSVHmassDown, "CSVHmassDown/F");
	TMVA_tree->Branch("CSVHmassFDown",&CSVHmassFDown, "CSVHmassFDown/F");
	
	
	
	TTree *BDT_tree = new TTree("BDT_tree","Tree for BDT output");
	BDT_tree->Branch("nJets",&nJets, "nJets/I");
	BDT_tree->Branch("naJets",&naJets, "naJets/I");
	BDT_tree->Branch("Naj",&Naj, "Naj/I");
	BDT_tree->Branch("Nab",&Nab, "Nab/I");
	BDT_tree->Branch("eventFlavor",&eventFlavor, "eventFlavor/I");
	BDT_tree->Branch("CSV0",&CSV0, "CSV0/F");
	BDT_tree->Branch("CSV1",&CSV1, "CSV1/F");
	BDT_tree->Branch("Emumass",&Emumass, "Emumass/F");
	BDT_tree->Branch("Hmass",&Hmass, "Hmass/F");
	BDT_tree->Branch("oldHmass",&oldHmass, "oldHmass/F");
	//BDT_tree->Branch("RegHmass",&RegHmass, "RegHmass/F");
	BDT_tree->Branch("DeltaPhiHV",&DeltaPhiHV, "DeltaPhiHV/F");
	BDT_tree->Branch("Hpt",&Hpt, "Hpt/F");
	BDT_tree->Branch("Zpt",&Zpt, "Zpt/F");
	BDT_tree->Branch("lep0pt",&lep0pt, "lep0pt/F");
	BDT_tree->Branch("ScalarSumPt",&ScalarSumPt, "ScalarSumPt/F");
	BDT_tree->Branch("ScalarSumJetPt",&ScalarSumJetPt, "ScalarSumJetPt/F");
	BDT_tree->Branch("ScalarSumHiggsJetPt",&ScalarSumHiggsJetPt, "ScalarSumHiggsJetPt/F");
	BDT_tree->Branch("Ht",&Ht, "Ht/F");
	BDT_tree->Branch("EtaStandDev",&EtaStandDev, "EtaStandDev/F");
	BDT_tree->Branch("UnweightedEta",&UnweightedEta, "UnweightedEta/F");
	BDT_tree->Branch("EvntShpCircularity",&EvntShpCircularity, "EvntShpCircularity/F");
	BDT_tree->Branch("alpha_j",&alpha_j, "alpha_j/F");
	BDT_tree->Branch("qtb1",&qtb1, "qtb1/F");
	BDT_tree->Branch("nSV",&nSV, "nSV/I");
	BDT_tree->Branch("Trigweight",&Trigweight, "Trigweight/F");
	BDT_tree->Branch("PUweight2012",&PUweight2012, "PUweight2012/F");
	BDT_tree->Branch("btag2CSF",&btag2CSF, "btag2CSF/F");
	BDT_tree->Branch("DetaJJ",&DetaJJ, "DetaJJ/F");
	BDT_tree->Branch("jetCHF0",&jetCHF[0], "jetCHF0/F");
	BDT_tree->Branch("jetCHF1",&jetCHF[1], "jetCHF1/F");
	BDT_tree->Branch("jetPt0",&jetPt[0], "jetPt0/F");
	BDT_tree->Branch("jetPt1",&jetPt[1], "jetPt1/F");
	BDT_tree->Branch("jetEta0",&jetEta[0], "jetEta0/F");
	BDT_tree->Branch("jetEta1",&jetEta[1], "jetEta1/F");
	BDT_tree->Branch("CSVNewShape0",&CSVNewShape[0], "CSVNewShape0/F");
	BDT_tree->Branch("CSVNewShape1",&CSVNewShape[1], "CSVNewShape1/F");
	BDT_tree->Branch("lep1pt",&leptonPt[1], "lep1pt/F");
	BDT_tree->Branch("lep_pfCombRelIso0",&lep_pfCombRelIso[0], "lep_pfCombRelIso0/F");
	BDT_tree->Branch("lep_pfCombRelIso1",&lep_pfCombRelIso[1], "lep_pfCombRelIso1/F");
	BDT_tree->Branch("DphiJJ",&DphiJJ, "DphiJJ/F");
	BDT_tree->Branch("RMS_eta",&RMS_eta, "RMS_eta/F");
	BDT_tree->Branch("PtbalZH",&PtbalZH, "PtbalZH/F");
	BDT_tree->Branch("EventPt",&EventPt, "EventPt/F");
	BDT_tree->Branch("EventMass",&EventMass, "EventMass/F");
	BDT_tree->Branch("AngleHemu",&AngleHemu, "AngleHemu/F");
	BDT_tree->Branch("Centrality",&Centrality, "Centrality/F");
	BDT_tree->Branch("MET",&MET, "MET/F");
	BDT_tree->Branch("EvntShpAplanarity",&EvntShpAplanarity, "EvntShpAplanarity/F");
	BDT_tree->Branch("EvntShpSphericity",&EvntShpSphericity, "EvntShpSphericity/F");
	BDT_tree->Branch("EvntShpIsotropy",&EvntShpIsotropy, "EvntShpIsotropy/F");
	BDT_tree->Branch("Zphi",&Zphi, "Zphi/F");
	BDT_tree->Branch("Hphi",&Hphi, "Hphi/F");
	BDT_tree->Branch("SV_mass",&SV_mass, "SV_mass/F");
	BDT_tree->Branch("Mte",&Mte, "Mte/F");
	BDT_tree->Branch("Mtmu",&Mtmu, "Mtmu/F");
	BDT_tree->Branch("delPullAngle",&delPullAngle, "delPullAngle/F");
	BDT_tree->Branch("delPullAngle2",&delPullAngle2, "delPullAngle2/F");
	BDT_tree->Branch("Mt",&Mt, "Mt/F");
	BDT_tree->Branch("dPhiHMET",&dPhiHMET, "dPhiHMET/F");
	BDT_tree->Branch("DeltaPhijetMETmin",&DeltaPhijetMETmin, "DeltaPhijetMETmin/F");
	BDT_tree->Branch("DeltaPhijetMETZtaumin",&DeltaPhijetMETZtaumin, "DeltaPhijetMETZtaumin/F");
	BDT_tree->Branch("AngleEMU",&AngleEMU, "AngleEMU/F");
	BDT_tree->Branch("AaronEleMissE",&AaronEleMissE, "AaronEleMissE/F");
	BDT_tree->Branch("AaronMuMissE",&AaronMuMissE, "AaronMuMissE/F");
	BDT_tree->Branch("Dphiemu",&Dphiemu, "Dphiemu/F");
	BDT_tree->Branch("delRjj",&delRjj, "delRjj/F");
	BDT_tree->Branch("Detaemu",&Detaemu, "Detaemu/F");
	BDT_tree->Branch("DphiEleMET",&DphiEleMET, "DphiEleMET/F");
	BDT_tree->Branch("dphiMuMET",&dphiMuMET, "dphiMuMET/F");
	BDT_tree->Branch("PtbalMETH",&PtbalMETH, "PtbalMETH/F");
	BDT_tree->Branch("topPt",&topPt, "topPt/F");
	BDT_tree->Branch("MassEleb0",&MassEleb0, "MassEleb0/F");
	BDT_tree->Branch("MassMub0",&MassMub0, "MassMub0/F");
	BDT_tree->Branch("MassEleb1",&MassEleb1, "MassEleb1/F");
	BDT_tree->Branch("MassMub1",&MassMub1, "MassMub1/F");
	BDT_tree->Branch("METsig",&METsig, "METsig/F");
	BDT_tree->Branch("delRemu",&delRemu, "delRemu/F");
	BDT_tree->Branch("PtbalZMET",&PtbalZMET, "PtbalZMET/F");
	BDT_tree->Branch("DphiZMET",&DphiZMET, "DphiZMET/F");
	BDT_tree->Branch("Zmass",&Zmass, "Zmass/F");
	BDT_tree->Branch("ZmassSVD",&ZmassSVD, "ZmassSVD/F");
	BDT_tree->Branch("ZmassSVDnegSol",&ZmassSVDnegSol, "ZmassSVDnegSol/F");
	BDT_tree->Branch("ZmassNegInclu",&ZmassNegInclu, "ZmassNegInclu/F");
	BDT_tree->Branch("DphiSecondMET",&DphiSecondMET, "DphiSecondMET/F");
	BDT_tree->Branch("DphiLeadMET",&DphiLeadMET, "DphiLeadMET/F");
	BDT_tree->Branch("topMass",&topMass, "topMass/F");
	BDT_tree->Branch("ProjVisT",&ProjVisT, "ProjVisT/F");
	BDT_tree->Branch("ProjMissT",&ProjMissT, "ProjMissT/F");
	BDT_tree->Branch("DitauMass",&DitauMass, "DitauMass/F");
	BDT_tree->Branch("pZeta25",&pZeta25, "pZeta25/F");
	BDT_tree->Branch("pZeta45",&pZeta45, "pZeta45/F");
	BDT_tree->Branch("pZeta65",&pZeta65, "pZeta65/F");
	BDT_tree->Branch("pZeta85",&pZeta85, "pZeta85/F");	
	BDT_tree->Branch("CSV0up",&CSVup[0], "CSV0up/F");
	BDT_tree->Branch("CSV1up",&CSVup[1], "CSV1up/F");
	BDT_tree->Branch("CSV0down",&CSVdown[0], "CSV0down/F");
	BDT_tree->Branch("CSV1down",&CSVdown[1], "CSV1down/F");
	BDT_tree->Branch("csv0Fup",&csvFup[0], "csv0Fup/F");
	BDT_tree->Branch("csv1Fup",&csvFup[1], "csv1Fup/F");
	BDT_tree->Branch("csv0Fdown",&csvFdown[0], "csv0Fdown/F");
	BDT_tree->Branch("csv1Fdown",&csvFdown[1], "csv1Fdown/F");
	BDT_tree->Branch("JER0Eup",&JER_e_up[0], "JER0Eup/F");
	BDT_tree->Branch("JER1Eup",&JER_e_up[1], "JER1Eup/F");
	BDT_tree->Branch("JER0Edown",&JER_e_down[0], "JER0Edown/F");
	BDT_tree->Branch("JER1Edown",&JER_e_down[1], "JER1Edown/F");
	BDT_tree->Branch("JER0PtUp",&JER_pt_up[0], "JER0PtUp/F");
	BDT_tree->Branch("JER1PtUp",&JER_pt_up[1], "JER1PtUp/F");
	BDT_tree->Branch("JER0PtDown",&JER_pt_down[0], "JER0PtDown/F");
	BDT_tree->Branch("JER1PtDown",&JER_pt_down[1], "JER1PtDown/F");
	BDT_tree->Branch("JES0Eup",&JES_e_up[0], "JES0Eup/F");
	BDT_tree->Branch("JES1Eup",&JES_e_up[1], "JES1Eup/F");
	BDT_tree->Branch("JES0Edown",&JES_e_down[0], "JES0Edown/F");
	BDT_tree->Branch("JES1Edown",&JES_e_down[1], "JES1Edown/F");
	BDT_tree->Branch("JES0PtUp",&JES_pt_up[0], "JES0PtUp/F");
	BDT_tree->Branch("JES1PtUp",&JES_pt_up[1], "JES1PtUp/F");
	BDT_tree->Branch("JES0PtDown",&JES_pt_down[0], "JES0PtDown/F");
	BDT_tree->Branch("JES1PtDown",&JES_pt_down[1], "JES1PtDown/F");
	BDT_tree->Branch("JERHmassUp",&JERHmassUp, "JERHmassUp/F");
	BDT_tree->Branch("JERHmassDown",&JERHmassDown, "JERHmassDown/F");
	BDT_tree->Branch("JESHmassUp",&JESHmassUp, "JESHmassUp/F");
	BDT_tree->Branch("JESHmassDown",&JESHmassDown, "JESHmassDown/F");
	BDT_tree->Branch("CSVHmassUp",&CSVHmassUp, "CSVHmassUp/F");
	BDT_tree->Branch("CSVHmassFUp",&CSVHmassFUp, "CSVHmassFUp/F");
	BDT_tree->Branch("CSVHmassDown",&CSVHmassDown, "CSVHmassDown/F");
	BDT_tree->Branch("CSVHmassFDown",&CSVHmassFDown, "CSVHmassFDown/F");
	
	
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
	int event =0, Ntree =0;
	float N_TopCR = 0.0, N_SingleTopCR = 0.0, N_LFCR = 0.0, N_HFCR =0.0;
	float FailedJetID=0.0, NTMVAtree =0.0, NBDTtree =0.0, NBDTbtree = 0.0;
	int  NNegEleMissE = 0, NPosEleMissE = 0, NPosMuMissE = 0, NNegMuMissE = 0;
	int NBothMissENeg = 0, NMixedEleMissENeg = 0, NMixedMuonMissENeg =0;
	float SFUnc_TOPCR = 0.0, SFUnc_SingleTOP = 0.0, SFUnc_LFCR = 0.0, SFUnc_HFCR = 0.0;
	bool firstevent = true;
	float Nweighted_Vtype =0.0, Nweighted_trigger = 0.0, Nweighted_preselect =0.0, Nweighted_Mjj =0.0, Nweighted_DphiZMET =0.0, Nweighted_Memu =0.0;
	float Nweighted_Naj =0.0, Nweighted_CSV0 =0.0, Nweighted_EfakeCuts =0.0, Nweighted_jetCHF0 = 0.0;
	float Nweighted_Mt = 0.0, Nweighted_Pzeta =0.0, Nweighted_DeltaPhiHV = 0.0, event_weighted = 0.0;
	float Nweighted_TopCR = 0.0, Nweighted_SingleTopCR = 0.0, Nweighted_LFCR = 0.0, Nweighted_HFCR =0.0;	
	
    do {
		
		if (debug)cout << "begining of event loop " << event << endl;
		
        weight = 1.0;
		float LFScaleFactor = 1.100305483;
		float TTScaleFactor = 1.073681624;
		event++;
		
		EleTrigWeight = 1.0, MuonTrigWeight = 1.0, WP95weight = 1.0, MuIDweight = 1.0;
		PUweight2012 = sample.PUweight;
		weight = LumiWeight*PUweight2012;
		Trigweight = 1.0;
		if (debug)cout << "weight is " << weight << " LumiWeight is " << LumiWeight << " PUweight2012 " << PUweight2012 << endl;		
		
		event_weighted = event_weighted + weight;
		
		if (((isM50sample && sample.genZpt < 100) || !isM50sample) && ((HFfile && sample.eventFlav == 5)||!HFfile) && ((LFfile && sample.eventFlav!=5)||!LFfile)){
			if (debug)cout << "HFfile is " << HFfile << " sample.eventFlav is " << sample.eventFlav << " isDATA " << isDATA << endl;
			if (sample.Vtype == 5 && sample.EVENT_json == 1 ){
				if (!(event%500))  std::cout << "entered event loop " << event << std::endl;
				if (firstevent)  std::cout << "entered event loop " << event << std::endl;
				
				
				
				if (debug)cout << "Electron pt "<<sample.vLepton_pt[0] << " Muon pt " <<sample.vLepton_pt[1] << endl;
				if (sample.vLepton_pt[1]>10 && sample.vLepton_pt[0]>10){
					//Trigger weight named after leading lepton
					EleTrigWeight = AssignWeight(EleTrigWeightTreehigh, sample.vLepton_pt[0], fabs(sample.vLepton_eta[0]));
					MuonTrigWeight =AssignWeight(MuonTrigWeightTreehigh, sample.vLepton_pt[1], fabs(sample.vLepton_eta[1]));
					EleTrigWeight = EleTrigWeight*AssignWeight(MuonTrigWeightTreelow, sample.vLepton_pt[0], fabs(sample.vLepton_eta[0]));
					MuonTrigWeight = MuonTrigWeight*AssignWeight(EleTrigWeightTreelow, sample.vLepton_pt[1], fabs(sample.vLepton_eta[1]));
					if (debug)cout << "Electron Trigger Weight is "<<EleTrigWeight<< endl;
					if (debug)cout << "Muon trigger weight is "<<MuonTrigWeight << endl;
					//Which trigger more likely to fire if both leptons above upper leg threshold
					if( sample.vLepton_pt[1]>20 && sample.vLepton_pt[0]>20){
						if (EleTrigWeight>MuonTrigWeight){MuonTrigWeight = 1.0;
						}else{ EleTrigWeight = 1.0;}
						if (debug)cout << "Both leptons above 20 GeV "<<MuonTrigWeight << " " <<  EleTrigWeight<< endl;
					}
					if( sample.vLepton_pt[1]>10 && sample.vLepton_pt[1]<20){ // Muon bottom leg
						EleTrigWeight = 1.0;
						if (debug)cout << "Inside the Muon bottom Leg "<<MuonTrigWeight << " " <<  EleTrigWeight<< endl;
					}
					if( sample.vLepton_pt[0]>10 && sample.vLepton_pt[0]<20){ // Electron bottom leg
						MuonTrigWeight = 1.0;
						if (debug)cout << "Inside the Electron bottom Leg "<< MuonTrigWeight << " " <<  EleTrigWeight<< endl;
					}
				}
				if (sample.vLepton_pt[1]>20 || sample.vLepton_pt[0]>20){
					if (fabs(MuonTrigWeight-1.0)>0.001 && fabs(EleTrigWeight-1.0)>0.001) {
						cout << "WE HAVE A BUG WITH THE TRIGGER WEIGHT" << endl;
						cout << "Electron pt "<<sample.vLepton_pt[0] << " Muon pt " <<sample.vLepton_pt[1] << endl;
						cout << "Electron trig weight "<<  EleTrigWeight << " Muon trig weight " << MuonTrigWeight << endl;
					}}
				// V21 VHbb electron weights
				if(sample.vLepton_pt[0]>10) {
					WP95weight = AssignWeight(lowptEleIDWeightTree, sample.vLepton_pt[0], sample.vLepton_eta[0]);
					if (debug)cout << "Electron WP95ID low pt Weight is "<<WP95weight<< endl;
					WP95weight = 1.0;
				}			
				if(sample.vLepton_pt[0]>20) {
					WP95weight = AssignWeight(EleIDWeightTree, sample.vLepton_pt[0], sample.vLepton_eta[0]);
					if (debug)cout << "Electron ID 95 Weight is "<<WP95weight<< endl;
				}
				// V21 VHbb muon weights
				if(sample.vLepton_pt[1]>10){
					MuIDweight = AssignWeight(lowptMuonIDWeightTree, sample.vLepton_pt[1], sample.vLepton_eta[1]);
					if (debug)cout << "Muon ID low pt Weight is "<<MuIDweight<< endl;
				}			
				if(sample.vLepton_pt[1]>10){
					MuIDweight = AssignWeight(MuonIDWeightTree, sample.vLepton_pt[1], sample.vLepton_eta[1]);
					if (debug)cout << "Muon ID Weight is "<<MuIDweight<< endl;
				}
				Trigweight = MuonTrigWeight*EleTrigWeight*WP95weight*MuIDweight;
				weight = weight*MuonTrigWeight*EleTrigWeight*WP95weight*MuIDweight;
				
				if (LFfile)LFScaleFactor = 1.00;
				if (TTfile)TTScaleFactor = 1.00;
				weight = weight*LFScaleFactor*TTScaleFactor;
				
				if (isDATA) {weight = 1.0; Trigweight=1.0; }
				if (debug)cout << "Trigweight is " << Trigweight << " weight is " << weight << " isDATA " << isDATA << endl;
				
				isdata= false;
				
				std::vector< std::pair<size_t,double> > indexedJetPt;
				std::vector<size_t> PtSortedJetIndex;
				std::vector< std::pair<size_t,double> > indexedJetPtJESup;
				std::vector<size_t> PtJESupSortedJetIndex;
				std::vector< std::pair<size_t,double> > indexedJetPtJESdown;
				std::vector<size_t> PtJESdownSortedJetIndex;
				std::vector< std::pair<size_t,double> > indexedJetPtJERup;
				std::vector<size_t> PtJERupSortedJetIndex;
				std::vector< std::pair<size_t,double> > indexedJetPtJERdown;
				std::vector<size_t> PtJERdownSortedJetIndex;
				std::vector< std::pair<size_t,double> > indexedJetCSV;
				std::vector<size_t> CSVSortedJetIndex;
				std::vector< std::pair<size_t,double> > indexedJetCSVup;
				std::vector<size_t> CSVupSortedJetIndex;
				std::vector< std::pair<size_t,double> > indexedJetCSVdown;
				std::vector<size_t> CSVdownSortedJetIndex;
				std::vector< std::pair<size_t,double> > indexedJetCSVFup;
				std::vector<size_t> CSVFupSortedJetIndex;
				std::vector< std::pair<size_t,double> > indexedJetCSVFdown;
				std::vector<size_t> CSVFdownSortedJetIndex;
				std::vector< std::pair<size_t,double> > indexedPt;
				std::vector<size_t> PtSortedIndex;
								
				// Analysis loop.
				// One can access the ntuple leaves directly from sample object
				
				if (debug)std::cout << "Before initilize variables" << std::endl;
				float jetPttmp[10], jetEtatmp[10], jetPhitmp[10], jetCSVtmp[10], jetCHFtmp[10], CSVNewShapetmp[10];
				float jetPtRawtmp[10], jetEtmp[10], jetVtx3dLtmp[10], jetVtx3deLtmp[10], jetVtxPttmp[10], jetVtxMasstmp[10], jetPtLeadTracktmp[10];
				float jetNconstintuentstmp[10], jetCEFtmp[10], jetNCHtmp[10], jetJECUnctmp[10], jetEttmp[10], jetMttmp[10], jetPtRawJERtmp[10], jetGenPttmp[10]; 
				float jetFlavor[10];
				
				for (int j=0; j < 15; ++j) {
					CSVNewShapetmp[j] = -99.99;
					jetPttmp[j] = -99.99;
					jetEtatmp[j] = -99.99;
					jetPhitmp[j] = -99.99;
					jetCSVtmp[j] = -99.99;
					jetCHFtmp[j] = -99.99;
					//RegjetPttmp[j] = -99.99;
					jetPtRawtmp[j] = -99.99;
					jetEtmp[j] = -99.99;
					jetVtx3dLtmp[j] = -99.99;
					jetVtx3deLtmp[j] = -99.99;
					jetVtxPttmp[j] = -99.99;
					jetVtxMasstmp[j] = -99.99;
					jetPtLeadTracktmp[j] = -99.99;
					jetNconstintuentstmp[j] = -99.99;
					jetCEFtmp[j] = -99.99;
					jetNCHtmp[j] = -99.99;
					jetJECUnctmp[j] = -99.99;
					jetEttmp[j] = -99.99;
					jetMttmp[j] = -99.99;
					jetPtRawJERtmp[j] = -99.99;
					jetGenPttmp[j] = -99.99;
					jetFlavor[j] = -99.99;
				}					
				for (int i=0; i < 5; ++i) {
					CSVNewShape[i] = -99.99;
					jetPt[i] = -99.99;
					jetEta[i] = -99.99;
					jetPhi[i] = -99.99;
					jetCSV[i] = -99.99;
					jetCHF[i] = -99.99;
					//RegjetPt[i] = -99.99;
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
				}
				if (debug)cout << "Number of Jets " <<sample.nhJets << " " <<sample.naJets << endl;
				if(debug)cout << "Number of Leptons " << sample.nvlep << " " << sample.nalep << endl;
				for (int j=0; j < 5; ++j) {
					CSVup[j] = -99.99;
					CSVdown[j] = -99.99;
					csvFup[j] = -99.99;
					csvFdown[j] = -99.99;
					JER_e_up[j] = -99.99;
					JER_e_down[j] = -99.99;
					JES_e_up[j] = -99.99;
					JES_e_down[j] = -99.99;
					JER_pt_up[j] = -99.99;
					JER_pt_down[j] = -99.99;
					JES_pt_up[j] = -99.99;
					JES_pt_down[j] = -99.99;
				}
				
				
				nJets =0, nSV =-99, nMuons = 0,  nElectrons = 0,  nLeptons = 0, Na_lep= 0, nPV= -99, MET= -99.99, Naj= 0, Nab =0, eventFlavor = -99;
				CSV0 = -1.0, CSV1 = -1.0, Emumass = -99.99, Hmass = -99.99, oldHmass = -99.99, DeltaPhiHV = -99.99, Hpt = -99.99, Zpt = -99.99;
				lep0pt = -99.99, ScalarSumPt = -99.99, EtaStandDev = -99.99, UnweightedEta = -99.99, EvntShpCircularity = -99.99;
				alpha_j = -99.99, qtb1 = 0.0, DphiJJ = -99.99,  btag2CSF = 1.0;
				RMS_eta = -99.99, PtbalZH= -99.99, EventPt= -99.99, AngleHemu= -99.99, Centrality = -99.99;
				qtlep1 = 0.0, alpha_lep = -99.99;
				//RegHmass = -99.99;
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
				JERHmassUp = -99.99, JERHmassDown = -99.99, JESHmassUp = -99.99, JESHmassDown = -99.99;
				CSVHmassUp = -99.99, CSVHmassFUp = -99.99, CSVHmassDown = -99.99, CSVHmassFDown = -99.99;

				
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
				double pllep1 = 0.0, pllep2 = 0.0;
				double plb1 = 0.0, plb2 = 0.0;
				
				AngleEMU = -99.99, CosThetaEle = -99.99, CosThetaMu = -99.99;
				EleMissE = -99.99, MuonMissE =-99.99, Dphiemu = -99.99;
				Zmass = -99.99, ZmassSVD = -99.99, ZmassSVDnegSol = -99.99, ZmassNegInclu = -99.99;
				AaronEleMissE = -599.99, AaronMuMissE = -599.99;
				
				double CSVshapeNew = -99.99;
				
				if(debug)std::cout << "variables initialized " << std::endl;
				
				for (int k=0;k<sample.nhJets;k++){
					if (debug) cout << "for "<< k << "th pass" << endl;
					indexedJetPt.push_back(std::pair<size_t,double>(k,(double) sample.hJet_pt[k]));
					hallhJet_pt.Fill(sample.hJet_pt[k]); 
					nJets++;
					//indexedJetCSV.push_back(std::pair<size_t,double>(k,(double) sample.hJet_csv[k]));
					jetPttmp[k] = sample.hJet_pt[k];
					jetEtatmp[k] = sample.hJet_eta[k];
					jetPhitmp[k] = sample.hJet_phi[k];
					jetCSVtmp[k] = sample.hJet_csv[k];
					jetCHFtmp[k] = sample.hJet_chf[k];
					jetPtRawtmp[k] = sample.hJet_ptRaw[k];
					jetEtmp[k] = sample.hJet_e[k];
					jetVtx3dLtmp[k] = sample.hJet_vtx3dL[k];
					jetVtx3deLtmp[k] = sample.hJet_vtx3deL[k];
					jetVtxPttmp[k] = sample.hJet_vtxPt[k];
					jetVtxMasstmp[k] = sample.hJet_vtxMass[k];
					jetPtLeadTracktmp[k] = sample.hJet_ptLeadTrack[k];
					jetNconstintuentstmp[k] = sample.hJet_nconstituents[k];
					jetCEFtmp[k] = sample.hJet_cef[k];
					jetNCHtmp[k] = sample.hJet_nch[k];
					jetJECUnctmp[k] = sample.hJet_JECUnc[k];
					jetEttmp[k] = evalEt(jetPt[k],jetEta[k],jetPhi[k],jetE[k]);
					jetMttmp[k] = evalMt(jetPt[k],jetEta[k],jetPhi[k],jetE[k]);
					jetPtRawJERtmp[k] = evalJERBias(jetPtRaw[k], sample.hJet_genPt[k], jetEta[k]);
					jetGenPttmp[k]=sample.hJet_genPt[k];
					jetFlavor[k] = sample.hJet_flavour[k];
					//if (indexedJetCSV[0].first == 0) RegjetPt[0] = sample.hJet_genPtReg0;
					if (debug)	cout << "CSV discriminator: " << sample.hJet_csv[k] << endl;
					CSVshapeNew = corrCSV(btagNew, sample.hJet_csv[k], sample.hJet_flavour[k]);
					//if(sample.hJet_csv[k]<=0 || sample.hJet_csv[k]>=1) CSVshapeNew=sample.hJet_csv[k];
					//else if(sample.hJet_flavour[k]==0) CSVshapeNew=sample.hJet_csv[k];
					//else if(fabs(sample.hJet_flavour[k])==5) CSVshapeNew=btagNew->ib->Eval(sample.hJet_csv[k]);
					//else if(fabs(sample.hJet_flavour[k])==4) CSVshapeNew=btagNew->ic->Eval(sample.hJet_csv[k]);
					//else if(fabs(sample.hJet_flavour[k])!=5 && fabs(sample.hJet_flavour[k])!=4)  CSVshapeNew=btagNew->il->Eval(sample.hJet_csv[k]);
					if (!(event%500)&&debug) cout << "The orignial CSV value was " << sample.hJet_csv[k] << "  the corrected value is " << CSVshapeNew << endl;
					indexedJetCSV.push_back(std::pair<size_t,double>(k,(double) CSVshapeNew));
					CSVNewShapetmp[k] = CSVshapeNew;
					CSVshapeNew = -99.99;
					indexedJetCSVup.push_back(std::pair<size_t,double>(k,(double) corrCSV(btagUp, sample.hJet_csv[k],sample.hJet_flavour[k])));
					indexedJetCSVdown.push_back(std::pair<size_t,double>(k,(double) corrCSV(btagDown, sample.hJet_csv[k],sample.hJet_flavour[k])));
					indexedJetCSVFup.push_back(std::pair<size_t,double>(k,(double) corrCSV(btagFUp, sample.hJet_csv[k],sample.hJet_flavour[k])));
					indexedJetCSVFdown.push_back(std::pair<size_t,double>(k,(double) corrCSV(btagFDown, sample.hJet_csv[k],sample.hJet_flavour[k])));
					indexedJetPtJESup.push_back(std::pair<size_t,double>(k,(double) sample.hJet_pt[k]*(1+sample.hJet_JECUnc[k])));
					indexedJetPtJESdown.push_back(std::pair<size_t,double>(k,(double) sample.hJet_pt[k]*(1-sample.hJet_JECUnc[k])));
					indexedJetPtJERup.push_back(std::pair<size_t,double>(k,(double) JERSys(true,sample.hJet_eta[k],sample.hJet_pt[k],sample.hJet_genPt[k])));
					indexedJetPtJERdown.push_back(std::pair<size_t,double>(k,(double) JERSys(false,sample.hJet_eta[k],sample.hJet_pt[k],sample.hJet_genPt[k])));
					CSVshapeNew = -99.99;
				}// end k for loop
				
				float MinDphiaJet = 99.99;
				if(debug)std::cout << "Number additional Jets " << sample.naJets << std::endl;
				for (int a=0;a<sample.naJets;a++){
					hallhJet_pt.Fill(sample.aJet_pt[a]);
					nJets++; 
					//if(sample.aJet_csv[a]<=0 || sample.aJet_csv[a]>=1) CSVshapeNew=sample.aJet_csv[a];
					//else if(sample.aJet_flavour[a]==0) CSVshapeNew=sample.aJet_csv[a];
					//else if(fabs(sample.aJet_flavour[a])==5) CSVshapeNew=btagNew->ib->Eval(sample.aJet_csv[a]);
					//else if(fabs(sample.aJet_flavour[a])==4) CSVshapeNew=btagNew->ic->Eval(sample.aJet_csv[a]);
					//else if(fabs(sample.aJet_flavour[a])!=5 && fabs(sample.aJet_flavour[a])!=4)  CSVshapeNew=btagNew->il->Eval(sample.aJet_csv[a]);	
					if (a<14){			
					jetPttmp[sample.nhJets+a] = sample.aJet_pt[a];
					jetEtatmp[sample.nhJets+a] = sample.aJet_eta[a];
					jetPhitmp[sample.nhJets+a] = sample.aJet_phi[a];
					jetCSVtmp[sample.nhJets+a] = sample.aJet_csv[a];
					jetCHFtmp[sample.nhJets+a] = sample.aJet_chf[a];
					//jetPtRawtmp[sample.nhJets+a] = sample.aJet_ptRaw[a];
					jetEtmp[sample.nhJets+a] = sample.aJet_e[a];
					jetVtx3dLtmp[sample.nhJets+a] = sample.aJet_vtx3dL[a];
					jetVtx3deLtmp[sample.nhJets+a] = sample.aJet_vtx3deL[a];
					//jetVtxPttmp[sample.nhJets+a] = sample.aJet_vtxPt[a];
					jetVtxMasstmp[sample.nhJets+a] = sample.aJet_vtxMass[a];
					//jetPtLeadTracktmp[sample.nhJets+a] = sample.aJet_ptLeadTrack[a];
					jetNconstintuentstmp[sample.nhJets+a] = sample.aJet_nconstituents[a];
					jetCEFtmp[sample.nhJets+a] = sample.aJet_cef[a];
					jetNCHtmp[sample.nhJets+a] = sample.aJet_nch[a];
					jetJECUnctmp[sample.nhJets+a] = sample.aJet_JECUnc[a];
					jetEttmp[sample.nhJets+a] = evalEt(sample.aJet_pt[a],sample.aJet_eta[a],sample.aJet_phi[a],sample.aJet_e[a]);
					jetMttmp[sample.nhJets+a] = evalMt(sample.aJet_pt[a],sample.aJet_eta[a],sample.aJet_phi[a],sample.aJet_e[a]);
					//jetPtRawJERtmp[sample.nhJets+a] = evalJERBias(sample.aJet_ptRaw[a], sample.aJet_genPt[a], sample.aJet_eta[a]);
					jetGenPttmp[sample.nhJets+a]= sample.aJet_genPt[a];
					jetFlavor[sample.nhJets+a]= sample.aJet_flavour[a];
					indexedJetPt.push_back(std::pair<size_t,double>(sample.nhJets+a,(double) sample.aJet_pt[a]));	
					CSVshapeNew = corrCSV(btagNew, sample.hJet_csv[a], sample.aJet_flavour[a]);
					CSVNewShapetmp[sample.nhJets+a] = CSVshapeNew;
						indexedJetCSVup.push_back(std::pair<size_t,double>(sample.nhJets+a,(double) corrCSV(btagUp, sample.aJet_csv[a],sample.aJet_flavour[a])));
						indexedJetCSVdown.push_back(std::pair<size_t,double>(sample.nhJets+a,(double) corrCSV(btagDown, sample.aJet_csv[a],sample.aJet_flavour[a])));
						indexedJetCSVFup.push_back(std::pair<size_t,double>(sample.nhJets+a,(double) corrCSV(btagFUp, sample.aJet_csv[a],sample.aJet_flavour[a])));
						indexedJetCSVFdown.push_back(std::pair<size_t,double>(sample.nhJets+a,(double) corrCSV(btagFDown, sample.aJet_csv[a],sample.aJet_flavour[a])));
						indexedJetPtJESup.push_back(std::pair<size_t,double>(sample.nhJets+a,(double) sample.aJet_pt[a]*(1+sample.aJet_JECUnc[a])));
						indexedJetPtJESdown.push_back(std::pair<size_t,double>(sample.nhJets+a,(double) sample.aJet_pt[a]*(1-sample.aJet_JECUnc[a])));
						indexedJetPtJERup.push_back(std::pair<size_t,double>(sample.nhJets+a,(double) JERSys(true,sample.aJet_eta[a],sample.aJet_pt[a],sample.aJet_genPt[a])));
						indexedJetPtJERdown.push_back(std::pair<size_t,double>(sample.nhJets+a,(double) JERSys(false,sample.aJet_eta[a],sample.aJet_pt[a],sample.aJet_genPt[a])));						
					}
					CSVshapeNew = -99.99;
					ScalarSumJetPt = ScalarSumJetPt + sample.aJet_pt[a];
				}
				indexedPt.push_back(std::pair<size_t,double>(2,(double) sample.vLepton_pt[0]));
				indexedPt.push_back(std::pair<size_t,double>(3,(double) sample.vLepton_pt[1]));
				
				
				
				std::sort(indexedJetPt.begin(),indexedJetPt.end(),::IndexedQuantityGreaterThan<double>);
				std::sort(indexedJetCSV.begin(),indexedJetCSV.end(),::IndexedQuantityGreaterThan<double>);
				std::sort(indexedJetCSVup.begin(),indexedJetCSVup.end(),::IndexedQuantityGreaterThan<double>);
				std::sort(indexedJetCSVdown.begin(),indexedJetCSVdown.end(),::IndexedQuantityGreaterThan<double>);
				std::sort(indexedJetCSVFup.begin(),indexedJetCSVFup.end(),::IndexedQuantityGreaterThan<double>);
				std::sort(indexedJetCSVFdown.begin(),indexedJetCSVFdown.end(),::IndexedQuantityGreaterThan<double>);
				std::sort(indexedJetPtJESup.begin(),indexedJetPtJESup.end(),::IndexedQuantityGreaterThan<double>);
				std::sort(indexedJetPtJESdown.begin(),indexedJetPtJESdown.end(),::IndexedQuantityGreaterThan<double>);
				std::sort(indexedJetPtJERup.begin(),indexedJetPtJERup.end(),::IndexedQuantityGreaterThan<double>);
				std::sort(indexedJetPtJERdown.begin(),indexedJetPtJERdown.end(),::IndexedQuantityGreaterThan<double>);
				
				for (size_t i = 0 ; (i != indexedJetPt.size()) ; ++i) {   PtSortedJetIndex.push_back(indexedJetPt[i].first);        }
				for (size_t i = 0 ; (i != indexedJetCSV.size()) ; ++i) {   CSVSortedJetIndex.push_back(indexedJetCSV[i].first);        }
				for (size_t i = 0 ; (i != indexedJetCSVup.size()) ; ++i) {   CSVupSortedJetIndex.push_back(indexedJetCSVup[i].first);        }
				for (size_t i = 0 ; (i != indexedJetCSVdown.size()) ; ++i) {   CSVdownSortedJetIndex.push_back(indexedJetCSVdown[i].first);        }
				for (size_t i = 0 ; (i != indexedJetCSVFup.size()) ; ++i) {   CSVFupSortedJetIndex.push_back(indexedJetCSVFup[i].first);        }
				for (size_t i = 0 ; (i != indexedJetCSVFdown.size()) ; ++i) {   CSVFdownSortedJetIndex.push_back(indexedJetCSVFdown[i].first);        }
				for (size_t i = 0 ; (i != indexedJetPtJESup.size()) ; ++i) {   PtJESupSortedJetIndex.push_back(indexedJetPtJESup[i].first);        }
				for (size_t i = 0 ; (i != indexedJetPtJESdown.size()) ; ++i) {   PtJESdownSortedJetIndex.push_back(indexedJetPtJESdown[i].first);        }
				for (size_t i = 0 ; (i != indexedJetPtJERup.size()) ; ++i) {   PtJERupSortedJetIndex.push_back(indexedJetPtJERup[i].first);        }
				for (size_t i = 0 ; (i != indexedJetPtJERdown.size()) ; ++i) {   PtJERdownSortedJetIndex.push_back(indexedJetPtJERdown[i].first);        }
				
				firstevent = false;
				
				//choose Higgs Jets with new Algo
				int CSV44mid1 = 99, CSV44mid2 = 99;
				HiggsCandBuilder(&CSV44mid1, &CSV44mid2, indexedJetCSV[0].second, indexedJetCSV[1].second, indexedJetCSV[0].first, indexedJetCSV[1].first, indexedJetPt[0].first, indexedJetPt[1].first);
				if (debug) cout << "HiggsAlgo candidate from jet " << CSV44mid1<< " and " <<CSV44mid2 << endl;
					
				//Jets with CVSup Systematic
				int CSVup44mid1 = 99, CSVup44mid2 = 99;
				HiggsCandBuilder(&CSVup44mid1, &CSVup44mid2, indexedJetCSVup[0].second, indexedJetCSVup[1].second, indexedJetCSVup[0].first, indexedJetCSVup[1].first, indexedJetPt[0].first, indexedJetPt[1].first);
				if (CSV44mid1 != CSVup44mid1) cout << "CVSup changes Higgs jet selection jet0 " << endl;
				if (CSV44mid2 != CSVup44mid2) cout << "CVSup changes Higgs jet selection jet2 " << endl;
				int CSVdown44mid1 = 99, CSVdown44mid2 = 99;
				HiggsCandBuilder(&CSVdown44mid1, &CSVdown44mid2, indexedJetCSVdown[0].second, indexedJetCSVdown[1].second, indexedJetCSVdown[0].first, indexedJetCSVdown[1].first, indexedJetPt[0].first, indexedJetPt[1].first);
				int CSVFup44mid1 = 99, CSVFup44mid2 = 99;
				HiggsCandBuilder(&CSVFup44mid1, &CSVFup44mid2, indexedJetCSVFup[0].second, indexedJetCSVFup[1].second, indexedJetCSVFup[0].first, indexedJetCSVFup[1].first, indexedJetPt[0].first, indexedJetPt[1].first);
				int CSVFdown44mid1 = 99, CSVFdown44mid2 = 99;
				HiggsCandBuilder(&CSVFdown44mid1, &CSVFdown44mid2, indexedJetCSVFdown[0].second, indexedJetCSVFdown[1].second, indexedJetCSVFdown[0].first, indexedJetCSVFdown[1].first, indexedJetPt[0].first, indexedJetPt[1].first);
				int PtJESup44mid1 = 99, PtJESup44mid2 = 99;
				HiggsCandBuilder(&PtJESup44mid1, &PtJESup44mid2, indexedJetCSV[0].second, indexedJetCSV[1].second, indexedJetCSV[0].first, indexedJetCSV[1].first, indexedJetPtJESup[0].first, indexedJetPtJESup[1].first);
				int PtJESdown44mid1 = 99, PtJESdown44mid2 = 99;
				HiggsCandBuilder(&PtJESdown44mid1, &PtJESdown44mid2, indexedJetCSV[0].second, indexedJetCSV[1].second, indexedJetCSV[0].first, indexedJetCSV[1].first, indexedJetPtJESdown[0].first, indexedJetPtJESdown[1].first);
				int PtJERup44mid1 = 99, PtJERup44mid2 = 99;
				HiggsCandBuilder(&PtJERup44mid1, &PtJERup44mid2, indexedJetCSV[0].second, indexedJetCSV[1].second, indexedJetCSV[0].first, indexedJetCSV[1].first, indexedJetPtJERup[0].first, indexedJetPtJERup[1].first);
				int PtJERdown44mid1 = 99, PtJERdown44mid2 = 99;
				HiggsCandBuilder(&PtJERdown44mid1, &PtJERdown44mid2, indexedJetCSV[0].second, indexedJetCSV[1].second, indexedJetCSV[0].first, indexedJetCSV[1].first, indexedJetPtJERdown[0].first, indexedJetPtJERdown[1].first);
				
				
				//Fill Pt for shape variables
				if (debug) for (size_t i = 0 ; (i != indexedPt.size()) ; ++i) {cout << "Unsorted pt of objects: " << indexedPt[i].second << endl; }
				
				indexedPt.push_back(std::pair<size_t,double>(0,(double) jetPt[CSV44mid1]));
				indexedPt.push_back(std::pair<size_t,double>(1,(double) jetPt[CSV44mid2]));
				std::sort(indexedPt.begin(),indexedPt.end(),::IndexedQuantityGreaterThan<double>);
				for (size_t i = 0 ; (i != indexedPt.size()) ; ++i) {   PtSortedIndex.push_back(indexedPt[i].first);        }
				if (debug) for (size_t i = 0 ; (i != indexedPt.size()) ; ++i) {cout << "Sorted pt of objects: " << indexedPt[i].second << endl; }
				//Fill Jets
				
				//if (CSV44mid1 == 1) RegjetPt[0] = sample.hJet_genPtReg1;
				CSV0 = max(CSVNewShapetmp[CSV44mid1],CSVNewShapetmp[CSV44mid2]);
				CSVNewShape[0] = CSVNewShapetmp[CSV44mid1];	
				jetCSV[0] = jetCSVtmp[CSV44mid1];
				jetPt[0] = jetPttmp[CSV44mid1];
				jetEta[0] = jetEtatmp[CSV44mid1];
				jetPhi[0] = jetPhitmp[CSV44mid1];
				jetCHF[0] = jetCHFtmp[CSV44mid1];
				jetPtRaw[0] = jetPtRawtmp[CSV44mid1];
				jetE[0] = jetEtmp[CSV44mid1];
				jetVtx3dL[0] = jetVtx3dLtmp[CSV44mid1];
				jetVtx3deL[0] = jetVtx3deLtmp[CSV44mid1];
				jetVtxPt[0] = jetVtxPttmp[CSV44mid1];
				jetVtxMass[0] = jetVtxMasstmp[CSV44mid1];
				jetPtLeadTrack[0] = jetPtLeadTracktmp[CSV44mid1];
				jetNconstintuents[0] = jetNconstintuentstmp[CSV44mid1];
				jetCEF[0] = jetCEFtmp[CSV44mid1];
				jetNCH[0] = jetNCHtmp[CSV44mid1];
				jetJECUnc[0] = jetJECUnctmp[CSV44mid1];
				jetEt[0] = jetEttmp[CSV44mid1];
				jetMt[0] = jetMttmp[CSV44mid1];
				jetPtRawJER[0] = jetPtRawJERtmp[CSV44mid1];
				jetGenPt[0]=jetGenPttmp[CSV44mid1];
				
				CSV1 = min(CSVNewShapetmp[CSV44mid1],CSVNewShapetmp[CSV44mid2]);
				CSVNewShape[1] = CSVNewShapetmp[CSV44mid2];	
				jetCSV[1] = jetCSVtmp[CSV44mid2];
				jetPt[1] = jetPttmp[CSV44mid2];
				jetEta[1] = jetEtatmp[CSV44mid2];
				jetPhi[1] = jetPhitmp[CSV44mid2];
				jetCHF[1] = jetCHFtmp[CSV44mid2];
				jetPtRaw[1] = jetPtRawtmp[CSV44mid2];
				jetE[1] = jetEtmp[CSV44mid2];
				jetVtx3dL[1] = jetVtx3dLtmp[CSV44mid2];
				jetVtx3deL[1] = jetVtx3deLtmp[CSV44mid2];
				jetVtxPt[1] = jetVtxPttmp[CSV44mid2];
				jetVtxMass[1] = jetVtxMasstmp[CSV44mid2];
				jetPtLeadTrack[1] = jetPtLeadTracktmp[CSV44mid2];
				jetNconstintuents[1] = jetNconstintuentstmp[CSV44mid2];
				jetCEF[1] = jetCEFtmp[CSV44mid2];
				jetNCH[1] = jetNCHtmp[CSV44mid2];
				jetJECUnc[1] = jetJECUnctmp[CSV44mid2];
				jetEt[1] = jetEttmp[CSV44mid2];
				jetMt[1] = jetMttmp[CSV44mid2];
				jetPtRawJER[1] = jetPtRawJERtmp[CSV44mid2];
				jetGenPt[1]=jetGenPttmp[CSV44mid2];
				
				int b = 0;
				for (int kk=0;kk<nJets;kk++){
					if (kk == CSV44mid2) continue;
					if (kk == CSV44mid1) continue;
					if (b<3){
					CSVNewShape[b+2] = CSVNewShapetmp[b];	
					jetCSV[b+2] = jetCSVtmp[b];
					jetPt[b+2] = jetPttmp[b];
					jetEta[b+2] = jetEtatmp[b];
					jetPhi[b+2] = jetPhitmp[b];
					jetCHF[b+2] = jetCHFtmp[b];
					jetPtRaw[b+2] = jetPtRawtmp[b];
					jetE[b+2] = jetEtmp[b];
					jetVtx3dL[b+2] = jetVtx3dLtmp[b];
					jetVtx3deL[b+2] = jetVtx3deLtmp[b];
					jetVtxPt[b+2] = jetVtxPttmp[b];
					jetVtxMass[b+2] = jetVtxMasstmp[b];
					jetPtLeadTrack[b+2] = jetPtLeadTracktmp[b];
					jetNconstintuents[b+2] = jetNconstintuentstmp[b];
					jetCEF[b+2] = jetCEFtmp[b];
					jetNCH[b+2] = jetNCHtmp[b];
					jetJECUnc[b+2] = jetJECUnctmp[b];
					jetEt[b+2] = jetEttmp[b];
					jetMt[b+2] = jetMttmp[b];
					jetPtRawJER[b+2] = jetPtRawJERtmp[b];
					jetGenPt[b+2]=jetGenPttmp[b];				
					CSVup[b+2] =  corrCSV(btagUp, jetCSVtmp[b],jetFlavor[b]);
					CSVdown[b+2] =  corrCSV(btagDown, jetCSVtmp[b],jetFlavor[b]);
					csvFup[b+2] =  corrCSV(btagFUp, jetCSVtmp[b],jetFlavor[b]);
					csvFdown[b+2] =  corrCSV(btagFDown, jetCSVtmp[b],jetFlavor[b]);
						JER_pt_up[b+2] = JERSys(true,jetEta[b+2],jetPt[b+2],jetGenPt[b+2]);
						JER_pt_down[b+2] = JERSys(false,jetEta[b+2],jetPt[b+2],jetGenPt[b+2]);
						JES_pt_up[b+2] = jetPt[b+2]*(1+jetJECUnc[b+2]);
						JES_pt_down[b+2] = jetPt[b+2]*(1-jetJECUnc[b+2]);
						JER_e_up[b+2] = jetE[b+2]*(JER_pt_up[b+2]/jetPt[b+2]);
						JER_e_down[b+2] =  jetE[b+2]*(JER_pt_down[b+2]/jetPt[b+2]);
						JES_e_up[b+2] = jetE[b+2]*(1+jetJECUnc[b+2]);
						JES_e_down[b+2] = jetE[b+2]*(1-jetJECUnc[b+2]);
					}
					b++;
					if(jetPttmp[kk]>30&&abs(jetEtatmp[kk])<2.5) {
						if(MinDphiaJet>deltaPhi(sample.MET_phi,jetPhitmp[kk])) MinDphiaJet =deltaPhi(sample.MET_phi,jetPhitmp[kk]);
					}
				}
				for (int c=2;c<nJets;c++){
					if ( fabs(jetEta[c]) < 2.4 && (jetPt[c] > 20)) Naj++;
					if (CSVNewShape[c] > 0.4 ) Nab++;
				}
				
				
				//btagging systematics
				float FirstCSVup, SecondCSVup, FirstFakeUp, SecondFakeUp;
				float FirstCSVdown, SecondCSVdown, FirstFakedown, SecondFakedown;
				SecondCSVup = corrCSV(btagUp, jetCSVtmp[CSVup44mid2],jetFlavor[CSVup44mid2]);
				SecondCSVdown = corrCSV(btagDown, jetCSVtmp[CSVdown44mid2],jetFlavor[CSVdown44mid2]);
				SecondFakeUp = corrCSV(btagFUp, jetCSVtmp[CSVFup44mid2],jetFlavor[CSVFup44mid2]);
				SecondFakedown = corrCSV(btagFDown, jetCSVtmp[CSVFdown44mid2],jetFlavor[CSVFdown44mid2]);
				FirstCSVup = corrCSV(btagUp, jetCSVtmp[CSVup44mid1],jetFlavor[CSVup44mid1]);
				FirstCSVdown = corrCSV(btagDown, jetCSVtmp[CSVdown44mid1],jetFlavor[CSVdown44mid1]);
				FirstFakeUp = corrCSV(btagFUp, jetCSVtmp[CSVFup44mid1],jetFlavor[CSVFup44mid1]);
				FirstFakedown = corrCSV(btagFDown, jetCSVtmp[CSVFdown44mid1],jetFlavor[CSVFdown44mid1]);
				CSVup[1] = min(FirstCSVup,SecondCSVup);
				CSVup[0] = max(FirstCSVup,SecondCSVup);
				FirstJet.SetPtEtaPhiE(jetPttmp[CSVup44mid1],jetEtatmp[CSVup44mid1],jetPhitmp[CSVup44mid1],jetEtmp[CSVup44mid1]);
				SecondJet.SetPtEtaPhiE(jetPttmp[CSVup44mid2],jetEtatmp[CSVup44mid2],jetPhitmp[CSVup44mid2],jetEtmp[CSVup44mid2]);
				Higgs = FirstJet+SecondJet;
				CSVHmassUp = Higgs.M();				
				csvFup[0] = max(FirstFakeUp,SecondFakeUp);
				csvFup[1] = min(FirstFakeUp,SecondFakeUp);
				FirstJet.SetPtEtaPhiE(jetPttmp[CSVFup44mid1],jetEtatmp[CSVFup44mid1],jetPhitmp[CSVFup44mid1],jetEtmp[CSVFup44mid1]);
				SecondJet.SetPtEtaPhiE(jetPttmp[CSVFup44mid2],jetEtatmp[CSVFup44mid2],jetPhitmp[CSVFup44mid2],jetEtmp[CSVFup44mid2]);
				Higgs = FirstJet+SecondJet;
				CSVHmassFUp = Higgs.M();								
				CSVdown[1] = min(FirstCSVdown,SecondCSVdown);
				CSVdown[0] = max(FirstCSVdown,SecondCSVdown);
				FirstJet.SetPtEtaPhiE(jetPttmp[CSVdown44mid1],jetEtatmp[CSVdown44mid1],jetPhitmp[CSVdown44mid1],jetEtmp[CSVdown44mid1]);
				SecondJet.SetPtEtaPhiE(jetPttmp[CSVdown44mid2],jetEtatmp[CSVdown44mid2],jetPhitmp[CSVdown44mid2],jetEtmp[CSVdown44mid2]);
				Higgs = FirstJet+SecondJet;
				CSVHmassDown = Higgs.M();				
				csvFdown[0] = max(FirstFakedown,SecondFakedown);
				csvFdown[1] = min(FirstFakedown,SecondFakedown);
				FirstJet.SetPtEtaPhiE(jetPttmp[CSVFdown44mid1],jetEtatmp[CSVFdown44mid1],jetPhitmp[CSVFdown44mid1],jetEtmp[CSVFdown44mid1]);
				SecondJet.SetPtEtaPhiE(jetPttmp[CSVFdown44mid2],jetEtatmp[CSVFdown44mid2],jetPhitmp[CSVFdown44mid2],jetEtmp[CSVFdown44mid2]);
				Higgs = FirstJet+SecondJet;
				CSVHmassFDown = Higgs.M();								
								
				//Jet energy scale systematics
				JES_pt_up[0] = jetPttmp[PtJESup44mid1]*(1+jetJECUnctmp[PtJESup44mid1]);
				JES_e_up[0] = jetEtmp[PtJESup44mid1]*(1+jetJECUnctmp[PtJESup44mid1]);
				JES_pt_down[0] = jetPttmp[PtJESdown44mid1]*(1-jetJECUnctmp[PtJESdown44mid1]);
				JES_e_down[0] = jetEtmp[PtJESdown44mid1]*(1-jetJECUnctmp[PtJESdown44mid1]);
				JES_pt_up[1] = jetPttmp[PtJESup44mid2]*(1+jetJECUnctmp[PtJESup44mid2]);
				JES_e_up[1] = jetEtmp[PtJESup44mid2]*(1+jetJECUnctmp[PtJESup44mid2]);
				JES_pt_down[1] = jetPttmp[PtJESdown44mid2]*(1-jetJECUnctmp[PtJESdown44mid2]);
				JES_e_down[1] = jetEtmp[PtJESdown44mid2]*(1-jetJECUnctmp[PtJESdown44mid2]);
				FirstJet.SetPtEtaPhiE(JES_pt_up[0],jetEtatmp[PtJESup44mid1],jetPhitmp[PtJESup44mid1],JES_e_up[0]);
				SecondJet.SetPtEtaPhiE(JES_pt_up[1],jetEtatmp[PtJESup44mid2],jetPhitmp[PtJESup44mid2],JES_e_up[1]);
				Higgs = FirstJet+SecondJet;
				JESHmassUp = Higgs.M();
				FirstJet.SetPtEtaPhiE(JES_pt_down[0],jetEtatmp[PtJESdown44mid1],jetPhitmp[PtJESdown44mid1],JES_e_down[0]);
				SecondJet.SetPtEtaPhiE(JES_pt_down[1],jetEtatmp[PtJESdown44mid2],jetPhitmp[PtJESdown44mid2],JES_e_down[1]);
				Higgs = FirstJet+SecondJet;
				JESHmassDown = Higgs.M();
				
				//jet energy resoltuion systematics
				JER_pt_up[0] = JERSys(true,jetEtatmp[PtJERup44mid1],jetPttmp[PtJERup44mid1],jetGenPttmp[PtJERup44mid1]);
				JER_pt_down[0] = JERSys(false,jetEtatmp[PtJERdown44mid1],jetPttmp[PtJERdown44mid1],jetGenPttmp[PtJERdown44mid1]);
				JER_e_up[0] = jetEtmp[PtJERup44mid1]*(JER_pt_up[0]/jetPttmp[PtJERup44mid1]);
				JER_e_down[0] = jetEtmp[PtJERdown44mid1]*(JER_pt_down[0]/jetPttmp[PtJERdown44mid1]);
				JER_pt_up[1] = JERSys(true,jetEtatmp[PtJERup44mid2],jetPttmp[PtJERup44mid2],jetGenPttmp[PtJERup44mid2]);
				JER_pt_down[1] = JERSys(false,jetEtatmp[PtJERdown44mid2],jetPttmp[PtJERdown44mid2],jetGenPttmp[PtJERdown44mid2]);
				JER_e_up[1] = jetEtmp[PtJERup44mid2]*(JER_pt_up[1]/jetPttmp[PtJERup44mid2]);
				JER_e_down[1] = jetEtmp[PtJERdown44mid2]*(JER_pt_down[1]/jetPttmp[PtJERdown44mid2]);
				FirstJet.SetPtEtaPhiE(JER_pt_up[0],jetEtatmp[PtJERup44mid1],jetPhitmp[PtJERup44mid1],JER_e_up[0]);
				SecondJet.SetPtEtaPhiE(JER_pt_up[1],jetEtatmp[PtJERup44mid2],jetPhitmp[PtJERup44mid2],JER_e_up[1]);
				Higgs = FirstJet+SecondJet;
				JERHmassUp = Higgs.M();
				FirstJet.SetPtEtaPhiE(JER_pt_down[0],jetEtatmp[PtJERdown44mid1],jetPhitmp[PtJERdown44mid1],JER_e_down[0]);
				SecondJet.SetPtEtaPhiE(JER_pt_down[1],jetEtatmp[PtJERdown44mid2],jetPhitmp[PtJERdown44mid2],JER_e_down[1]);
				Higgs = FirstJet+SecondJet;
				JERHmassDown = Higgs.M();
				
				//if (CSV44mid2 == 0) RegjetPt[1] = sample.hJet_genPtReg0;
				//if (CSV44mid2 == 1) RegjetPt[1] = sample.hJet_genPtReg1;					
				
				DetaJJ = jetEta[0]-jetEta[1];	
				FirstJet.SetPtEtaPhiM(jetPt[0],jetEta[0],jetPhi[0],4.2/1000.0);
				SecondJet.SetPtEtaPhiM(jetPt[1],jetEta[1],jetPhi[1],4.2/1000.0);
				Higgs = FirstJet+SecondJet;
				if (debug&&Higgs.M()<20) {
					cout << "pt " << jetPt[0] << " eta " << jetEta[0] << " phi " << jetPhi[0] << " csv " << CSVNewShapetmp[0]<< endl;
					cout << "pt " << jetPt[1] << " eta " << jetEta[1] << " phi " << jetPhi[1] << " csv " << CSVNewShapetmp[1]<< endl;
				}

				
				Hmass = Higgs.M();
				oldHmass = sample.H_mass;
				//RegHmass = sample.newHiggsMass;
				//Hpt = sample.newHiggsPt;
				Hpt = Higgs.Pt();
				ScalarSumHiggsJetPt = jetPt[1] + jetPt[0];
				ScalarSumJetPt = ScalarSumJetPt+ jetPt[1] + jetPt[0];
				Rho25 = sample.rho25;
				MindPhiMEThJetOR30aJet = min(sample.minDeltaPhijetMET,MinDphiaJet);
				MinMET = min(sample.METtype1corr_et,sample.METnoPUCh_et);
				if(debug) cout << "weight of the Jet Histograms: " << weight << endl;
				JetDistributions("allEvts", weight);
				
				if(debug)std::cout << "NJet information filled " << std::endl;
				
				nSV = sample.nSvs;
				nPV = sample.nPVs;
				MET = sample.METtype1corr_et;
				METsig = sample.METtype1corr_sig;
				eventFlavor = sample.eventFlav;
				SV_mass = sample.Sv_massSv[0];
				naJets = sample.naJets;
				
				btag2CSF = sample.btag2CSF;
				if (debug) cout << "halfway through filling histos" << endl;
				
				for (int z = 0; z<sample.nvlep;z++){
					leptonPt[z] = sample.vLepton_pt[z];
					leptonEta[z] = sample.vLepton_eta[z];
					leptonPhi[z] = sample.vLepton_phi[z];
					lep_pfCombRelIso[z] = sample.vLepton_pfCombRelIso[z];
					lep_id95[z] = sample.vLepton_id95[z];
					if (abs(sample.vLepton_type[z] == 13)) nMuons++;
					if (abs(sample.vLepton_type[z] == 11)) nElectrons++;
					nLeptons++;
				}
				for (int zz = 0; zz< sample.nalep;zz++){
					if (zz<3){
						if ((sample.aLepton_pt[zz] == sample.vLepton_pt[0]) || (sample.aLepton_pt[zz] == sample.vLepton_pt[1])){
							if(debug)cout << "Additional Letpon " << zz+sample.nvlep << " pt is the same as vector lepton" << endl;
						} else {
							leptonPt[zz+sample.nvlep] = sample.aLepton_pt[zz];
							nLeptons++;
						}
						leptonPhi[zz+sample.nvlep] = sample.aLepton_phi[zz];
						leptonEta[zz+sample.nvlep] = sample.aLepton_eta[zz];
						lep_pfCombRelIso[zz+sample.nvlep] = sample.aLepton_pfCombRelIso[zz];
					}
					if (abs(sample.aLepton_type[zz] == 13)) nMuons++;
					if (abs(sample.aLepton_type[zz] == 11)) nElectrons++;
					Na_lep++;
					
				}
				if(sample.svmass<999)DitauMass = sample.svmass;
				pZeta25=sample.Pzeta25;
				pZeta45=sample.Pzeta45;
				pZeta65=sample.Pzeta65;
				pZeta85=sample.Pzeta85;
				if (sample.nvlep > 1) {
					Emumass = sample.V_mass;
					Zpt = sample.V_pt;
					DphiZMET = deltaPhi(sample.METtype1corr_phi,sample.V_phi);
				}
				
				LeptonDistributions("allEvts", weight);
				//dPhiHMET = sample.HMETdPhi;
				//DeltaPhiHV = sample.HVdPhi;
				dPhiHMET = deltaPhi(sample.METtype1corr_phi,Higgs.Phi());
				DeltaPhiHV = deltaPhi(sample.V_phi,Higgs.Phi());
				Mt = sample.VMt;
				//DeltaPhijetMETmin = sample.minDeltaPhijetMET; 
				DeltaPhijetMETmin = min(deltaPhi(jetPhi[0],Higgs.Phi()),deltaPhi(jetPhi[1],Higgs.Phi()));
				DeltaPhijetMETZtaumin = sample.minDeltaPhijetMETZtau;
				//Mte = sample.VMte;
				//Mtmu = sample.VMtmu;
				Ht = sample.MHT_mht;
				delPullAngle = sample.deltaPullAngle;
				delPullAngle2 = sample.deltaPullAngle2;
				topMass = sample.top_mass;
				if (debug)cout << "top mass " << topMass << endl;
				topPt = sample.top_pt;
				topWmass = sample.top_wMass;
				
				if( (sample.nvlep > 1) && (sample.nhJets > 1)) {
					RMS_eta = sample.vLepton_eta[1]*sample.vLepton_eta[1]+sample.vLepton_eta[0]*sample.vLepton_eta[0]+jetEta[0]*jetEta[0]+jetEta[1]*jetEta[1];
					RMS_eta = sqrt(RMS_eta/4);
					StandDevEta[0] =sample.vLepton_eta[1];
					StandDevEta[1] =jetEta[0];
					StandDevEta[2] =sample.vLepton_eta[0];
					StandDevEta[3] =jetEta[1];
					EtaStandDev = TMath::RMS(4,StandDevEta);
					AverageEta = sample.vLepton_eta[0]+sample.vLepton_eta[1]+jetEta[0];
					AverageEta = (AverageEta+jetEta[1])/4;
					UnweightedEta = (sample.vLepton_eta[0]-AverageEta)*(sample.vLepton_eta[0]-AverageEta);
					UnweightedEta = UnweightedEta + (sample.vLepton_eta[1]-AverageEta)*(sample.vLepton_eta[1]-AverageEta);
					UnweightedEta = UnweightedEta + (jetEta[0]-AverageEta)*(jetEta[0]-AverageEta);
					UnweightedEta = UnweightedEta + (jetEta[1]-AverageEta)*(jetEta[1]-AverageEta);
					
					PtbalZH = (Hpt-Zpt);
					PtbalMETH = Hpt-MET;
					PtbalZMET = Zpt - MET;
					lep0pt = sample.vLepton_pt[0];
					
					FirstJet.SetPtEtaPhiM(jetPt[0],jetEta[0],jetPhi[0],4.2/1000.0);
					SecondJet.SetPtEtaPhiM(jetPt[1],jetEta[1],jetPhi[1],4.2/1000.0);
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
					
					Zphi = sample.V_phi;
					Hphi = sample.H_phi;
					
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
					if(debug)cout << "Missing Energy via SVD: Electron " << c_svd(0)<< " Muon  " << c_svd(1)<< endl;
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
					if(debug)cout << "Missing Energy via Aaron linear algebra: Electron " << AaronEleMissE << " Muon " << AaronMuMissE << endl;
					
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
				Nweighted_Vtype = Nweighted_Vtype + weight;
				
				float HmassMaxCut = 170, HmassMinCut = 85, CSV0cut = 0.4, DphiZMETcut = 0.85, dphiHVcut = 1.6;
				float ZmassMinCut = 35, ZmassMaxCut = 170, CHFb0cut = 0.2, PzetaCut = 10;
				float MtmuCut =115, MteCut = 95;
				
				float Pzeta = ProjMissT-(0.25*ProjVisT);
				if ((sample.vLepton_pt[1] > 20 || sample.vLepton_pt[0] > 20) && (sample.vLepton_pt[1]>10 && sample.vLepton_pt[0]> 10)){
					//	if (sample.triggerFlags[51]||sample.triggerFlags[52]||sample.triggerFlags[53]||sample.triggerFlags[54]){
					Ntrigger++;
					Nweighted_trigger = Nweighted_trigger + weight;
					if (debug) cout << "PUweight2012 inside trigger requirement " << PUweight2012 << endl;
					JetDistributions("HLT", weight);
					LeptonDistributions("HLT", weight);
					TH2FDistributions("HLT", weight);
					EventShapeDistributions("HLT", weight);
					EventDistributions("HLT", weight);
					JetDistributions("NoWeightHLT", LumiWeight);
					LeptonDistributions("NoWeightHLT", LumiWeight);
					TH2FDistributions("NoWeightHLT", LumiWeight);
					EventShapeDistributions("NoWeightHLT", LumiWeight);
					EventDistributions("NoWeightHLT", LumiWeight);
					JetDistributions("PUWeightHLT", LumiWeight*PUweight2012);
					LeptonDistributions("PUWeightHLT", LumiWeight*PUweight2012);
					TH2FDistributions("PUWeightHLT", LumiWeight*PUweight2012);
					EventShapeDistributions("PUWeightHLT", LumiWeight*PUweight2012);
					EventDistributions("PUWeightHLT", LumiWeight*PUweight2012);	
					JetDistributions("TrigWeightHLT", LumiWeight*EleTrigWeight*MuonTrigWeight);
					LeptonDistributions("TrigWeightHLT", LumiWeight*EleTrigWeight*MuonTrigWeight);
					TH2FDistributions("TrigWeightHLT", LumiWeight*EleTrigWeight*MuonTrigWeight);
					EventShapeDistributions("TrigWeightHLT", LumiWeight*EleTrigWeight*MuonTrigWeight);
					EventDistributions("TrigWeightHLT", LumiWeight*EleTrigWeight*MuonTrigWeight);
					JetDistributions("IDWeightHLT", LumiWeight*MuIDweight*WP95weight);
					LeptonDistributions("IDWeightHLT", LumiWeight*MuIDweight*WP95weight);
					TH2FDistributions("IDWeightHLT", LumiWeight*MuIDweight*WP95weight);
					EventShapeDistributions("IDWeightHLT", LumiWeight*MuIDweight*WP95weight);
					EventDistributions("IDWeightHLT", LumiWeight*MuIDweight*WP95weight);					
					if ( jetPt[0] > 20 && jetPt[1] > 20  &&
						fabs(sample.vLepton_eta[0]) < 2.5 && fabs(sample.vLepton_eta[1]) < 2.4 && fabs(jetEta[0]) < 2.5 &&
						fabs(jetEta[1]) < 2.5 && sample.hJet_id[0]==1 && sample.hJet_id[1]==1 && sample.hbhe){
						Npreselect++;
						Nweighted_preselect = Nweighted_preselect + weight;
						isdata = false;
						JetDistributions("PreSelect", weight);
						LeptonDistributions("PreSelect", weight);
						TH2FDistributions("PreSelect", weight);
						EventShapeDistributions("PreSelect", weight);
						EventDistributions("PreSelect", weight);
						if (isDATA) isdata = true;
						if (fabs(DphiZMET) < DphiZMETcut && (oldHmass>=HmassMinCut)&&(oldHmass<=HmassMaxCut) && CSV0>CSV0cut && DeltaPhiHV > dphiHVcut && Pzeta > PzetaCut && Emumass< 68 && DitauMass > ZmassMinCut && DitauMass < ZmassMaxCut && jetCHF[0]>CHFb0cut) hdelRemu_NotThisCut.Fill(delRemu, weight);
						if(delRemu>0.3){
							N_EfakeCuts++;
							Nweighted_EfakeCuts = Nweighted_EfakeCuts + weight;
							isdata = false;
							JetDistributions("delRemu", weight);
							LeptonDistributions("delRemu", weight);
							TH2FDistributions("delRemu", weight);
							EventShapeDistributions("delRemu", weight);
							EventDistributions("delRemu", weight);
							JetDistributions("NoWeightdelRemu", LumiWeight);
							LeptonDistributions("NoWeightdelRemu", LumiWeight);
							TH2FDistributions("NoWeightdelRemu", LumiWeight);
							EventShapeDistributions("NoWeightdelRemu", LumiWeight);
							EventDistributions("NoWeightdelRemu", LumiWeight);
							JetDistributions("PUWeightdelRemu", LumiWeight*PUweight2012);
							LeptonDistributions("PUWeightdelRemu", LumiWeight*PUweight2012);
							TH2FDistributions("PUWeightdelRemu", LumiWeight*PUweight2012);
							EventShapeDistributions("PUWeightdelRemu", LumiWeight*PUweight2012);
							EventDistributions("PUWeightdelRemu", LumiWeight*PUweight2012);	
							JetDistributions("TrigWeightdelRemu", LumiWeight*EleTrigWeight*MuonTrigWeight);
							LeptonDistributions("TrigWeightdelRemu", LumiWeight*EleTrigWeight*MuonTrigWeight);
							TH2FDistributions("TrigWeightdelRemu", LumiWeight*EleTrigWeight*MuonTrigWeight);
							EventShapeDistributions("TrigWeightdelRemu", LumiWeight*EleTrigWeight*MuonTrigWeight);
							EventDistributions("TrigWeightdelRemu", LumiWeight*EleTrigWeight*MuonTrigWeight);
							JetDistributions("IDWeightdelRemu", LumiWeight*MuIDweight*WP95weight);
							LeptonDistributions("IDWeightdelRemu", LumiWeight*MuIDweight*WP95weight);
							TH2FDistributions("IDWeightdelRemu", LumiWeight*MuIDweight*WP95weight);
							EventShapeDistributions("IDWeightdelRemu", LumiWeight*MuIDweight*WP95weight);
							EventDistributions("IDWeightdelRemu", LumiWeight*MuIDweight*WP95weight);												
							FOM_tree->Fill();
							if (fabs(DphiZMET) < DphiZMETcut && (oldHmass>=HmassMinCut)&&(oldHmass<=HmassMaxCut) && CSV0>CSV0cut && DeltaPhiHV > dphiHVcut && Pzeta > PzetaCut && jetCHF[0]>CHFb0cut) hMemu_NotThisCut.Fill(DitauMass, weight);
//								if (fabs(DphiZMET) < DphiZMETcut && fabs(DphiSecondMET) < 1.5 && (oldHmass>HmassMinCut)&&(oldHmass<=HmassMaxCut) && DeltaPhiHV > dphiHVcut && CSV0<CSV0cut && Naj < 2 && nSV ==0){	
							if (CSV0<CSV0cut && ProjMissT>-50 && DeltaPhiHV > dphiHVcut && fabs(DphiZMET) < DphiZMETcut ){	
									N_LFCR++;
									Nweighted_LFCR = Nweighted_LFCR + (LumiWeight*PUweight2012*Trigweight);
									isdata = false;
									JetDistributions("LFCR", weight);
									LeptonDistributions("LFCR", weight);
									TH2FDistributions("LFCR", weight);
									EventShapeDistributions("LFCR", weight);
									EventDistributions("LFCR", weight);
									SFUnc_LFCR = SFUnc_LFCR + (LumiWeight*PUweight2012*Trigweight)*(LumiWeight*PUweight2012*Trigweight);
									JetDistributions("NoWeightLFCR", LumiWeight);
									LeptonDistributions("NoWeightLFCR", LumiWeight);
									TH2FDistributions("NoWeightLFCR", LumiWeight);
									EventShapeDistributions("NoWeightLFCR", LumiWeight);
									EventDistributions("NoWeightLFCR", LumiWeight);	
									JetDistributions("PUWeightLFCR", LumiWeight*PUweight2012);
									LeptonDistributions("PUWeightLFCR", LumiWeight*PUweight2012);
									TH2FDistributions("PUWeightLFCR", LumiWeight*PUweight2012);
									EventShapeDistributions("PUWeightLFCR", LumiWeight*PUweight2012);
									EventDistributions("PUWeightLFCR", LumiWeight*PUweight2012);	
									JetDistributions("TrigWeightLFCR", LumiWeight*EleTrigWeight*MuonTrigWeight);
									LeptonDistributions("TrigWeightLFCR", LumiWeight*EleTrigWeight*MuonTrigWeight);
									TH2FDistributions("TrigWeightLFCR", LumiWeight*EleTrigWeight*MuonTrigWeight);
									EventShapeDistributions("TrigWeightLFCR", LumiWeight*EleTrigWeight*MuonTrigWeight);
									EventDistributions("TrigWeightLFCR", LumiWeight*EleTrigWeight*MuonTrigWeight);
									JetDistributions("IDWeightLFCR", LumiWeight*MuIDweight*WP95weight);
									LeptonDistributions("IDWeightLFCR", LumiWeight*MuIDweight*WP95weight);
									TH2FDistributions("IDWeightLFCR", LumiWeight*MuIDweight*WP95weight);
									EventShapeDistributions("IDWeightLFCR", LumiWeight*MuIDweight*WP95weight);
									EventDistributions("IDWeightLFCR", LumiWeight*MuIDweight*WP95weight);									
								}// LF Control Region	
								if (fabs(DphiZMET) < DphiZMETcut && (oldHmass>=HmassMinCut)&&(oldHmass<=HmassMaxCut) && DeltaPhiHV > dphiHVcut && jetCHF[0]>CHFb0cut && Pzeta > PzetaCut) hCSV0_NotThisCut.Fill(CSV0, weight);							
								if(CSV0>CSV0cut){
									N_CSV0++;
									Nweighted_CSV0 = Nweighted_CSV0 + weight;
									if(isDATA)isdata=true;
									JetDistributions("CSV0", weight);
									LeptonDistributions("CSV0", weight);
									TH2FDistributions("CSV0", weight);
									EventShapeDistributions("CSV0", weight);
									EventDistributions("CSV0", weight);
									TMVA_tree->Fill();							
									//if (fabs(DphiZMET) > DphiZMETcut && (oldHmass>HmassMinCut)&&(oldHmass<=HmassMaxCut) && DeltaPhiHV > dphiHVcut && jetCHF[0]>CHFb0cut ){
									if (ProjMissT<-50 && DeltaPhiHV > dphiHVcut ){
										if (CSV1 > 0.5 ){
											N_TopCR++;
											Nweighted_TopCR = Nweighted_TopCR + (LumiWeight*PUweight2012*Trigweight);
											isdata = false;
											JetDistributions("TopCR", weight);
											LeptonDistributions("TopCR", weight);
											TH2FDistributions("TopCR", weight);
											EventShapeDistributions("TopCR", weight);
											EventDistributions("TopCR", weight);
											SFUnc_TOPCR = SFUnc_TOPCR + (LumiWeight*PUweight2012*Trigweight)*(LumiWeight*PUweight2012*Trigweight);
											JetDistributions("NoWeightTopCR", LumiWeight);
											LeptonDistributions("NoWeightTopCR", LumiWeight);
											TH2FDistributions("NoWeightTopCR", LumiWeight);
											EventShapeDistributions("NoWeightTopCR", LumiWeight);
											EventDistributions("NoWeightTopCR", LumiWeight);									
											JetDistributions("PUWeightTopCR", LumiWeight*PUweight2012);
											LeptonDistributions("PUWeightTopCR", LumiWeight*PUweight2012);
											TH2FDistributions("PUWeightTopCR", LumiWeight*PUweight2012);
											EventShapeDistributions("PUWeightTopCR", LumiWeight*PUweight2012);
											EventDistributions("PUWeightTopCR", LumiWeight*PUweight2012);	
											JetDistributions("TrigWeightTopCR", LumiWeight*EleTrigWeight*MuonTrigWeight);
											LeptonDistributions("TrigWeightTopCR", LumiWeight*EleTrigWeight*MuonTrigWeight);
											TH2FDistributions("TrigWeightTopCR", LumiWeight*EleTrigWeight*MuonTrigWeight);
											EventShapeDistributions("TrigWeightTopCR", LumiWeight*EleTrigWeight*MuonTrigWeight);
											EventDistributions("TrigWeightTopCR", LumiWeight*EleTrigWeight*MuonTrigWeight);	
											JetDistributions("IDWeightTopCR", LumiWeight*MuIDweight*WP95weight);
											LeptonDistributions("IDWeightTopCR", LumiWeight*MuIDweight*WP95weight);
											TH2FDistributions("IDWeightTopCR", LumiWeight*MuIDweight*WP95weight);
											EventShapeDistributions("IDWeightTopCR", LumiWeight*MuIDweight*WP95weight);
											EventDistributions("IDWeightTopCR", LumiWeight*MuIDweight*WP95weight);	
										}//Top CR
										if (debug)std::cout << "Looking for bug in event " << event << std::endl;
										if (CSV1 < 0.5){
											N_SingleTopCR++;
											Nweighted_SingleTopCR = Nweighted_SingleTopCR + (LumiWeight*PUweight2012*Trigweight);
											isdata = false;
											JetDistributions("SingleTopCR", weight);
											LeptonDistributions("SingleTopCR", weight);
											TH2FDistributions("SingleTopCR", weight);
											EventShapeDistributions("SingleTopCR", weight);
											EventDistributions("SingleTopCR", weight);
											SFUnc_SingleTOP = SFUnc_SingleTOP + (LumiWeight*PUweight2012*Trigweight)*(LumiWeight*PUweight2012*Trigweight);
										}//Single Top CR
									}// Top and Single Top Orthogonal to Signal, reverse DphiZMET cut	
									if (Emumass< 68 && DitauMass > ZmassMinCut && DitauMass < ZmassMaxCut && oldHmass<300 ) { 
										NMemu++;
										Nweighted_Memu = Nweighted_Memu + weight;
										if(isDATA)isdata=true;
										JetDistributions("MemuCut", weight);
										LeptonDistributions("MemuCut", weight);
										TH2FDistributions("MemuCut", weight);
										EventShapeDistributions("MemuCut", weight);
										EventDistributions("MemuCut", weight); 										
									if ((oldHmass>=HmassMinCut)&&(oldHmass<=HmassMaxCut) && DeltaPhiHV > dphiHVcut && Pzeta > PzetaCut && jetCHF[0]>CHFb0cut) hDphiZMET_NotThisCut.Fill(DphiZMET, weight);
									if (fabs(DphiZMET) < DphiZMETcut){
										N_DphiZMET++;
										Nweighted_DphiZMET = Nweighted_DphiZMET + weight;
										if(isDATA)isdata=true;
										JetDistributions("DphiZMET", weight);
										LeptonDistributions("DphiZMET", weight);
										TH2FDistributions("DphiZMET", weight);
										EventShapeDistributions("DphiZMET", weight);
										EventDistributions("DphiZMET", weight);
										if ((oldHmass>=HmassMinCut)&&(oldHmass<=HmassMaxCut) && Pzeta > PzetaCut && jetCHF[0]>CHFb0cut) hDeltaPhiHV_NotThisCut.Fill(DeltaPhiHV, weight);
										if (DeltaPhiHV > dphiHVcut){
											N_DeltaPhiHV++;
											Nweighted_DeltaPhiHV = Nweighted_DeltaPhiHV + weight;
											if(isDATA)isdata=true;
											JetDistributions("DeltaPhiHV", weight);
											LeptonDistributions("DeltaPhiHV", weight);
											TH2FDistributions("DeltaPhiHV", weight);
											EventShapeDistributions("DeltaPhiHV", weight);
											EventDistributions("DeltaPhiHV", weight);
											if ((oldHmass>=HmassMinCut)&&(oldHmass<=HmassMaxCut) && jetCHF[0]>CHFb0cut) hPzeta_NotThisCut.Fill(Pzeta, weight);
											if (Pzeta > PzetaCut){
												N_Pzeta++;
												Nweighted_Pzeta = Nweighted_Pzeta + weight;
												if(isDATA)isdata=true;
												JetDistributions("Pzeta", weight);
												LeptonDistributions("Pzeta", weight);
												TH2FDistributions("Pzeta", weight);
												EventShapeDistributions("Pzeta", weight);
												EventDistributions("Pzeta", weight);												
												if (Naj<2) hMjj_NotThisCut.Fill(oldHmass, weight);
												if((oldHmass>=HmassMinCut)&&(oldHmass<=HmassMaxCut)){
													N_Mjj++;
													Nweighted_Mjj = Nweighted_Mjj + weight; 
													if(isDATA)isdata=true;
													JetDistributions("Mjj", weight);
													LeptonDistributions("Mjj", weight);
													TH2FDistributions("Mjj", weight);
													EventShapeDistributions("Mjj", weight);
													EventDistributions("Mjj", weight);
													if (nJets<4){
														//Naj
														N_Naj++;
														Nweighted_Naj = Nweighted_Naj + weight;
														if(isDATA)isdata=true;
														JetDistributions("Naj", weight);
														LeptonDistributions("Naj", weight);
														TH2FDistributions("Naj", weight);
														EventShapeDistributions("Naj", weight);
														EventDistributions("Naj", weight);
													}	//Naj
													if (jetCHF[0]>CHFb0cut){
														N_jetCHF0++;
														Nweighted_jetCHF0 = Nweighted_jetCHF0 + weight;
														if(isDATA)isdata=true;
														JetDistributions("CHF0", weight);
														LeptonDistributions("CHF0", weight);
														TH2FDistributions("CHF0", weight);
														EventShapeDistributions("CHF0", weight);
														EventDistributions("CHF0", weight);
														if (Mte < MteCut && Mtmu < MtmuCut){
															if(isDATA)isdata=true;
															N_Mt++;
															Nweighted_Mt = Nweighted_Mt + weight;
															JetDistributions("Mt", weight);
															LeptonDistributions("Mt", weight);
															TH2FDistributions("Mt", weight);
															EventShapeDistributions("Mt", weight);
															EventDistributions("Mt", weight);
															/*Ntree++;
															 if (Ntree%2){
															 TMVA_tree->Fill();
															 NTMVAtree = NTMVAtree + 1;
															 }
															 BDT_tree->Fill();*/
														}	//Mte and Mtmu
													}	//CHF0
												}//Mjj							
											} //PzetaCut
										}//DeltaPhiHV
									}// Delta Phi Z, MET (Z is emu vectors only)
								}//CSV
							}//emu mass requirement
							//if ((oldHmass<HmassMinCut || oldHmass>HmassMaxCut) && oldHmass < 300 && Mte < MteCut && Mtmu < MtmuCut && jetCHF[0]>CHFb0cut){
							if (CSV0>CSV0cut && ProjMissT>-50 ){	
								BDT_tree->Fill();
								//if (Pzeta> -20 && ProjMissT > -20 && CSV1 < 0.679 && jetCEF[1]>0.5 && fabs(Dphiemu)>1.0){
									N_HFCR++;
									Nweighted_HFCR = Nweighted_HFCR + (LumiWeight*PUweight2012*Trigweight);
									isdata = false;
									JetDistributions("HFCR", weight);
									LeptonDistributions("HFCR", weight);
									TH2FDistributions("HFCR", weight);
									EventShapeDistributions("HFCR", weight);
									EventDistributions("HFCR", weight);	
									SFUnc_HFCR = SFUnc_HFCR + (LumiWeight*PUweight2012*Trigweight)*(LumiWeight*PUweight2012*Trigweight);
									if(isDATA)isdata=true;
								//}// Purity cuts									
							}// HF Control Region																
						}//EleFakeCuts
					} else {
						FailedJetID = FailedJetID + 1;
					}//Jet ID and eta requirement
				}//end trigger emulation
				isdata = false;
				EventDistributions("allEvts", weight);
			}//end requirement Zemu event
			if (debug) std::cout << "Vtype or json requirement " << event << std::endl;
			
		}//end isM50sample genZpt cut
		if(isDATA)isdata=true;
		if (debug) std::cout << "Zpt and Flavor isDATA is " << isDATA << std::endl;
		
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
	std::cout << "CHF0: " << N_jetCHF0 << endl;
	std::cout << "Mte&Mtmu: " << N_Mt << endl;
	std::cout << "Naj: " << N_Naj << endl;
	
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
	
	ofstream myfile;
	myfile.open(TString::Format("%sDataEntry.txt",directory.c_str()).Data());	
	myfile<<TString::Format("\t %s",directory.c_str()).Data()<<endl;
	myfile<<TString::Format("Number of events \t %0.5f \t %0.5f",(double)event,event_weighted).Data()<<endl;
	myfile<<TString::Format("Vtype 5 \t %0.5f \t %0.5f",(double)N_Vtype,Nweighted_Vtype).Data()<<endl;
	myfile<<TString::Format("emu trigger \t %0.5f \t %0.5f",Ntrigger,Nweighted_trigger).Data()<<endl;
	myfile<<TString::Format("PreSelection \t %0.5f \t %0.5f",Npreselect,Nweighted_preselect).Data()<<endl;
	myfile<<TString::Format("EleFakeCuts \t %0.5f \t %0.5f",N_EfakeCuts,Nweighted_EfakeCuts).Data()<<endl;
	myfile<<TString::Format("Memu \t %0.5f \t %0.5f",NMemu,Nweighted_Memu).Data()<<endl;
	myfile<<TString::Format("CSV0 \t %0.5f \t %0.5f",N_CSV0,Nweighted_CSV0).Data()<<endl;
	myfile<<TString::Format("DphiZMET \t %0.5f \t %0.5f",N_DphiZMET,Nweighted_DphiZMET).Data()<<endl;
	myfile<<TString::Format("DphiZH \t %0.5f \t %0.5f",N_DeltaPhiHV,Nweighted_DeltaPhiHV).Data()<<endl;
	myfile<<TString::Format("Pzeta \t %0.5f \t %0.5f",N_Pzeta,Nweighted_Pzeta).Data()<<endl;
	myfile<<TString::Format("Mjj \t %0.5f \t %0.5f",N_Mjj,Nweighted_Mjj).Data()<<endl;
	myfile<<TString::Format("CHF0 \t %0.5f \t %0.5f",N_jetCHF0,Nweighted_jetCHF0).Data()<<endl;
	myfile<<TString::Format("Mte&Mtmu \t %0.5f \t %0.5f",N_Mt,Nweighted_Mt).Data()<<endl;
	myfile<<TString::Format("Naj \t %0.5f \t %0.5f",N_Naj,Nweighted_Naj).Data()<<endl;
	myfile<<endl;
	myfile<<TString::Format("Top CR \t %0.5f \t %0.5f \t %0.5f \t %0.5f",N_TopCR,Nweighted_TopCR,SFUnc_TOPCR*SFUnc_TOPCR,SFUnc_TOPCR).Data()<<endl;
	myfile<<TString::Format("Single Top CR \t %0.5f \t %0.5f \t %0.5f \t %0.5f",N_SingleTopCR,Nweighted_SingleTopCR,
							SFUnc_SingleTOP*SFUnc_SingleTOP,SFUnc_SingleTOP).Data()<<endl;
	myfile<<TString::Format("LF CR \t %0.5f \t %0.5f \t %0.5f \t %0.5f",N_LFCR,Nweighted_LFCR,SFUnc_LFCR*SFUnc_LFCR,SFUnc_LFCR).Data()<<endl;
	myfile<<TString::Format("HF CR \t %0.5f \t %0.5f \t %0.5f \t %0.5f",N_HFCR,Nweighted_HFCR,SFUnc_HFCR*SFUnc_HFCR,SFUnc_HFCR).Data()<<endl;
	myfile<<endl;
	myfile<<endl;
	myfile<<TString::Format("Events \t %s \t %0.5f \t %0.5f \t %0.5f \t %0.5f",directory.c_str(),Nweighted_TopCR,Nweighted_SingleTopCR,
							Nweighted_LFCR,Nweighted_HFCR).Data()<<endl;
	myfile<<TString::Format("Errors \t %s \t %0.5f \t %0.5f \t %0.5f \t %0.5f",directory.c_str(),SFUnc_TOPCR,
							SFUnc_SingleTOP,SFUnc_LFCR,SFUnc_HFCR).Data()<<endl;
	myfile<<TString::Format("Table \t %s \t %0.2fat%0.2f \t %0.2fat%0.2f \t %0.2fat%0.2f \t %0.2fat%0.2f",
							directory.c_str(),Nweighted_TopCR,SFUnc_TOPCR,Nweighted_SingleTopCR,SFUnc_SingleTOP,Nweighted_LFCR,SFUnc_LFCR,Nweighted_HFCR,SFUnc_HFCR).Data()<<endl;
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
	
	TMVA_tree->Write();
	BDT_tree->Write();
	FOM_tree->Write();
	
	////////
	// write histos to file:
	////////
	for (map<std::string,TH1*>::iterator it=histmap.begin(); it!=histmap.end();it++) {
		(*it).second->Write();
		c1.Clear();
		(*it).second->Draw();
		c1.Print((directory+"/"+(*it).first+suffixps).c_str());
		if (debug) cout << "This is the histmap string: " << (*it).first << endl;
	}
	for (map<string,TH2*>::iterator it2=bidimhistmap.begin(); it2!=bidimhistmap.end();it2++) {
		(*it2).second->Write();
		c1.Clear();
		(*it2).second->Draw();
		//	c1.Print((directory+"/"+(*it2).first+suffixps).c_str());
		delete (*it2).second;
	}
	
	cout << "macro done, now clean up" << endl;
	// Write and Close the output file.
	//ofile.Write();
	
	ofile.Close();
	
	/*	MuonTrigWeightFilelow->Close();
	 EleTrigWeightFilelow->Close();
	 MuonTrigWeightFilehigh->Close();
	 EleTrigWeightFilehigh->Close();
	 MuonIDWeightFile->Close();
	 EleIDWeightFile->Close();
	 lowptMuonIDWeightFile->Close();
	 lowptEleIDWeightFile->Close();
	 cout << "files closed" << endl;
	 */
	cout << "all I have to do is return 0" << endl;	
	return 0;
}


double SetWeight( std::string filename){
	double SampleWeight = 1.0;
	if (findString(filename, "ZH_ZToLL_125HToBB")){ 
		SampleWeight = lumi/lumiZH125;
	cout << "found ZH_ZToLL_HToBB_M-125 string" << endl;}
	if (findString(filename, "ZH_ZToLL_115HToBB")){ SampleWeight = lumi/lumiZH115;}
	if (findString(filename, "ZH_ZToLL_135HToBB")){ SampleWeight = lumi/lumiZH135;}
	if (findString(filename, "ZH_ZToLL_110HToBB")){ SampleWeight = lumi/lumiZH110;}
	if (findString(filename, "ZH_ZToLL_130HToBB")){ SampleWeight = lumi/lumiZH130;}
	if (findString(filename, "ZH_ZToLL_120HToBB")){ SampleWeight = lumi/lumiZH120;}
	if (findString(filename, "DYPtZ")){ SampleWeight = lumi/lumiZJH;}
	if (findString(filename, "DYM50")){ SampleWeight = lumi/lumiZJL;}
	if (findString(filename, "120to170")){ SampleWeight = 1.0;}
	if (findString(filename, "170to300")){ SampleWeight = 1.0;}
	if (findString(filename, "300to470")){ SampleWeight = 1.0;}
	if (findString(filename, "470to600")){ SampleWeight = 1.0;}
	if (findString(filename, "80to120")){ SampleWeight = 1.0;}
	if (findString(filename, "TT")){ SampleWeight = lumi/(lumiTT);}
	if (findString(filename, "T_TuneZ2_s")){ SampleWeight = lumi/lumiTs;}
	if (findString(filename, "T_TuneZ2_t-channel")){ SampleWeight = lumi/lumiTt;}
	if (findString(filename, "T_tW")){ SampleWeight = lumi/lumiTtW;}
	if (findString(filename, "Tbar_TuneZ2_s")){ SampleWeight = lumi/lumiTsb;}
	if (findString(filename, "Tbar_TuneZ2_t-channel")){ SampleWeight = lumi/lumiTtb;}
	if (findString(filename, "Tbar_tW")){ SampleWeight = lumi/lumiTtWb;}
	if (findString(filename, "WJetsToLNuPtW100")){ SampleWeight = lumi/lumiWJ100Mad;}
	if (findString(filename, "WJets.root")){ SampleWeight = lumi/lumiWJ;}
	if (findString(filename, "WW")){ SampleWeight = lumi/lumiWW;}
	if (findString(filename, "WZ")){ SampleWeight = lumi/lumiWZ;}
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
		
	}
	
	string mjj = "hMjj_" + cut;
	//string mjjReg = "hMjjReg_" + cut;
	string mjjOld = "hMjjOld_" + cut;
	string mjjAC = "hMjjAC_" + cut;
	string ptjj = "hPtjj_" +cut;
	string detajj = "hdetaJJ_" +cut;
	string scalarSumHiggsJetPt = "hScalarSumHiggsJetPt_" +cut;
	string scalarSumJetPt = "hScalarSumJetPt_" +cut;
	//string Regressedjetpt0 = "hRegJetPt0_"+cut;
	//string Regressedjetpt1 = "hRegJetPt1_"+cut;
	string rho25 = "hRho25_"+cut;
	string hminDeltaPhijetMET = "hminDeltaPhijetMET_" +cut;
	string hminMET = "hMinMET_"+cut;
	string hminDeltaPhijetMETOR30aJET = "hMindPhiMEThJetOR30aJet_"+cut;
	//string RegRes = "hRegressionResolutionMjj_"+cut;
	string RecoRes = "hRecoResolutionMjj_"+cut;
	string Recob0Res = "hRecoResolutionb0_"+cut;
	string Recob1Res = "hRecoResolutionb1_"+cut;
	//string Regb0Res = "hRegResolutionb0_"+cut;
	//string Regb1Res = "hRegResolutionb1_"+cut;
	
	//if (!isdata || (Hmass < 90 || Hmass > 150)) fillhisto(mjjReg, RegHmass, ph_weight, "Regressed Mass", 30, 0, 300);
	if (!isdata || (oldHmass < 90 || oldHmass > 150)) fillhisto(mjjOld, oldHmass, ph_weight, "Old Higgs Mass", 30, 0, 300);
	if (!isdata || (Hmass < 90 || Hmass > 150)) fillhisto(mjj, Hmass, ph_weight, "Invariant Mass of two Jets", 30, 0, 300);
	if (!isdata || (Hmass < 90 || Hmass > 150)) fillhisto(mjjAC, Hmass, ph_weight, "Invariant Mass of two Jets", 20, 75, 175);
	fillhisto(ptjj, Hpt, ph_weight, "Pt of two b jets with highest CSV", 25, 0, 250);
	fillhisto(detajj, fabs(DetaJJ), ph_weight, "Delta eta between two jets", 10, 0, 5);
	fillhisto(scalarSumHiggsJetPt, ScalarSumHiggsJetPt, ph_weight, "scalar sum higgs jet pt", 25, 25, 375);
	fillhisto(scalarSumJetPt, ScalarSumJetPt, ph_weight, "scalar sum alljet pt", 25, 30, 300);
	//fillhisto(Regressedjetpt0, RegjetPt[0], ph_weight, "Regressed jet0 pT", 10, 0.0, 200);
	//fillhisto(Regressedjetpt1, RegjetPt[1], ph_weight, "Regressed jet1 pT", 10, 0.0, 200);
	fillhisto(rho25, Rho25, ph_weight, "Energy density eta<2.5", 10, 0.0, 20);	
	fillhisto(hminDeltaPhijetMET,DeltaPhijetMETmin, ph_weight, "Delta phi between MET and nearest jet", 16, 0, 3.14159265);
	fillhisto(hminMET, MinMET, ph_weight, "min(METtype1corr,METnoPUCh)",		13, 0.0, 260);
	fillhisto(hminDeltaPhijetMETOR30aJET,MindPhiMEThJetOR30aJet, ph_weight, "min(DeltaPhijetMET,Dphiajetpt30MET)", 16, 0, 3.14159265);
	//if (!isdata) fillhisto(RegRes, GenHiggsMass-RegHmass, ph_weight, "Resoltuion Higgs Regressed", 100, -100, 100);
	if (!isdata) fillhisto(RecoRes, GenHiggsMass-oldHmass, ph_weight, "Resolution Higgs Reco", 100, -100, 100);
	//if (!isdata) fillhisto(Regb0Res, jetGenPt[0]-RegjetPt[0], ph_weight, "Resoltuion jet 0 Regressed", 100, -100, 100);
	if (!isdata) fillhisto(Recob0Res, jetGenPt[0]-jetPt[0], ph_weight, "Resolution jet 0 Reco", 100, -100, 100);
	//if (!isdata) fillhisto(Regb1Res, jetGenPt[1]-RegjetPt[1], ph_weight, "Resoltuion jet 1 Regressed", 100, -100, 100);
	if (!isdata) fillhisto(Recob1Res, jetGenPt[1]-jetPt[1], ph_weight, "Resolution jet 1 Reco", 100, -100, 100);
}// JetDistributions


void LeptonDistributions(string cut, double ph_weight){
	
	for (int i=0; i != nLeptons && i < 10; ++i) {
		string leptonpT = Form("hPtlep%i_", i) + cut;
		string leptoneta = Form("hEtalep%i_", i) + cut;
		string leptonphi = Form("hPhilep%i_", i) + cut;
		string leppfCombRelIso = Form("hPFRelIso%i_", i) + cut;
		
		fillhisto(leptonpT, leptonPt[i], ph_weight, "Lepton pT", 10, 0.0, 100);
		fillhisto(leptoneta, leptonEta[i], ph_weight, "Lepton #eta", 13, -3, 3.5);
		fillhisto(leptonphi, leptonPhi[i], ph_weight, "Lepton #phi", 20, -3.14159265, 4.7123889);
		if (lep_pfCombRelIso[i] > 0 ) fillhisto(leppfCombRelIso, lep_pfCombRelIso[i], ph_weight, "PF Rel Iso of Lepton", 25, 0, 0.5);
		
	}
	
	string memu = "hMemu_" + cut;
	string ditaumass = "hDiTauMass_" + cut;
	string memuAC = "hMemuAC_" + cut;
	string ditaumassData = "hDiTauMassData_" + cut;
	string ptemu = "hPtemu_" +cut;
	string LeadLepPt = "hLeadLepPt_" +cut;
	string SecLepPt = "hSecLepPt_" +cut;
	string DiTauRes = "hResolutionDiTau_"+cut;
	string VisRes = "hResolutionVis_"+cut;
	
	fillhisto(memu, Emumass, ph_weight, "Invariant Mass of two Leptons", 50, 0, 150);
	fillhisto(ditaumass, DitauMass, ph_weight, "Di Tau Mass", 50, 0, 150);
	fillhisto(memuAC, Emumass, ph_weight, "Invariant Mass of two Leptons ", 30, 0, 150);
	fillhisto(ditaumassData, DitauMass, ph_weight, "Di Tau Mass", 30, 0, 150);
	fillhisto(ptemu, Zpt, ph_weight, "Pt of two Leptons", 10, 0, 150);
	if(leptonPt[0]> leptonPt[1]){
		fillhisto(LeadLepPt, leptonPt[0], ph_weight, "Leading Lepton Pt", 10, 0.0, 100);
		fillhisto(SecLepPt, leptonPt[1], ph_weight, "Second Lepton Pt", 10, 0.0, 100);
	}else{
		fillhisto(LeadLepPt, leptonPt[1], ph_weight, "Leading Lepton Pt", 10, 0.0, 100);
		fillhisto(SecLepPt, leptonPt[0], ph_weight, "Second Lepton Pt", 10, 0.0, 100);
	}
	//if (!isdata) fillhisto(DiTauRes, GenZMass-RegHmass, ph_weight, "Resoltuion Ditau Mass", 100, -100, 100);
	if (!isdata) fillhisto(VisRes, GenZMass-Emumass, ph_weight, "Resolution Visible Mass", 100, -100, 100);
	
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
	string hRecoPt0vsGenPt = "hRecoPt0vsGenPt_"+cut;
	string hRecoPt1vsGenPt = "hRecoPt1vsGenPt_"+cut;
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
	if (!isdata){
		fill2Dhisto(hRecoPt0vsGenPt, jetGenPt[0], jetPt[0], ph_weight, "RecoPt vs GenPt jet0", 100, 0, 200, 100, 0.0, 200);
		fill2Dhisto(hRecoPt1vsGenPt, jetGenPt[1], jetPt[1], ph_weight, "RecoPt vs GenPt jet1", 100, 0, 200, 100, 0.0, 200);
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
	}
}// TH2FDistributions		

void EventShapeDistributions(string cut, double ph_weight){
	
	string hPtbalance = "hPtbalance_"+cut;
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
	string NegIncRes  = "hResolutionZmassSVDIncNegSol_"+cut;
	string SVDRes = "hResolutionZmassSVD_"+cut;
	
	fillhisto(hPtbalance, Hpt/Zpt, ph_weight, "Pt balance ratio of  Z and H", 20, 0, 10);
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
	if (!isdata && ph_weight>-98) fillhisto(NegIncRes, GenZMass-ZmassNegInclu, ph_weight, "Resoltuion Zmass SVD Including Negative Solutions", 100, -100, 100);
	if (!isdata && ZmassSVD>-98) fillhisto(SVDRes, GenZMass-ZmassSVD, ph_weight, "Resolution Zmass SVD", 100, -100, 100);
	
	
	
	
	
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
