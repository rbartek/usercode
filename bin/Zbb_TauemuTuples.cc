/**
 @brief Example of an analysis code to read tree ntuple
 and create histogram
 Usage:
 
 \verbatim
 mkdir [directoryName]
 Zbb_TauemuTuples <inputfile> [outputfilename] [directoryName]
 \endverbatim
 
 @param inputfile Either a ROOT file or an ASCII file containing list of
 ROOT files.
 
 @param outputfile Name of the ROOT file which contains the histogram.
 Defaulted to 'output.root'
 
 @author Rachel Wilken <rachel.wilken@cern.ch>
 
 @date Fri Nov 9 2011
 
 */

#include "Zbb_TauemuTuples.h"
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
    tautreeReader sample(ifilename, std::string("tautree"));
	
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
	TMVA_tree->Branch("DeltaPhiHV",&DeltaPhiHV, "DeltaPhiHV/F");
	TMVA_tree->Branch("Hpt",&Hpt, "Hpt/F");
	TMVA_tree->Branch("Zpt",&Zpt, "Zpt/F");
	TMVA_tree->Branch("lep0pt",&lep0pt, "lep0pt/F");
	TMVA_tree->Branch("ScalarSumPt",&ScalarSumPt, "ScalarSumPt/F");
	TMVA_tree->Branch("Ht",&Ht, "Ht/F");
	TMVA_tree->Branch("EtaStandDev",&EtaStandDev, "EtaStandDev/F");
	TMVA_tree->Branch("UnweightedEta",&UnweightedEta, "UnweightedEta/F");
	TMVA_tree->Branch("EvntShpCircularity",&EvntShpCircularity, "EvntShpCircularity/F");
	TMVA_tree->Branch("alpha_j",&alpha_j, "alpha_j/F");
	TMVA_tree->Branch("qtb1",&qtb1, "qtb1/F");
	TMVA_tree->Branch("nSV",&nSV, "nSV/I");
	TMVA_tree->Branch("Trigweight",&Trigweight, "Trigweight/F");
	TMVA_tree->Branch("B2011PUweight",&B2011PUweight, "B2011PUweight/F");
	TMVA_tree->Branch("A2011PUweight",&A2011PUweight, "A2011PUweight/F");
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
	BDT_tree->Branch("DeltaPhiHV",&DeltaPhiHV, "DeltaPhiHV/F");
	BDT_tree->Branch("Hpt",&Hpt, "Hpt/F");
	BDT_tree->Branch("Zpt",&Zpt, "Zpt/F");
	BDT_tree->Branch("lep0pt",&lep0pt, "lep0pt/F");
	BDT_tree->Branch("ScalarSumPt",&ScalarSumPt, "ScalarSumPt/F");
	BDT_tree->Branch("Ht",&Ht, "Ht/F");
	BDT_tree->Branch("EtaStandDev",&EtaStandDev, "EtaStandDev/F");
	BDT_tree->Branch("UnweightedEta",&UnweightedEta, "UnweightedEta/F");
	BDT_tree->Branch("EvntShpCircularity",&EvntShpCircularity, "EvntShpCircularity/F");
	BDT_tree->Branch("alpha_j",&alpha_j, "alpha_j/F");
	BDT_tree->Branch("qtb1",&qtb1, "qtb1/F");
	BDT_tree->Branch("nSV",&nSV, "nSV/I");
	BDT_tree->Branch("Trigweight",&Trigweight, "Trigweight/F");
	BDT_tree->Branch("B2011PUweight",&B2011PUweight, "B2011PUweight/F");
	BDT_tree->Branch("A2011PUweight",&A2011PUweight, "A2011PUweight/F");
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
	
	
	
	
	TTree *BDT_btree = new TTree("BDT_btree","Tree of b jets for BDT output");
	if (isZjets){
		BDT_btree->Branch("nJets",&nJets, "nJets/I");
		BDT_btree->Branch("Naj",&Naj, "Naj/I");
		BDT_btree->Branch("Nab",&Nab, "Nab/I");
		BDT_btree->Branch("naJets",&naJets, "naJets/I");
		BDT_btree->Branch("eventFlavor",&eventFlavor, "eventFlavor/I");
		BDT_btree->Branch("CSV0",&CSV0, "CSV0/F");
		BDT_btree->Branch("CSV1",&CSV1, "CSV1/F");
		BDT_btree->Branch("Emumass",&Emumass, "Emumass/F");
		BDT_btree->Branch("Hmass",&Hmass, "Hmass/F");
		BDT_btree->Branch("DeltaPhiHV",&DeltaPhiHV, "DeltaPhiHV/F");
		BDT_btree->Branch("Hpt",&Hpt, "Hpt/F");
		BDT_btree->Branch("Zpt",&Zpt, "Zpt/F");
		BDT_btree->Branch("lep0pt",&lep0pt, "lep0pt/F");
		BDT_btree->Branch("ScalarSumPt",&ScalarSumPt, "ScalarSumPt/F");
		BDT_btree->Branch("Ht",&Ht, "Ht/F");
		BDT_btree->Branch("EtaStandDev",&EtaStandDev, "EtaStandDev/F");
		BDT_btree->Branch("UnweightedEta",&UnweightedEta, "UnweightedEta/F");
		BDT_btree->Branch("EvntShpCircularity",&EvntShpCircularity, "EvntShpCircularity/F");
		BDT_btree->Branch("alpha_j",&alpha_j, "alpha_j/F");
		BDT_btree->Branch("qtb1",&qtb1, "qtb1/F");
		BDT_btree->Branch("nSV",&nSV, "nSV/I");
		BDT_btree->Branch("Trigweight",&Trigweight, "Trigweight/F");
		BDT_btree->Branch("B2011PUweight",&B2011PUweight, "B2011PUweight/F");
		BDT_btree->Branch("A2011PUweight",&A2011PUweight, "A2011PUweight/F");
		BDT_btree->Branch("btag2CSF",&btag2CSF, "btag2CSF/F");
		BDT_btree->Branch("DetaJJ",&DetaJJ, "DetaJJ/F");
		BDT_btree->Branch("jetCHF0",&jetCHF[0], "jetCHF0/F");
		BDT_btree->Branch("jetCHF1",&jetCHF[1], "jetCHF1/F");
		BDT_btree->Branch("jetPt0",&jetPt[0], "jetPt0/F");
		BDT_btree->Branch("jetPt1",&jetPt[1], "jetPt1/F");
		BDT_btree->Branch("jetEta0",&jetEta[0], "jetEta0/F");
		BDT_btree->Branch("jetEta1",&jetEta[1], "jetEta1/F");
		BDT_btree->Branch("CSVNewShape0",&CSVNewShape[0], "CSVNewShape0/F");
		BDT_btree->Branch("CSVNewShape1",&CSVNewShape[1], "CSVNewShape1/F");
		BDT_btree->Branch("lep1pt",&leptonPt[1], "lep1pt/F");
		BDT_btree->Branch("lep_pfCombRelIso0",&lep_pfCombRelIso[0], "lep_pfCombRelIso0/F");
		BDT_btree->Branch("lep_pfCombRelIso1",&lep_pfCombRelIso[1], "lep_pfCombRelIso1/F");
		BDT_btree->Branch("DphiJJ",&DphiJJ, "DphiJJ/F");
		BDT_btree->Branch("RMS_eta",&RMS_eta, "RMS_eta/F");
		BDT_btree->Branch("PtbalZH",&PtbalZH, "PtbalZH/F");
		BDT_btree->Branch("EventPt",&EventPt, "EventPt/F");
		BDT_btree->Branch("EventMass",&EventMass, "EventMass/F");
		BDT_btree->Branch("AngleHemu",&AngleHemu, "AngleHemu/F");
		BDT_btree->Branch("Centrality",&Centrality, "Centrality/F");
		BDT_btree->Branch("MET",&MET, "MET/F");
		BDT_btree->Branch("EvntShpAplanarity",&EvntShpAplanarity, "EvntShpAplanarity/F");
		BDT_btree->Branch("EvntShpSphericity",&EvntShpSphericity, "EvntShpSphericity/F");
		BDT_btree->Branch("EvntShpIsotropy",&EvntShpIsotropy, "EvntShpIsotropy/F");
		BDT_btree->Branch("Zphi",&Zphi, "Zphi/F");
		BDT_btree->Branch("Hphi",&Hphi, "Hphi/F");
		BDT_btree->Branch("SV_mass",&SV_mass, "SV_mass/F");
		BDT_btree->Branch("Mte",&Mte, "Mte/F");
		BDT_btree->Branch("Mtmu",&Mtmu, "Mtmu/F");
		BDT_btree->Branch("delPullAngle",&delPullAngle, "delPullAngle/F");
		BDT_btree->Branch("delPullAngle2",&delPullAngle2, "delPullAngle2/F");
		BDT_btree->Branch("Mt",&Mt, "Mt/F");
		BDT_btree->Branch("dPhiHMET",&dPhiHMET, "dPhiHMET/F");
		BDT_btree->Branch("DeltaPhijetMETmin",&DeltaPhijetMETmin, "DeltaPhijetMETmin/F");
		BDT_btree->Branch("DeltaPhijetMETZtaumin",&DeltaPhijetMETZtaumin, "DeltaPhijetMETZtaumin/F");
		BDT_btree->Branch("AngleEMU",&AngleEMU, "AngleEMU/F");
		BDT_btree->Branch("AaronEleMissE",&AaronEleMissE, "AaronEleMissE/F");
		BDT_btree->Branch("AaronMuMissE",&AaronMuMissE, "AaronMuMissE/F");
		BDT_btree->Branch("Dphiemu",&Dphiemu, "Dphiemu/F");
		BDT_btree->Branch("delRjj",&delRjj, "delRjj/F");
		BDT_btree->Branch("Detaemu",&Detaemu, "Detaemu/F");
		BDT_btree->Branch("DphiEleMET",&DphiEleMET, "DphiEleMET/F");
		BDT_btree->Branch("dphiMuMET",&dphiMuMET, "dphiMuMET/F");
		BDT_btree->Branch("PtbalMETH",&PtbalMETH, "PtbalMETH/F");
		BDT_btree->Branch("topPt",&topPt, "topPt/F");
		BDT_btree->Branch("MassEleb0",&MassEleb0, "MassEleb0/F");
		BDT_btree->Branch("MassMub0",&MassMub0, "MassMub0/F");
		BDT_btree->Branch("MassEleb1",&MassEleb1, "MassEleb1/F");
		BDT_btree->Branch("MassMub1",&MassMub1, "MassMub1/F");
		BDT_btree->Branch("METsig",&METsig, "METsig/F");
		BDT_btree->Branch("delRemu",&delRemu, "delRemu/F");
		BDT_btree->Branch("PtbalZMET",&PtbalZMET, "PtbalZMET/F");
		BDT_btree->Branch("DphiZMET",&DphiZMET, "DphiZMET/F");
		BDT_btree->Branch("Zmass",&Zmass, "Zmass/F");
		BDT_btree->Branch("ZmassSVD",&ZmassSVD, "ZmassSVD/F");
		BDT_btree->Branch("ZmassSVDnegSol",&ZmassSVDnegSol, "ZmassSVDnegSol/F");
		BDT_btree->Branch("ZmassNegInclu",&ZmassNegInclu, "ZmassNegInclu/F");
		BDT_btree->Branch("DphiSecondMET",&DphiSecondMET, "DphiSecondMET/F");
		BDT_btree->Branch("DphiLeadMET",&DphiLeadMET, "DphiLeadMET/F");
		BDT_btree->Branch("topMass",&topMass, "topMass/F");
		BDT_btree->Branch("ProjVisT",&ProjVisT, "ProjVisT/F");
		BDT_btree->Branch("ProjMissT",&ProjMissT, "ProjMissT/F");
		
		
	}
	
	
	
    // Here one can declare histograms
    // In compiled C++, it possible not to use pointer type (TH1F*) and
    // use the non-pointer type (TH1F) of ROOT object.
	
	//Declare Histograms
	
    TH1F hallhJet_pt  ("hallhJet_pt","Pt of all jets in event",		150, 0.0, 150);
    TH1F hEleIDReco  ("hEleIDReco","Electron ID/Reco weight after EleFakeCuts ",		101, -0.5, 2.0);
	TH1F hCutFlow	("hCutFlow",  "Selection",					11, 0.5, 11.5);
	
	//NotThisCut
	TH1F hMemu_NotThisCut		("hMemu_NotThisCut",  "Invariant Mass of two Leptons Not This Cut", 100, 0, 200);
	TH1F hDphiZMET_NotThisCut		("hDphiZMET_NotThisCut","Delta phi between Z and MET Not This Cut",  30, 0, 4);
	TH1F hCSV0_NotThisCut	("hCSV0_NotThisCut","CSV0 Not This Cut", 30, 0, 1.5);
	TH1F hPtjj_NotThisCut		("hPtjj_NotThisCut","Pt Higgs Not This Cut", 100, 0, 200);
	TH1F hAngleEMU_NotThisCut		("hAngleEMU_NotThisCut", "Angle(e,#mu) Not This Cut", 30, 0, 3.5);
	TH1F hCHFb0_NotThisCut ("hCHFb0_NotThisCut","Charged Hadron Fraction b0 Not This Cut", 60, 0.0, 1.2);
	
	if (debug) std::cout << "all histograms declared " << std::endl;
	
    // Loop over all events.
    // For other methods to access event/navigate through the sample,
    // see the documentation of RootTreeReader class.
	float  N_Vtype =0.0, Ntrigger = 0.0, Npreselect =0.0, N_Mjj =0.0, N_ptH =0.0, N_DphiZMET =0.0, NMemu =0.0, N_Naj =0.0, N_CSV0 =0.0, N_EfakeCuts =0.0, N_jetCHF0 = 0.0;
	int event =0, Ntree =0, N_isomu = 0, N_tightDoubleele = 0, N_loosedoubleEle = 0, N_sigleEleWP80 = 0, N_singleEleDiJet = 0;
	float N_TopCR = 0.0, N_SingleTopCR = 0.0, N_LFCR = 0.0, N_HFCR =0.0, N_LFCR_HF = 0.0, N_LFCR_LF =0.0, N_HFCR_LF = 0.0, N_HFCR_HF= 0.0;
	float FailedJetID=0.0, NtreePU = 0.0, NTMVAtree =0.0, NBDTtree =0.0, NBDTbtree = 0.0;
	int  NNegEleMissE = 0, NPosEleMissE = 0, NPosMuMissE = 0, NNegMuMissE = 0;
	int NBothMissENeg = 0, NMixedEleMissENeg = 0, NMixedMuonMissENeg =0;
	bool firstevent = true;
	bool eleposclose = true, eleposfar = true, elenegclose= true, elenegfar= true;
	
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
					if(AaronMuMissE < 0 || AaronEleMissE < 0){
						if (AaronMuMissE < 0 && AaronEleMissE < 0) NBothMissENeg++;
						if (AaronMuMissE > 0 && AaronEleMissE < 0) NMixedEleMissENeg++;
						if (AaronMuMissE < 0 && AaronEleMissE > 0) NMixedMuonMissENeg++;
					}				
					if (AaronEleMissE < 0) NNegEleMissE++;
					if (AaronEleMissE > 0) NPosEleMissE++;
					if (AaronMuMissE > 0) NPosMuMissE++;
					if (AaronMuMissE < 0) NNegMuMissE++;
					
					
					/*			if (sample.vLepton_eta[0] > 0 && phielectron > 0  && AngleEMU < 0.5 && eleposclose) {
					 cout << "Electron has eta of " << sample.vLepton_eta[0]<< " phi of " <<phielectron << " costheta is " << CosThetaEle<< " Angle from Muon " <<AngleEMU<< endl;
					 cout << "Muon has eta of " <<sample.vLepton_eta[1] << " phi of " <<phimuon << " costheta is " << CosThetaMu<<" Det is " << DetA << endl;
					 cout << "sin phi ele " << sin(phielectron) << " cos phi ele " << cos(phielectron)<< endl;
					 cout << "sin phi mu " << sin(phimuon) << " cos phi mu " << cos(phimuon)<< endl;
					 cout << "MET is " << MET << " METx " << METx << " METy " << METy << endl;
					 eleposclose = false;
					 }
					 if (sample.vLepton_eta[0] > 0 && phielectron > 0  && AngleEMU > 2.9 && eleposfar) {
					 cout << "Electron has eta of " << sample.vLepton_eta[0]<< " phi of " <<phielectron << " costheta is " << CosThetaEle<< " Angle from Muon " <<AngleEMU<< endl;
					 cout << "Muon has eta of " <<sample.vLepton_eta[1] << " phi of " <<phimuon << " costheta is " << CosThetaMu<<" Det is " << DetA << endl;
					 cout << "sin phi ele " << sin(phielectron) << " cos phi ele " << cos(phielectron)<< endl;
					 cout << "sin phi mu " << sin(phimuon) << " cos phi mu " << cos(phimuon)<< endl;
					 cout << "MET is " << MET << " METx " << METx << " METy " << METy << endl;
					 eleposfar = false;
					 }
					 if (sample.vLepton_eta[0] < 0 && phielectron < 0  && AngleEMU < 0.5 && elenegclose) {
					 cout << "Electron has eta of " << sample.vLepton_eta[0]<< " phi of " <<phielectron << " costheta is " << CosThetaEle<< " Angle from Muon " <<AngleEMU<< endl;
					 cout << "Muon has eta of " <<sample.vLepton_eta[1] << " phi of " <<phimuon << " costheta is " << CosThetaMu<<" Det is " << DetA << endl;
					 cout << "sin phi ele " << sin(phielectron) << " cos phi ele " << cos(phielectron)<< endl;
					 cout << "sin phi mu " << sin(phimuon) << " cos phi mu " << cos(phimuon)<< endl;
					 cout << "MET is " << MET << " METx " << METx << " METy " << METy << endl;
					 elenegclose = false;
					 }
					 if (sample.vLepton_eta[0] < 0 && phielectron < 0  && AngleEMU >2.9 && elenegfar) {
					 cout << "Electron has eta of " << sample.vLepton_eta[0]<< " phi of " <<phielectron << " costheta is " << CosThetaEle<< " Angle from Muon " <<AngleEMU<< endl;
					 cout << "Muon has eta of " <<sample.vLepton_eta[1] << " phi of " <<phimuon << " costheta is " << CosThetaMu<<" Det is " << DetA << endl;
					 cout << "sin phi ele " << sin(phielectron) << " cos phi ele " << cos(phielectron)<< endl;
					 cout << "sin phi mu " << sin(phimuon) << " cos phi mu " << cos(phimuon)<< endl;
					 cout << "MET is " << MET << " METx " << METx << " METy " << METy  << endl;
					 elenegfar = false;
					 }
					 
					 if (fabs(Detaemu)<0.1 && fabs(DphiZMET)< 0.1 ) {
					 cout << "Electron has eta of " << sample.vLepton_eta[0]<< " phi of " <<phielectron << " costheta is " << CosThetaEle<< " Angle from Muon " <<AngleEMU<< endl;
					 cout << "Muon has eta of " <<sample.vLepton_eta[1] << " phi of " <<phimuon << " costheta is " << CosThetaMu<<" Det is " << DetA << endl;
					 cout << "sin phi ele " << sin(phielectron) << " cos phi ele " << cos(phielectron)<< endl;
					 cout << "sin phi mu " << sin(phimuon) << " cos phi mu " << cos(phimuon)<< endl;
					 cout << "MET is " << MET << " METx " << METx << " METy " << METy << " phi of MET is " << sample.MET_phi << endl;
					 }
					 */			
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
				
				
				
				
				if ((sample.vLepton_pt[1] > 20 || sample.vLepton_pt[0] > 20) && (sample.vLepton_pt[1]>10 && sample.vLepton_pt[0]> 10)){
					//	if (sample.triggerFlags[39]||sample.triggerFlags[40]||sample.triggerFlags[41]){
					Ntrigger = Ntrigger + 1*PUweight2011;
					//cout << "PUweight2011 inside trigger requirement" << PUweight2011 << endl;
					JetDistributions("HLT", weight);
					LeptonDistributions("HLT", weight);
					TH2FDistributions("HLT", weight);
					EventShapeDistributions("HLT", weight);
					EventDistributions("HLT", weight);				
					if ( jetPt[0] > 20 && jetPt[1] > 20 && 
						fabs(sample.vLepton_eta[0]) < 2.5 && fabs(sample.vLepton_eta[1]) < 2.4 && fabs(sample.hJet_eta[0]) < 2.5 &&
						fabs(sample.hJet_eta[1]) < 2.5 && sample.hJet_id[0]==1 && sample.hJet_id[1]==1 && sample.hbhe){
						Npreselect = Npreselect + 1*PUweight2011;
						FOM_tree->Fill();
						/*						if (Hmass<250){
						// Ntree = Ntree + 1*PUweight2011;
						 Ntree++;
						 if (Ntree%2){
						 TMVA_tree->Fill();
						 NTMVAtree =NTMVAtree + 1*PUweight2011;}else{
						 if (isZjets){
						 if (sample.eventFlav == 5){
						 NBDTbtree = NBDTbtree + 1*PUweight2011;
						 BDT_btree->Fill();
						 }else{
						 BDT_tree->Fill();
						 NBDTtree = NBDTtree + 1*PUweight2011;
						 }
						 }
						 else{
						 BDT_tree->Fill();
						 NBDTtree = NBDTtree + 1*PUweight2011;
						 }
						 }	
						 }//if higgs mass reasonable */
						JetDistributions("PreSelect", weight);
						LeptonDistributions("PreSelect", weight);
						TH2FDistributions("PreSelect", weight);
						EventShapeDistributions("PreSelect", weight);
						EventDistributions("PreSelect", weight);
						if (isDATA) isdata = true;
						if ((AngleEMU>0.2&&delRemu>0.3) &&fabs(DphiZMET) < 1.250 && (Hmass>=80)&&(Hmass<=150) && CSV0>0.244  && Emumass< 85 && Emumass>10) hPtjj_NotThisCut.Fill(Hpt, weight);
						if (fabs(DphiZMET) < 1.250 && (Hmass>=80)&&(Hmass<=150) && CSV0>0.244 && Emumass< 85 && Emumass>10 && Hpt>20) hAngleEMU_NotThisCut.Fill(AngleEMU, weight);
						if(AngleEMU>0.25&&delRemu>0.4){
							N_EfakeCuts = N_EfakeCuts + 1*PUweight2011;
							//+ 1*PUweight2011;
							JetDistributions("EleFakeCuts", weight);
							LeptonDistributions("EleFakeCuts", weight);
							TH2FDistributions("EleFakeCuts", weight);
							EventShapeDistributions("EleFakeCuts", weight);
							EventDistributions("EleFakeCuts", weight);									
							if (Hpt>20) { 
								N_ptH = N_ptH + 1*PUweight2011;
								JetDistributions("PtH", weight);
								LeptonDistributions("PtH", weight);
								TH2FDistributions("PtH", weight);
								EventShapeDistributions("PtH", weight);
								EventDistributions("PtH", weight);						
								if (fabs(DphiZMET) < 1.250 && (Hmass>=80)&&(Hmass<=150) && CSV0>0.244 ) hMemu_NotThisCut.Fill(Emumass, weight);
								if (Emumass> 85 && fabs(DphiZMET) < 1.250 && (Hmass>=80)&&(Hmass<=150) && CSV0>0.244 ){
									if (jetCHF[0] > 0.35 && CSV0 > 0.898 && CSV1< 0.5 && METsig < 50 && Nab < 1 && EvntShpAplanarity < 0.15
										&& EventPt<125 && topPt < 125 && nJets < 6 & AngleEMU < 2.25 && jetCHF[0] < 0.95 && Centrality > 0.2
										&& Centrality < 0.45 && EvntShpCircularity > 0.2 && fabs(Detaemu)<1.9 && fabs(DphiEleMET) > 1
										&& Dphiemu>2 && EventMass < 350 && EvntShpIsotropy > 0.3 && ProjMissT < 25 && ProjMissT > -75
										&& Hpt < 130 && ScalarSumJetPt < 225 && ScalarSumPt < 275 && Mt>20 && delRemu < 2.5 && delRjj<4
										&& fabs(dphiMuMET)>1 && DeltaPhiHV>0.5 && topMass< 250 ){
										if ((sample.vLepton_pt[1] > 15 || sample.vLepton_pt[0] > 15) && (sample.vLepton_pt[1]>20 && sample.vLepton_pt[0]> 20)){
											N_SingleTopCR = N_SingleTopCR + 1*PUweight2011;
											JetDistributions("SingleTopCR", weight);
											LeptonDistributions("SingleTopCR", weight);
											TH2FDistributions("SingleTopCR", weight);
											EventShapeDistributions("SingleTopCR", weight);
											EventDistributions("SingleTopCR", weight);
										}//HLT emulation											
									}//CHF cut
								}// Single Top CR								
								if (Emumass< 85 && Emumass>10) { 
									 NMemu = NMemu + 1*PUweight2011;
									JetDistributions("MemuCut", weight);
									LeptonDistributions("MemuCut", weight);
									TH2FDistributions("MemuCut", weight);
									EventShapeDistributions("MemuCut", weight);
									EventDistributions("MemuCut", weight);
									if (fabs(DphiZMET) < 1.250 && (Hmass>=80)&&(Hmass<=150) && CSV0<0.244 && fabs(Dphiemu) > 1.5 && Naj < 2 
									&& EventMass < 300 && ((jetCHF[2] < 0)||(jetCHF[2] > 0.1 && jetCHF[2] < 0.7))
									&& fabs(dPhiHMET) > 0.5 && MET < 75 && Emumass < 70 && ProjMissT < 75 && Mt<125){	
										N_LFCR = N_LFCR + 1*PUweight2011;
										JetDistributions("LFCR", weight);
										LeptonDistributions("LFCR", weight);
										TH2FDistributions("LFCR", weight);
										EventShapeDistributions("LFCR", weight);
										EventDistributions("LFCR", weight);	
										if (isZjets	) {
											if (eventFlavor == 5){
												N_LFCR_HF = N_LFCR_HF + 1*PUweight2011;
												JetDistributions("LFCR_HF", weight);
												LeptonDistributions("LFCR_HF", weight);
												TH2FDistributions("LFCR_HF", weight);
												EventShapeDistributions("LFCR_HF", weight);
												EventDistributions("LFCR_HF", weight);	
											}else{
												N_LFCR_LF = N_LFCR_LF + 1*PUweight2011;
												JetDistributions("LFCR_LF", weight);
												LeptonDistributions("LFCR_LF", weight);
												TH2FDistributions("LFCR_LF", weight);
												EventShapeDistributions("LFCR_LF", weight);
												EventDistributions("LFCR_LF", weight);												
											}}}// LF Control Region								
									if (fabs(DphiZMET) < 1.250 && (Hmass>=80)&&(Hmass<=150) && DeltaPhiHV>=2.9 ) hCSV0_NotThisCut.Fill(CSV0, weight);
									if(CSV0>0.244){
										N_CSV0 = N_CSV0 + 1*PUweight2011;
										JetDistributions("CSV0", weight);
										LeptonDistributions("CSV0", weight);
										TH2FDistributions("CSV0", weight);
										EventShapeDistributions("CSV0", weight);
										EventDistributions("CSV0", weight);							
										if (fabs(DphiZMET) > 1.250 && Hmass>=80&&Hmass<=150){
											if (jetCHF[0] > 0.15 && CSV0>0.5 && METsig > 10 && EvntShpAplanarity > 0.1 && CSV1>0.244){
												if ((sample.vLepton_pt[1] > 25 || sample.vLepton_pt[0] > 25) && (sample.vLepton_pt[1]>15 && sample.vLepton_pt[0]> 15)){
													N_TopCR = N_TopCR + 1*PUweight2011;
													JetDistributions("TopCR", weight);
													LeptonDistributions("TopCR", weight);
													TH2FDistributions("TopCR", weight);
													EventShapeDistributions("TopCR", weight);
													EventDistributions("TopCR", weight);
												}//HLT emulation											
											}//CHF cut
										}// Top CR
										if ((Hmass>=80)&&(Hmass<=150)) hDphiZMET_NotThisCut.Fill(DphiZMET, weight);
										if (fabs(DphiZMET) < 1.250){
											N_DphiZMET = N_DphiZMET + 1*PUweight2011;
											JetDistributions("DphiZMET", weight);
											LeptonDistributions("DphiZMET", weight);
											TH2FDistributions("DphiZMET", weight);
											EventShapeDistributions("DphiZMET", weight);
											EventDistributions("DphiZMET", weight);
											if((Hmass>=80)&&(Hmass<=150)){
												N_Mjj = N_Mjj + 1*PUweight2011; 
												JetDistributions("Mjj", weight);
												LeptonDistributions("Mjj", weight);
												TH2FDistributions("Mjj", weight);
												EventShapeDistributions("Mjj", weight);
												EventDistributions("Mjj", weight);
												//NtreePU = NtreePU + 1*PUweight2011;
												Ntree++;
												if (Ntree%2){
													TMVA_tree->Fill();
												NTMVAtree = NTMVAtree + 1*PUweight2011;}else{
													if (isZjets){
														if (sample.eventFlav == 5){
															NBDTbtree = NBDTbtree + 1*PUweight2011;
															BDT_btree->Fill();
														}else{
															BDT_tree->Fill();
															NBDTtree = NBDTtree + 1*PUweight2011;
														}
													}
													else{
														BDT_tree->Fill();
														NBDTtree = NBDTtree + 1*PUweight2011;
													}		
												}											
											}//Mjj			
											if ((Hmass<80 || Hmass>150) && CSV0>0.244 && CSV1<0.5 && (ProjMissT-0.25*ProjVisT)>-1 && fabs(Dphiemu) > 1 
											&& EventMass < 150 && EventPt < 65 && MET < 60 && METsig < 25 && Nab < 1 && (ProjMissT-0.25*ProjVisT)<50
											&& Emumass < 70 && Emumass > 30 && jetCHF[0] > 0.3 && ScalarSumJetPt<120){	
												N_HFCR = N_HFCR + 1*PUweight2011;
												JetDistributions("HFCR", weight);
												LeptonDistributions("HFCR", weight);
												TH2FDistributions("HFCR", weight);
												EventShapeDistributions("HFCR", weight);
												EventDistributions("HFCR", weight);											
												if (isZjets	) {
													if (eventFlavor == 5){
														N_HFCR_HF = N_HFCR_HF + 1*PUweight2011;
														JetDistributions("HFCR_HF", weight);
														LeptonDistributions("HFCR_HF", weight);
														TH2FDistributions("LFCR_HF", weight);
														EventShapeDistributions("HFCR_HF", weight);
														EventDistributions("HFCR_HF", weight);	
													}else{
														N_HFCR_LF = N_HFCR_LF + 1*PUweight2011;
														JetDistributions("HFCR_LF", weight);
														LeptonDistributions("HFCR_LF", weight);
														TH2FDistributions("HFCR_LF", weight);
														EventShapeDistributions("HFCR_LF", weight);
														EventDistributions("HFCR_LF", weight);												
													}}}// LF Control Region
											if((Hmass>=80)&&(Hmass<=150)){
												if (Naj <2) hCHFb0_NotThisCut.Fill(jetCHF[0], weight);
												if(jetCHF[0] > 0.15) {
													N_jetCHF0 = N_jetCHF0 + 1*PUweight2011;
													JetDistributions("jetCHF0", weight);
													LeptonDistributions("jetCHF0", weight);
													TH2FDistributions("jetCHF0", weight);
													EventShapeDistributions("jetCHF0", weight);
													EventDistributions("jetCHF0", weight);														
													if (Naj<2){
														N_Naj = N_Naj + 1*PUweight2011;
														JetDistributions("Naj", weight);
														LeptonDistributions("Naj", weight);
														TH2FDistributions("Naj", weight);
														EventShapeDistributions("Naj", weight);
														EventDistributions("Naj", weight);
													}	//Naj
												}//jetCHF[0]
											}//Mjj							
										}// Delta Phi Z, MET (Z is emu vectors only)
									}//CSV
								}//emu mass requirement
							}//pt of higgs
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
	std::cout << "Hpt: " << N_ptH << endl;
	std::cout << "Memu: " << NMemu << endl;	
	std::cout << "CSV0: " << N_CSV0 << endl;
	std::cout << "DphiZMET: " << N_DphiZMET << endl;
	std::cout << "Mjj: " << N_Mjj << endl;
	std::cout << "CHFb0: " << N_jetCHF0 << endl;
	std::cout << "Naj: " << N_Naj << endl;
	
	
	
	std::cout << endl << endl;
	std::cout << "Top CR: " << N_TopCR << endl;
	std::cout << "Single Top CR: " << N_SingleTopCR << endl;
	std::cout << "LF CR: " << N_LFCR << endl;
	std::cout << "HF CR: " << N_HFCR << endl;
	std::cout << "For DY LF CR light jet: " << N_LFCR_LF << endl;
	std::cout << "For DY HF CR light jet: " << N_HFCR_LF << endl;
	std::cout << "For DY LF CR b jet: " << N_LFCR_HF << endl;
	std::cout << "For DY HF CR b jet: " << N_HFCR_HF << endl;
	
	
	std::cout << endl << endl;
	std::cout << "Number of Events that failed JetID " << FailedJetID << endl;
	std::cout << "BDT tree: " << NBDTtree << endl;
	std::cout << "BDT  bjet tree: " << NBDTbtree << endl;
	std::cout << "TMVA tree: " << NTMVAtree << endl;
	
ofstream myfile;
myfile.open(TString::Format("%sDataEntry.txt",directory.c_str()).Data());	
myfile<<TString::Format("\t %s",directory.c_str()).Data()<<endl;
//	myfile<<TString::Format("Number of events \t %f \t %0.5f",event,event*LumiWeight).Data()<<endl;
//	myfile<<TString::Format("Vtype 5 \t %f \t %0.5f",N_Vtype,N_Vtype*LumiWeight).Data()<<endl;
	myfile<<TString::Format("emu trigger \t %0.5f \t %0.5f",Ntrigger,Ntrigger*LumiWeight).Data()<<endl;
	myfile<<TString::Format("PreSelection \t %0.5f \t %0.5f",Npreselect,Npreselect*LumiWeight).Data()<<endl;
	myfile<<TString::Format("EleFakeCuts \t %0.5f \t %0.5f",N_EfakeCuts,N_EfakeCuts*LumiWeight).Data()<<endl;
	myfile<<TString::Format("Hpt \t %0.5f \t %0.5f",N_ptH,N_ptH*LumiWeight).Data()<<endl;
	myfile<<TString::Format("Memu \t %0.5f \t %0.5f",NMemu,NMemu*LumiWeight).Data()<<endl;
	myfile<<TString::Format("CSV0 \t %0.5f \t %0.5f",N_CSV0,N_CSV0*LumiWeight).Data()<<endl;
	myfile<<TString::Format("DphiZMET \t %0.5f \t %0.5f",N_DphiZMET,N_DphiZMET*LumiWeight).Data()<<endl;
	myfile<<TString::Format("Mjj\t %0.5f \t %0.5f",N_Mjj,N_Mjj*LumiWeight).Data()<<endl;
	myfile<<TString::Format("CHFb0 \t %0.5f \t %0.5f",N_jetCHF0,N_jetCHF0*LumiWeight).Data()<<endl;
	myfile<<TString::Format("Naj \t %0.5f \t %0.5f",N_Naj,N_Naj*LumiWeight).Data()<<endl;
	myfile<<endl;
	myfile<<TString::Format("Top CR	\t %0.5f \t %0.5f",N_TopCR,N_TopCR*LumiWeight).Data()<<endl;
	myfile<<TString::Format("Single Top CR \t %0.5f \t %0.5f",N_SingleTopCR,N_SingleTopCR*LumiWeight).Data()<<endl;
	myfile<<TString::Format("LF CR \t %0.5f \t %0.5f",N_LFCR,N_LFCR*LumiWeight).Data()<<endl;
	myfile<<TString::Format("HF CR \t %0.5f \t %0.5f",N_HFCR,N_HFCR*LumiWeight).Data()<<endl;
	
myfile.close();


/*	
	std::cout << endl << endl;
	std::cout << "Number of Events with postive Electron Missing Energy " << NPosEleMissE << endl;
	std::cout << "Number of Events with negative Electron Missing Energy " << NNegEleMissE << endl;
	std::cout << "Number of Events with postive Muon Missing Energy " << NPosMuMissE << endl;
	std::cout << "Number of Events with negative Muon Missing Energy " << NNegMuMissE << endl;
	std::cout << "Number of Events with both leptons having negative missing enegery " << NBothMissENeg << endl;
	std::cout << "Number of Events with negative electron Missing Energy postitive muon Emiss " << NMixedEleMissENeg << endl;
	std::cout << "Number of Events with negative Muon Missing Energy postitive electron Emiss " << NMixedMuonMissENeg << endl;
	
	std::cout << endl << endl;
	std::cout << "Number of Events passing single muon trigger " << N_isomu << endl;
	std::cout << "Number of Events passing double electron tight: " << N_tightDoubleele << " loose: " <<  N_loosedoubleEle<< endl;
	std::cout << "Number of Events passing single electron WP80 " << N_sigleEleWP80 << endl;
	std::cout << "Number of Events passing single electron diJet " << N_singleEleDiJet << endl;
	
	
	*/
	
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
	hDphiZMET_NotThisCut.Draw();
	c1.Print((directory+"/DphiZMET_NotThisCut"+suffixps).c_str());
	
	c1.Clear(); // don't create a new canvas
	hCSV0_NotThisCut.Draw();
	c1.Print((directory+"/CSV0_NotThisCut"+suffixps).c_str());
	
	c1.Clear(); // don't create a new canvas
	hAngleEMU_NotThisCut.Draw();
	c1.Print((directory+"/AngleEMU_NotThisCut"+suffixps).c_str());
	
	c1.Clear(); // don't create a new canvas
	hCHFb0_NotThisCut.Draw();
	c1.Print((directory+"/CHFb0_NotThisCut"+suffixps).c_str());
	
	c1.Clear(); // don't create a new canvas
	hPtjj_NotThisCut.Draw();
	c1.Print((directory+"/PtH_NotThisCut"+suffixps).c_str());		
	
	hCutFlow.SetBinContent(1,N_Vtype );
	hCutFlow.SetBinContent(2,Ntrigger );
	hCutFlow.SetBinContent(3,Npreselect );
	hCutFlow.SetBinContent(4,N_EfakeCuts );
	hCutFlow.SetBinContent(5,N_ptH );
	hCutFlow.SetBinContent(6,NMemu );
	hCutFlow.SetBinContent(7,N_CSV0 );
	hCutFlow.SetBinContent(8,N_DphiZMET );
	hCutFlow.SetBinContent(9,N_Mjj );
	hCutFlow.SetBinContent(10,N_jetCHF0 );
	hCutFlow.SetBinContent(11,N_Naj );
	hCutFlow.Draw();
	c1.Print((directory+"/CutFlow"+suffixps).c_str());
	
	TMVA_tree->Write();
	BDT_tree->Write();
	
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
	fillhisto(scalarSumHiggsJetPt, ScalarSumHiggsJetPt, ph_weight, "scalar sum higgs jet pt", 25, 0, 300);
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
	fillhisto(hScalarSumPt, ScalarSumPt, ph_weight, "scalar sum of pt of four particles", 80, 0.0, 400);
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
	fillhisto(hPzetaCut, ProjMissT-0.25*ProjVisT, ph_weight, "P_#zeta Cut", 75, -200, 100);
	
	
	
	
	
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
