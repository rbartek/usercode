/**
 @brief Example of an analysis code to read tree ntuple
 and create histogram
 Usage:
 
 \verbatim
 mkdir [directoryName]
 Zbb_JetStudy <inputfile> [outputfilename] [directoryName]
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
		BDT_btree->Branch("ScalarSumJetPt",&ScalarSumJetPt, "ScalarSumJetPt/F");
		BDT_btree->Branch("ScalarSumHiggsJetPt",&ScalarSumHiggsJetPt, "ScalarSumHiggsJetPt/F");
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
	
	TH1F hLeadingPtHiggs  ("hLeadingPtHiggs","Pt of Leading Higgs Jet",		150, 0.0, 300);
    TH1F hLeadingPtJets  ("hLeadingPtJets","Pt of Leading Jet",		150, 0.0, 300);
    TH1F hSecondPtHiggs  ("hSecondPtHiggs","Pt of Second Higgs Jet",		150, 0.0, 300);
	TH1F hSecondPtJets  ("hSecondPtJets","Pt of Second Jet",		150, 0.0, 300);
	TH1F hClosestMatch	("hClosestMatch",  "Delta R closet match",		40, 0, .2);
    TH1F hLeadingCSVHiggs  ("hLeadingCSVHiggs","CSV of Leading Higgs Jet",		44, -1.1, 1.1);
	TH1F hLeadingCSVJets  ("hLeadingCSVJets","CSV of Leading Jet",		44, -1.1, 1.1);
    TH1F hSecondCSVHiggs  ("hSecondCSVHiggs","CSV of Second Higgs Jet",		44, -1.1, 1.1);
	TH1F hSecondCSVJets  ("hSecondCSVJets","CSV of Second Jet",		44, -1.1, 1.1);
	TH1F hHiggsMass103Matched  ("hHiggsMass103Matched","Mass of two matched jets mH=103", 70, 0, 350);
	TH1F hHiggsMass125Matched  ("hHiggsMass125Matched","Mass of two matched jets mH=125", 70, 0, 350);
	TH1F hHiggsMass114Matched  ("hHiggsMass114Matched","Mass of two matched jets mH=114", 70, 0, 350);
	TH1F hHiggsMass109Matched  ("hHiggsMass109Matched","Mass of two matched jets mH=109", 70, 0, 350);
	TH1F hHiggsMassPt  ("hHiggsMassPt","Mass of two higest pt jets", 70, 0, 350);
	TH1F hHiggsMassCSV  ("hHiggsMassCSV","Mass of two jets with highest CSV", 70, 0, 350);
	TH1F hHiggsMassgenInfo  ("hHiggsMassgenInfo","Mass of two gen jets", 70, 0, 350);
	TH1F hJetHiggsMass  ("hJetHiggsMass","Mass of two hjets", 70, 0, 350);
	TH1F hHiggsMassMatched  ("hHiggsMassMatched","Mass of two matched jets", 70, 0, 350);
	TH1F hHiggsMassCSVthenPtT  ("hHiggsMassCSVthenPtT","Mass of two jets N_CSVthenPtT", 70, 0, 350);
	TH1F hHiggsMassCSVthenPtM  ("hHiggsMassCSVthenPtM","Mass of two jets N_CSVthenPtM", 70, 0, 350);
	TH1F hHiggsMassCSVthenPtTM  ("hHiggsMassCSVthenPtTM","Mass of two jets N_CSVthenPtTM", 70, 0, 350);
	TH1F hHiggsMassCSVthenPtL  ("hHiggsMassCSVthenPtL","Mass of two jets N_CSVthenPtL", 70, 0, 350);
	TH1F hHiggsMassCSVthenPtTL  ("hHiggsMassCSVthenPtTL","Mass of two jets N_CSVthenPtTL", 70, 0, 350);
	TH1F hHiggsMassCSVthenPtML  ("hHiggsMassCSVthenPtML","Mass of two jets N_CSVthenPtML", 70, 0, 350);
	
	if (debug) std::cout << "all histograms declared " << std::endl;
	
    // Loop over all events.
    // For other methods to access event/navigate through the sample,
    // see the documentation of RootTreeReader class.
	float  N_Vtype =0.0, Ntrigger = 0.0, Npreselect =0.0, N_Mjj =0.0, N_DphiZMET =0.0, NMemu =0.0, N_Naj =0.0, N_CSV0 =0.0, N_EfakeCuts =0.0, N_jetCHF0 = 0.0;
	int event =0, Ntree =0;
	float N_TopCR = 0.0, N_SingleTopCR = 0.0, N_LFCR = 0.0, N_HFCR =0.0;
	float FailedJetID=0.0, NTMVAtree =0.0, NBDTtree =0.0, NBDTbtree = 0.0;
	int  NNegEleMissE = 0, NPosEleMissE = 0, NPosMuMissE = 0, NNegMuMissE = 0;
	int NBothMissENeg = 0, NMixedEleMissENeg = 0, NMixedMuonMissENeg =0;
	bool firstevent = true;
	
	int NVtype =0, N_NonHiggsPt = 0, N_NonHiggsCSV = 0, N_NonhJetMatched = 0, N_count = 0, highPtCSV = 0, secondPtCSV = 0;
	int ReasonbleMass = 0;
	int N_HighPtcorrect = 0, N_HighCSVcorrect = 0, N_MassMethod = 0, N_CSVthenPt = 0;
	int N_CSVthenPtM = 0, N_CSVthenPtTM = 0, N_CSVthenPtL = 0, N_CSVthenPtTL = 0, N_CSVthenPtML =0;
	int N_114MassMethod = 0, N_125MassMethod = 0, N_109MassMethod = 0;
	int N_CSVthenPt5 = 0, N_CSVthenPtTT = 0, N_CSVthenPtMM = 0, N_CSVthenPtLL = 0, N_CSVthenPt55 =0;
	int N_CSVthenPt4 = 0, N_CSVthenPt5L = 0, N_CSVthenPtM5 = 0, N_CSVthenPtT4 = 0, N_CSVthenPtM4 =0;
	int N_CSVthenPt44 = 0, N_CSVthenPt4L = 0, N_CSVthenPt54 = 0,N_CSVthenPtT5 = 0 ;
	int N_CSVclosest = 0, N_PtClosestJet = 0, N_CSVHighPtT = 0, N_CSVHighPtM = 0, N_CSVHighPtL = 0, N_CSVHighPt5 = 0, N_CSVHighPt4 = 0;
	int N_complicated = 0;
	
	//	double Emumass = 91.1976;
	
    do {
		//event = event + 1*PUweight2011;
		event++;
		cout << "   Event " << event << endl;
		isdata = false;
		//	if (isM50sample && sample.genZpt < 100) cout << "event cut out to remove overlap with DY_PtZ" << endl;
		if (((isM50sample && sample.genZpt < 100) || !isM50sample) && ((HFfile && sample.eventFlav == 5)||!HFfile) && ((LFfile && sample.eventFlav!=5)||!LFfile)){
			//cout << "Vtype is " << sample.Vtype << endl;
			if (sample.Vtype == 5 && sample.EVENT_json == 1 ){
				N_Vtype++;
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
				//One can access the ntuple leaves directly from sample object
				float jetgenEta[10], jetgenPt[10], jetgenPhi[10];
				for (int i=0; i < 10; ++i) {
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
					jetgenEta[i] = -99.99;
					jetgenPt[i] = -99.99;
					jetgenPhi[i] = -99.99;
				}
				
				int nhJets= 0;
				
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
					if ( sample.hJet_pt[k] > 20 && fabs(sample.hJet_eta[k]) < 2.5 && sample.hJet_id[k]==1 && sample.hbhe){
						indexedJetPt.push_back(std::pair<size_t,double>(k,(double) sample.hJet_pt[k]));
						indexedPt.push_back(std::pair<size_t,double>(k,(double) sample.hJet_pt[k]));
						nJets++;
						nhJets++;
						if (debug)	cout << "CSV discriminator: " << sample.hJet_csv[k] << endl;
						if(sample.hJet_csv[k]<=0 || sample.hJet_csv[k]>=1) CSVshapeNew=sample.hJet_csv[k];
						else if(sample.hJet_flavour[k]==0) CSVshapeNew=sample.hJet_csv[k];
						else if(fabs(sample.hJet_flavour[k])==5) CSVshapeNew=btagNew->ib->Eval(sample.hJet_csv[k]);
						else if(fabs(sample.hJet_flavour[k])==4) CSVshapeNew=btagNew->ic->Eval(sample.hJet_csv[k]);
						else if(fabs(sample.hJet_flavour[k])!=5 && fabs(sample.hJet_flavour[k])!=4)  CSVshapeNew=btagNew->il->Eval(sample.hJet_csv[k]);
						if (!(event%500)) cout << "The orignial CSV value was " << sample.hJet_csv[k] << "  the corrected value is " << CSVshapeNew << endl;
						indexedJetCSV.push_back(std::pair<size_t,double>(k,(double) CSVshapeNew));
						jetPt[k] = sample.hJet_pt[k];
						jetEta[k] = sample.hJet_eta[k];
						jetPhi[k] = sample.hJet_phi[k];
						jetCSV[k] = sample.hJet_csv[k];
						jetCHF[k] = sample.hJet_chf[k];	
						jetgenEta[k] = sample.hJet_genEta[k];
						jetgenPt[k] = sample.hJet_genPt[k];
						jetgenPhi[k] = sample.hJet_genPhi[k];						
						CSVNewShape[k] = CSVshapeNew;
					}
					CSVshapeNew = -99.99;
				}// end k for loop
				
				for (int a=0;a<sample.naJets;a++){
					if ( sample.aJet_pt[a]> 20 && fabs(sample.aJet_eta[a]) < 2.5 && sample.aJet_id[a]==1 && sample.hbhe){
					bool isduplicate = false;
						if(sample.aJet_csv[a]<=0 || sample.aJet_csv[a]>=1) CSVshapeNew=sample.aJet_csv[a];
						else if(sample.aJet_flavour[a]==0) CSVshapeNew=sample.aJet_csv[a];
						else if(fabs(sample.aJet_flavour[a])==5) CSVshapeNew=btagNew->ib->Eval(sample.aJet_csv[a]);
						else if(fabs(sample.aJet_flavour[a])==4) CSVshapeNew=btagNew->ic->Eval(sample.aJet_csv[a]);
						else if(fabs(sample.aJet_flavour[a])!=5 && fabs(sample.aJet_flavour[a])!=4)  CSVshapeNew=btagNew->il->Eval(sample.aJet_csv[a]);				
						if ((sample.aJet_genPt[a]-sample.hJet_genPt[0])<0.1) isduplicate = true;
						if ((sample.aJet_genPt[a]-sample.hJet_genPt[1])<0.1) isduplicate = true;
						if (sample.aJet_genPt[a]==-99.99) isduplicate = true;
						if (a<10 && !isduplicate){
							nJets++; 
							jetPt[nhJets+a] = sample.aJet_pt[a];
							jetEta[nhJets+a] = sample.aJet_eta[a];
							jetPhi[nhJets+a] = sample.aJet_phi[a];
							jetCSV[nhJets+a] = sample.aJet_csv[a];
							jetCHF[nhJets+a] = sample.aJet_chf[a];
							jetgenEta[nhJets+a] = sample.aJet_genEta[a];
							jetgenPt[nhJets+a] = sample.aJet_genPt[a];
							jetgenPhi[nhJets+a] = sample.aJet_genPhi[a];													
							CSVNewShape[nhJets+a] = CSVshapeNew;
							indexedJetPt.push_back(std::pair<size_t,double>(nhJets+a,(double) sample.aJet_pt[a]));	
							indexedJetCSV.push_back(std::pair<size_t,double>(nhJets+a,(double) CSVshapeNew));
						}
						if ( fabs(sample.aJet_eta[a]) < 2.4 && (sample.aJet_pt[a] > 20)) Naj++;
						if (CSVshapeNew > 0.5 ) Nab++;
						ScalarSumJetPt = ScalarSumJetPt + sample.aJet_pt[a];
					}
					CSVshapeNew = -99.99;
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
				
				
				//CSVNewShape[0] = indexedJetCSV[0].second;
				if (nJets > 1) { 
					hJetHiggsMass.Fill(sample.H_mass);
					unsigned int ClosestToMass1 = 99, ClosestToMass2 = 99;
					float ClosestToHiggsMass = 9999.999;
					for (int kk = 0; kk<nJets; kk++){
						for (int jj = 1; jj<nJets; jj++){
							if (jj==kk) continue;
							FirstJet.SetPtEtaPhiM(jetPt[kk],jetEta[kk],jetPhi[kk],4.2/1000.0);
							SecondJet.SetPtEtaPhiM(jetPt[jj],jetEta[jj],jetPhi[jj],4.2/1000.0);
							
							Higgs = FirstJet+SecondJet;
							if (abs(Higgs.M()-103)<ClosestToHiggsMass) {
								ClosestToHiggsMass = abs(Higgs.M()-103);
								ClosestToMass1 = jj;
								ClosestToMass2 = kk;
							}
						}
					}
					FirstJet.SetPtEtaPhiM(jetPt[ClosestToMass1],jetEta[ClosestToMass1],jetPhi[ClosestToMass1],4.2/1000.0);
					SecondJet.SetPtEtaPhiM(jetPt[ClosestToMass2],jetEta[ClosestToMass2],jetPhi[ClosestToMass2],4.2/1000.0);
					Higgs = FirstJet+SecondJet;
					hHiggsMass103Matched.Fill(Higgs.M());
					unsigned int ClosestTo125Mass1 = 99, ClosestTo125Mass2 = 99;
					ClosestToHiggsMass = 9999.999;
					for (int kk = 0; kk<nJets; kk++){
						for (int jj = 1; jj<nJets; jj++){
							if (jj==kk) continue;
							FirstJet.SetPtEtaPhiM(jetPt[kk],jetEta[kk],jetPhi[kk],4.2/1000.0);
							SecondJet.SetPtEtaPhiM(jetPt[jj],jetEta[jj],jetPhi[jj],4.2/1000.0);
							
							Higgs = FirstJet+SecondJet;
							if (abs(Higgs.M()-125)<ClosestToHiggsMass) {
								ClosestToHiggsMass = abs(Higgs.M()-125);
								ClosestTo125Mass1 = jj;
								ClosestTo125Mass2 = kk;
							}
						}
					}
					FirstJet.SetPtEtaPhiM(jetPt[ClosestTo125Mass1],jetEta[ClosestTo125Mass1],jetPhi[ClosestTo125Mass1],4.2/1000.0);
					SecondJet.SetPtEtaPhiM(jetPt[ClosestTo125Mass2],jetEta[ClosestTo125Mass2],jetPhi[ClosestTo125Mass2],4.2/1000.0);
					Higgs = FirstJet+SecondJet;
					hHiggsMass125Matched.Fill(Higgs.M());
					unsigned int ClosestTo114Mass1 = 99, ClosestTo114Mass2 = 99;
					ClosestToHiggsMass = 9999.999;
					for (int kk = 0; kk<nJets; kk++){
						for (int jj = 1; jj<nJets; jj++){
							if (jj==kk) continue;
							FirstJet.SetPtEtaPhiM(jetPt[kk],jetEta[kk],jetPhi[kk],4.2/1000.0);
							SecondJet.SetPtEtaPhiM(jetPt[jj],jetEta[jj],jetPhi[jj],4.2/1000.0);
							
							Higgs = FirstJet+SecondJet;
							if (abs(Higgs.M()-114)<ClosestToHiggsMass) {
								ClosestToHiggsMass = abs(Higgs.M()-114);
								ClosestTo114Mass1 = jj;
								ClosestTo114Mass2 = kk;
							}
						}
					}
					FirstJet.SetPtEtaPhiM(jetPt[ClosestTo114Mass1],jetEta[ClosestTo114Mass1],jetPhi[ClosestTo114Mass1],4.2/1000.0);
					SecondJet.SetPtEtaPhiM(jetPt[ClosestTo114Mass2],jetEta[ClosestTo114Mass2],jetPhi[ClosestTo114Mass2],4.2/1000.0);
					Higgs = FirstJet+SecondJet;
					hHiggsMass114Matched.Fill(Higgs.M());
					unsigned int ClosestTo109Mass1 = 99, ClosestTo109Mass2 = 99;
					ClosestToHiggsMass = 9999.999;
					for (int kk = 0; kk<nJets; kk++){
						for (int jj = 1; jj<nJets; jj++){
							if (jj==kk) continue;
							FirstJet.SetPtEtaPhiM(jetPt[kk],jetEta[kk],jetPhi[kk],4.2/1000.0);
							SecondJet.SetPtEtaPhiM(jetPt[jj],jetEta[jj],jetPhi[jj],4.2/1000.0);
							
							Higgs = FirstJet+SecondJet;
							if (abs(Higgs.M()-109)<ClosestToHiggsMass) {
								ClosestToHiggsMass = abs(Higgs.M()-109);
								ClosestTo109Mass1 = jj;
								ClosestTo109Mass2 = kk;
							}
						}
					}
					FirstJet.SetPtEtaPhiM(jetPt[ClosestTo109Mass1],jetEta[ClosestTo109Mass1],jetPhi[ClosestTo109Mass1],4.2/1000.0);
					SecondJet.SetPtEtaPhiM(jetPt[ClosestTo109Mass2],jetEta[ClosestTo109Mass2],jetPhi[ClosestTo109Mass2],4.2/1000.0);
					Higgs = FirstJet+SecondJet;
					hHiggsMass109Matched.Fill(Higgs.M());
					
					unsigned int CSVthenPt1 = 99, CSVthenPt2 = 99;
					if(indexedJetCSV[0].second>0.898){
						CSVthenPt1 = indexedJetCSV[0].first;
						CSVthenPt2 = indexedJetCSV[1].first;							
					}else{
						CSVthenPt1 = indexedJetPt[0].first;
						CSVthenPt2 = indexedJetPt[1].first;
					}
					FirstJet.SetPtEtaPhiM(jetPt[CSVthenPt1],jetEta[CSVthenPt1],jetPhi[CSVthenPt1],4.2/1000.0);
					SecondJet.SetPtEtaPhiM(jetPt[CSVthenPt2],jetEta[CSVthenPt2],jetPhi[CSVthenPt2],4.2/1000.0);
					Higgs = FirstJet+SecondJet;
					hHiggsMassCSVthenPtT.Fill(Higgs.M());
					
					unsigned int CSVthenPtM1 = 99, CSVthenPtM2 = 99;
					if(indexedJetCSV[0].second>0.679){
						CSVthenPtM1 = indexedJetCSV[0].first;
						CSVthenPtM2 = indexedJetCSV[1].first;							
					}else{
						CSVthenPtM1 = indexedJetPt[0].first;
						CSVthenPtM2 = indexedJetPt[1].first;
					}
					FirstJet.SetPtEtaPhiM(jetPt[CSVthenPtM1],jetEta[CSVthenPtM1],jetPhi[CSVthenPtM1],4.2/1000.0);
					SecondJet.SetPtEtaPhiM(jetPt[CSVthenPtM2],jetEta[CSVthenPtM2],jetPhi[CSVthenPtM2],4.2/1000.0);
					Higgs = FirstJet+SecondJet;
					hHiggsMassCSVthenPtM.Fill(Higgs.M());
					
					unsigned int CSVthenPtTM1 = 99, CSVthenPtTM2 = 99;
					if(indexedJetCSV[0].second>0.898&&indexedJetCSV[1].second>0.679){
						CSVthenPtTM1 = indexedJetCSV[0].first;
						CSVthenPtTM2 = indexedJetCSV[1].first;							
					}else{
						CSVthenPtTM1 = indexedJetPt[0].first;
						CSVthenPtTM2 = indexedJetPt[1].first;
					}
					FirstJet.SetPtEtaPhiM(jetPt[CSVthenPtTM1],jetEta[CSVthenPtTM1],jetPhi[CSVthenPtTM1],4.2/1000.0);
					SecondJet.SetPtEtaPhiM(jetPt[CSVthenPtTM2],jetEta[CSVthenPtTM2],jetPhi[CSVthenPtTM2],4.2/1000.0);
					Higgs = FirstJet+SecondJet;
					hHiggsMassCSVthenPtTM.Fill(Higgs.M());
					
					unsigned int CSVthenPtL1 = 99, CSVthenPtL2 = 99;
					if(indexedJetCSV[0].second>0.244){
						CSVthenPtL1 = indexedJetCSV[0].first;
						CSVthenPtL2 = indexedJetCSV[1].first;							
					}else{
						CSVthenPtL1 = indexedJetPt[0].first;
						CSVthenPtL2 = indexedJetPt[1].first;
					}
					FirstJet.SetPtEtaPhiM(jetPt[CSVthenPtL1],jetEta[CSVthenPtL1],jetPhi[CSVthenPtL1],4.2/1000.0);
					SecondJet.SetPtEtaPhiM(jetPt[CSVthenPtL2],jetEta[CSVthenPtL2],jetPhi[CSVthenPtL2],4.2/1000.0);
					Higgs = FirstJet+SecondJet;
					hHiggsMassCSVthenPtL.Fill(Higgs.M());
					
					unsigned int CSVthenPtTL1 = 99, CSVthenPtTL2 = 99;
					if(indexedJetCSV[0].second>0.898&&indexedJetCSV[1].second>0.244){
						CSVthenPtTL1 = indexedJetCSV[0].first;
						CSVthenPtTL2 = indexedJetCSV[1].first;							
					}else{
						CSVthenPtTL1 = indexedJetPt[0].first;
						CSVthenPtTL2 = indexedJetPt[1].first;
					}
					FirstJet.SetPtEtaPhiM(jetPt[CSVthenPtTL1],jetEta[CSVthenPtTL1],jetPhi[CSVthenPtTL1],4.2/1000.0);
					SecondJet.SetPtEtaPhiM(jetPt[CSVthenPtTL2],jetEta[CSVthenPtTL2],jetPhi[CSVthenPtTL2],4.2/1000.0);
					Higgs = FirstJet+SecondJet;
					hHiggsMassCSVthenPtTL.Fill(Higgs.M());
					
					unsigned int CSVthenPtML1 = 99, CSVthenPtML2 = 99;
					if(indexedJetCSV[0].second>0.679&&indexedJetCSV[1].second>0.244){
						CSVthenPtML1 = indexedJetCSV[0].first;
						CSVthenPtML2 = indexedJetCSV[1].first;							
					}else{
						CSVthenPtML1 = indexedJetPt[0].first;
						CSVthenPtML2 = indexedJetPt[1].first;
					}
					FirstJet.SetPtEtaPhiM(jetPt[CSVthenPtML1],jetEta[CSVthenPtML1],jetPhi[CSVthenPtML1],4.2/1000.0);
					SecondJet.SetPtEtaPhiM(jetPt[CSVthenPtML2],jetEta[CSVthenPtML2],jetPhi[CSVthenPtML2],4.2/1000.0);
					Higgs = FirstJet+SecondJet;
					hHiggsMassCSVthenPtML.Fill(Higgs.M());
					unsigned int CSVthenPt51 = 99, CSVthenPt52 = 99;
					if(indexedJetCSV[0].second>0.5){
						CSVthenPt51 = indexedJetCSV[0].first;
						CSVthenPt52 = indexedJetCSV[1].first;							
					}else{
						CSVthenPt51 = indexedJetPt[0].first;
						CSVthenPt52 = indexedJetPt[1].first;
					}
					unsigned int CSVthenPtTT1 = 99, CSVthenPtTT2 = 99;
					if(indexedJetCSV[0].second>0.898&&indexedJetCSV[1].second>0.898){
						CSVthenPtTT1 = indexedJetCSV[0].first;
						CSVthenPtTT2 = indexedJetCSV[1].first;							
					}else{
						CSVthenPtTT1 = indexedJetPt[0].first;
						CSVthenPtTT2 = indexedJetPt[1].first;
					}
					unsigned int CSVthenPtMM1 = 99, CSVthenPtMM2 = 99;
					if(indexedJetCSV[0].second>0.679&&indexedJetCSV[1].second>0.679){
						CSVthenPtMM1 = indexedJetCSV[0].first;
						CSVthenPtMM2 = indexedJetCSV[1].first;							
					}else{
						CSVthenPtMM1 = indexedJetPt[0].first;
						CSVthenPtMM2 = indexedJetPt[1].first;
					}
					unsigned int CSVthenPtLL1 = 99, CSVthenPtLL2 = 99;
					if(indexedJetCSV[0].second>0.244&&indexedJetCSV[1].second>0.244){
						CSVthenPtLL1 = indexedJetCSV[0].first;
						CSVthenPtLL2 = indexedJetCSV[1].first;							
					}else{
						CSVthenPtLL1 = indexedJetPt[0].first;
						CSVthenPtLL2 = indexedJetPt[1].first;
					}
					unsigned int CSVthenPt551 = 99, CSVthenPt552 = 99;
					if(indexedJetCSV[0].second>0.5&&indexedJetCSV[1].second>0.5){
						CSVthenPt551 = indexedJetCSV[0].first;
						CSVthenPt552 = indexedJetCSV[1].first;							
					}else{
						CSVthenPt551 = indexedJetPt[0].first;
						CSVthenPt552 = indexedJetPt[1].first;
					}
					unsigned int CSVthenPt41 = 99, CSVthenPt42 = 99;
					if(indexedJetCSV[0].second>0.4){
						CSVthenPt41 = indexedJetCSV[0].first;
						CSVthenPt42 = indexedJetCSV[1].first;							
					}else{
						CSVthenPt41 = indexedJetPt[0].first;
						CSVthenPt42 = indexedJetPt[1].first;
					}
					unsigned int CSVthenPt5L1 = 99, CSVthenPt5L2 = 99;
					if(indexedJetCSV[0].second>0.5&&indexedJetCSV[1].second>0.244){
						CSVthenPt5L1 = indexedJetCSV[0].first;
						CSVthenPt5L2 = indexedJetCSV[1].first;							
					}else{
						CSVthenPt5L1 = indexedJetPt[0].first;
						CSVthenPt5L2 = indexedJetPt[1].first;
					}
					unsigned int CSVthenPtM51 = 99, CSVthenPtM52 = 99;
					if(indexedJetCSV[0].second>0.679&&indexedJetCSV[1].second>0.5){
						CSVthenPtM51 = indexedJetCSV[0].first;
						CSVthenPtM52 = indexedJetCSV[1].first;							
					}else{
						CSVthenPtM51 = indexedJetPt[0].first;
						CSVthenPtM52 = indexedJetPt[1].first;
					}
					unsigned int CSVthenPtT41 = 99, CSVthenPtT42 = 99;
					if(indexedJetCSV[0].second>0.898&&indexedJetCSV[1].second>0.4){
						CSVthenPtT41 = indexedJetCSV[0].first;
						CSVthenPtT42 = indexedJetCSV[1].first;							
					}else{
						CSVthenPtT41 = indexedJetPt[0].first;
						CSVthenPtT42 = indexedJetPt[1].first;
					}
					unsigned int CSVthenPtM41 = 99, CSVthenPtM42 = 99;
					if(indexedJetCSV[0].second>0.679&&indexedJetCSV[1].second>0.4){
						CSVthenPtM41 = indexedJetCSV[0].first;
						CSVthenPtM42 = indexedJetCSV[1].first;							
					}else{
						CSVthenPtM41 = indexedJetPt[0].first;
						CSVthenPtM42 = indexedJetPt[1].first;
					}
					unsigned int CSVthenPt441 = 99, CSVthenPt442 = 99;
					if(indexedJetCSV[0].second>0.4&&indexedJetCSV[1].second>0.4){
						CSVthenPt441 = indexedJetCSV[0].first;
						CSVthenPt442 = indexedJetCSV[1].first;							
					}else{
						CSVthenPt441 = indexedJetPt[0].first;
						CSVthenPt442 = indexedJetPt[1].first;
					}
					unsigned int CSVthenPt4L1 = 99, CSVthenPt4L2 = 99;
					if(indexedJetCSV[0].second>0.4&&indexedJetCSV[1].second>0.244){
						CSVthenPt4L1 = indexedJetCSV[0].first;
						CSVthenPt4L2 = indexedJetCSV[1].first;							
					}else{
						CSVthenPt4L1 = indexedJetPt[0].first;
						CSVthenPt4L2 = indexedJetPt[1].first;
					}
					unsigned int CSVthenPt541 = 99, CSVthenPt542 = 99;
					if(indexedJetCSV[0].second>0.5&&indexedJetCSV[1].second>0.4){
						CSVthenPt541 = indexedJetCSV[0].first;
						CSVthenPt542 = indexedJetCSV[1].first;							
					}else{
						CSVthenPt541 = indexedJetPt[0].first;
						CSVthenPt542 = indexedJetPt[1].first;
					}
					unsigned int CSVthenPtT51 = 99, CSVthenPtT52 = 99;
					if(indexedJetCSV[0].second>0.898&&indexedJetCSV[1].second>0.5){
						CSVthenPtT51 = indexedJetCSV[0].first;
						CSVthenPtT52 = indexedJetCSV[1].first;							
					}else{
						CSVthenPtT51 = indexedJetPt[0].first;
						CSVthenPtT52 = indexedJetPt[1].first;
					}
					unsigned int CSVclosest1 = 99, CSVclosest2 = 99;
					float DELR = 999.99;
					CSVclosest1 = indexedJetCSV[0].first;
					for (int iter = 0; iter<nJets; iter++){
						if(iter==(int)CSVclosest1)continue;
						float diffR = sqrt((jetEta[iter]-jetEta[CSVclosest1])*(jetEta[iter]-jetEta[CSVclosest1])+(jetPhi[iter]-jetPhi[CSVclosest1])*(jetPhi[iter]-jetPhi[CSVclosest1]));
						if (diffR<DELR)CSVclosest2 = iter;
					}				
					unsigned int PtClosestJet1 = 99, PtClosestJet2 = 99;
					DELR = 999.99;
					PtClosestJet1 = indexedJetCSV[0].first;
					for (int iter = 0; iter<nJets; iter++){
						if(iter==(int)PtClosestJet1)continue;
						float diffR = sqrt((jetEta[iter]-jetEta[PtClosestJet1])*(jetEta[iter]-jetEta[PtClosestJet1])+(jetPhi[iter]-jetPhi[PtClosestJet1])*(jetPhi[iter]-jetPhi[PtClosestJet1]));
						if (diffR<DELR)PtClosestJet2 = iter;
					}				
					unsigned int CSVHighPtT1 = 99, CSVHighPtT2=99;
					CSVHighPtT1 = indexedJetPt[0].first;
					for (int iter = 0; iter<nJets; iter++){
						if ((indexedJetCSV[iter].first!=CSVHighPtT1)&&indexedJetCSV[iter].second>0.898) CSVHighPtT2 = indexedJetCSV[iter].first;
					}
					if(CSVHighPtT2==99)CSVHighPtT2=indexedJetPt[1].first;
					unsigned int CSVHighPtM1 = 99, CSVHighPtM2=99;
					CSVHighPtM1 = indexedJetPt[0].first;
					for (int iter = 0; iter<nJets; iter++){
						if ((indexedJetCSV[iter].first!=CSVHighPtM1)&&indexedJetCSV[iter].second>0.679) CSVHighPtM2 = indexedJetCSV[iter].first;
					}
					if(CSVHighPtM2==99)CSVHighPtM2=indexedJetPt[1].first;
					unsigned int CSVHighPtL1 = 99, CSVHighPtL2=99;
					CSVHighPtL1 = indexedJetPt[0].first;
					for (int iter = 0; iter<nJets; iter++){
						if ((indexedJetCSV[iter].first!=CSVHighPtL1)&&indexedJetCSV[iter].second>0.244) CSVHighPtL2 = indexedJetCSV[iter].first;
					}
					if(CSVHighPtM2==99)CSVHighPtM2=indexedJetPt[1].first;
					unsigned int CSVHighPt51 = 99, CSVHighPt52=99;
					CSVHighPt51 = indexedJetPt[0].first;
					for (int iter = 0; iter<nJets; iter++){
						if ((indexedJetCSV[iter].first!=CSVHighPt51)&&indexedJetCSV[iter].second>0.5) CSVHighPt52 = indexedJetCSV[iter].first;
					}
					if(CSVHighPt52==99)CSVHighPt52=indexedJetPt[1].first;
					unsigned int CSVHighPt41 = 99, CSVHighPt42=99;
					CSVHighPt41 = indexedJetPt[0].first;
					for (int iter = 0; iter<nJets; iter++){
						if ((indexedJetCSV[iter].first!=CSVHighPt41)&&indexedJetCSV[iter].second>0.4) CSVHighPt42 = indexedJetCSV[iter].first;
					}
					if(CSVHighPt42==99)CSVHighPt42=indexedJetPt[1].first;
					
					unsigned int complicated1 = 99, complicated2 = 99;
					if(indexedJetCSV[0].second>0.679&&indexedJetCSV[1].second>0.244){
						complicated1 = indexedJetCSV[0].first;
						complicated2 = indexedJetCSV[1].first;							
					}else if (indexedJetCSV[0].second>0.679&&indexedJetCSV[1].second<0.244){
						complicated1 = indexedJetCSV[0].first;
						if(indexedJetCSV[0].first!=indexedJetPt[0].first)complicated2 = indexedJetPt[0].first;
						if(indexedJetCSV[0].first==indexedJetPt[0].first)complicated2 = indexedJetPt[1].first;	
					}else if(indexedJetCSV[0].second<0.679&&indexedJetCSV[1].second>0.244){
						complicated1 = indexedJetCSV[0].first;
						complicated2 = indexedJetCSV[1].first;							
					}else{
						complicated1 = indexedJetPt[0].first;
						complicated2 = indexedJetPt[1].first;
					}
					
					
					if ( indexedJetPt[0].first == indexedJetCSV[0].first) highPtCSV++;
					if ( indexedJetPt[1].first == indexedJetCSV[1].first) secondPtCSV++;
					
					hLeadingPtJets.Fill(indexedJetPt[0].second);
					hSecondPtJets.Fill(indexedJetPt[1].second);
					hLeadingCSVJets.Fill(indexedJetCSV[0].second);
					hSecondCSVJets.Fill(indexedJetCSV[1].second);
					
					FirstJet.SetPtEtaPhiM(jetPt[indexedJetPt[0].first],jetEta[indexedJetPt[0].first],jetPhi[indexedJetPt[0].first],4.2/1000.0);
					SecondJet.SetPtEtaPhiM(jetPt[indexedJetPt[1].first],jetEta[indexedJetPt[1].first],jetPhi[indexedJetPt[1].first],4.2/1000.0);
					
					Higgs = FirstJet+SecondJet;
					hHiggsMassPt.Fill(Higgs.M());
					
					FirstJet.SetPtEtaPhiM(jetPt[indexedJetCSV[0].first],jetEta[indexedJetCSV[0].first],jetPhi[indexedJetCSV[0].first],4.2/1000.0);
					SecondJet.SetPtEtaPhiM(jetPt[indexedJetCSV[1].first],jetEta[indexedJetCSV[1].first],jetPhi[indexedJetCSV[1].first],4.2/1000.0);
					
					Higgs = FirstJet+SecondJet;
					hHiggsMassCSV.Fill(Higgs.M());
					
					
					CSV1 = indexedJetCSV[1].second;	
					CSV0 = indexedJetCSV[1].second;	
					
					//match to gen object
					float delpt = 999.99;
					float deleta = 99.99;
					float delphi = 99.99;
					float delR = 99.99;
					float EventMinDelR = 99.99;
					unsigned int bestmatchB = 99;
					unsigned int bestmatchBbar = 99;
					bool goodGenInfo = true;
					NVtype++;
					
					if ((sample.genB_pt-sample.genB_eta)<0.1 || (sample.genBbar_pt-sample.genBbar_eta)<0.1) goodGenInfo = false;
					if (goodGenInfo){ 
						N_count++;
						
						FirstJet.SetPtEtaPhiM(sample.genB_pt,sample.genB_eta,sample.genB_phi,4.2/1000.0);
						SecondJet.SetPtEtaPhiM(sample.genBbar_pt,sample.genBbar_eta,sample.genBbar_phi,4.2/1000.0);
						
						Higgs = FirstJet+SecondJet;
						if (Higgs.M()>20){
						hHiggsMassgenInfo.Fill(Higgs.M());
						ReasonbleMass++;
						for (int k=0;k<nJets;k++){
							if (nhJets>1){
								if((CSVNewShape[0]<CSVNewShape[k]) && (CSVNewShape[1]<CSVNewShape[k]) ) N_NonHiggsCSV++;
								if((jetPt[0]<jetPt[k]) && (jetPt[1]<jetPt[k])) N_NonHiggsPt++;
							}							
							//float diffpt = abs(jetPt[k]-sample.genB_pt);
							//float diffeta = abs(jetEta[k]-sample.genB_eta);
							//float diffphi	= abs(jetPhi[k]-sample.genB_phi);
							float diffR = sqrt((jetgenEta[k]-sample.genB_eta)*(jetgenEta[k]-sample.genB_eta)+(jetgenPhi[k]-sample.genB_phi)*(jetgenPhi[k]-sample.genB_phi));
							if (diffR<EventMinDelR) EventMinDelR = diffR;
							if (diffR<delR) {
								//								if (diffR<delR && diffphi<delphi && diffeta<deleta &&diffpt<delpt) {
								bestmatchB = k;
								//	delpt = diffpt;
								//	deleta = diffeta;
								//	delphi = diffphi;
								delR = diffR;
							}
							cout << "jet pt " << jetgenPt[k] << " eta " << jetgenEta[k] << " phi " << jetgenPhi[k] << endl;
						}
							cout << "gen B info pt " << sample.genB_pt << " eta " << sample.genB_eta << " phi " << sample.genB_phi <<  endl;						
							cout << "gen Bbar info pt " << sample.genBbar_pt << " eta " << sample.genBbar_eta << " phi " << sample.genBbar_phi <<  endl;						
						if (bestmatchB == 99) {
							cout << "not a match to b quark from Jets collection" << endl;
							for (int k=0;k<nJets;k++){
								float diffpt = abs(jetgenPt[k]-sample.genB_pt);
								float diffeta = abs(jetgenEta[k]-sample.genB_eta);
								float diffphi	= abs(jetgenPhi[k]-sample.genB_phi);
								float diffR = sqrt((jetgenEta[k]-sample.genB_eta)*(jetgenEta[k]-sample.genB_eta)+(jetgenPhi[k]-sample.genB_phi)*(jetgenPhi[k]-sample.genB_phi));
								cout << "diffpt " << diffpt << " diffeta " << diffeta << " diffphi " << diffphi << " diffR " << diffR << endl;
								cout << "jet pt " << jetgenPt[k] << " eta " << jetgenEta[k] << " phi " << jetgenPhi[k] << endl;
							}	
							cout << "gen info pt " << sample.genB_pt << " eta " << sample.genB_eta << " phi " << sample.genB_phi <<  endl;						
						}
						delpt = 999.99, deleta = 99.99, delphi = 99.99, delR = 99.99;
						for (int kbar=0;kbar<nJets;kbar++){
							float diffR = sqrt((jetgenEta[kbar]-sample.genBbar_eta)*(jetgenEta[kbar]-sample.genBbar_eta)+(jetgenPhi[kbar]-sample.genBbar_phi)*(jetgenPhi[kbar]-sample.genBbar_phi));
							if (diffR<EventMinDelR) EventMinDelR = diffR;
							if (diffR<delR ) {
								bestmatchBbar = kbar;
								delR = diffR;
							}					
						}
						hClosestMatch.Fill(EventMinDelR);
						if (bestmatchBbar == 99){
							cout << "not a match to Bbar from hJets collection" << endl;
							for (int kbar=0;kbar<nJets;kbar++){
								float diffpt = abs(jetgenPt[kbar]-sample.genBbar_pt);
								float diffeta = abs(jetgenEta[kbar]-sample.genBbar_eta);
								float diffphi	= abs(jetgenPhi[kbar]-sample.genBbar_phi);
								float diffR = sqrt((jetgenEta[kbar]-sample.genBbar_eta)*(jetgenEta[kbar]-sample.genBbar_eta)+(jetgenPhi[kbar]-sample.genBbar_phi)*(jetgenPhi[kbar]-sample.genBbar_phi));
								cout << "diffpt " << diffpt << " diffeta " << diffeta << " diffphi " << diffphi << " diffR " << diffR << endl;
								cout << "Jet pt " << jetgenPt[kbar] << " eta " << jetgenEta[kbar] << " phi " << jetgenPhi[kbar] <<  endl;
							}
							cout << "Gen Info pt " << sample.genBbar_pt << " eta " << sample.genBbar_eta << " phi " << sample.genBbar_phi << endl;
						}
						bool isBquark = false, isBbarquark = false;
						if ((bestmatchBbar == bestmatchB)&&bestmatchBbar<9){
							cout << "    same jet" << endl;
							float diffR = sqrt((jetgenEta[bestmatchB]-sample.genB_eta)*(jetgenEta[bestmatchB]-sample.genB_eta)+(jetgenPhi[bestmatchB]-sample.genB_phi)*(jetgenPhi[bestmatchB]-sample.genB_phi));
							float diffRbar = sqrt((jetgenEta[bestmatchBbar]-sample.genBbar_eta)*(jetgenEta[bestmatchBbar]-sample.genBbar_eta)+(jetgenPhi[bestmatchBbar]-sample.genBbar_phi)*(jetgenPhi[bestmatchBbar]-sample.genBbar_phi));
							if ( diffR<diffRbar) {
								isBquark = true;
								bestmatchBbar = 99;
								delpt = 999.99, deleta = 99.99, delphi = 99.99, delR = 99.99;
								for (int kbar=0;kbar<nJets;kbar++){
									if ((unsigned int)kbar==bestmatchB) continue;
									float bdiffR = sqrt((jetgenEta[kbar]-sample.genBbar_eta)*(jetgenEta[kbar]-sample.genBbar_eta)+(jetgenPhi[kbar]-sample.genBbar_phi)*(jetgenPhi[kbar]-sample.genBbar_phi));
									if (bdiffR<delR ) {
										bestmatchBbar = kbar;
										delR = bdiffR;
									}					
								}		
								if (bestmatchBbar == 99){
									cout << "not a match to Bbar from hJets collection" << endl;
									for (int kbar=0;kbar<nJets;kbar++){
										float diffpt = abs(jetgenPt[kbar]-sample.genBbar_pt);
										float diffeta = abs(jetgenEta[kbar]-sample.genBbar_eta);
										float diffphi	= abs(jetgenPhi[kbar]-sample.genBbar_phi);
										float diffR = sqrt((jetgenEta[kbar]-sample.genBbar_eta)*(jetgenEta[kbar]-sample.genBbar_eta)+(jetgenPhi[kbar]-sample.genBbar_phi)*(jetgenPhi[kbar]-sample.genBbar_phi));
										cout << "diffpt " << diffpt << " diffeta " << diffeta << " diffphi " << diffphi << " diffR " << diffR << endl;
										cout << "Jet pt " << jetgenPt[kbar] << " eta " << jetgenEta[kbar] << " phi " << jetgenPhi[kbar] <<  endl;
									}
									cout << "Gen Info pt " << sample.genBbar_pt << " eta " << sample.genBbar_eta << " phi " << sample.genBbar_phi << endl;
								}
								
							}
							if (diffR>diffRbar) {
								isBbarquark = true;		
								bestmatchB = 99;
								delpt = 999.99, deleta = 99.99, delphi = 99.99, delR = 99.99;
								for (int k=0;k<nJets;k++){
									if((unsigned int)k==bestmatchBbar) continue;
									float adiffR = sqrt((jetgenEta[k]-sample.genB_eta)*(jetgenEta[k]-sample.genB_eta)+(jetgenPhi[k]-sample.genB_phi)*(jetgenPhi[k]-sample.genB_phi));
									if (adiffR<delR ) {
										bestmatchB = k;
										delR = adiffR;
									}
								}
								if (bestmatchB == 99) {
									cout << "not a match to b quark from Jets collection" << endl;
									for (int k=0;k<nJets;k++){
										float diffpt = abs(jetgenPt[k]-sample.genB_pt);
										float diffeta = abs(jetgenEta[k]-sample.genB_eta);
										float diffphi	= abs(jetgenPhi[k]-sample.genB_phi);
										float diffR = sqrt((jetgenEta[k]-sample.genB_eta)*(jetgenEta[k]-sample.genB_eta)+(jetgenPhi[k]-sample.genB_phi)*(jetgenPhi[k]-sample.genB_phi));
										cout << "diffpt " << diffpt << " diffeta " << diffeta << " diffphi " << diffphi << " diffR " << diffR << endl;
										cout << "jet pt " << jetgenPt[k] << " eta " << jetgenEta[k] << " phi " << jetgenPhi[k] << endl;
									}	
									cout << "gen info pt " << sample.genB_pt << " eta " << sample.genB_eta << " phi " << sample.genB_phi <<  endl;						
								}
							}
						}
						
							float diffeta = abs(jetgenEta[bestmatchB]-sample.genB_eta);
							float diffphi	= abs(jetgenPhi[bestmatchB]-sample.genB_phi);
							float diffR = sqrt(diffeta*diffeta+diffphi*diffphi);
							if(diffR>1)bestmatchB = 99;
							 diffeta = abs(jetgenEta[bestmatchBbar]-sample.genBbar_eta);
							 diffphi	= abs(jetgenPhi[bestmatchBbar]-sample.genBbar_phi);
							 diffR = sqrt(diffeta*diffeta+diffphi*diffphi);
							if(diffR>1)bestmatchBbar = 99;
							if (nhJets>1 && (bestmatchBbar>1||bestmatchB>1)) N_NonhJetMatched++;
							//if (nhJets == 1 && (bestmatchBbar>0||bestmatchB>0) )N_NonhJetMatched++;
							//if (nhJets == 0 && (bestmatchBbar!=99||bestmatchB!=99) )N_NonhJetMatched++;
							
							
							cout << "b quark matched to jet" << 		bestmatchB << " Bbar matched to jet" << 		bestmatchBbar << endl;

						if(bestmatchBbar != bestmatchB && bestmatchB!=99 &&bestmatchBbar!=99){
							if (jetPt[bestmatchBbar]>jetPt[bestmatchB]){
								hLeadingPtHiggs.Fill(jetPt[bestmatchBbar]);
								hSecondPtHiggs.Fill(jetPt[bestmatchB]);
							}else{
								hLeadingPtHiggs.Fill(jetPt[bestmatchB]);
								hSecondPtHiggs.Fill(jetPt[bestmatchBbar]);
							}
							if (CSVNewShape[bestmatchBbar]>CSVNewShape[bestmatchB]){
								hLeadingCSVHiggs.Fill(CSVNewShape[bestmatchBbar]);
								hSecondCSVHiggs.Fill(CSVNewShape[bestmatchB]);
							}else{
								hLeadingCSVHiggs.Fill(CSVNewShape[bestmatchB]);
								hSecondCSVHiggs.Fill(CSVNewShape[bestmatchBbar]);
							}
						
						
						FirstJet.SetPtEtaPhiM(jetPt[bestmatchB],jetEta[bestmatchB],jetPhi[bestmatchB],4.2/1000.0);
						SecondJet.SetPtEtaPhiM(jetPt[bestmatchBbar],jetEta[bestmatchBbar],jetPhi[bestmatchBbar],4.2/1000.0);
						
						Higgs = FirstJet+SecondJet;
						hHiggsMassMatched.Fill(Higgs.M());
						
						
						if (((indexedJetPt[0].first == bestmatchBbar) || (indexedJetPt[0].first == bestmatchB)) &&
							((indexedJetPt[1].first == bestmatchBbar) || (indexedJetPt[1].first == bestmatchB))) N_HighPtcorrect++;
						if (((indexedJetCSV[0].first == bestmatchBbar) || (indexedJetCSV[0].first == bestmatchB)) &&
							((indexedJetCSV[1].first == bestmatchBbar) || (indexedJetCSV[1].first == bestmatchB))) N_HighCSVcorrect++;
						if (((ClosestToMass1 == bestmatchBbar) || (ClosestToMass1 == bestmatchB)) &&
							((ClosestToMass2 == bestmatchBbar) || (ClosestToMass2 == bestmatchB))) N_MassMethod++;
						if (((ClosestTo125Mass1 == bestmatchBbar) || (ClosestTo125Mass1 == bestmatchB)) &&
							((ClosestTo125Mass2 == bestmatchBbar) || (ClosestTo125Mass2 == bestmatchB))) N_125MassMethod++;
						if (((ClosestTo114Mass1 == bestmatchBbar) || (ClosestTo114Mass1 == bestmatchB)) &&
							((ClosestTo114Mass2 == bestmatchBbar) || (ClosestTo114Mass2 == bestmatchB))) N_114MassMethod++;
						if (((ClosestTo109Mass1 == bestmatchBbar) || (ClosestTo109Mass1 == bestmatchB)) &&
							((ClosestTo109Mass2 == bestmatchBbar) || (ClosestTo109Mass2 == bestmatchB))) N_109MassMethod++;
						if (((CSVthenPt1 == bestmatchBbar) || (CSVthenPt1 == bestmatchB)) &&
							((CSVthenPt2 == bestmatchBbar) || (CSVthenPt2 == bestmatchB))) N_CSVthenPt++;
						if (((CSVthenPtM1 == bestmatchBbar) || (CSVthenPtM1 == bestmatchB)) &&
							((CSVthenPtM2 == bestmatchBbar) || (CSVthenPtM2 == bestmatchB))) N_CSVthenPtM++;
						if (((CSVthenPtTM1 == bestmatchBbar) || (CSVthenPtTM1 == bestmatchB)) &&
							((CSVthenPtTM2 == bestmatchBbar) || (CSVthenPtTM2 == bestmatchB))) N_CSVthenPtTM++;
						if (((CSVthenPtL1 == bestmatchBbar) || (CSVthenPtL1 == bestmatchB)) &&
							((CSVthenPtL2 == bestmatchBbar) || (CSVthenPtL2 == bestmatchB))) N_CSVthenPtL++;
						if (((CSVthenPtTL1 == bestmatchBbar) || (CSVthenPtTL1 == bestmatchB)) &&
							((CSVthenPtTL2 == bestmatchBbar) || (CSVthenPtTL2 == bestmatchB))) N_CSVthenPtTL++;
						if (((CSVthenPtML1 == bestmatchBbar) || (CSVthenPtML1 == bestmatchB)) &&
							((CSVthenPtML2 == bestmatchBbar) || (CSVthenPtML2 == bestmatchB))) N_CSVthenPtML++;
						if (((CSVthenPt51 == bestmatchBbar) || (CSVthenPt51 == bestmatchB)) &&
							((CSVthenPt52 == bestmatchBbar) || (CSVthenPt52 == bestmatchB))) N_CSVthenPt5++;
						if (((CSVthenPtTT1 == bestmatchBbar) || (CSVthenPtTT1 == bestmatchB)) &&
							((CSVthenPtTT2 == bestmatchBbar) || (CSVthenPtTT2 == bestmatchB))) N_CSVthenPtTT++;
						if (((CSVthenPtMM1 == bestmatchBbar) || (CSVthenPtMM1 == bestmatchB)) &&
							((CSVthenPtMM2 == bestmatchBbar) || (CSVthenPtMM2 == bestmatchB))) N_CSVthenPtMM++;
						if (((CSVthenPtLL1 == bestmatchBbar) || (CSVthenPtLL1 == bestmatchB)) &&
							((CSVthenPtLL2 == bestmatchBbar) || (CSVthenPtLL2 == bestmatchB))) N_CSVthenPtLL++;
						if (((CSVthenPt551 == bestmatchBbar) || (CSVthenPt551 == bestmatchB)) &&
							((CSVthenPt552 == bestmatchBbar) || (CSVthenPt552 == bestmatchB))) N_CSVthenPt55++;
						if (((CSVthenPt41 == bestmatchBbar) || (CSVthenPt41 == bestmatchB)) &&
							((CSVthenPt42 == bestmatchBbar) || (CSVthenPt42 == bestmatchB))) N_CSVthenPt4++;
						if (((CSVthenPt5L1 == bestmatchBbar) || (CSVthenPt5L1 == bestmatchB)) &&
							((CSVthenPt5L2 == bestmatchBbar) || (CSVthenPt5L2 == bestmatchB))) N_CSVthenPt5L++;
						if (((CSVthenPtM51 == bestmatchBbar) || (CSVthenPtM51 == bestmatchB)) &&
							((CSVthenPtM52 == bestmatchBbar) || (CSVthenPtM52 == bestmatchB))) N_CSVthenPtM5++;
						if (((CSVthenPtT41 == bestmatchBbar) || (CSVthenPtT41 == bestmatchB)) &&
							((CSVthenPtT42 == bestmatchBbar) || (CSVthenPtT42 == bestmatchB))) N_CSVthenPtT4++;
						if (((CSVthenPtM41 == bestmatchBbar) || (CSVthenPtM41 == bestmatchB)) &&
							((CSVthenPtM42 == bestmatchBbar) || (CSVthenPtM42 == bestmatchB))) N_CSVthenPtM4++;
						if (((CSVthenPt441 == bestmatchBbar) || (CSVthenPt441 == bestmatchB)) &&
							((CSVthenPt442 == bestmatchBbar) || (CSVthenPt442 == bestmatchB))) N_CSVthenPt44++;
						if (((CSVthenPt4L1 == bestmatchBbar) || (CSVthenPt4L1 == bestmatchB)) &&
							((CSVthenPt4L2 == bestmatchBbar) || (CSVthenPt4L2 == bestmatchB))) N_CSVthenPt4L++;
						if (((CSVthenPt541 == bestmatchBbar) || (CSVthenPt541 == bestmatchB)) &&
							((CSVthenPt542 == bestmatchBbar) || (CSVthenPt542 == bestmatchB))) N_CSVthenPt54++;
						if (((CSVthenPtT51 == bestmatchBbar) || (CSVthenPtT51 == bestmatchB)) &&
							((CSVthenPtT52 == bestmatchBbar) || (CSVthenPtT52 == bestmatchB))) N_CSVthenPtT5++;
						if (((CSVclosest1 == bestmatchBbar) || (CSVclosest1 == bestmatchB)) &&
							((CSVclosest2 == bestmatchBbar) || (CSVclosest2 == bestmatchB))) N_CSVclosest++;
						if (((PtClosestJet1 == bestmatchBbar) || (PtClosestJet1 == bestmatchB)) &&
							((PtClosestJet2 == bestmatchBbar) || (PtClosestJet2 == bestmatchB))) N_PtClosestJet++;
						if (((CSVHighPtT1 == bestmatchBbar) || (CSVHighPtT1 == bestmatchB)) &&
							((CSVHighPtT2 == bestmatchBbar) || (CSVHighPtT2 == bestmatchB))) N_CSVHighPtT++;
						if (((CSVHighPtM1 == bestmatchBbar) || (CSVHighPtM1 == bestmatchB)) &&
							((CSVHighPtM2 == bestmatchBbar) || (CSVHighPtM2 == bestmatchB))) N_CSVHighPtM++;
						if (((CSVHighPtL1 == bestmatchBbar) || (CSVHighPtL1 == bestmatchB)) &&
							((CSVHighPtL2 == bestmatchBbar) || (CSVHighPtL2 == bestmatchB))) N_CSVHighPtL++;
						if (((CSVHighPt51 == bestmatchBbar) || (CSVHighPt51 == bestmatchB)) &&
							((CSVHighPt52 == bestmatchBbar) || (CSVHighPt52 == bestmatchB))) N_CSVHighPt5++;
						if (((CSVHighPt41 == bestmatchBbar) || (CSVHighPt41 == bestmatchB)) &&
							((CSVHighPt42 == bestmatchBbar) || (CSVHighPt42 == bestmatchB))) N_CSVHighPt4++;
						if (((complicated1 == bestmatchBbar) || (complicated1 == bestmatchB)) &&
							((complicated2 == bestmatchBbar) || (complicated2 == bestmatchB))) N_complicated++;
						}//found a match
						
					}//gen Higgs mass requirement
					}//good genInfo requirement
					
					
					
					
					DetaJJ = sample.hJet_eta[indexedJetCSV[0].first]-sample.hJet_eta[indexedJetCSV[1].first];				
					Hmass = sample.H_mass;
					Hpt = sample.H_pt;
					ScalarSumHiggsJetPt = sample.hJet_pt[indexedJetCSV[1].first] + sample.hJet_pt[indexedJetCSV[0].first];
					ScalarSumJetPt = ScalarSumJetPt+ sample.hJet_pt[indexedJetCSV[1].first] + sample.hJet_pt[indexedJetCSV[0].first];
					//cout << "weight of the Jet Histograms: " << weight << endl;
					//JetDistributions("allEvts", weight);
					
					
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
					
					//	LeptonDistributions("allEvts", weight);
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
						//TH2FDistributions("allEvts", weight);
						
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
						
						//	EventShapeDistributions("allEvts", weight);
					}// two leptons and two jets requirement
					
					
					
					if ((sample.vLepton_pt[1] > 20 || sample.vLepton_pt[0] > 20) && (sample.vLepton_pt[1]>10 && sample.vLepton_pt[0]> 10)){
						//	if (sample.triggerFlags[39]||sample.triggerFlags[40]||sample.triggerFlags[41]){
						Ntrigger = Ntrigger + 1*PUweight2011;
						if ( jetPt[0] > 20 && jetPt[1] > 20 && 
							fabs(sample.vLepton_eta[0]) < 2.5 && fabs(sample.vLepton_eta[1]) < 2.4 && fabs(sample.hJet_eta[0]) < 2.5 &&
							fabs(sample.hJet_eta[1]) < 2.5 && sample.hJet_id[0]==1 && sample.hJet_id[1]==1 && sample.hbhe){
							Npreselect = Npreselect + 1*PUweight2011;
							if (isDATA) isdata = true;
							if(delRemu>0.4){
								N_EfakeCuts = N_EfakeCuts + 1*PUweight2011;
								if(CSV0 > 0.244 && DeltaPhiHV> 2 && Hmass>45 && Hmass < 150 && Emumass > 10 && Emumass < 70)FOM_tree->Fill();
								if (Emumass< 85 && Emumass>10) { 
									NMemu = NMemu + 1*PUweight2011;
									if (fabs(DphiZMET) < 1.250 && (Hmass>=80)&&(Hmass<=150) && CSV0<0.244 && Naj < 2 ){	
										if ((sample.vLepton_pt[1] > 25 || sample.vLepton_pt[0] > 25) && (sample.vLepton_pt[1]>15 && sample.vLepton_pt[0]> 15)){
											N_LFCR = N_LFCR + 1*PUweight2011;
										}//HLT plateau cuts
									}// LF Control Region								
									if(CSV0>0.244){
										N_CSV0 = N_CSV0 + 1*PUweight2011;
										if (fabs(DphiZMET) > 1.250 && Hmass>=80&&Hmass<=150){
											if (CSV1>0.244 && EvntShpAplanarity > 0.1 ){
												if ((sample.vLepton_pt[1] > 25 || sample.vLepton_pt[0] > 25) && (sample.vLepton_pt[1]>15 && sample.vLepton_pt[0]> 15)){
													N_TopCR = N_TopCR + 1*PUweight2011;
												}//HLT emulation											
											}//CHF cut
										}// Top CR
										if (fabs(DphiZMET) > 1.250 && (Hmass>=80)&&(Hmass<=150) ){
											if (CSV0 > 0.5 && EvntShpAplanarity < 0.1 && CSV1 < 0.898  && Nab < 1){
												if ((sample.vLepton_pt[1] > 25 || sample.vLepton_pt[0] > 25) && (sample.vLepton_pt[1]>15 && sample.vLepton_pt[0]> 15)){
													N_SingleTopCR = N_SingleTopCR + 1*PUweight2011;
												}//HLT emulation											
											}//CHF cut
										}// Single Top CR																		
										if (fabs(DphiZMET) < 1.250){
											N_DphiZMET = N_DphiZMET + 1*PUweight2011;
											if((Hmass>=80)&&(Hmass<=150)){
												N_Mjj = N_Mjj + 1*PUweight2011; 
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
											if ((Hmass<80 || Hmass>150) && Hmass < 250 && CSV0>0.898 && Nab < 1 ){
												N_HFCR = N_HFCR + 1*PUweight2011;
											}// HF Control Region
											if((Hmass>=80)&&(Hmass<=150)){
												if(jetCHF[0] > 0.15) {
													N_jetCHF0 = N_jetCHF0 + 1*PUweight2011;
													if (Naj<2){
														N_Naj = N_Naj + 1*PUweight2011;
													}	//Naj
												}//jetCHF[0]
											}//Mjj							
										}// Delta Phi Z, MET (Z is emu vectors only)
									}//CSV
								}//emu mass requirement
							}//EleFakeCuts
						} else {
							FailedJetID = FailedJetID + 1*PUweight2011;
						}//Jet ID and eta requirement
					}//end trigger emulation
					isdata = false;
					//				EventDistributions("allEvts", weight);
				}//end requirement Zemu event
			}//end isM50sample genZpt cut
			
		}// end njet requirement
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
	std::cout << "Mjj: " << N_Mjj << endl;
	std::cout << "CHFb0: " << N_jetCHF0 << endl;
	std::cout << "Naj: " << N_Naj << endl;
	
	
	
	std::cout << endl << endl;
	std::cout << "Top CR: " << N_TopCR << endl;
	std::cout << "Single Top CR: " << N_SingleTopCR << endl;
	std::cout << "LF CR: " << N_LFCR << endl;
	std::cout << "HF CR: " << N_HFCR << endl;
	
	
	std::cout << endl << endl;
	std::cout << "Number of Events that failed JetID " << FailedJetID << endl;
	std::cout << "Number of Events Vtype 5 " << NVtype << endl;
	std::cout << "Number of events with good genInfo " << N_count << endl;
	std::cout << "Number of events with reasonble Mass" << ReasonbleMass << endl;
	
	std::cout << "Number of Events not Higgs  with highest Pt " << N_NonHiggsPt << endl;
	std::cout << "Number of Events not Higgs  with highest CSV " << N_NonHiggsCSV << endl;
	std::cout << "Number of Events non hJet matched " << N_NonhJetMatched << endl;
	
	std::cout << "same jet highest CSV and Pt " << highPtCSV << endl;
	std::cout << "same jet second CSV and Pt " << secondPtCSV << endl;
	std::cout << "Two Higest Pt jets match genInfo " << N_HighPtcorrect << endl;
	std::cout << "Two Higest CSV jets match genInfo " << N_HighCSVcorrect << endl;
	std::cout << "Mass 103 matched jets match genInfo " << N_MassMethod << endl;
	std::cout << "Mass 125 matched jets match genInfo " << N_125MassMethod << endl;
	std::cout << "Mass 114 matched jets match genInfo " << N_114MassMethod << endl;
	std::cout << "Mass 109 matched jets match genInfo " << N_109MassMethod << endl;
	std::cout << "CSVT then Pt jets match genInfo " << N_CSVthenPt << endl;
	std::cout << "CSVM then Pt jets match genInfo " << N_CSVthenPtM << endl;
	std::cout << "CSVL then Pt jets match genInfo " << N_CSVthenPtL << endl;
	std::cout << "CSV4 then Pt jets match genInfo " << N_CSVthenPt4 << endl;
	std::cout << "CSV5 then Pt jets match genInfo " << N_CSVthenPt5 << endl;
	std::cout << "CSVTT then Pt jets match genInfo " << N_CSVthenPtTT << endl;
	std::cout << "CSVMM then Pt jets match genInfo " << N_CSVthenPtMM << endl;
	std::cout << "CSVLL then Pt jets match genInfo " << N_CSVthenPtLL << endl;
	std::cout << "CSV55 then Pt jets match genInfo " << N_CSVthenPt55 << endl;
	std::cout << "CSV44 then Pt jets match genInfo " << N_CSVthenPt44 << endl;
	std::cout << "CSVTM then Pt jets match genInfo " << N_CSVthenPtTM << endl;
	std::cout << "CSVTL then Pt jets match genInfo " << N_CSVthenPtTL << endl;
	std::cout << "CSVT4 then Pt jets match genInfo " << N_CSVthenPtT4 << endl;
	std::cout << "CSVT5 then Pt jets match genInfo " << N_CSVthenPtT5 << endl;
	std::cout << "CSVML then Pt jets match genInfo " << N_CSVthenPtML << endl;
	std::cout << "CSVM4 then Pt jets match genInfo " << N_CSVthenPtM4 << endl;
	std::cout << "CSVM5 then Pt jets match genInfo " << N_CSVthenPtM5 << endl;
	std::cout << "CSV54 then Pt jets match genInfo " << N_CSVthenPt54 << endl;
	std::cout << "CSV5L then Pt jets match genInfo " << N_CSVthenPt5L << endl;
	std::cout << "CSV4L then Pt jets match genInfo " << N_CSVthenPt4L << endl;
	std::cout << "High Pt CSVT match genInfo " << N_CSVHighPtT << endl;
	std::cout << "High Pt CSVM match genInfo " << N_CSVHighPtM << endl;
	std::cout << "High Pt CSVL match genInfo " << N_CSVHighPtL << endl;
	std::cout << "High Pt CSV4 match genInfo " << N_CSVHighPt4 << endl;
	std::cout << "High Pt CSV5 match genInfo " << N_CSVHighPt5 << endl;
	std::cout << "High CSV then closest jet match genInfo " << N_CSVclosest << endl;
	std::cout << "High Pt then closest jet match genInfo " << N_PtClosestJet << endl;
	std::cout << "complicated matches genInfo " << N_complicated << endl;
	
	
	
	// Here one can create canvas, draw and print the histogram.
	TCanvas c1("c1","c1");
	c1.cd();
	c1.SetFillColor(kWhite);
	string suffixps = ".gif";
	
	c1.Clear(); // don't create a new canvas
	hSecondCSVHiggs.Draw();
	c1.Print((directory+"/SecondCSVHiggs"+suffixps).c_str());
	
	c1.Clear(); // don't create a new canvas
	hLeadingCSVHiggs.Draw();
	c1.Print((directory+"/LeadingCSVHiggs"+suffixps).c_str());
	
	c1.Clear(); // don't create a new canvas
	hClosestMatch.Draw();
	c1.Print((directory+"/ClosestMatch"+suffixps).c_str());
	
	c1.Clear(); // don't create a new canvas
	hSecondPtHiggs.Draw();
	c1.Print((directory+"/SecondPtHiggs"+suffixps).c_str());
	
	c1.Clear(); // don't create a new canvas
	hLeadingPtHiggs.Draw();
	c1.Print((directory+"/LeadingPtHiggs"+suffixps).c_str());
	
	
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
	
	for (int i=0; i != nJets && i < 5; ++i) {
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
	
	for (int i=0; i != nLeptons && i < 5; ++i) {
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
