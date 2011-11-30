/**
    @brief Example of an analysis code to read tree ntuple
    and create histogram
    Usage:

    \verbatim
    mkdir [directoryName]
    Zbb_tree <inputfile> [outputfilename] [directoryName]
    \endverbatim

    @param inputfile Either a ROOT file or an ASCII file containing list of
    ROOT files.

    @param outputfile Name of the ROOT file which contains the histogram.
    Defaulted to 'output.root'

    @author Rachel Wilken <rachel.wilken@cern.ch>

    @date Thu Nov 17 2011

 */


#include <UserCode/wilken/interface/treeReader.h>
#include <UserCode/wilken/interface/GenericTool.h>
#include "PhysicsTools/CandUtils/interface/EventShapeVariables.h"

#include <string>
#include <iostream>
#include <vector>
#include <algorithm>
#include <map>

#include "TLorentzVector.h"
#include "TVector3.h"
#include "TMath.h"

#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TTree.h>
#include <TCanvas.h>
#include "TString.h"
#include <TStyle.h>

using namespace std;

double deltaphi(double phi1, double phi2) {
	double result = phi1 - phi2;
	const double mypi  = 3.141592653589;
    while (result > mypi) result -= 2*mypi;
    while (result <= -mypi) result += 2*mypi;
    return result;	
}// end deltaphi function


int main(int argc, char** argv) {

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
    treeReader sample(ifilename,
                                     std::string("tree"));

    // If the input file(s) doesn't contain any event, exit.
    if (!(sample.readEvent(0))) {
        return 0;
    }

    // Set the sumw2 option for all histogram
    TH1::SetDefaultSumw2(true);

    // Create the output ROOT file
    TFile ofile(ofilename.c_str(),"RECREATE");
    ofile.cd();

    // Check how many events and file(s) are analyzed
    long int nEvent     = sample.nEvent();
    size_t   nFile      = sample.nFile();

    std::cout << "There are " << nEvent << " events in the chain."
              << std::endl ;
    std::cout << "There are " << nFile  << " files in the chain."
              << std::endl ;
    std::cout << "The tree name  is: " << sample.treeName()  << std::endl;
    std::cout << "The tree title is: " << sample.treeTitle() << std::endl;
	
TTree *TMVA_tree = new TTree("TMVA_tree","Tree for TMVA input");
	int nJets, nSV;
	float CSV1, CSV2, Zmass, Hmass, DeltaPhiHV, Hpt, Zpt; 
	float mu1pt, Ht, EtaStandDev, UnweightedEta, EvntShpCircularity;
	float alpha_j, qtb1,DphiJJ, Trigweight; 
	TMVA_tree->Branch("nJets",&nJets, "nJets/I");
	TMVA_tree->Branch("CSV1",&CSV1, "CSV1/F");
	TMVA_tree->Branch("CSV2",&CSV2, "CSV2/F");
	TMVA_tree->Branch("Zmass",&Zmass, "Zmass/F");
	TMVA_tree->Branch("Hmass",&Hmass, "Hmass/F");
	TMVA_tree->Branch("DeltaPhiHV",&DeltaPhiHV, "DeltaPhiHV/F");
	TMVA_tree->Branch("Hpt",&Hpt, "Hpt/F");
	TMVA_tree->Branch("Zpt",&Zpt, "Zpt/F");
	TMVA_tree->Branch("mu1pt",&mu1pt, "mu1pt/F");
	TMVA_tree->Branch("Ht",&Ht, "Ht/F");
	TMVA_tree->Branch("EtaStandDev",&EtaStandDev, "EtaStandDev/F");
	TMVA_tree->Branch("UnweightedEta",&UnweightedEta, "UnweightedEta/F");
	TMVA_tree->Branch("EvntShpCircularity",&EvntShpCircularity, "EvntShpCircularity/F");
	TMVA_tree->Branch("alpha_j",&alpha_j, "alpha_j/F");
	TMVA_tree->Branch("qtb1",&qtb1, "qtb1/F");
	TMVA_tree->Branch("nSV",&nSV, "nSV/I");
	TMVA_tree->Branch("Trigweight",&Trigweight, "Trigweight/F");

						  

bool debug = false;	


    // Here one can declare histograms
    // In compiled C++, it possible not to use pointer type (TH1F*) and
    // use the non-pointer type (TH1F) of ROOT object.

	//Declare Histograms

    TH1F hnJets ("hnJets",  "Number of Good Jets", 11, -0.5, 10.5);
    TH1F hnMuons  ("hnMuons","Number of Good Muons", 6, -0.5, 5.5);
    TH1F hnSV  ("hnSV","Number of Secondary Verticies", 6, -0.5, 5.5);
    TH1F hnPV  ("hnPV","Number of Primary Verticies", 21, -0.5, 20.5);
    TH1F hallhJet_pt  ("hallhJet_pt","Pt of all jets in event",		150, 0.0, 150);
    TH1F hMET		("hMET","Missing Et",		150, 0.0, 1200);
	TH1F hEtab1		("hEtab1","Eta of b jet with highest CSV", 75, -3, 3.5);
	TH1F hEtab2		("hEtab2","Eta of b jet with second highest CSV", 75, -3, 3.5);
	TH1F hPhib1		("hPhib1","Phi of b jet with highest CSV", 75, -3.25, 4.5);
	TH1F hPhib2		("hPhib2","Phi of b jet with second highest CSV", 75, -3.25, 4.5);
    TH1F hPtb1		("hPtb1","Pt b jet with highest CSV",		150, 0.0, 150);
    TH1F hPtb2		("hPtb2","Pt b jet with second highest CSV", 150, 0.0, 150);
	TH1F hCHFb1	        ("hCHFb1","charged Hadron Energy Fraction b1", 101, 0.0, 1.2);
	TH1F hCHFb2	        ("hCHFb2","charged Hadron Energy Fraction b2", 101, 0.0, 1.2);
	TH1F hPtjj		("hPtjj","Pt of two b jets with highest CSV", 150, 0.0, 300);
	TH1F hPtmumu	("hPtmumu","Pt of two muons with highest pt", 150, 0.0, 300);
	TH1F hPtbalZH	("hPtbalZH","Pt balance of Z and H", 101, -75.0, 75);
	TH1F hPtmu1		("hPtmu1","Pt of muon with highest pt", 75, 0.0, 200);
	TH1F hPtmu2		("hPtmu2","Pt of muon with second highest pt", 75, 0.0, 150);
	TH1F hEtamu1		("hEtamu1","Eta of muon with highest pt", 75, -3, 3.5);
	TH1F hEtamu2		("hEtamu2","Eta of muon with second highest pt", 75, -3, 3.5);
	TH1F hPhimu1		("hPhimu1","Phi of muon with highest pt", 75, -3.25, 4.5);
	TH1F hPhimu2		("hPhimu2","Phi of muon with second highest pt", 75, -3.25, 4.5);
	TH1F hPFRelIsomu1		("hPFRelIsomu1","PF Rel Iso of muon with highest pt", 101, 0.0001, 0.2);
	TH1F hPFRelIsomu2		("hPFRelIsomu2","Pf Rel Iso of muon with second highest pt", 101, 0.0001, 0.2);	
	TH1F hCSV1		("hCSV1","Jet with highest CSV",			50, -0.1, 1.5);
	TH1F hCSV2		("hCSV2","Jet with second highest CSV",		50, -0.1, 1.5);
	TH1F hdphiVH	("hdphiVH","Delta phi between Z and Higgs", 50, -0.1, 4.5);
	TH1F hdetaJJ	("hdetaJJ","Delta eta between two jets", 101, -3, 3);
	TH1F hNaj		("hNaj",  "Number of Additional Jets",		13, -2.5, 10.5);
	TH1F hMjj		("hMjj",  "Invariant Mass of two Jets",		75, 0, 200);
	TH1F hMmumu		("hMmumu",  "Invariant Mass of two muons",	75, 0, 200);
	TH1F hCutFlow	("hCutFlow",  "Selection",					8, 0, 8);
	TH1F hRMSeta	("hRMSeta",  "RMS Eta",		71, 0, 3);	
	TH1F hStaDeveta	("hStaDeveta",  "Standard Deviation Eta",		71, 0, 3);	
	TH1F hUnweightedEta	("hUnweightedEta",  "Unweighted Eta",		101, 0, 10);		
	TH2F hDphiDetajj("hDphiDetajj", "#Delta#phi vs #Delta#eta JJ", 75, -5, 5, 75, -3.5, 3.5);
	
//	Lorentz vectors:
	TH1F hdphiVH_vect	("hdphiVH_vect","Delta phi between Z and Higgs from Lorentz Vector", 50, -3.5, 4.5);
	TH1F hdphiJJ_vect	("hdphiJJ_vect","Delta phi between two jets from Lorentz Vector", 71, -3.5, 4);
	TH1F hPtjj_vect		("hPtjj_vect","Pt of two b jets with highest CSV from Lorentz Vector", 150, 0.0, 300);
	TH1F hPtmumu_vect	("hPtmumu_vect","Pt of two muons with highest pt from Lorentz Vector", 150, 0.0, 300);
	TH1F hHt	        ("hHt","scalar sum of pt of four particles", 101, 0.0, 500);
	TH1F hCentrality ("hCentrality", "Centrality", 71, 0.0, 1.0);
	TH1F hEventPt ("hEventPt", "Pt of HV system", 50, 0.0, 100);
	TH1F hAngle ("hAngle", "Angle between H and Z", 101, -0.1, 4);
	TH2F hqtvsalphaJJ ("hqtvsalphaJJ", "Armenteros-Podolansky Plot Jets", 70, -2, 2, 75, 0, 75);
	TH2F hqtvsalphaZ ("hqtvsalphaZ", "Armenteros-Podolansky Plot Z", 70, -2, 2, 75, 0, 75);
	
//EventShapeVariables
	TH1F hSphericity_EvtShp ("hSphericity_EvtShp", "EventShapeVariables sphericity", 71, 0.0, 1);
	TH1F hAplanarity_EvtShp ("hAplanarity_EvtShp", "EventShapeVariables Aplanarity", 50, 0.0, .5);
	TH1F hCircularity_EvtShp ("hCircularity_EvtShp", "EventShapeVariables circularity", 45, 0.0, 1.2);
	TH1F hIsotropy_EvtShp ("hIsotropy_EvtShp", "EventShapeVariables isotropy", 50, 0.0, 1.2);
	
//AfterCuts
	TH1F hMjj_aftercuts		("hMjj_aftercuts",  "Invariant Mass of two Jets After Cuts",		55, 0, 200);
	TH1F hEtab1_aftercuts		("hEtab1_aftercuts","Eta of b jet with highest CSV After Cuts", 55, -3, 3.5);
	TH1F hEtab2_aftercuts		("hEtab2_aftercuts","Eta of b jet with second highest CSV After Cuts", 55, -3, 3.5);
	TH1F hPhib1_aftercuts		("hPhib1_aftercuts","Phi of b jet with highest CSV After Cuts", 55, -3.25, 4.5);
	TH1F hPhib2_aftercuts		("hPhib2_aftercuts","Phi of b jet with second highest CSV After Cuts", 55, -3.25, 4.5);
    TH1F hPtb1_aftercuts		("hPtb1_aftercuts","Pt b jet with highest CSV After Cuts",		75, 0.0, 150);
    TH1F hPtb2_aftercuts		("hPtb2_aftercuts","Pt b jet with second highest CSV After Cuts", 75, 0.0, 150);
	TH1F hCHFb1_aftercuts	        ("hCHFb1_aftercuts","charged Hadron Energy Fraction b1 After Cuts", 101, 0.0, 1.2);
	TH1F hCHFb2_aftercuts	        ("hCHFb2_aftercuts","charged Hadron Energy Fraction b2 After Cuts", 101, 0.0, 1.2);
	TH1F hPtjj_aftercuts		("hPtjj_aftercuts","Pt of two b jets with highest CSV After Cuts", 150, 0.0, 300);
	TH1F hPtmumu_aftercuts	("hPtmumu_aftercuts","Pt of two muons with highest pt After Cuts", 150, 0.0, 300);
	TH1F hPtbalZH_aftercuts	("hPtbalZH_aftercuts","Pt balance of Z and H After Cuts", 71, -75.0, 75);
	TH1F hPtmu1_aftercuts		("hPtmu1_aftercuts","Pt of muon with highest pt After Cuts", 75, 0.0, 200);
	TH1F hPtmu2_aftercuts		("hPtmu2_aftercuts","Pt of muon with second highest pt After Cuts", 75, 0.0, 150);
	TH1F hEtamu1_aftercuts		("hEtamu1_aftercuts","Eta of muon with highest pt After Cuts", 55, -3, 3.5);
	TH1F hEtamu2_aftercuts		("hEtamu2_aftercuts","Eta of muon with second highest pt After Cuts", 55, -3, 3.5);
	TH1F hPhimu1_aftercuts		("hPhimu1_aftercuts","Phi of muon with highest pt After Cuts", 55, -3.25, 4.5);
	TH1F hPhimu2_aftercuts		("hPhimu2_aftercuts","Phi of muon with second highest pt After Cuts", 55, -3.25, 4.5);
	TH1F hPFRelIsomu1_aftercuts		("hPFRelIsomu1_aftercuts","PF Rel Iso of muon with highest pt After Cuts", 71, 0.0001, 0.2);
	TH1F hPFRelIsomu2_aftercuts		("hPFRelIsomu2_aftercuts","Pf Rel Iso of muon with second highest pt After Cuts", 71, 0.0001, 0.2);	
	TH1F hCSV1_aftercuts		("hCSV1_aftercuts","Jet with highest CSV After Cuts",			50, -0.1, 1.5);
	TH1F hCSV2_aftercuts		("hCSV2_aftercuts","Jet with second highest CSV After Cuts",		50, -0.1, 1.5);
	TH1F hdphiVH_aftercuts	("hdphiVH_aftercuts","Delta phi between Z and Higgs After Cuts", 50, -0.1, 4.5);
	TH1F hdetaJJ_aftercuts	("hdetaJJ_aftercuts","Delta eta between two jets After Cuts", 71, -3, 3);
	TH1F hNaj_aftercuts		("hNaj_aftercuts",  "Number of Additional Jets After Cuts",		13, -2.5, 10.5);
	TH1F hMmumu_aftercuts		("hMmumu_aftercuts",  "Invariant Mass of two muons After Cuts",	75, 0, 200);
	TH1F hRMSeta_aftercuts	("hRMSeta_aftercuts",  "RMS Eta After Cuts",		71, 0, 3);	
	TH1F hStaDeveta_aftercuts	("hStaDeveta_aftercuts",  "Standard Deviation Eta After Cuts",		71, 0, 3);	
	TH1F hUnweightedEta_aftercuts	("hUnweightedEta_aftercuts",  "Unweighted Eta After Cuts",		101, 0, 10);		
	TH2F hDphiDetajj_aftercuts ("hDphiDetajj_aftercuts", "#Delta#phi vs #Delta#eta JJ After Cuts", 55, -5, 5, 55, -3.5, 3.5);
	TH1F hdphiVH_vect_aftercuts	("hdphiVH_vect_aftercuts","Delta phi between Z and Higgs from Lorentz Vector After Cuts", 50, -3.5, 4.5);
	TH1F hdphiJJ_vect_aftercuts	("hdphiJJ_vect_aftercuts","Delta phi between two jets from Lorentz Vector After Cuts", 71, -3.5, 4);
	TH1F hPtjj_vect_aftercuts		("hPtjj_vect_aftercuts","Pt of two b jets with highest CSV from Lorentz Vector After Cuts", 150, 0.0, 300);
	TH1F hPtmumu_vect_aftercuts	("hPtmumu_vect_aftercuts","Pt of two muons with highest pt from Lorentz Vector After Cuts", 150, 0.0, 300);
	TH1F hHt_aftercuts	        ("hHt_aftercuts","scalar sum of pt of four particles After Cuts", 71, 0.0, 500);
	TH1F hCentrality_aftercuts ("hCentrality_aftercuts", "Centrality After Cuts", 71, 0.0, 1.0);
	TH1F hEventPt_aftercuts ("hEventPt_aftercuts", "Pt of HV system After Cuts", 50, 0.0, 100);
	TH1F hAngle_aftercuts ("hAngle_aftercuts", "Angle between H and Z After Cuts", 71, -0.1, 4);
	TH2F hqtvsalphaJJ_aftercuts ("hqtvsalphaJJ_aftercuts", "Armenteros-Podolansky Plot Jets After Cuts", 70, -2, 2, 75, 0, 75);
	TH2F hqtvsalphaZ_aftercuts ("hqtvsalphaZ_aftercuts", "Armenteros-Podolansky Plot Z After Cuts", 70, -2, 2, 75, 0, 75);
	TH1F hSphericity_EvtShp_aftercuts ("hSphericity_EvtShp_aftercuts", "EventShapeVariables sphericity After Cuts", 71, 0.0, 1);
	TH1F hAplanarity_EvtShp_aftercuts ("hAplanarity_EvtShp_aftercuts", "EventShapeVariables Aplanarity After Cuts", 50, 0.0, .5);
	TH1F hCircularity_EvtShp_aftercuts ("hCircularity_EvtShp_aftercuts", "EventShapeVariables circularity After Cuts", 45, 0.0, 1.2);
	TH1F hIsotropy_EvtShp_aftercuts ("hIsotropy_EvtShp_aftercuts", "EventShapeVariables isotropy After Cuts", 50, 0.0, 1.2);
	
	//not pt cuts
	TH1F hMjj_NoPtCut		("hMjj_NoPtCut",  "Invariant Mass of two Jets After Cuts No Pt Cuts",		55, 0, 200);
	TH1F hMmumu_NoPtCut		("hMmumu_NoPtCut",  "Invariant Mass of two muonsAfter Cuts No Pt Cuts",	75, 0, 200);
	TH1F hPtmumu_NoPtCut	("hPtmumu_NoPtCut","Pt of two muons with highest ptAfter Cuts No Pt Cuts", 150, 0.0, 300);
	TH1F hPtjj_NoPtCut		("hPtjj_NoPtCut","Pt of two b jets with highest CSVAfter Cuts No Pt Cuts", 150, 0.0, 300);
	TH1F hPtbalZH_NoPtCut	("hPtbalZH_NoPtCut","Pt balance of Z and HAfter Cuts No Pt Cuts", 71, -75.0, 75);
	TH2F hqtvsalphaJJ_NoPtCut ("hqtvsalphaJJ_NoPtCut", "Armenteros-Podolansky Plot JetsAfter Cuts No Pt Cuts", 70, -2, 2, 75, 0, 75);
	TH2F hqtvsalphaZ_NoPtCut ("hqtvsalphaZ_NoPtCut", "Armenteros-Podolansky Plot ZAfter Cuts No Pt Cuts", 70, -2, 2, 75, 0, 75);
	TH1F hHt_NoPtCut	        ("hHt_NoPtCut","scalar sum of pt of four particlesAfter Cuts No Pt Cuts", 71, 0.0, 500);
	TH1F hCentrality_NoPtCut ("hCentrality_NoPtCut", "CentralityAfter Cuts No Pt Cuts", 71, 0.0, 1.0);
	TH1F hEventPt_NoPtCut ("hEventPt_NoPtCut", "Pt of HV system After Cuts No Pt Cuts", 50, 0.0, 100);
	TH1F hPtjj_vect_NoPtCut		("hPtjj_vect_NoPtCut","Pt of two b jets with highest CSV from Lorentz VectorAfter Cuts No Pt Cuts", 150, 0.0, 300);
	TH1F hPtmumu_vect_NoPtCut	("hPtmumu_vect_NoPtCut","Pt of two muons with highest pt from Lorentz VectorAfter Cuts No Pt Cuts", 150, 0.0, 300);
	

//no angle cuts	
	TH1F hAngle_NoDphiCut ("hAngle_NoDphiCut", "Angle between H and Z After Cuts No DeltaPhi Cut", 71, -0.1, 4);
	TH1F hPhimu1_NoDphiCut		("hPhimu1_NoDphiCut","Phi of muon with highest pt After Cuts No DeltaPhi Cut", 75, -3.25, 4.5);
	TH1F hPhimu2_NoDphiCut		("hPhimu2_NoDphiCut","Phi of muon with second highest pt After Cuts No DeltaPhi Cut", 75, -3.25, 4.5);
	TH1F hPhib1_NoDphiCut		("hPhib1_NoDphiCut","Phi of b jet with highest CSV After Cuts No DeltaPhi Cut", 75, -3.25, 4.5);
	TH1F hPhib2_NoDphiCut		("hPhib2_NoDphiCut","Phi of b jet with second highest CSV After Cuts No DeltaPhi Cut", 75, -3.25, 4.5);
	TH1F hdphiJJ_vect_NoDphiCut	("hdphiJJ_vect_NoDphiCut","Delta phi between two jets from Lorentz Vector After Cuts No DeltaPhi Cut", 71, -3.5, 4);
	TH2F hDphiDetajj_NoDphiCut ("hDphiDetajj_NoDphiCut", "#Delta#phi vs #Delta#eta JJ After Cuts No DeltaPhi Cut", 55, -5, 5, 55, -3.5, 3.5);

	
	
if (debug) std::cout << "all histograms declared " << std::endl;


    // Loop over all events.
    // For other methods to access event/navigate through the sample,
    // see the documentation of RootTreeReader class.
	int event =0, Npreselect =0, NpT_jj =0, NpT_Z =0, NCSV1 =0, NCSV2 =0;
	int NdPhi =0, N_Naj =0, N_Mjj =0, Ntree =0;
//	double Zmass = 91.1976;


    do {
event++;
		if (sample.Vtype == 0){
	if (debug) std::cout << "entered event loop " << event << std::endl;
	std::vector< std::pair<size_t,double> > indexedJetPt;
        std::vector<size_t> PtSortedJetIndex;
	std::vector< std::pair<size_t,double> > indexedJetCSV;
        std::vector<size_t> CSVSortedJetIndex;
  std::vector< std::pair<size_t,double> > indexedPt;
        std::vector<size_t> PtSortedIndex;

        // Analysis loop.
        // One can access the ntuple leaves directly from sample object

nJets =0, nSV =0;
CSV1 = -1.0, CSV2 = -1.0, Zmass = -99.99, Hmass = -99.99, DeltaPhiHV = -99.99, Hpt = -99.99, Zpt = -99.99;
mu1pt = -99.99, Ht = -99.99, EtaStandDev = -99.99, UnweightedEta = -99.99, EvntShpCircularity = -99.99;
alpha_j = -99.99, qtb1 = 0.0, DphiJJ = -99.99, Trigweight = 1.0;
		double RMS_eta = -99.99;
		double StandDevEta[4];
		double AverageEta = -99.99;
		TLorentzVector FirstJet;
		TLorentzVector SecondJet;
		TLorentzVector FirstMuon;
		TLorentzVector SecondMuon;
		TLorentzVector Higgs;
		TLorentzVector ZBoson;
		TLorentzVector HVsystem;
		double alpha_mu = -99.99;
		double plmu1 = 0.0, plmu2 = 0.0 , qtmu1 = 0.0, qtmu2 = 0.0;
		double plb1 = 0.0, plb2 = 0.0 , qtb2 = 0.0;
		double EvntShpAplanarity = -99.99, EvntShpSphericity = -99.99, EvntShpIsotropy = -99.99;



for (int k=0;k<sample.nhJets;k++){
	if (debug) cout << "for "<< k << "th pass" << endl;
	indexedJetPt.push_back(std::pair<size_t,double>(k,(double) sample.hJet_pt[k]));
	indexedPt.push_back(std::pair<size_t,double>(k,(double) sample.hJet_pt[k]));
	hallhJet_pt.Fill(sample.hJet_pt[k]); 
	nJets++;
	indexedJetCSV.push_back(std::pair<size_t,double>(k,(double) sample.hJet_csv[k]));
	if (debug)	cout << "CSV discriminator: " << sample.hJet_csv[k] << endl;

}// end k for loop
		for (int a=0;a<sample.naJets;a++){
			hallhJet_pt.Fill(sample.aJet_pt[a]);
			nJets++; 
		}
		indexedPt.push_back(std::pair<size_t,double>(2,(double) sample.vLepton_pt[0]));
		indexedPt.push_back(std::pair<size_t,double>(3,(double) sample.vLepton_pt[1]));

		
if (event == 10) for (size_t i = 0 ; (i != indexedPt.size()) ; ++i) {cout << "Unsorted pt of objects: " << indexedPt[i].second << endl; }

std::sort(indexedJetPt.begin(),indexedJetPt.end(),::IndexedQuantityGreaterThan<double>);
		std::sort(indexedJetCSV.begin(),indexedJetCSV.end(),::IndexedQuantityGreaterThan<double>);
		std::sort(indexedPt.begin(),indexedPt.end(),::IndexedQuantityGreaterThan<double>);

for (size_t i = 0 ; (i != indexedJetPt.size()) ; ++i) {   PtSortedJetIndex.push_back(indexedJetPt[i].first);        }
		for (size_t i = 0 ; (i != indexedJetCSV.size()) ; ++i) {   CSVSortedJetIndex.push_back(indexedJetCSV[i].first);        }
		for (size_t i = 0 ; (i != indexedPt.size()) ; ++i) {   PtSortedIndex.push_back(indexedPt[i].first);        }

if (event == 10) for (size_t i = 0 ; (i != indexedPt.size()) ; ++i) {cout << "Sorted pt of objects: " << indexedPt[i].second << endl; }
			CSV1 = sample.hJet_csv[indexedJetCSV[0].first];
			CSV2 = sample.hJet_csv[indexedJetCSV[1].first];		
double dummy = -99.99;		
		if (sample.nhJets > 1) { 
		if (CSV2 > -1) { 	
		hPtb1.Fill(dummy);
			hPtb1.Fill(sample.hJet_pt[indexedJetCSV[0].first]);
			hPtb2.Fill(sample.hJet_pt[indexedJetCSV[1].first]);
			hCSV1.Fill(CSV1);
			hCSV2.Fill(CSV2);
			hEtab1.Fill(sample.hJet_eta[indexedJetCSV[0].first]);
			hEtab2.Fill(sample.hJet_eta[indexedJetCSV[1].first]);
			hPhib1.Fill(sample.hJet_phi[indexedJetCSV[0].first]);
			hPhib2.Fill(sample.hJet_phi[indexedJetCSV[1].first]);
			hCHFb1.Fill(sample.hJet_chf[indexedJetCSV[0].first]);
			hCHFb2.Fill(sample.hJet_chf[indexedJetCSV[1].first]);
			hdetaJJ.Fill(sample.hJet_eta[indexedJetCSV[0].first]-sample.hJet_eta[indexedJetCSV[1].first]);
Hmass = sample.H_mass;
			hMjj.Fill(Hmass);
Hpt = sample.H_pt;
		    hPtjj.Fill(Hpt);
		}//end csv non -1 requirement
		}// end njet requirement			
Zmass = sample.V_mass;
Zpt = sample.V_pt;
mu1pt = sample.vLepton_pt[0];
nSV = sample.nSvs;
Trigweight = sample.weightTrig;
	if (debug) cout << "halfway through filling histos" << endl;
		if (sample.nhJets > -1) { hnJets.Fill(nJets); }
		if (sample.nvlep > -1) { hnMuons.Fill(sample.nvlep); }
		hnSV.Fill(nSV);
		hnPV.Fill(sample.nPVs);
		hNaj.Fill(sample.naJets);
		hMET.Fill(sample.MET_sumet);
if ((sample.vLepton_type[0] == 13) && (sample.vLepton_type[1] == 13)){
		if (sample.nvlep > 1) {
		hMmumu.Fill(Zmass);
		hPtmumu.Fill(Zpt);
			hPtmu1.Fill(mu1pt);
			hPtmu2.Fill(sample.vLepton_pt[1]);		
			hEtamu1.Fill(sample.vLepton_eta[0]);
			hEtamu2.Fill(sample.vLepton_eta[1]);
			hPhimu1.Fill(sample.vLepton_phi[0]);
			hPhimu2.Fill(sample.vLepton_phi[1]);
			hPFRelIsomu1.Fill(sample.vLepton_pfCombRelIso[0]);
			hPFRelIsomu2.Fill(sample.vLepton_pfCombRelIso[1]);
		}
			if( (sample.nvlep > 1) && (sample.nhJets > 1)) {
			RMS_eta = sample.vLepton_eta[1]*sample.vLepton_eta[1]+sample.vLepton_eta[0]*sample.vLepton_eta[0]+sample.hJet_eta[indexedJetCSV[0].first]*sample.hJet_eta[indexedJetCSV[0].first]+sample.hJet_eta[indexedJetCSV[1].first]*sample.hJet_eta[indexedJetCSV[1].first];
			RMS_eta = sqrt(RMS_eta/4);
			hRMSeta.Fill(RMS_eta);
			StandDevEta[0] =sample.vLepton_eta[1];
			StandDevEta[1] =sample.hJet_eta[indexedJetCSV[0].first];
			StandDevEta[2] =sample.vLepton_eta[0];
			StandDevEta[3] =sample.hJet_eta[indexedJetCSV[1].first];
			EtaStandDev = TMath::RMS(4,StandDevEta);
			hStaDeveta.Fill(EtaStandDev);
			AverageEta = sample.vLepton_eta[0]+sample.vLepton_eta[1]+sample.hJet_eta[indexedJetCSV[0].first];
			AverageEta = (AverageEta+sample.hJet_eta[indexedJetCSV[1].first])/4;
			UnweightedEta = (sample.vLepton_eta[0]-AverageEta)*(sample.vLepton_eta[0]-AverageEta);
			UnweightedEta = UnweightedEta + (sample.vLepton_eta[1]-AverageEta)*(sample.vLepton_eta[1]-AverageEta);
			UnweightedEta = UnweightedEta + (sample.hJet_eta[indexedJetCSV[0].first]-AverageEta)*(sample.hJet_eta[indexedJetCSV[0].first]-AverageEta);
			UnweightedEta = UnweightedEta + (sample.hJet_eta[indexedJetCSV[1].first]-AverageEta)*(sample.hJet_eta[indexedJetCSV[1].first]-AverageEta);
			hUnweightedEta.Fill(UnweightedEta);
			

DeltaPhiHV = sample.HVdPhi;			
				hdphiVH.Fill(DeltaPhiHV);	
				hPtbalZH.Fill(Hpt-Zpt);
				
FirstJet.SetPtEtaPhiM(sample.hJet_pt[indexedJetCSV[0].first],sample.hJet_eta[indexedJetCSV[0].first],sample.hJet_phi[indexedJetCSV[0].first],4.2/1000.0);
SecondJet.SetPtEtaPhiM(sample.hJet_pt[indexedJetCSV[1].first],sample.hJet_eta[indexedJetCSV[1].first],sample.hJet_phi[indexedJetCSV[1].first],4.2/1000.0);
FirstMuon.SetPtEtaPhiM(mu1pt,sample.vLepton_eta[0],sample.vLepton_phi[0],105.65836668/1000.0);
SecondMuon.SetPtEtaPhiM(sample.vLepton_pt[1],sample.vLepton_eta[1],sample.vLepton_phi[1],105.65836668/1000.0);
Higgs = FirstJet+SecondJet;
ZBoson = FirstMuon+SecondMuon;

Ht = FirstJet.Pt()+SecondJet.Pt()+FirstMuon.Pt()+SecondMuon.Pt();
DphiJJ = FirstJet.DeltaPhi(SecondJet);

hdphiVH_vect.Fill(Higgs.DeltaPhi(ZBoson));
hdphiJJ_vect.Fill(DphiJJ);
hDphiDetajj.Fill((sample.hJet_eta[indexedJetCSV[0].first]-sample.hJet_eta[indexedJetCSV[1].first]),(EtaStandDev));
hPtjj_vect.Fill(Higgs.Pt());
hPtmumu_vect.Fill(ZBoson.Pt());
hHt.Fill(Ht);
hCentrality.Fill((indexedPt[2].second+indexedPt[3].second)/Ht);
hAngle.Fill(Higgs.Angle(ZBoson.Vect()));

HVsystem = 	Higgs+ZBoson;
hEventPt.Fill(HVsystem.Pt());


//Armenteros-Podolansky Plots
TVector3 VectZ = ZBoson.Vect();
TVector3 Unit_Z = VectZ.Unit();
				plmu1 = Unit_Z.Dot(FirstMuon.Vect());
				plmu2 = Unit_Z.Dot(SecondMuon.Vect());
alpha_mu = (plmu1-plmu2)/(plmu1+plmu2);
qtmu1= FirstMuon.Perp(ZBoson.Vect()); 
qtmu2 = SecondMuon.Perp(ZBoson.Vect());
hqtvsalphaZ.Fill(alpha_mu,qtmu1);

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
				hqtvsalphaJJ.Fill(alpha_j,qtb1);
				
//EventShapeVariables
std::vector<math::XYZVector> particles;
 TVector3 TVect_particle = FirstMuon.Vect();
 math::XYZVector Vparticle1(TVect_particle.X(), TVect_particle.Y(), TVect_particle.Z());
				particles.push_back(Vparticle1);
				TVect_particle = SecondMuon.Vect();
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
	 hAplanarity_EvtShp.Fill(EvntShpAplanarity);
				hSphericity_EvtShp.Fill(EvntShpSphericity);
				hCircularity_EvtShp.Fill(EvntShpCircularity);
				hIsotropy_EvtShp.Fill(EvntShpIsotropy);
				

			}// two muons and two jets requirement
			
		
if (debug) cout << "done filling histograms for event: " << event << endl;

//Calculate cut efficiencies
	if ((sample.nhJets > 1) && (sample.nvlep > 1)){
//all samples are of typle Zmumu
		if (CSV2 > -1 ){
			if (Zmass >75 && Zmass<105) {
			TMVA_tree->Fill();
				Ntree++;
				cout << "njets: " << nJets << endl;
			}
		if ( (CSV1 > 0.898) && (CSV2>0.5) && (sample.naJets < 1) && ((DeltaPhiHV>=2.90)||(DeltaPhiHV<-2.90)) ){
			if (sample.hJet_pt[indexedJetCSV[0].first] >=20 && sample.hJet_pt[indexedJetCSV[1].first] >=20 ) {
				hMmumu_NoPtCut.Fill(Zmass);
				hPtmumu_NoPtCut.Fill(Zpt);
				hPtmumu_vect_NoPtCut.Fill(ZBoson.Pt());
				hqtvsalphaZ_NoPtCut.Fill(alpha_mu,qtmu1);
				hPtmu1_aftercuts.Fill(mu1pt);
				if (mu1pt > 20) hPtmu2_aftercuts.Fill(sample.vLepton_pt[1]);		
			}// end jet pt cuts
			if (sample.vLepton_pt[1] > 20 && (Zmass >75 && Zmass<105) ) {
			hMjj_NoPtCut.Fill(Hmass);
			hPtjj_NoPtCut.Fill(Hpt);
			hPtjj_vect_NoPtCut.Fill(Higgs.Pt());
			hqtvsalphaJJ_NoPtCut.Fill(alpha_j,qtb1);
		if (sample.hJet_pt[indexedJetCSV[1].first] >=20)hPtb1_aftercuts.Fill(sample.hJet_pt[indexedJetCSV[0].first]);
		if (sample.hJet_pt[indexedJetCSV[0].first] >=20)hPtb2_aftercuts.Fill(sample.hJet_pt[indexedJetCSV[1].first]);
			}//end muon pt cuts
			hPtbalZH_NoPtCut.Fill(Hpt-Zpt);
			hHt_NoPtCut.Fill(Ht);
			hCentrality_NoPtCut.Fill((indexedPt[2].second+indexedPt[3].second)/Ht);
			hEventPt_NoPtCut.Fill(HVsystem.Pt());
		}// end non pt cuts
	if (sample.hJet_pt[indexedJetCSV[0].first] >=20 && sample.hJet_pt[indexedJetCSV[1].first] >=20 && sample.vLepton_pt[1] > 20) {
			if ((Hpt > 100) && (Zpt>100) && (CSV1 > 0.898) && (CSV2>0.5) && ((DeltaPhiHV>=2.90)||(DeltaPhiHV<-2.90)) && (sample.naJets < 1) ) hMmumu_aftercuts.Fill(Zmass);
			if (Zmass >75 && Zmass<105){
	Npreselect++;
				if ((Zpt>100) && (CSV1 > 0.898) && (CSV2>0.5) && ((DeltaPhiHV>=2.90)||(DeltaPhiHV<-2.90)) && (sample.naJets < 1) ) {
					hPtjj_aftercuts.Fill(Hpt);
					hPtjj_vect_aftercuts.Fill(Higgs.Pt());
					}
		if(Hpt > 100) {
		NpT_jj++;
			if ((CSV1 > 0.898) && (CSV2>0.5) && ((DeltaPhiHV>=2.90)||(DeltaPhiHV<-2.90)) && (sample.naJets < 1) ) {
			hPtmumu_aftercuts.Fill(Zpt);				   
			hPtmumu_vect_aftercuts.Fill(ZBoson.Pt());
		}
			if(Zpt > 100) {
			NpT_Z++;
				if ((CSV2>0.5) && ((DeltaPhiHV>=2.90)||(DeltaPhiHV<-2.90)) && (sample.naJets < 1) ) hCSV1_aftercuts.Fill(CSV1);
			   if(CSV1 > 0.898){
			   NCSV1++;
				   if (((DeltaPhiHV>=2.90)||(DeltaPhiHV<-2.90)) && (sample.naJets < 1) ) hCSV2_aftercuts.Fill(CSV2);
			   if(CSV2>0.5){
			   NCSV2++;
				   if (sample.naJets < 1) {
		hdphiVH_vect_aftercuts.Fill(Higgs.DeltaPhi(ZBoson));
		hdphiVH_aftercuts.Fill(DeltaPhiHV);	
		hAngle_NoDphiCut.Fill(Higgs.Angle(ZBoson.Vect()));
		hPhimu1_NoDphiCut.Fill(sample.vLepton_phi[0]);
		hPhimu2_NoDphiCut.Fill(sample.vLepton_phi[1]);
		hPhib1_NoDphiCut.Fill(sample.hJet_phi[indexedJetCSV[0].first]);
		hPhib2_NoDphiCut.Fill(sample.hJet_phi[indexedJetCSV[1].first]);	
		hdphiJJ_vect_NoDphiCut.Fill(EtaStandDev);
		hDphiDetajj_NoDphiCut.Fill((sample.hJet_eta[indexedJetCSV[0].first]-sample.hJet_eta[indexedJetCSV[1].first]),EtaStandDev);		
				   }
			   if((DeltaPhiHV>=2.90)||(DeltaPhiHV<-2.90)){
			   NdPhi++;
				   hNaj_aftercuts.Fill(sample.naJets);
			   if(sample.naJets < 2){
			   N_Naj++;
				   hEtab1_aftercuts.Fill(sample.hJet_eta[indexedJetCSV[0].first]);
				   hEtab2_aftercuts.Fill(sample.hJet_eta[indexedJetCSV[1].first]);
				   hCHFb1_aftercuts.Fill(sample.hJet_chf[indexedJetCSV[0].first]);
				   hCHFb2_aftercuts.Fill(sample.hJet_chf[indexedJetCSV[1].first]);	
				   hdetaJJ_aftercuts.Fill(sample.hJet_eta[indexedJetCSV[0].first]-sample.hJet_eta[indexedJetCSV[1].first]);		   
			   hMjj_aftercuts.Fill(Hmass);
				   hEtamu1_aftercuts.Fill(sample.vLepton_eta[0]);
				   hEtamu2_aftercuts.Fill(sample.vLepton_eta[1]);
				   hPhimu1_aftercuts.Fill(sample.vLepton_phi[0]);
				   hPhimu2_aftercuts.Fill(sample.vLepton_phi[1]);
				   hPFRelIsomu1_aftercuts.Fill(sample.vLepton_pfCombRelIso[0]);
				   hPFRelIsomu2_aftercuts.Fill(sample.vLepton_pfCombRelIso[1]);	
				   hRMSeta_aftercuts.Fill(RMS_eta);
			   hStaDeveta_aftercuts.Fill(EtaStandDev);
				   hUnweightedEta_aftercuts.Fill(UnweightedEta);
				   hPtbalZH_aftercuts.Fill(Hpt-Zpt);
		hPhib1_aftercuts.Fill(sample.hJet_phi[indexedJetCSV[0].first]);
		hPhib2_aftercuts.Fill(sample.hJet_phi[indexedJetCSV[1].first]);	
		hDphiDetajj_aftercuts.Fill((sample.hJet_eta[indexedJetCSV[0].first]-sample.hJet_eta[indexedJetCSV[1].first]),(sample.hJet_phi[indexedJetCSV[0].first]-sample.hJet_phi[indexedJetCSV[1].first]));		
				   hdphiJJ_vect_aftercuts.Fill(EtaStandDev);
				   hHt_aftercuts.Fill(Ht);
				   hCentrality_aftercuts.Fill((indexedPt[2].second+indexedPt[3].second)/Ht);
				   hAngle_aftercuts.Fill(Higgs.Angle(ZBoson.Vect()));
				   hEventPt_aftercuts.Fill(HVsystem.Pt());
				   hqtvsalphaZ_aftercuts.Fill(alpha_mu,qtmu1);
				   hqtvsalphaJJ_aftercuts.Fill(alpha_j,qtb1);
				   hAplanarity_EvtShp_aftercuts.Fill(EvntShpAplanarity);
				   hSphericity_EvtShp_aftercuts.Fill(EvntShpSphericity);
				   hCircularity_EvtShp_aftercuts.Fill(EvntShpCircularity);
				   hIsotropy_EvtShp_aftercuts.Fill(EvntShpIsotropy);				   
			   if((Hmass>=95)&&(Hmass<=125)){N_Mjj++; }//higgs mass window
			   }//number of additional Jets cut
			   }//Delta Phi Cut
			   }//CSV2 cut
			   }//CVS1 cut
			   }// Z pt cut
		}// dijet pt cut
	}// preselection z mass requirement
	}//preselection CSV not -1
	}// preselection pt requirement
}// preselection number of muons/jets requirement
		}//end requirement muon not electron
		}//end requirement Zmumu event
    } while (sample.nextEvent());
if (debug){	std::cout << "Number of events " << event << endl;}
	std::cout << "tree: " << Ntree << endl;
	std::cout << "PreSelection: " << Npreselect << endl;
std::cout << "pT_jj: " << NpT_jj << endl;
std::cout << "pT_Z: " << NpT_Z << endl;
std::cout << "CSV1: " << NCSV1 << endl;
std::cout << "CSV2: " << NCSV2 << endl;
std::cout << "DPhi: " << NdPhi << endl;
std::cout << "Naj: " << N_Naj << endl;
std::cout << "Mjj: " << N_Mjj << endl;
std::cout << "Number of Events Passing the Selection: " << N_Mjj << endl;


    // Here one can create canvas, draw and print the histogram.
    TCanvas c1("c1","c1");
    c1.cd();
    c1.SetFillColor(kWhite);
    string suffixps = ".gif";
	
	hnJets.Draw();
    c1.Print((directory+"/nhJets"+suffixps).c_str());
	
    c1.Clear(); // don't create a new canvas
    hnMuons.Draw();
    c1.Print((directory+"/nMuons"+suffixps).c_str());
	
	c1.Clear(); // don't create a new canvas
    hnSV.Draw();
    c1.Print((directory+"/nSV"+suffixps).c_str());
	
	c1.Clear(); // don't create a new canvas
    hnPV.Draw();
    c1.Print((directory+"/nPV"+suffixps).c_str());
	
    c1.Clear(); // don't create a new canvas
    hallhJet_pt.Draw();
    c1.Print((directory+"/PtAllJets"+suffixps).c_str());
	
	c1.Clear(); // don't create a new canvas
    hMET.Draw();
    c1.Print((directory+"/MET"+suffixps).c_str());	

    c1.Clear(); // don't create a new canvas
    hPtb1.Draw();
    c1.Print((directory+"/Ptb1"+suffixps).c_str());
	
    c1.Clear(); // don't create a new canvas
    hPtb2.Draw();
    c1.Print((directory+"/Ptb2"+suffixps).c_str());
	
    c1.Clear(); // don't create a new canvas
    hEtab1.Draw();
    c1.Print((directory+"/Etab1"+suffixps).c_str());
	
    c1.Clear(); // don't create a new canvas
    hEtab2.Draw();
    c1.Print((directory+"/Etab2"+suffixps).c_str());
	
    c1.Clear(); // don't create a new canvas
    hPhib1.Draw();
    c1.Print((directory+"/Phib1"+suffixps).c_str());
	
    c1.Clear(); // don't create a new canvas
    hPhib2.Draw();
    c1.Print((directory+"/Phib2"+suffixps).c_str());

    c1.Clear(); // don't create a new canvas
    hCHFb1.Draw();
    c1.Print((directory+"/CHFb1"+suffixps).c_str());

    c1.Clear(); // don't create a new canvas
    hCHFb2.Draw();
    c1.Print((directory+"/CHFb2"+suffixps).c_str());
	
	c1.Clear(); // don't create a new canvas
    hPtjj.Draw();
    c1.Print((directory+"/Ptjj"+suffixps).c_str());
	
	c1.Clear(); // don't create a new canvas
    hPtmumu.Draw();
    c1.Print((directory+"/Ptmumu"+suffixps).c_str());

	c1.Clear(); // don't create a new canvas
    hPtmu1.Draw();
    c1.Print((directory+"/Ptmu1"+suffixps).c_str());
	
	c1.Clear(); // don't create a new canvas
    hPtmu2.Draw();
    c1.Print((directory+"/Ptmu2"+suffixps).c_str());

	c1.Clear(); // don't create a new canvas
    hEtamu1.Draw();
    c1.Print((directory+"/Etamu1"+suffixps).c_str());
	
	c1.Clear(); // don't create a new canvas
    hEtamu2.Draw();
    c1.Print((directory+"/Etamu2"+suffixps).c_str());
	
	c1.Clear(); // don't create a new canvas
    hPhimu1.Draw();
    c1.Print((directory+"/Phimu1"+suffixps).c_str());
	
	c1.Clear(); // don't create a new canvas
    hPhimu2.Draw();
    c1.Print((directory+"/Phimu2"+suffixps).c_str());
	
	c1.Clear(); // don't create a new canvas
    hPFRelIsomu1.Draw();
    c1.Print((directory+"/PFRelIsomu1"+suffixps).c_str());

	c1.Clear(); // don't create a new canvas
    hPFRelIsomu2.Draw();
    c1.Print((directory+"/PFRelIsomu2"+suffixps).c_str());
	
	c1.Clear(); // don't create a new canvas
    hCSV1.Draw();
    c1.Print((directory+"/CSV1"+suffixps).c_str());
	
	c1.Clear(); // don't create a new canvas
    hCSV2.Draw();
    c1.Print((directory+"/CSV2"+suffixps).c_str());
	
	c1.Clear(); // don't create a new canvas
	hNaj.Draw();
	c1.Print((directory+"/Naj"+suffixps).c_str());
		
	c1.Clear(); // don't create a new canvas
	hMjj.Draw();
	c1.Print((directory+"/Mjj"+suffixps).c_str());

	c1.Clear(); // don't create a new canvas
	hMmumu.Draw();
	c1.Print((directory+"/Mmumu"+suffixps).c_str());

	c1.Clear(); // don't create a new canvas
	hdetaJJ.Draw();
	c1.Print((directory+"/detaJJ"+suffixps).c_str());

	c1.Clear(); // don't create a new canvas
	hdphiVH.Draw();
	c1.Print((directory+"/dphiVH"+suffixps).c_str());	
	
	c1.Clear(); // don't create a new canvas
	hPtbalZH.Draw();
	c1.Print((directory+"/PtbalZH"+suffixps).c_str());	
	
	c1.Clear(); // don't create a new canvas
	hDphiDetajj.Draw();
	c1.Print((directory+"/DphivsDeta_JJ"+suffixps).c_str());	

	c1.Clear(); // don't create a new canvas
	hdphiVH_vect.Draw();
	c1.Print((directory+"/dphiVH_vect"+suffixps).c_str());

	c1.Clear(); // don't create a new canvas
	hdphiJJ_vect.Draw();
	c1.Print((directory+"/dphiJJ_vect"+suffixps).c_str());
	
	c1.Clear(); // don't create a new canvas
	hPtjj_vect.Draw();
	c1.Print((directory+"/Ptjj_vect"+suffixps).c_str());

	c1.Clear(); // don't create a new canvas
	hPtmumu_vect.Draw();
	c1.Print((directory+"/Ptmumu_vect"+suffixps).c_str());

	c1.Clear(); // don't create a new canvas
	hHt.Draw();
	c1.Print((directory+"/Ht"+suffixps).c_str());

	c1.Clear(); // don't create a new canvas
	hCentrality.Draw();
	c1.Print((directory+"/Centrality"+suffixps).c_str());

	c1.Clear(); // don't create a new canvas
	hEventPt.Draw();
	c1.Print((directory+"/EventPt"+suffixps).c_str());
	
	c1.Clear(); // don't create a new canvas
	hAngle.Draw();
	c1.Print((directory+"/Angle"+suffixps).c_str());

	c1.Clear(); // don't create a new canvas
	hqtvsalphaJJ.Draw();
	c1.Print((directory+"/qtvsalphaJJ"+suffixps).c_str());

	c1.Clear(); // don't create a new canvas
	hqtvsalphaZ.Draw();
	c1.Print((directory+"/qtvsalphaZ"+suffixps).c_str());

	c1.Clear(); // don't create a new canvas
	hRMSeta.Draw();
	c1.Print((directory+"/RMSeta"+suffixps).c_str());

	c1.Clear(); // don't create a new canvas
	hStaDeveta.Draw();
	c1.Print((directory+"/StaDeveta"+suffixps).c_str());

	c1.Clear(); // don't create a new canvas
	hUnweightedEta.Draw();
	c1.Print((directory+"/UnweightedEta"+suffixps).c_str());

	c1.Clear(); // don't create a new canvas
	hSphericity_EvtShp.Draw();
	c1.Print((directory+"/Sphericity_EvtShp"+suffixps).c_str());

	c1.Clear(); // don't create a new canvas
	hAplanarity_EvtShp.Draw();
	c1.Print((directory+"/Aplanarity_EvtShp"+suffixps).c_str());

	c1.Clear(); // don't create a new canvas
	hCircularity_EvtShp.Draw();
	c1.Print((directory+"/Circularity_EvtShp"+suffixps).c_str());

	c1.Clear(); // don't create a new canvas
	hIsotropy_EvtShp.Draw();
	c1.Print((directory+"/Isotropy_EvtShp"+suffixps).c_str());
		
	
	c1.Clear(); // don't create a new canvas
/*	hCutFlow.GetXaxis().SetBinLabel(1,"PreSelection");
	hCutFlow.GetXaxis().SetBinLabel(2,"p_{T}(jj)");
	hCutFlow.GetXaxis().SetBinLabel(3,"p_{T}(Z)");
	hCutFlow.GetXaxis().SetBinLabel(4,"CSV1");
	hCutFlow.GetXaxis().SetBinLabel(5,"CSV2");
	hCutFlow.GetXaxis().SetBinLabel(6,"#Delta#phi");
	hCutFlow.GetXaxis().SetBinLabel(7,"N_{aj}");
	hCutFlow.GetXaxis().SetBinLabel(8,"M_{jj}");*/
	hCutFlow.SetBinContent(1,Npreselect );
	hCutFlow.SetBinContent(2,NpT_jj );
	hCutFlow.SetBinContent(3,NpT_Z );
	hCutFlow.SetBinContent(4,NCSV1 );
	hCutFlow.SetBinContent(5,NCSV2 );
	hCutFlow.SetBinContent(6,NdPhi );
	hCutFlow.SetBinContent(7,N_Naj );
	hCutFlow.SetBinContent(8,N_Mjj );
	hCutFlow.Draw();
	c1.Print((directory+"/CutFlow"+suffixps).c_str());
	
	c1.Clear(); // don't create a new canvas
    hPtb1_aftercuts.Draw();
    c1.Print((directory+"/Ptb1_aftercuts"+suffixps).c_str());
	
    c1.Clear(); // don't create a new canvas
    hPtb2_aftercuts.Draw();
    c1.Print((directory+"/Ptb2_aftercuts"+suffixps).c_str());
	
    c1.Clear(); // don't create a new canvas
    hEtab1_aftercuts.Draw();
    c1.Print((directory+"/Etab1_aftercuts"+suffixps).c_str());
	
    c1.Clear(); // don't create a new canvas
    hEtab2_aftercuts.Draw();
    c1.Print((directory+"/Etab2_aftercuts"+suffixps).c_str());
	
    c1.Clear(); // don't create a new canvas
    hPhib1_aftercuts.Draw();
    c1.Print((directory+"/Phib1_aftercuts"+suffixps).c_str());
	
    c1.Clear(); // don't create a new canvas
    hPhib2_aftercuts.Draw();
    c1.Print((directory+"/Phib2_aftercuts"+suffixps).c_str());
	
    c1.Clear(); // don't create a new canvas
    hCHFb1_aftercuts.Draw();
    c1.Print((directory+"/CHFb1_aftercuts"+suffixps).c_str());
	
    c1.Clear(); // don't create a new canvas
    hCHFb2_aftercuts.Draw();
    c1.Print((directory+"/CHFb2_aftercuts"+suffixps).c_str());
	
	c1.Clear(); // don't create a new canvas
    hPtjj_aftercuts.Draw();
    c1.Print((directory+"/Ptjj_aftercuts"+suffixps).c_str());
	
	c1.Clear(); // don't create a new canvas
    hPtmumu_aftercuts.Draw();
    c1.Print((directory+"/Ptmumu_aftercuts"+suffixps).c_str());
	
	c1.Clear(); // don't create a new canvas
    hPtmu1_aftercuts.Draw();
    c1.Print((directory+"/Ptmu1_aftercuts"+suffixps).c_str());
	
	c1.Clear(); // don't create a new canvas
    hPtmu2_aftercuts.Draw();
    c1.Print((directory+"/Ptmu2_aftercuts"+suffixps).c_str());
	
	c1.Clear(); // don't create a new canvas
    hEtamu1_aftercuts.Draw();
    c1.Print((directory+"/Etamu1_aftercuts"+suffixps).c_str());
	
	c1.Clear(); // don't create a new canvas
    hEtamu2_aftercuts.Draw();
    c1.Print((directory+"/Etamu2_aftercuts"+suffixps).c_str());
	
	c1.Clear(); // don't create a new canvas
    hPhimu1_aftercuts.Draw();
    c1.Print((directory+"/Phimu1_aftercuts"+suffixps).c_str());
	
	c1.Clear(); // don't create a new canvas
    hPhimu2_aftercuts.Draw();
    c1.Print((directory+"/Phimu2_aftercuts"+suffixps).c_str());
	
	c1.Clear(); // don't create a new canvas
    hPFRelIsomu1_aftercuts.Draw();
    c1.Print((directory+"/PFRelIsomu1_aftercuts"+suffixps).c_str());
	
	c1.Clear(); // don't create a new canvas
    hPFRelIsomu2_aftercuts.Draw();
    c1.Print((directory+"/PFRelIsomu2_aftercuts"+suffixps).c_str());
	
	c1.Clear(); // don't create a new canvas
    hCSV1_aftercuts.Draw();
    c1.Print((directory+"/CSV1_aftercuts"+suffixps).c_str());
	
	c1.Clear(); // don't create a new canvas
    hCSV2_aftercuts.Draw();
    c1.Print((directory+"/CSV2_aftercuts"+suffixps).c_str());
	
	c1.Clear(); // don't create a new canvas
	hNaj_aftercuts.Draw();
	c1.Print((directory+"/Naj_aftercuts"+suffixps).c_str());
	
	c1.Clear(); // don't create a new canvas
	hMjj_aftercuts.Draw();
	c1.Print((directory+"/Mjj_aftercuts"+suffixps).c_str());
	
	c1.Clear(); // don't create a new canvas
	hMmumu_aftercuts.Draw();
	c1.Print((directory+"/Mmumu_aftercuts"+suffixps).c_str());
	
	c1.Clear(); // don't create a new canvas
	hdetaJJ_aftercuts.Draw();
	c1.Print((directory+"/detaJJ_aftercuts"+suffixps).c_str());
	
	c1.Clear(); // don't create a new canvas
	hdphiVH_aftercuts.Draw();
	c1.Print((directory+"/dphiVH_aftercuts"+suffixps).c_str());	
	
	c1.Clear(); // don't create a new canvas
	hPtbalZH_aftercuts.Draw();
	c1.Print((directory+"/PtbalZH_aftercuts"+suffixps).c_str());	
		
	c1.Clear(); // don't create a new canvas
	hDphiDetajj_aftercuts.Draw();
	c1.Print((directory+"/DphivsDeta_JJ_aftercuts"+suffixps).c_str());	
	
	c1.Clear(); // don't create a new canvas
	hdphiVH_vect_aftercuts.Draw();
	c1.Print((directory+"/dphiVH_vect_aftercuts"+suffixps).c_str());
	
	c1.Clear(); // don't create a new canvas
	hdphiJJ_vect_aftercuts.Draw();
	c1.Print((directory+"/dphiJJ_vect_aftercuts"+suffixps).c_str());
	
	c1.Clear(); // don't create a new canvas
	hPtjj_vect_aftercuts.Draw();
	c1.Print((directory+"/Ptjj_vect_aftercuts"+suffixps).c_str());
	
	c1.Clear(); // don't create a new canvas
	hPtmumu_vect_aftercuts.Draw();
	c1.Print((directory+"/Ptmumu_vect_aftercuts"+suffixps).c_str());
	
	c1.Clear(); // don't create a new canvas
	hHt_aftercuts.Draw();
	c1.Print((directory+"/Ht_aftercuts"+suffixps).c_str());
	
	c1.Clear(); // don't create a new canvas
	hCentrality_aftercuts.Draw();
	c1.Print((directory+"/Centrality_aftercuts"+suffixps).c_str());
	
	c1.Clear(); // don't create a new canvas
	hEventPt_aftercuts.Draw();
	c1.Print((directory+"/EventPt_aftercuts"+suffixps).c_str());
	
	c1.Clear(); // don't create a new canvas
	hAngle_aftercuts.Draw();
	c1.Print((directory+"/Angle_aftercuts"+suffixps).c_str());
	
	c1.Clear(); // don't create a new canvas
	hqtvsalphaJJ_aftercuts.Draw();
	c1.Print((directory+"/qtvsalphaJJ_aftercuts"+suffixps).c_str());
	
	c1.Clear(); // don't create a new canvas
	hqtvsalphaZ_aftercuts.Draw();
	c1.Print((directory+"/qtvsalphaZ_aftercuts"+suffixps).c_str());
	
	c1.Clear(); // don't create a new canvas
	hRMSeta_aftercuts.Draw();
	c1.Print((directory+"/RMSeta_aftercuts"+suffixps).c_str());
	
	c1.Clear(); // don't create a new canvas
	hStaDeveta_aftercuts.Draw();
	c1.Print((directory+"/StaDeveta_aftercuts"+suffixps).c_str());
	
	c1.Clear(); // don't create a new canvas
	hUnweightedEta_aftercuts.Draw();
	c1.Print((directory+"/UnweightedEta_aftercuts"+suffixps).c_str());
	
	c1.Clear(); // don't create a new canvas
	hSphericity_EvtShp_aftercuts.Draw();
	c1.Print((directory+"/Sphericity_EvtShp_aftercuts"+suffixps).c_str());
	
	c1.Clear(); // don't create a new canvas
	hAplanarity_EvtShp_aftercuts.Draw();
	c1.Print((directory+"/Aplanarity_EvtShp_aftercuts"+suffixps).c_str());
	
	c1.Clear(); // don't create a new canvas
	hCircularity_EvtShp_aftercuts.Draw();
	c1.Print((directory+"/Circularity_EvtShp_aftercuts"+suffixps).c_str());
	
	c1.Clear(); // don't create a new canvas
	hIsotropy_EvtShp_aftercuts.Draw();
	c1.Print((directory+"/Isotropy_EvtShp_aftercuts"+suffixps).c_str());
	
    c1.Clear(); // don't create a new canvas
    hPhib1_NoDphiCut.Draw();
    c1.Print((directory+"/Phib1_NoDphiCut"+suffixps).c_str());
	
    c1.Clear(); // don't create a new canvas
    hPhib2_NoDphiCut.Draw();
    c1.Print((directory+"/Phib2_NoDphiCut"+suffixps).c_str());

	c1.Clear(); // don't create a new canvas
    hPtjj_NoPtCut.Draw();
    c1.Print((directory+"/Ptjj_NoPtCut"+suffixps).c_str());
	
	c1.Clear(); // don't create a new canvas
    hPtmumu_NoPtCut.Draw();
    c1.Print((directory+"/Ptmumu_NoPtCut"+suffixps).c_str());
	
	c1.Clear(); // don't create a new canvas
    hPhimu1.Draw();
    c1.Print((directory+"/Phimu1_NoDphiCut"+suffixps).c_str());
	
	c1.Clear(); // don't create a new canvas
    hPhimu2.Draw();
    c1.Print((directory+"/Phimu2_NoDphiCut"+suffixps).c_str());
	
	c1.Clear(); // don't create a new canvas
	hMjj_NoPtCut.Draw();
	c1.Print((directory+"/Mjj_NoPtCut"+suffixps).c_str());
	
	c1.Clear(); // don't create a new canvas
	hMmumu_NoPtCut.Draw();
	c1.Print((directory+"/Mmumu_NoPtCut"+suffixps).c_str());
	
	c1.Clear(); // don't create a new canvas
	hPtbalZH_NoPtCut.Draw();
	c1.Print((directory+"/PtbalZH_NoPtCut"+suffixps).c_str());	
		
	c1.Clear(); // don't create a new canvas
	hDphiDetajj_NoDphiCut.Draw();
	c1.Print((directory+"/DphiDetajj_NoDphiCut"+suffixps).c_str());	
	
	c1.Clear(); // don't create a new canvas
	hdphiJJ_vect_NoDphiCut.Draw();
	c1.Print((directory+"/dphiJJ_vect_NoDphiCut"+suffixps).c_str());
	
	c1.Clear(); // don't create a new canvas
	hPtjj_vect_NoPtCut.Draw();
	c1.Print((directory+"/Ptjj_vect_NoPtCut"+suffixps).c_str());
	
	c1.Clear(); // don't create a new canvas
	hPtmumu_vect_NoPtCut.Draw();
	c1.Print((directory+"/Ptmumu_vect_NoPtCut"+suffixps).c_str());
	
	c1.Clear(); // don't create a new canvas
	hHt_NoPtCut.Draw();
	c1.Print((directory+"/Ht_NoPtCut"+suffixps).c_str());
	
	c1.Clear(); // don't create a new canvas
	hCentrality_NoPtCut.Draw();
	c1.Print((directory+"/Centrality_NoPtCut"+suffixps).c_str());
	
	c1.Clear(); // don't create a new canvas
	hEventPt_NoPtCut.Draw();
	c1.Print((directory+"/EventPt_NoPtCut"+suffixps).c_str());
	
	c1.Clear(); // don't create a new canvas
	hAngle_NoDphiCut.Draw();
	c1.Print((directory+"/Angle_NoDphiCut"+suffixps).c_str());
	
	c1.Clear(); // don't create a new canvas
	hqtvsalphaJJ_NoPtCut.Draw();
	c1.Print((directory+"/qtvsalphaJJ_NoPtCut"+suffixps).c_str());
	
	c1.Clear(); // don't create a new canvas
	hqtvsalphaZ_NoPtCut.Draw();
	c1.Print((directory+"/qtvsalphaZ_NoPtCut"+suffixps).c_str());

	TMVA_tree->Write();


    // Write and Close the output file.
    ofile.Write();
    ofile.Close();

    return 0;
}
