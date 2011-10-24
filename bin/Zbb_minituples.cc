/**
    @file exampleAnalysis.cc

    @brief Example of an analysis code to read hbb ntuple
    and create histogram
    Usage:

    \verbatim
    mkdir [directoryName]
    exampleAnalysis <inputfile> [outputfilename] [directoryName]
    \endverbatim

    @param inputfile Either a ROOT file or an ASCII file containing list of
    ROOT files.

    @param outputfile Name of the ROOT file which contains the histogram.
    Defaulted to 'output.root'

    @author Rachel Wilken <rachel.wilken@cern.ch>

    @date Thu Oct 20 2011

    @version $Id: exampleSort.cc,v 1.2 2011/10/20 19:33:58 wilken Exp $
 */


#include <UserCode/wilken/interface/hbbReader.h>
#include <UserCode/wilken/interface/GenericTool.h>
#include <string>
#include <iostream>
#include <vector>
#include <algorithm>
#include <map>

#include "TLorentzVector.h"

#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
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


    // Create the hbbReader object.
    hbbReader sample(ifilename,
                                     std::string("hbb"));

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

bool debug = false;	


    // Here one can declare histograms
    // In compiled C++, it possible not to use pointer type (TH1F*) and
    // use the non-pointer type (TH1F) of ROOT object.

	//Declare Histograms
	TH1F *hjetPt[9];
	TH1F *hjetCSV[9];

	string Num[10];
	Num[0] = "0";
	Num[1] = "1";
	Num[2] = "2";
	Num[3] = "3";
	Num[4] = "4";
	Num[5] = "5";
	Num[6] = "6";
	Num[7] = "7";
	Num[8] = "8"; 
	Num[9] = "9";

	for (int i=0; i<9; i++){
		hjetPt[i] = new TH1F(("jetPt"+Num[i]).c_str(), ("p_{T} of Jet " +Num[i]+"(ordered in p_{T})").c_str(), 100, 0., 200.);
		hjetCSV[i] = new TH1F(("jetCSV"+Num[i]).c_str(), ("CSV of Jet " +Num[i]+"(ordered in CSV)").c_str(), 50, -0.5, 1.25);
	}// end i loop

    TH1F hnJets ("hnJets",  "Number of Good Jets", 11, -0.5, 10.5);
    TH1F hnMuons  ("hnMuons","Number of Good Muons", 6, -0.5, 5.5);
    TH1F halljetPt  ("halljetPt","Pt of all jets in event",		150, 0.0, 150);
    TH1F hPtb1		("hPtb1","Pt b jet with highest CSV",		150, 0.0, 150);
    TH1F hPtb2		("hPtb2","Pt b jet with second highest CSV", 150, 0.0, 150);
	TH1F hPtjj		("hPtjj","Pt of two b jets with highest CSV", 300, 0.0, 300);
	TH1F hPtmumu	("hPtmumu","Pt of two muons with highest pt", 300, 0.0, 300);
	TH1F hPtmu1		("hPtmu1","Pt of muon with highest pt", 75, 0.0, 150);
	TH1F hPtmu2		("hPtmu2","Pt of muon with second highest pt", 75, 0.0, 150);
	TH1F hCSV1		("hCSV1","Jet with highest CSV",			50, -0.5, 1.25);
	TH1F hCSV2		("hCSV2","Jet with second highest CSV",		50, -0.5, 1.25);
	TH1F hdphiVH	("hdphiVH","Delta phi between Z and Higgs", 50, -3.5, 3.5);
	TH1F hNaj		("hNaj",  "Number of Additional Jets",		13, -2.5, 10.5);
	TH1F hMjj		("hMjj",  "Invariant Mass of two Jets",		75, 0, 200);
	TH1F hMmumu		("hMmumu",  "Invariant Mass of two muons",	75, 0, 200);
	TH1F hCutFlow	("hCutFlow",  "Selection",					9, 0, 9);
	
	
if (debug) std::cout << "all histograms declared " << std::endl;


    // Loop over all events.
    // For other methods to access event/navigate through the sample,
    // see the documentation of RootTreeReader class.
	int event =0, Npreselect =0, NpT_jj =0, NpT_Z =0, NCSV1 =0, NCSV2 =0;
	int NdPhi =0, N_Naj =0, N_Mjj =0, N_HLT =0;
//	double Zmass = 91.1976;


    do {
event++;
	if (debug) std::cout << "entered event loop " << event << std::endl;
	std::vector< std::pair<size_t,double> > indexedJetPt;
        std::vector<size_t> PtSortedJetIndex;
	std::vector< std::pair<size_t,double> > indexedJetCSV;
        std::vector<size_t> CSVSortedJetIndex;

        // Analysis loop.
        // One can access the ntuple leaves directly from sample object
int AdditionalJets = 0;

for (int k=0;k<sample.nJets;k++){
	if (debug) cout << "for "<< k << "th pass" << endl;
	indexedJetPt.push_back(std::pair<size_t,double>(k,(double) sample.jetPt[k]));
	halljetPt.Fill(sample.jetPt[k]); 
	indexedJetCSV.push_back(std::pair<size_t,double>(k,(double) sample.bDisc_CSV[k]));
	if ( 	(sample.jetPt[k] > 20) && (sample.jetEta[k] > -2.5 ) && (sample.jetEta[k] < 2.5) ) AdditionalJets++;
	if (debug)	cout << "CSV discriminator: " << sample.bDisc_CSV[k] << endl;

}// end k for loop

if (debug) for (size_t i = 0 ; (i != indexedJetCSV.size()) ; ++i) {cout << "Unsorted jet CSV: " << indexedJetCSV[i].second << endl; }

std::sort(indexedJetPt.begin(),indexedJetPt.end(),::IndexedQuantityGreaterThan<double>);
std::sort(indexedJetCSV.begin(),indexedJetCSV.end(),::IndexedQuantityGreaterThan<double>);

for (size_t i = 0 ; (i != indexedJetPt.size()) ; ++i) {   PtSortedJetIndex.push_back(indexedJetPt[i].first);        }
for (size_t i = 0 ; (i != indexedJetCSV.size()) ; ++i) {   CSVSortedJetIndex.push_back(indexedJetCSV[i].first);        }

if (debug) for (size_t i = 0 ; (i != indexedJetCSV.size()) ; ++i) {cout << "Sorted jet CSV: " << indexedJetCSV[i].second << endl; }

double Phi_Z = -99.99, Phi_H = -99.99;
double Mass_jj = -99.99, Pt_jj = -99.99;
		if (indexedJetCSV.size() > 1) {
TLorentzVector FirstJet = TLorentzVector(sample.jetPx[indexedJetCSV[0].first],sample.jetPy[indexedJetCSV[0].first],sample.jetPz[indexedJetCSV[0].first],sample.jetP[indexedJetCSV[0].first]); 
TLorentzVector SecondJet = TLorentzVector(sample.jetPx[indexedJetCSV[1].first],sample.jetPy[indexedJetCSV[1].first],sample.jetPz[indexedJetCSV[1].first],sample.jetP[indexedJetCSV[1].first]); 
		
TLorentzVector jj = FirstJet + SecondJet;
Mass_jj= jj.M();
Pt_jj = jj.Pt();
Phi_H = jj.Phi();
		
		} //if have two b jets

double Mass_mumu = -99.99, Pt_mumu = -99.99;
		if (sample.nMuons >1){
TLorentzVector FirstMuon = TLorentzVector(sample.muonPx[0], sample.muonPy[0], sample.muonPz[0], sample.muonP[0]);
TLorentzVector SecondMuon = TLorentzVector(sample.muonPx[1], sample.muonPy[1], sample.muonPz[1], sample.muonP[1]);

TLorentzVector mumu = FirstMuon + SecondMuon;

Mass_mumu = mumu.M();
Phi_Z = mumu.Phi();
Pt_mumu = mumu.Pt();
		}// if have two muons
		

double dPhiVH = deltaphi(Phi_H, Phi_Z);
		
for (size_t i = 0 ; i != indexedJetCSV.size() && i != indexedJetPt.size() && (i < 8); ++i) {
				hjetPt[i+1]->Fill(sample.jetPt[indexedJetPt[i].first]);
				hjetCSV[i+1]->Fill(sample.bDisc_CSV[indexedJetCSV[i].first]);
			}
	if (debug) cout << "Array of Histograms Filled " << endl;
		
		if (sample.nJets > 1) { 
		if (sample.bDisc_CSV[1] > -1) { 	
			hPtb1.Fill(sample.jetPt[indexedJetCSV[0].first]);
			hPtb2.Fill(sample.jetPt[indexedJetCSV[1].first]);
			hCSV1.Fill(sample.bDisc_CSV[indexedJetCSV[0].first]);
			hCSV2.Fill(sample.bDisc_CSV[indexedJetCSV[1].first]);
		}
			hMjj.Fill(Mass_jj);
		    hPtjj.Fill(Pt_jj);
			if (sample.nMuons > 1) hdphiVH.Fill(dPhiVH);			
		}
	if (debug) cout << "halfway through filling histos" << endl;
		if (sample.nJets > -1) { hnJets.Fill(sample.nJets); }
		if (sample.nMuons > -1) { hnMuons.Fill(sample.nMuons); }
		hNaj.Fill(AdditionalJets-2);
		if (sample.nMuons > 1) {
		hMmumu.Fill(Mass_mumu);
		hPtmumu.Fill(Pt_mumu);
		hPtmu1.Fill(sample.muonPt[0]);
		hPtmu2.Fill(sample.muonPt[1]);		
		}

if (debug) cout << "done filling histograms for event: " << event << endl;

//Calculate cut efficiencies
//		if ( (sample.HLT_IsoMu17_v5_trig) || (sample.HLT_IsoMu13_v4_trig) || (sample.HLT_IsoMu15_v4_trig) || (sample.HLT_IsoMu17_v4_trig) || (sample.HLT_IsoMu17_v6_trig) || (sample.HLT_CentralJet80_MET80_trig) || (sample.HLT_pfMHT150_trig) || (sample.HLT_L1_ETM30_trig) || (sample.HLT_MET100_trig) || (sample.HLT_IsoMu17_v8_trig) || (sample.HLT_IsoMu17_v7_trig)  ){
			if ( (sample.HLT_IsoMu17_v5_trig) ){
	N_HLT++;
	if ((sample.nJets > 1) && (sample.nMuons > 1)){
	if (sample.jetPt[indexedJetCSV[0].first] >=20 && sample.jetPt[indexedJetCSV[1].first] >= 20 && sample.muonPt[1] >= 20) {
	Npreselect++;
		if(Pt_jj >= 100) {
		NpT_jj++;
			if(Pt_mumu >= 100) {
			NpT_Z++;
			   if(sample.bDisc_CSV[indexedJetCSV[0].first] >= 0.898){
			   NCSV1++;
			   if(sample.bDisc_CSV[indexedJetCSV[1].first]>=0.5){
			   NCSV2++;
			   if((dPhiVH>=2.95)||(dPhiVH<=-2.95)){
			   NdPhi++;
			   if(sample.nJets < 5){
			   N_Naj++;
			   if((Mass_jj>=95)&&(Mass_jj<=125)){N_Mjj++; }//higgs mass window
			   }//number of additional Jets cut
			   }//Delta Phi Cut
			   }//CSV2 cut
			   }//CVS1 cut
			   }// Z pt cut
		}// dijet pt cut
	}// preselection Pt requirement
}// preselection number of muons/jets requirement
	}//trigger requirement


    } while (sample.nextEvent());
if (debug){	std::cout << "Number of events " << event << endl;}
std::cout << "HLT: " << N_HLT << endl;
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
    c1.Print((directory+"/nJets"+suffixps).c_str());
	
    c1.Clear(); // don't create a new canvas
    hnMuons.Draw();
    c1.Print((directory+"/nMuons"+suffixps).c_str());
	
    
	for (int a=1; a<9; a++){
		hjetPt[a]->Draw();
		hjetPt[a]->Write();
		c1.Print((directory+"/JetPt"+Num[a]+suffixps).c_str());
		c1.Clear();

		hjetCSV[a]->Draw();
		hjetCSV[a]->Write();
		c1.Print((directory+"/JetCSV"+Num[a]+suffixps).c_str());
		c1.Clear();

	}//end loop a


    c1.Clear(); // don't create a new canvas
    halljetPt.Draw();
    c1.Print((directory+"/PtAllJets"+suffixps).c_str());

    c1.Clear(); // don't create a new canvas
    hPtb1.Draw();
    c1.Print((directory+"/Ptb1"+suffixps).c_str());
	
    c1.Clear(); // don't create a new canvas
    hPtb2.Draw();
    c1.Print((directory+"/Ptb2"+suffixps).c_str());
	
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
	hdphiVH.Draw();
	c1.Print((directory+"/dphiVH"+suffixps).c_str());

	c1.Clear(); // don't create a new canvas
/*	hCutFlow.GetXaxis().SetBinLabel(1,"HLT");
	hCutFlow.GetXaxis().SetBinLabel(2,"PreSelection");
	hCutFlow.GetXaxis().SetBinLabel(3,"p_{T}(jj)");
	hCutFlow.GetXaxis().SetBinLabel(4,"p_{T}(Z)");
	hCutFlow.GetXaxis().SetBinLabel(5,"CSV1");
	hCutFlow.GetXaxis().SetBinLabel(6,"CSV2");
	hCutFlow.GetXaxis().SetBinLabel(7,"#Delta#phi");
	hCutFlow.GetXaxis().SetBinLabel(8,"N_{aj}");
	hCutFlow.GetXaxis().SetBinLabel(9,"M_{jj}");*/
	hCutFlow.SetBinContent(1,N_HLT );
	hCutFlow.SetBinContent(2,Npreselect );
	hCutFlow.SetBinContent(3,NpT_jj );
	hCutFlow.SetBinContent(4,NpT_Z );
	hCutFlow.SetBinContent(5,NCSV1 );
	hCutFlow.SetBinContent(6,NCSV2 );
	hCutFlow.SetBinContent(7,NdPhi );
	hCutFlow.SetBinContent(8,N_Naj );
	hCutFlow.SetBinContent(9,N_Mjj );
	hCutFlow.Draw();
	c1.Print((directory+"/CutFlow"+suffixps).c_str());
	
	
    // Write and Close the output file.
    ofile.Write();
    ofile.Close();

    return 0;
}
