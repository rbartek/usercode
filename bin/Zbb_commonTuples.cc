/**
    @brief Example of an analysis code to read tree ntuple
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

 */


#include <UserCode/wilken/interface/treeReader.h>
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

bool debug = false;	


    // Here one can declare histograms
    // In compiled C++, it possible not to use pointer type (TH1F*) and
    // use the non-pointer type (TH1F) of ROOT object.

	//Declare Histograms

    TH1F hnJets ("hnJets",  "Number of Good Jets", 11, -0.5, 10.5);
    TH1F hnMuons  ("hnMuons","Number of Good Muons", 6, -0.5, 5.5);
    TH1F hnSV  ("hnSV","Number of Secondary Verticies", 6, -0.5, 5.5);
    TH1F hnPV  ("hnPV","Number of Primary Verticies", 6, -0.5, 5.5);
    TH1F hallhJet_pt  ("hallhJet_pt","Pt of all jets in event",		150, 0.0, 150);
    TH1F hMET		("hMET","Missing Et",		150, 0.0, 150);
    TH1F hPtb1		("hPtb1","Pt b jet with highest CSV",		150, 0.0, 150);
    TH1F hPtb2		("hPtb2","Pt b jet with second highest CSV", 150, 0.0, 150);
	TH1F hPtjj		("hPtjj","Pt of two b jets with highest CSV", 300, 0.0, 300);
	TH1F hPtmumu	("hPtmumu","Pt of two muons with highest pt", 300, 0.0, 300);
	TH1F hPtbalZH	("hPtbalZH","Pt balance of Z and H", 101, -75.0, 75);
	TH1F hPtmu1		("hPtmu1","Pt of muon with highest pt", 75, 0.0, 200);
	TH1F hPtmu2		("hPtmu2","Pt of muon with second highest pt", 75, 0.0, 150);
	TH1F hEtamu1		("hEtamu1","Eta of muon with highest pt", 75, -3, 3);
	TH1F hEtamu2		("hEtamu2","Eta of muon with second highest pt", 75, -3, 3);
	TH1F hPFRelIsomu1		("hPFRelIsomu1","PF Rel Iso of muon with highest pt", 101, 0.0001, 0.2);
	TH1F hPFRelIsomu2		("hPFRelIsomu2","Pf Rel Iso of muon with second highest pt", 101, 0.0001, 0.2);	
	TH1F hCSV1		("hCSV1","Jet with highest CSV",			50, -0.1, 1.2);
	TH1F hCSV2		("hCSV2","Jet with second highest CSV",		50, -0.1, 1.2);
	TH1F hdphiVH	("hdphiVH","Delta phi between Z and Higgs", 50, -3.5, 3.5);
	TH1F hdetaJJ	("hdetaJJ","Delta eta between two jets", 101, -3, 3);
	TH1F hdphiJJ	("hdphiJJ","Delta phi between two jets", 101, -3.5, 3.5);
	TH1F hNaj		("hNaj",  "Number of Additional Jets",		13, -2.5, 10.5);
	TH1F hMjj		("hMjj",  "Invariant Mass of two Jets",		75, 0, 200);
	TH1F hMjj_aftercuts		("hMjj_aftercuts",  "Invariant Mass of two Jets After Cuts",		75, 0, 200);
	TH1F hMmumu		("hMmumu",  "Invariant Mass of two muons",	75, 0, 200);
	TH1F hCutFlow	("hCutFlow",  "Selection",					8, 0, 8);
	TH2F hDphiDetajj("hDphiDetajj", "#Delta#phi vs #Delta#eta JJ", 101, -5, 5, 101, -5, 5);
	
if (debug) std::cout << "all histograms declared " << std::endl;


    // Loop over all events.
    // For other methods to access event/navigate through the sample,
    // see the documentation of RootTreeReader class.
	int event =0, Npreselect =0, NpT_jj =0, NpT_Z =0, NCSV1 =0, NCSV2 =0;
	int NdPhi =0, N_Naj =0, N_Mjj =0;
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

int njets =0;

for (int k=0;k<sample.nhJets;k++){
	if (debug) cout << "for "<< k << "th pass" << endl;
	indexedJetPt.push_back(std::pair<size_t,double>(k,(double) sample.hJet_pt[k]));
	hallhJet_pt.Fill(sample.hJet_pt[k]); 
	njets++;
	indexedJetCSV.push_back(std::pair<size_t,double>(k,(double) sample.hJet_csv[k]));
	if (debug)	cout << "CSV discriminator: " << sample.hJet_csv[k] << endl;

}// end k for loop
		for (int a=0;a<sample.naJets;a++){
			hallhJet_pt.Fill(sample.aJet_pt[a]);
			njets++; 
		}

		
if (event == 11) for (size_t i = 0 ; (i != indexedJetCSV.size()) ; ++i) {cout << "Unsorted jet CSV: " << indexedJetCSV[i].second << endl; }

std::sort(indexedJetPt.begin(),indexedJetPt.end(),::IndexedQuantityGreaterThan<double>);
std::sort(indexedJetCSV.begin(),indexedJetCSV.end(),::IndexedQuantityGreaterThan<double>);

for (size_t i = 0 ; (i != indexedJetPt.size()) ; ++i) {   PtSortedJetIndex.push_back(indexedJetPt[i].first);        }
for (size_t i = 0 ; (i != indexedJetCSV.size()) ; ++i) {   CSVSortedJetIndex.push_back(indexedJetCSV[i].first);        }

if (event == 11) for (size_t i = 0 ; (i != indexedJetCSV.size()) ; ++i) {cout << "Sorted jet CSV: " << indexedJetCSV[i].second << endl; }
		
//double dPhiVH = deltaphi(sample.H_phi, sample.V_phi);
				
double dummy = -99.99;		
		if (sample.nhJets > 1) { 
		if (sample.hJet_csv[1] > -1) { 	
		hPtb1.Fill(dummy);
			hPtb1.Fill(sample.hJet_pt[indexedJetCSV[0].first]);
			hPtb2.Fill(sample.hJet_pt[indexedJetCSV[1].first]);
			hCSV1.Fill(sample.hJet_csv[indexedJetCSV[0].first]);
			hCSV2.Fill(sample.hJet_csv[indexedJetCSV[1].first]);
			hdetaJJ.Fill(sample.hJet_eta[indexedJetCSV[0].first]-sample.hJet_eta[indexedJetCSV[1].first]);
			hdphiJJ.Fill(sample.hJet_phi[indexedJetCSV[0].first]-sample.hJet_phi[indexedJetCSV[1].first]);
			hDphiDetajj.Fill((sample.hJet_eta[indexedJetCSV[0].first]-sample.hJet_eta[indexedJetCSV[1].first]),(sample.hJet_phi[indexedJetCSV[0].first]-sample.hJet_phi[indexedJetCSV[1].first]));
		}
			hMjj.Fill(sample.H_mass);
		    hPtjj.Fill(sample.H_pt);
			
			if ((sample.vLepton_type[0] == 13) && (sample.vLepton_type[1] == 13)){
			if (sample.nvlep > 1) {
			hdphiVH.Fill(sample.HVdPhi);	
			hPtbalZH.Fill(sample.H_pt-sample.V_pt);
			}
		}
	if (debug) cout << "halfway through filling histos" << endl;
		if (sample.nhJets > -1) { hnJets.Fill(njets); }
		if (sample.nvlep > -1) { hnMuons.Fill(sample.nvlep); }
		hnSV.Fill(sample.nSvs);
		hnPV.Fill(sample.nPVs);
		hNaj.Fill(sample.naJets);
		hMET.Fill(sample.MET_sumet);
		if (sample.nvlep > 1) {
		hMmumu.Fill(sample.V_mass);
		hPtmumu.Fill(sample.V_pt);
			hPtmu1.Fill(sample.vLepton_pt[0]);
			hPtmu2.Fill(sample.vLepton_pt[1]);		
			hEtamu1.Fill(sample.vLepton_eta[0]);
			hEtamu2.Fill(sample.vLepton_eta[1]);
			hPFRelIsomu1.Fill(sample.vLepton_pfCombRelIso[0]);
			hPFRelIsomu2.Fill(sample.vLepton_pfCombRelIso[1]);
		}

if (debug) cout << "done filling histograms for event: " << event << endl;

//Calculate cut efficiencies
	if ((sample.nhJets > 1) && (sample.nvlep > 1)){
	if (sample.hJet_pt[indexedJetCSV[0].first] >=20 && sample.hJet_pt[indexedJetCSV[1].first] >=20 && sample.vLepton_pt[1] > 20) {
		if (sample.hJet_csv[indexedJetCSV[0].first] > -1 && sample.hJet_csv[indexedJetCSV[1].first] > -1 ){
			if (sample.V_mass >75 && sample.V_mass<105){
	Npreselect++;
		if(sample.H_pt > 100) {
		NpT_jj++;
			if(sample.V_pt > 100) {
			NpT_Z++;
			   if(sample.hJet_csv[indexedJetCSV[0].first] > 0.898){
			   NCSV1++;
			   if(sample.hJet_csv[indexedJetCSV[1].first]>0.5){
			   NCSV2++;
			   if((sample.HVdPhi>=2.90)||(sample.HVdPhi<-2.90)){
			   NdPhi++;
			   if(sample.nhJets < 4){
			   N_Naj++;
			   hMjj_aftercuts.Fill(sample.H_mass);
			   if((sample.H_mass>=95)&&(sample.H_mass<=125)){N_Mjj++; }//higgs mass window
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
    } while (sample.nextEvent());
if (debug){	std::cout << "Number of events " << event << endl;}
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
	hdphiJJ.Draw();
	c1.Print((directory+"/dphiJJ"+suffixps).c_str());	

	c1.Clear(); // don't create a new canvas
	hDphiDetajj.Draw();
	c1.Print((directory+"/DphivsDeta_JJ"+suffixps).c_str());	

	c1.Clear(); // don't create a new canvas
	hMjj_aftercuts.Draw();
	c1.Print((directory+"/Mjj_aftercuts"+suffixps).c_str());
	
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
	
	
    // Write and Close the output file.
    ofile.Write();
    ofile.Close();

    return 0;
}
