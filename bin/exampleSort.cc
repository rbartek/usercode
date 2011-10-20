/**
    @file exampleAnalysis.cc

    @brief Example of an analysis code to read hbb ntuple
    and create histogram
    Usage:

    \verbatim
    mkdir plots/[directoryName]
    exampleAnalysis <inputfile> [outputfilename] [directoryName]
    \endverbatim

    @param inputfile Either a ROOT file or an ASCII file containing list of
    ROOT files.

    @param outputfile Name of the ROOT file which contains the histogram.
    Defaulted to 'output.root'

    @author Rachel Wilken <rachel.wilken@cern.ch>

    @date Thu Oct 20 2011

    @version $Id: exampleAnalysis.cc,v 1.3 2011/10/20 16:49:30 wilken Exp $
 */


#include <UserCode/wilken/interface/hbbReader.h>
#include <UserCode/wilken/interface/GenericTool.h>
#include <string>
#include <iostream>
#include <vector>
#include <algorithm>
#include <map>

#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TCanvas.h>
#include "TString.h"
#include <TStyle.h>

using namespace std;



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
		hjetPt[i] = new TH1F(("jetPt"+Num[i]).c_str(), ("p_{T} of " +Num[i]+ " Jet").c_str(), 50, 0., 100.);
		hjetCSV[i] = new TH1F(("jetCSV"+Num[i]).c_str(), ("CSV of " +Num[i]+ " Jet").c_str(), 50, -1.5, 1.5);
	}// end i loop

    TH1F halljetPt  ("halljetPt","Pt of all jets in event", 125, 0.0, 150);

if (debug) std::cout << "all histograms declared " << std::endl;


    // Loop over all events.
    // For other methods to access event/navigate through the sample,
    // see the documentation of RootTreeReader class.
	int event =0;
	double Zmass = 91.1976;


    do {
event++;
	if (debug) std::cout << "entered event loop " << event << std::endl;
	std::vector< std::pair<size_t,double> > indexedJetPt;
        std::vector<size_t> PtSortedJetIndex;
	std::vector< std::pair<size_t,double> > indexedJetCSV;
        std::vector<size_t> CSVSortedJetIndex;

        // Analysis loop.
        // One can access the ntuple leaves directly from sample object

for (int k=0;k<9;k++){
	if (debug) cout << "for "<< k << "th pass" << endl;
	if (sample.jetPt[k] > -100) { 
	indexedJetPt.push_back(std::pair<size_t,double>(k,(double) sample.jetPt[k]));
	halljetPt.Fill(sample.jetPt[k]); 
	}// end if
	if (sample.bDisc_CSV[k] > -100) { 
	indexedJetCSV.push_back(std::pair<size_t,double>(k,(double) sample.bDisc_CSV[k]));
	}// end if
}// end k for loop

if (event < 10 )for (size_t i = 0 ; (i != indexedJetCSV.size()) ; ++i) {cout << "Unsorted jet CSV: " << indexedJetCSV[i].second << endl; }

std::sort(indexedJetPt.begin(),indexedJetPt.end(),::IndexedQuantityGreaterThan<double>);
std::sort(indexedJetCSV.begin(),indexedJetCSV.end(),::IndexedQuantityGreaterThan<double>);

for (size_t i = 0 ; (i != indexedJetPt.size()) ; ++i) {   PtSortedJetIndex.push_back(indexedJetPt[i].first);        }
for (size_t i = 0 ; (i != indexedJetCSV.size()) ; ++i) {   CSVSortedJetIndex.push_back(indexedJetCSV[i].first);        }

if (event < 10 )for (size_t i = 0 ; (i != indexedJetCSV.size()) ; ++i) {cout << "Sorted jet CSV: " << indexedJetCSV[i].second << endl; }

for (size_t i = 0 ; i != indexedJetPt.size() && (i < 8); ++i) {
				hjetPt[i+1]->Fill(sample.jetPt[indexedJetPt[i].first]);
				hjetCSV[i+1]->Fill(sample.bDisc_CSV[indexedJetCSV[i].first]);
			}


if (debug) cout << "done filling histograms for event: " << event << endl;

    } while (sample.nextEvent());
if (debug){	std::cout << "Number of events " << event << endl;}

    // Here one can create canvas, draw and print the histogram.
    TCanvas c1("c1","c1");
    c1.cd();
    c1.SetFillColor(kWhite);
    string suffixps = ".gif";
    
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


    // Write and Close the output file.
    ofile.Write();
    ofile.Close();

    return 0;
}
