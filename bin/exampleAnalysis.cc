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

    // Here one can declare histograms
    // In compiled C++, it possible not to use pointer type (TH1F*) and
    // use the non-pointer type (TH1F) of ROOT object.

    TH1F hnJets ("hnJets",  "Number of Good Jets", 11, -0.5, 10.5);
    TH1F hnMuons  ("hnMuons","Number of Good Muons", 6, -0.5, 5.5);
    TH1F hjetPt  ("hjetPt","Pt of all jets in event", 125, 0.0, 250);

    // Loop over all events.
    // For other methods to access event/navigate through the sample,
    // see the documentation of RootTreeReader class.
    do {
        // Analysis loop.

        // One can access the ntuple leaves directly from sample object
            if (sample.nJets > -1) { hnJets.Fill(sample.nJets); }
            if (sample.nMuons > -1) { hnMuons.Fill(sample.nMuons); }

for (int k=0;k<9;k++){
    if (sample.jetPt[k] > -100) { hjetPt.Fill(sample.jetPt[k]); }
}// end k for loop

    } while (sample.nextEvent());

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

    c1.Clear(); // don't create a new canvas
    hjetPt.Draw();
    c1.Print((directory+"/jetPt"+suffixps).c_str());


    // Write and Close the output file.
    ofile.Write();
    ofile.Close();

    return 0;
}
