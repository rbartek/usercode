/**
    @file templateAnalysis.cc

    @brief Template of analysis code to read hbb ntuple
    and create histogram.
    Usage:

    \verbatim
    templateAnalysis <inputfile> [outputfile]
    \endverbatim

    @param inputfile Either a ROOT file or an ASCII file containing list of
    ROOT files.

    @param outputfile Name of the ROOT file which contains the histogram.
    Defaulted to 'output.root'

    @author Rachel Wilken <rachel.wilken@cern.ch>

    @date Oct 19 2011

    @version $Id: templateAnalysis.cc,v 1.1 2011/10/19 wilken Exp $
 */


#include <UserCode/wilken/interface/hbbReader.h>
#include <UserCode/wilken/interface/GenericTool.h>
#include <string>
#include <iostream>
#include <vector>
#include <algorithm>
#include <map>

// Include the appropriate ROOT header file(s) here
#include <TFile.h>
#include <TH1.h>
#include <TCanvas.h>


int main(int argc, char** argv) {

    // There must be at least one argument, the input filename.
    if (argc < 2) {
        std::cout << "Input filename is not specified! Exiting!" << std::endl;
        exit(1);
    }

    // First argument: a ROOT file or an ASCII file containing list of
    // ROOT files.
    std::string ifilename(argv[1]);
    std::string ofilename("output.root");

    if (argc == 3) {
        // If there is a second argument, use it for the output filename.
        ofilename = std::string(argv[2]);
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

    // Here one can declare histograms

    // Loop over all events.
    // For other methods to access event/navigate through the sample,
    // see the documentation of RootTreeReader class.
    do {
        // One can access the ntuple leaves directly from sample object
        std::cout << "Run: " << sample.runNumber << std::endl;
        std::cout << "LumiSection: " << sample.lumiBlock << std::endl;
        std::cout << "Event: " << sample.evtNumber << std::endl;

    } while (sample.nextEvent());

    // Here one can create canvas, draw and print the histogram.
    TCanvas c1("c1","c1");
    c1.cd();

    // Write and Close the output file.
    ofile.Write();
    ofile.Close();

    return 0;
}

