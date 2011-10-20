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

    @version $Id: exampleAnalysis.cc,v 1.1 2011/10/20 15:36:46 wilken Exp $
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
#include <TCanvas.h>



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

    TH1F hRECOELE_PT ("hRECOELE_PT",
                      "RECO electron p_{T}",
                      200,0.0,100.0);
    TH1F hRECOMU_PT  ("hRECOMU_PT",
                      "RECO muon p_{T}",
                      200,0.0,100.0);

    // Loop over all events.
    // For other methods to access event/navigate through the sample,
    // see the documentation of RootTreeReader class.
    do {
        // Analysis loop.

        // One can access the ntuple leaves directly from sample object
        for (size_t i = 0 ; i != 100 ; ++i) {
            // 100 is the size of RECOELE_PT and RECOMU_PT array
            if (sample.RECOELE_PT[i] > 0.0) {
                hRECOELE_PT.Fill(sample.RECOELE_PT[i]);
            }
            if (sample.RECOMU_PT[i] > 0.0) {
                hRECOMU_PT.Fill(sample.RECOMU_PT[i]);
            }
        }

    } while (sample.nextEvent());

    // Here one can create canvas, draw and print the histogram.
    TCanvas c1("c1","c1");
    c1.cd();

    c1.SetFillColor(kWhite);

    hRECOELE_PT.Draw();
    c1.Print("hRECOELE_PT.eps");

    c1.Clear(); // don't create a new canvas
    hRECOMU_PT.Draw();
    c1.Print("hRECOMU_PT.eps");

    // Write and Close the output file.
    ofile.Write();
    ofile.Close();

    return 0;
}
