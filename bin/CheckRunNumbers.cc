/**
 @brief Example of an analysis code to read tree ntuple
 and create histogram
 Usage:
 
 \verbatim
 mkdir [directoryName]
 SkeletonReader <inputfile> [outputfilename] [directoryName]
 \endverbatim
 
 @param inputfile Either a ROOT file or an ASCII file containing list of
 ROOT files.
 
 @param outputfile Name of the ROOT file which contains the histogram.
 Defaulted to 'output.root'
 
 @author Rachel Wilken <rachel.wilken@cern.ch>
 
 @date Fri Nov 9 2011
 
 */

#include "Hbb_TauemuTuples.h"
#include <iostream>
#include <sstream>
#include <fstream>

using namespace std;


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
    std::string ofilename("RunNumbers.root");
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
	
	BTagShapeNew* btagNew;
	btagNew = new BTagShapeNew("csvdiscrNew.root");
	btagNew->computeFunctions();
 
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
	
	
	//open file
//    string inname = "Aug05json.txt";
    string inname = "Run2011Prompt.txt";
  //  string inname = "May10json.txt";
    ifstream jsonfile(inname);
	
    if (!jsonfile) {
        cout << "There was a problem opening file " << jsonfile << " for reading." << endl;
        return 0;
    }
	std::vector<int> ListOfRuns;
	std::vector< std::pair<int,int> > ListOfLumis;

int i;	
int TheRunIs = 0;
int FirstLumi = 0;
    while (jsonfile >> i) {
		// cout << "Value from file is " << i << endl;
		//PromptV4
		//if (i<163870 && i > 10000) i = 160000;
		//if (i>170248) i = 160000;
		//PromptV6
		//if (i<172620 && i > 10000) i = 160000;
		//if (i>175859) i = 160000;
		//Run2011B Prompt
		if (i< 175860 && i > 10000) i = 160000;
		if (i > 160404) {
			ListOfRuns.push_back(i);
			TheRunIs = i;
			FirstLumi = 0;
		}
		//cout << "The run is " << TheRunIs << " i is " << i << endl;
		if (i < 10000 && TheRunIs > 160404){
			if (FirstLumi==0){
				ListOfLumis.push_back(std::pair<int,int>(i,(int)TheRunIs));
				FirstLumi = i;
			} else {
				for (int k = FirstLumi+1; k<=i;k++){
					ListOfLumis.push_back(std::pair<int,int>(k,(int)TheRunIs));
				}
				FirstLumi = 0;
			}
		}
    }
	
	
   // Here one can declare histograms
    // In compiled C++, it possible not to use pointer type (TH1F*) and
    // use the non-pointer type (TH1F) of ROOT object.
	
	//Declare Histograms
	

	TH2F *hRunNumberVsnPU = new TH2F	("hRunNumberVsnPU",  "hRunNumberVsnPU",100,16404,180252,20,0,21);
	TH1F *hRunNumber = new TH1F	("hRunNumber",  "RunNumber",200,16404,180252);
		
	if (debug) std::cout << "all histograms declared " << std::endl;
	
    // Loop over all events.
    // For other methods to access event/navigate through the sample,
    // see the documentation of RootTreeReader class.
	int  JsonFalse =0, Duplicates =0, JsonRuns =0, Lumis = 0;
	int event =0;
	int MatchedRun=0;
	bool firstevent = true;
	std::vector<int> EventNumbers;
	std::vector<int> FileRunList;
	//	double Emumass = 91.1976;
	JsonRuns = ListOfRuns.size();
	Lumis = ListOfLumis.size();
    do {
		if(sample.Vtype==5){
			event++;
			if (sample.EVENT_json == 0) JsonFalse++;
			hRunNumberVsnPU->Fill(sample.EVENT_run,sample.nPVs);
			hRunNumber->Fill(sample.EVENT_run);
			EventNumbers.push_back(sample.EVENT_event);
			bool newRun = true;
			for(unsigned int lumi_iter = 0; lumi_iter<ListOfLumis.size();lumi_iter++){
				if (fabs(ListOfLumis[lumi_iter].second-sample.EVENT_run)<0.1){
					if (fabs(ListOfLumis[lumi_iter].first-sample.EVENT_lumi)<0.1){
						ListOfLumis.erase(ListOfLumis.begin()+lumi_iter);
					}
				}
			}
			for(unsigned int iter = 0; iter<FileRunList.size();iter++){
			if(fabs(FileRunList[iter]-sample.EVENT_run)<0.1) newRun = false;
			}
			if(newRun) {
				FileRunList.push_back(sample.EVENT_run);
				bool matched_run = false;
				for(unsigned int k = 0; k<ListOfRuns.size();k++){
				 if (fabs(ListOfRuns[k]-sample.EVENT_run)<0.1){ 
						MatchedRun++;
						matched_run = true;
						ListOfRuns.erase(ListOfRuns.begin()+k);
					}
				}
				if(!matched_run&&sample.EVENT_json == 1) cout << "Run in root file not matched " << sample.EVENT_run << endl;
			}	// if new run
		}//if Vtype
	} while (sample.nextEvent());
	
	
	//Check for duplicate events
	for(unsigned int i = 0; i<EventNumbers.size();i++){
		for(unsigned int j = 0; j<EventNumbers.size();j++){
			if (i==j) continue;
			if (i<j) continue;//only count duplicate once
			if (EventNumbers[i]==EventNumbers[j]) {
			Duplicates++;
			cout << "Duplicate Event: " << EventNumbers[i] << " at position " << i << " and " << j << endl;
			}
	}
}
	for(unsigned int a = 0; a<ListOfRuns.size();a++){
		cout << "Missing runs from json " << ListOfRuns[a] << endl;
	}
/*	for(unsigned int a = 0; a<ListOfLumis.size();a++){
		cout << "Missing lumis from json " << ListOfLumis[a].second<< " " << ListOfLumis[a].first << endl;
	}*/
	cout << "Number Missing lumis from json " << ListOfLumis.size() << endl;
	cout << "Number lumis from json " << Lumis << endl;

	
	std::cout << endl << endl;
	std::cout << "Number of events " << event << endl;
	std::cout << "JsonFalse " << JsonFalse << endl;
	std::cout << "Duplicate Events " << Duplicates << endl;
	std::cout << "Number of Runs in File " << FileRunList.size() << endl;
	std::cout << "Number of Runs in JSON " << JsonRuns << endl;
	std::cout << "Number of Runs Matched " << MatchedRun << endl;
	// Here one can create canvas, draw and print the histogram.
	TCanvas c1("c1","c1");
	c1.cd();
	c1.SetFillColor(kWhite);
	string suffixps = ".gif";
		
	c1.Clear(); // don't create a new canvas
	hRunNumberVsnPU->Draw();
	c1.Print((directory+"/RunNumberVsnPU"+suffixps).c_str());
	
	c1.Clear(); // don't create a new canvas
	hRunNumber->Draw();
	c1.Print((directory+"/RunNumber"+suffixps).c_str());
			
	// Write and Close the output file.
	ofile.Write();
	ofile.Close();
	
	return 0;
}

bool findString(std::string strToSearch, std::string strPattern){
	size_t found;
	
	bool foundStr = false;
	
	found=strToSearch.find(strPattern);
	if (found!=string::npos)
		foundStr = true;
	
	return foundStr;		
}





		





