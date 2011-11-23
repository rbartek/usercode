#ifndef Zbb_commonTuples_h
#define Zbb_commonTuples_h


#include <UserCode/wilken/interface/treeReader.h>
#include <UserCode/wilken/interface/GenericTool.h>
#include "PhysicsTools/CandUtils/interface/EventShapeVariables.h"

#include <string>
#include <iostream>
#include <vector>
#include <algorithm>
#include <map>
#include <cstring>
#include <string>

#include "TLorentzVector.h"
#include "TVector3.h"
#include "TMath.h"

#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TCanvas.h>
//#include "TString.h"
#include <TStyle.h>
#include <map>

// needed for use of C++ map //////////////////////////
#ifdef __MAKECINT__
#pragma link C++ class std::map<std::string, TH1F*>+;
#endif
//////////////////////////////////////////////////////

//virtual void     fillhisto(std::string histname,double val, double weight, std::string title, int nbins, double min, double max);
//virtual void     fill2Dhisto(string histname,double valx, double valy, double weight, string title, int nxbins, double xmin, double xmax, int nybins, double ymin, double ymax);
//virtual void     MuonDistributions(string cut, double ph_weight);
//virtual void     JetDistributions(string cut, double ph_weight);
//virtual void     EventDistributions(string cut, double ph_weight);
double   SetWeight( std::string filename);
bool     findString(std::string strToSearch, std::string strPattern);

// maps for histo function
//std::map<std::string,TH1*> histmap;
//std::map<std::string,TH2*> bidimhistmap;



/*void fillhisto(string histname,double val, double weight, string title, int nbins, double min, double max) {
	if (!histmap[histname]) {
		if (nbins==0) {
			cout << " ERROR: must define histogram " << histname << " properties on first call of fillhisto(...)" << endl;
			exit(1);
		}
		
		histmap[histname]=new TH1F(histname.c_str(),title.c_str(),nbins,min,max);
		
	}
	histmap[histname]->Fill(val, weight);
	
} // finish fillhisto
void fill2Dhisto(string histname,double valx, double valy, double weight, string title, int nxbins, double xmin, double xmax, int nybins, double ymin, double ymax) {
	if (!bidimhistmap[histname]) {
		if (nxbins==0 || nybins==0) {
			//cout << " ERROR: must define bi dim histogram " << histname << " properties on first call of fillbidimhisto(...)" << endl;
			exit(1);
		}
		fill2Dhisto[histname]=new TH2F(histname.c_str(),title.c_str(),nxbins,xmin,xmax,nybins,ymin,ymax);
		
	}
	fill2Dhisto[histname]->Fill(valx,valy,weight);
}// fill2Dhisto
*/




#endif // #ifdef Zbb_commonTuple_cc

