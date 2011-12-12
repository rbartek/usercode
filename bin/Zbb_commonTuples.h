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

void     fillhisto(std::string histname, double val, double weight, std::string title, int nbins, double min, double max);
void     fill2Dhisto(std::string histname, double valx, double valy, double weight, std::string title, int nxbins, double xmin, double xmax, int nybins, double ymin, double ymax);
void     MuonDistributions(std::string cut, double ph_weight);
void     JetDistributions(std::string cut, double ph_weight);
void     TH2FDistributions(std::string cut, double ph_weight);
void     EventShapeDistributions(std::string cut, double ph_weight);
void     EventDistributions(std::string cut, double ph_weight);
double   SetWeight( std::string filename);
bool     findString(std::string strToSearch, std::string strPattern);

double weight;
int nJets, nSV, nMuons, Na_mu, nPV, Naj;
float jetPt[10], jetEta[10], jetPhi[10], jetCSV[10], jetCHF[10];
float muonPt[10], muonEta[10], muonPhi[10], muonPFiso[10];
float DetaJJ;
float CSV0, CSV1, Zmass, Hmass, DeltaPhiHV, Hpt, Zpt; 
float mu0pt, Ht, EtaStandDev, UnweightedEta, EvntShpCircularity;
float alpha_j, qtb1, DphiJJ, Trigweight; 
double RMS_eta, PtbalZH, EventPt, Angle, Centrality, MET;
double alpha_mu, qtmu1;
double EvntShpAplanarity, EvntShpSphericity, EvntShpIsotropy;


// maps for histo function
std::map<std::string,TH1*> histmap;
std::map<std::string,TH2*> bidimhistmap;


#endif // #ifdef Zbb_commonTuple_cc



