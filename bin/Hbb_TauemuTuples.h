#ifndef Zbb_TauemuTuples_h
#define Zbb_TauemuTuples_h


#include <UserCode/wilken/interface/treeReader.h>
#include <UserCode/wilken/interface/GenericTool.h>
#include "UserCode/wilken/interface/btagshapeNew.h"
//#include "UserCode/wilken/interface/xsectauFall11.h"
#include "UserCode/wilken/interface/xsecPUweightedTauCand.h"
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
#include "TMatrixDSym.h"
#include "TDecompSVD.h"
#include "TF1.h"

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
void     LeptonDistributions(std::string cut, double ph_weight);
void     JetDistributions(std::string cut, double ph_weight);
void     TH2FDistributions(std::string cut, double ph_weight);
void     EventShapeDistributions(std::string cut, double ph_weight);
void     EventDistributions(std::string cut, double ph_weight);
double   SetWeight( std::string filename);
bool     findString(std::string strToSearch, std::string strPattern);

double lumi = 4.982;

double weight, LumiWeight;
int nJets, nSV, nMuons, nElectrons, nLeptons, Na_lep, nPV, Naj, eventFlavor, naJets, Nab;
float Zphi, Hphi;
float jetPt[5], jetEta[5], jetPhi[5], jetCSV[5], jetCHF[5], CSVNewShape[5];
float RegjetPt[5], jetPtRaw[5], jetE[5], jetVtx3dL[5], jetVtx3deL[5], jetVtxPt[5], jetVtxMass[5], jetPtLeadTrack[5];
float jetNconstintuents[5], jetCEF[5], jetNCH[5], jetJECUnc[5], jetEt[5], jetMt[5], jetPtRawJER[5], jetGenPt[5]; 
float SVnumTracks[5], jetNHF[5];
float leptonPt[5], leptonEta[5], leptonPhi[5], lep_pfCombRelIso[5], lep_id95[5], SortedMuonPt[5];
int lep_flavor[5];
float DetaJJ, btag2CSF;
float CSV0, CSV1, Emumass, Hmass, DeltaPhiHV, Hpt, Zpt, DitauMass; 
float pZeta25, pZeta45, pZeta65, pZeta85;
float oldHmass, RegHmass, GenHiggsMass;
float GenZMass;
float Rho25, MinMET, MindPhiMEThJetOR30aJet, SVmass[5];
float lep0pt, ScalarSumPt, EtaStandDev, UnweightedEta, EvntShpCircularity;
float alpha_j, qtb1, DphiJJ, Trigweight, B2011PUweight, A2011PUweight; 
float RMS_eta, PtbalZH, EventPt, AngleHemu, Centrality, MET;
float alpha_lep, qtlep1;
float EvntShpAplanarity, EvntShpSphericity, EvntShpIsotropy;
float dPhiHMET, Mt, DeltaPhijetMETmin, DeltaPhijetMETZtaumin;
float SV_mass, Mte, Mtmu, Ht, delPullAngle, delPullAngle2;
float AngleEMU, CosThetaEle, CosThetaMu, EleMissE, MuonMissE, Dphiemu, AaronEleMissE, AaronMuMissE;
float Zmass, ZmassSVD, ZmassSVDnegSol, ZmassNegInclu, ZmassMatrix;
float delRjj, Detaemu, DphiEleMET, dphiMuMET, PtbalMETH, topMass, topPt, topWmass;
float MassEleb0, MassMub0, MassEleb1, MassMub1, METsig;
float delRemu, PtbalZMET, DphiZMET;
float DphiLeadMET, DphiSecondMET;
float ProjVisT, ProjMissT;
float ScalarSumJetPt, ScalarSumHiggsJetPt, EventMass;
float PUweight2011, EleTrigWeight, MuonTrigWeight, WP95weight, EleRecoWeight, MuIDweight;
bool isdata, isDATA;

float CSVup[5], CSVdown[5], csvFup[5], csvFdown[5];
float JER_e_up[5], JER_e_down[5], JES_e_up[5], JES_e_down[5];
float JER_pt_up[5], JER_pt_down[5], JES_pt_up[5], JES_pt_down[5];
float JERHmassUp, JERHmassDown, JESHmassUp, JESHmassDown;

// maps for histo function
std::map<std::string,TH1*> histmap;
std::map<std::string,TH2*> bidimhistmap;


#endif // #ifdef Zbb_commonTuple_cc



