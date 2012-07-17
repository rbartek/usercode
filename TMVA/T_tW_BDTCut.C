/**********************************************************************************
 * Project   : TMVA - a Root-integrated toolkit for multivariate data analysis    *
 * Package   : TMVA                                                               *
 * Exectuable: T_tW_BDTCut                                     *
 *                                                                                *
 * This macro provides a simple example on how to use the trained classifiers     *
 * within an analysis module                                                      *
 **********************************************************************************/

#include <cstdlib>
#include <vector>
#include <iostream>
#include <map>
#include <string>
#include "/home/hep/wilken/taus/CMSSW_4_4_2_patch8/src/UserCode/wilken/interface/xsecV21.h" 
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TStopwatch.h"

#include "TMVAGui.C"

#if not defined(__CINT__) || defined(__MAKECINT__)
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/MethodCuts.h"
#endif

using namespace TMVA;

void T_tW_BDTCut( TString myMethodList = "" ) 
{   
#ifdef __CINT__
   gROOT->ProcessLine( ".O0" ); // turn off optimization in CINT
#endif

   //---------------------------------------------------------------

   // This loads the library
   TMVA::Tools::Instance();

   // Default MVA methods to be trained + tested
   std::map<std::string,int> Use;

   // --- Cut optimisation
   Use["Cuts"]            = 0;
   Use["CutsD"]           = 0;
   Use["CutsPCA"]         = 0;
   Use["CutsGA"]          = 0;
   Use["CutsSA"]          = 0;
   // 
   // --- 1-dimensional likelihood ("naive Bayes estimator")
   Use["Likelihood"]      = 0;
   Use["LikelihoodD"]     = 0; // the "D" extension indicates decorrelated input variables (see option strings)
   Use["LikelihoodPCA"]   = 0; // the "PCA" extension indicates PCA-transformed input variables (see option strings)
   Use["LikelihoodKDE"]   = 0;
   Use["LikelihoodMIX"]   = 0;
   //
   // --- Mutidimensional likelihood and Nearest-Neighbour methods
   Use["PDERS"]           = 0;
   Use["PDERSD"]          = 0;
   Use["PDERSPCA"]        = 0;
   Use["PDEFoam"]         = 0;
   Use["PDEFoamBoost"]    = 0; // uses generalised MVA method boosting
   Use["KNN"]             = 0; // k-nearest neighbour method
   //
   // --- Linear Discriminant Analysis
   Use["LD"]              = 0; // Linear Discriminant identical to Fisher
   Use["Fisher"]          = 0;
   Use["FisherG"]         = 0;
   Use["BoostedFisher"]   = 0; // uses generalised MVA method boosting
   Use["HMatrix"]         = 0;
   //
   // --- Function Discriminant analysis
   Use["FDA_GA"]          = 0; // minimisation of user-defined function using Genetics Algorithm
   Use["FDA_SA"]          = 0;
   Use["FDA_MC"]          = 0;
   Use["FDA_MT"]          = 0;
   Use["FDA_GAMT"]        = 0;
   Use["FDA_MCMT"]        = 0;
   //
   // --- Neural Networks (all are feed-forward Multilayer Perceptrons)
   Use["MLP"]             = 0; // Recommended ANN
   Use["MLPBFGS"]         = 0; // Recommended ANN with optional training method
   Use["MLPBNN"]          = 0; // Recommended ANN with BFGS training method and bayesian regulator
   Use["CFMlpANN"]        = 0; // Depreciated ANN from ALEPH
   Use["TMlpANN"]         = 0; // ROOT's own ANN
   //
   // --- Support Vector Machine 
   Use["SVM"]             = 0;
   // 
   // --- Boosted Decision Trees using this
   Use["BDT"]             = 1; // uses Adaptive Boost
   Use["BDTG"]            = 0; // uses Gradient Boost
   Use["BDTB"]            = 0; // uses Bagging
   Use["BDTD"]            = 0; // decorrelation + Adaptive Boost
   // 
   // --- Friedman's RuleFit method, ie, an optimised series of cuts ("rules")
   Use["RuleFit"]         = 0;
   // ---------------------------------------------------------------
   Use["Plugin"]          = 0;
   Use["Category"]        = 0;
   Use["SVM_Gauss"]       = 0;
   Use["SVM_Poly"]        = 0;
   Use["SVM_Lin"]         = 0;

   std::cout << std::endl;
   std::cout << "==> Start T_tW_BDTCut" << std::endl;

   // Select methods (don't look at this code - not of interest)
   if (myMethodList != "") {
      for (std::map<std::string,int>::iterator it = Use.begin(); it != Use.end(); it++) it->second = 0;

      std::vector<TString> mlist = gTools().SplitString( myMethodList, ',' );
      for (UInt_t i=0; i<mlist.size(); i++) {
         std::string regMethod(mlist[i]);

         if (Use.find(regMethod) == Use.end()) {
            std::cout << "Method \"" << regMethod 
                      << "\" not known in TMVA under this name. Choose among the following:" << std::endl;
            for (std::map<std::string,int>::iterator it = Use.begin(); it != Use.end(); it++) {
               std::cout << it->first << " ";
            }
            std::cout << std::endl;
            return;
         }
         Use[regMethod] = 1;
      }
   }

   // --------------------------------------------------------------------------------------------------

   // --- Create the Reader object

   TMVA::Reader *reader = new TMVA::Reader( "!Color:!Silent" );    

   // Create a set of variables and declare them to the reader
   // - the variable names MUST corresponds in name and type to those given in the weight file(s) used
   Float_t Hmass, Emumass;
	Float_t Hpt, Zpt;
	Float_t CSV0, CSV1;
	Float_t DeltaPhiHV, DetaJJ;
	Int_t  nJets, eventFlavor, Naj; 
 	Float_t BDTvalue, Trigweight, B2011PUweight, A2011PUweight, btag2CSF, MET;
	Float_t alpha_j, qtb1, jetPhi0, jetPhi1, jetEta0, jetEta1, Zphi, Hphi;
	Float_t Ht, EvntShpCircularity, jetCHF0, jetCHF1;
	Float_t EtaStandDev, lep_pfCombRelIso0, lep_pfCombRelIso1, EvntShpIsotropy;
	Float_t lep0pt, lep1pt, UnweightedEta, DphiJJ, RMS_eta, EvntShpSphericity;
	Float_t PtbalZH, EventPt, AngleEMU, Centrality, EvntShpAplanarity;
	 Float_t UnweightedEta, lep0pt;
	 	Int_t naJets, nSV;
	 Float_t Mte, Mtmu, dPhiHMET, Dphiemu, delRjj, delRemu, DphiZMET, DeltaPhijetMETmin;
	 Float_t MassEleb0, DphiEleMET, PtbalMETH,dphiEleMET, dEtaJJ, dphiZMET, ScalarSumPt;
	Float_t ZmassSVD, AngleHemu, ProjVisT, topMass, topPt, ZmassSVDnegSol, Mt, Zmass, ZmassNegInclu;
	
	Float_t dphiEMU, dphiZMET;
	
	reader->AddVariable( "Hmass", &Hmass );
	//reader->AddVariable( "Naj", &Naj );
	reader->AddVariable( "CSV0", &CSV0 );
	reader->AddVariable( "Emumass", &Emumass );
	reader->AddVariable( "DeltaPhiHV", &DeltaPhiHV );
	reader->AddVariable( "Mt", &Mt );
	reader->AddVariable( "dPhiHMET", &dPhiHMET );
	reader->AddVariable( "dphiEMU := abs(Dphiemu)", &dphiEMU );
	reader->AddVariable( "dphiZMET:=abs(DphiZMET)", &dphiZMET );
	reader->AddVariable( "PtbalMETH", &PtbalMETH );
	reader->AddVariable( "EtaStandDev", &EtaStandDev );
	reader->AddVariable( "jetCHF0", &jetCHF0 );
	reader->AddVariable( "ProjVisT", &ProjVisT );

	
   // Spectator variables declared in the training have to be added to the reader, too
 
//   reader->AddSpectator( "UnweightedEta",   &UnweightedEta );
  // reader->AddSpectator( "lep0pt",   &lep0pt );

/*   Float_t Category_cat1, Category_cat2, Category_cat3;
   if (Use["Category"]){
      // Add artificial spectators for distinguishing categories
      reader->AddSpectator( "Category_cat1 := var3<=0",             &Category_cat1 );
      reader->AddSpectator( "Category_cat2 := (var3>0)&&(var4<0)",  &Category_cat2 );
      reader->AddSpectator( "Category_cat3 := (var3>0)&&(var4>=0)", &Category_cat3 );
   }
*/
   // --- Book the MVA methods

   TString dir    = "weights/";
   TString prefix = "TMVAClassification";

   // Book method(s)
   for (std::map<std::string,int>::iterator it = Use.begin(); it != Use.end(); it++) {
      if (it->second) {
         TString methodName = TString(it->first) + TString(" method");
         TString weightfile = dir + prefix + TString("_") + TString(it->first) + TString(".weights.xml");
         reader->BookMVA( methodName, weightfile ); 
      }
   }
   
	// Book output histograms
   	TH1F* hCHFb0_OpenSelection= new TH1F		("hCHFb0_OpenSelection", "charged Hadron Energy Fraction b1", 40, 0.0, 1.2);
   	TH1F* hCHFb1_OpenSelection= new TH1F		("hCHFb1_OpenSelection", "charged Hadron Energy Fraction b2", 40, 0.0, 1.2);
 	TH1F* hPtjj_OpenSelection= new TH1F		("hPtjj_OpenSelection","Pt of two b jets with highest CSV ", 50, 0.0, 400);
 	TH1F* hPtmumu_OpenSelection= new TH1F	("hPtmumu_OpenSelection","Pt of two muons with highest pt ", 50, 0.0, 400);
	TH1F* hPtbalZH_OpenSelection= new TH1F	("hPtbalZH_OpenSelection", "Pt balance of Z and H", 40, -80, 80);
	TH1F* hPtmu0_OpenSelection= new TH1F		("hPtmu0_OpenSelection","Pt of muon with highest pt ", 30, 0.0, 300);
	TH1F* hPtmu1_OpenSelection= new TH1F		("hPtmu1_OpenSelection","Pt of muon with second highest pt ", 30, 0.0, 300);
	TH1F* hPFRelIsomu0_OpenSelection= new TH1F		("hPFRelIsomu0_OpenSelection", "PF Rel Iso of muon with highest Pt", 40, 0., 0.2);
	TH1F* hPFRelIsomu1_OpenSelection= new TH1F		("hPFRelIsomu1_OpenSelection", "PF Rel Iso of muon with second highest Pt", 40, 0., 0.2);
	TH1F* hCSV0_OpenSelection= new TH1F		("hCSV0_OpenSelection","Jet with highest CSV ",			40, 0, 1.5);
	TH1F* hCSV1_OpenSelection= new TH1F		("hCSV1_OpenSelection","Jet with second highest CSV ",		40, 0, 1.5);
	TH1F* hdphiVH_OpenSelection= new TH1F	("hdphiVH_OpenSelection","Delta phi between Z and Higgs ", 50, -0.1, 4.5);
	TH1F* hdetaJJ_OpenSelection= new TH1F	("hdetaJJ_OpenSelection","Delta eta between two jets ", 60, -4, 4);
	TH1F* hNjets_OpenSelection= new TH1F	("hNjets_OpenSelection", "Number of Jets",		13, -2.5, 10.5);
	TH1F* hMjj_OpenSelection	= new TH1F	("hMjj_OpenSelection",  "Invariant Mass of two Jets ",		50, 0, 300);
	TH1F* hMmumu_OpenSelection	= new TH1F	("hMmumu_OpenSelection",  "Invariant Mass of two muons ",	75, 0, 200);
    TH1F* hRMSeta_OpenSelection= new TH1F	("hRMSeta_OpenSelection", "RMS Eta",		30, 0, 3);
    TH1F* hStaDeveta_OpenSelection= new TH1F	("hStaDeveta_OpenSelection", "Standard Deviation Eta",		30, 0, 3);
	TH1F* hUnweightedEta_OpenSelection= new TH1F	("hUnweightedEta_OpenSelection",  "Unweighted Eta ",		50, 0, 15);		
	TH1F* hdphiJJ_vect_OpenSelection= new TH1F	("hdphiJJ_vect_OpenSelection", "Delta phi between two jets",  30, -3.5, 4);
	TH1F* hCircularity_OpenSelection= new TH1F("hCircularity_OpenSelection","EventShapeVariables circularity", 30, 0.0, 1.2);
	TH1F* hHt_OpenSelection= new TH1F("hHt_OpenSelection","scalar sum of pt of four particles", 50, 0.0, 500);
    TH1F* hCentrality_OpenSelection= new TH1F	("hCentrality_OpenSelection", "Centrality", 40, 0.0, 0.8);
    TH1F* hEventPt_OpenSelection= new TH1F	("hEventPt_OpenSelection", "Pt of HV system", 50, 0.0, 100);
    TH1F* hAngleEMU_OpenSelection= new TH1F	("hAngleEMU_OpenSelection", "AngleEMU between H and Z", 45, 0, 4.5);
    TH1F* hSphericity_OpenSelection= new TH1F	("hSphericity_OpenSelection", "EventShapeVariables sphericity", 50, 0.0, 1);
    TH1F* hAplanarity_OpenSelection= new TH1F	("hAplanarity_OpenSelection", "EventShapeVariables Aplanarity", 50, -0.1, .4);
    TH1F* hIsotropy_OpenSelection= new TH1F	("hIsotropy_OpenSelection",  "EventShapeVariables isotropy", 30, 0.0, 1.3);
    TH2F* hDphiDetajj_OpenSelection= new TH2F	("hDphiDetajj_OpenSelection", "#Delta#phi vs #Delta#eta JJ", 25, -5, 5, 25, -5, 5);
	TH1F* hMtmu_OpenSelection= new TH1F ("hMtmu_OpenSelection", "Mt with respect to Muon", 101, -0.1, 200);
	TH1F* hMte_OpenSelection= new TH1F ("hMte_OpenSelection", "Mt with respect to Electron", 101, -0.1, 200);
	TH1F* hdPhiHMET_OpenSelection= new TH1F ("hdPhiHMET_OpenSelection", "Delta phi between MET and Higgs", 50, -0.1, 4.5);
	TH1F* hDphiemu_OpenSelection= new TH1F ("hDphiemu_OpenSelection", "Delta phi between e and muon", 50, -3.5, 4.5);
	TH1F* hdelRjj_OpenSelection= new TH1F ("hdelRjj_OpenSelection", "Delta R jj", 55, 0, 5.5);
	TH1F* hdelRemu_OpenSelection= new TH1F ("hdelRemu_OpenSelection", "Delta R emu", 55, 0, 5.5);
	TH1F* hDphiZMET_OpenSelection= new TH1F ("hDphiZMET_OpenSelection", "Delta phi between Z and MET",  71, -3.5, 4);
	TH1F* hDeltaPhijetMETmin_OpenSelection= new TH1F ("hDeltaPhijetMETmin_OpenSelection", "Delta phi between MET and nearest jet", 50, -0.1, 4.5);
	TH1F* hAngleHemu_OpenSelection= new TH1F ("hAngleHemu_OpenSelection", "Angle between H and Z", 30, 0, 3.5);
	TH1F* hProjVisT_OpenSelection= new TH1F ("hProjVisT_OpenSelection", "Transverse componenet of Projection of Z onto bisector", 80, 0, 200);
	TH1F* htopMass_OpenSelection= new TH1F ("htopMass_OpenSelection", "Top Mass single lepton", 100, 75, 375);
	TH1F* htopPt_OpenSelection= new TH1F ("htopPt_OpenSelection", "Pt of Top", 125, 0, 250);
	TH1F* hVMt_OpenSelection= new TH1F ("hVMt_OpenSelection", "VMt", 75, 0, 150);
	TH1F* hZmassSVD_OpenSelection= new TH1F ("hZmassSVD_OpenSelection", "Invariant Mass of two Leptons corrected SVD", 75, 0, 150);
	TH1F* hZmassSVDnegSol_OpenSelection= new TH1F ("hZmassSVDnegSol_OpenSelection", "Invariant Mass of two Leptons corrected SVD", 100, -50, 200);
	TH1F* hZmass_OpenSelection= new TH1F ("hZmass_OpenSelection", "Zmass ", 5, 0, 150);
	TH1F* hZmassNegInclu_OpenSelection= new TH1F ("hZmassNegInclu_OpenSelection", "Invariant Mass of two Leptons corrected SVD", 100, -50, 200);
	
	
	TTree *treeWithBDT = new TTree("treeWithBDT","Tree wiht BDT output");
	treeWithBDT->SetDirectory(0);
	treeWithBDT->Branch("nJets",&nJets, "nJets/I");
	treeWithBDT->Branch("naJets",&naJets, "naJets/I");
	treeWithBDT->Branch("eventFlavor",&eventFlavor, "eventFlavor/I");
	treeWithBDT->Branch("CSV0",&CSV0, "CSV0/F");
	treeWithBDT->Branch("CSV1",&CSV1, "CSV1/F");
	treeWithBDT->Branch("Emumass",&Emumass, "Emumass/F");
	treeWithBDT->Branch("Hmass",&Hmass, "Hmass/F");
	treeWithBDT->Branch("DeltaPhiHV",&DeltaPhiHV, "DeltaPhiHV/F");
	treeWithBDT->Branch("Hpt",&Hpt, "Hpt/F");
	treeWithBDT->Branch("Zpt",&Zpt, "Zpt/F");
	treeWithBDT->Branch("lep0pt",&lep0pt, "lep0pt/F");
	treeWithBDT->Branch("Ht",&Ht, "Ht/F");
	treeWithBDT->Branch("EtaStandDev",&EtaStandDev, "EtaStandDev/F");
	treeWithBDT->Branch("UnweightedEta",&UnweightedEta, "UnweightedEta/F");
	treeWithBDT->Branch("EvntShpCircularity",&EvntShpCircularity, "EvntShpCircularity/F");
	treeWithBDT->Branch("alpha_j",&alpha_j, "alpha_j/F");
	treeWithBDT->Branch("qtb1",&qtb1, "qtb1/F");
	treeWithBDT->Branch("nSV",&nSV, "nSV/I");
	treeWithBDT->Branch("Trigweight",&Trigweight, "Trigweight/F");
	treeWithBDT->Branch("B2011PUweight",&B2011PUweight, "B2011PUweight/F");
	treeWithBDT->Branch("A2011PUweight",&A2011PUweight, "A2011PUweight/F");
	treeWithBDT->Branch("btag2CSF",&btag2CSF, "btag2CSF/F");
	treeWithBDT->Branch("DetaJJ",&DetaJJ, "DetaJJ/F");
	treeWithBDT->Branch("jetCHF0",&jetCHF0, "jetCHF0/F");
	treeWithBDT->Branch("jetCHF1",&jetCHF1, "jetCHF1/F");
	treeWithBDT->Branch("jetPhi0",&jetPhi0, "jetPhi0/F");
	treeWithBDT->Branch("jetPhi1",&jetPhi1, "jetPhi1/F");
	treeWithBDT->Branch("jetEta0",&jetEta0, "jetEta0/F");
	treeWithBDT->Branch("jetEta1",&jetEta1, "jetEta1/F");
	treeWithBDT->Branch("lep1pt",&lep1pt, "lep1pt/F");
	treeWithBDT->Branch("lep_pfCombRelIso0",&lep_pfCombRelIso0, "lep_pfCombRelIso0/F");
	treeWithBDT->Branch("lep_pfCombRelIso1",&lep_pfCombRelIso1, "lep_pfCombRelIso1/F");
	treeWithBDT->Branch("DphiJJ",&DphiJJ, "DphiJJ/F");
	treeWithBDT->Branch("RMS_eta",&RMS_eta, "RMS_eta/F");
	treeWithBDT->Branch("PtbalZH",&PtbalZH, "PtbalZH/F");
	treeWithBDT->Branch("EventPt",&EventPt, "EventPt/F");
	treeWithBDT->Branch("AngleEMU",&AngleEMU, "AngleEMU/F");
	treeWithBDT->Branch("Centrality",&Centrality, "Centrality/F");
	treeWithBDT->Branch("MET",&MET, "MET/F");
	treeWithBDT->Branch("EvntShpAplanarity",&EvntShpAplanarity, "EvntShpAplanarity/F");
	treeWithBDT->Branch("EvntShpSphericity",&EvntShpSphericity, "EvntShpSphericity/F");
	treeWithBDT->Branch("EvntShpIsotropy",&EvntShpIsotropy, "EvntShpIsotropy/F");
	treeWithBDT->Branch("Zphi",&Zphi, "Zphi/F");
	treeWithBDT->Branch("Hphi",&Hphi, "Hphi/F");
	treeWithBDT->Branch("Mte",&Mte, "Mte/F");
	treeWithBDT->Branch("Mtmu",&Mtmu, "Mtmu/F");
	treeWithBDT->Branch("dPhiHMET",&dPhiHMET, "dPhiHMET/F");
	treeWithBDT->Branch("Dphiemu",&Dphiemu, "Dphiemu/F");
	treeWithBDT->Branch("delRjj",&delRjj, "delRjj/F");
	treeWithBDT->Branch("delRemu",&delRemu, "delRemu/F");
	treeWithBDT->Branch("DphiZMET",&DphiZMET, "DphiZMET/F");
	treeWithBDT->Branch("DeltaPhijetMETmin",&DeltaPhijetMETmin, "DeltaPhijetMETmin/F");
	treeWithBDT->Branch("BDTvalue",&BDTvalue, "BDTvalue/F");
	treeWithBDT->Branch("AngleHemu",&AngleHemu, "AngleHemu/F");
	treeWithBDT->Branch("ProjVisT",&ProjVisT, "ProjVisT/F");
	treeWithBDT->Branch("topMass",&topMass, "topMass/F");
	treeWithBDT->Branch("topPt",&topPt, "topPt/F");
	treeWithBDT->Branch("ZmassSVD",&ZmassSVD, "ZmassSVD/F");
	treeWithBDT->Branch("ZmassSVDnegSol",&ZmassSVDnegSol, "ZmassSVDnegSol/F");
	treeWithBDT->Branch("Zmass",&Zmass, "Zmass/F");
	treeWithBDT->Branch("ZmassNegInclu",&ZmassNegInclu, "ZmassNegInclu/F");
	treeWithBDT->Branch("topMass",&topMass, "topMass/F");
	treeWithBDT->Branch("Mt",&Mt, "Mt/F");

	
   UInt_t nbin = 15;
   TH1F   *histLk(0), *histLkD(0), *histLkPCA(0), *histLkKDE(0), *histLkMIX(0), *histPD(0), *histPDD(0);
   TH1F   *histPDPCA(0), *histPDEFoam(0), *histPDEFoamErr(0), *histPDEFoamSig(0), *histKNN(0), *histHm(0);
   TH1F   *histFi(0), *histFiG(0), *histFiB(0), *histLD(0), *histNn(0),*histNnbfgs(0),*histNnbnn(0);
   TH1F   *histNnC(0), *histNnT(0), *histBdt(0), *histBdtG(0), *histBdtD(0), *histRf(0), *histSVMG(0);
   TH1F   *histSVMP(0), *histSVML(0), *histFDAMT(0), *histFDAGA(0), *histCat(0), *histPBdt(0);

   if (Use["Likelihood"])    histLk      = new TH1F( "MVA_Likelihood",    "MVA_Likelihood",    nbin, -1, 1 );
   if (Use["LikelihoodD"])   histLkD     = new TH1F( "MVA_LikelihoodD",   "MVA_LikelihoodD",   nbin, -1, 0.9999 );
   if (Use["LikelihoodPCA"]) histLkPCA   = new TH1F( "MVA_LikelihoodPCA", "MVA_LikelihoodPCA", nbin, -1, 1 );
   if (Use["LikelihoodKDE"]) histLkKDE   = new TH1F( "MVA_LikelihoodKDE", "MVA_LikelihoodKDE", nbin,  -0.00001, 0.99999 );
   if (Use["LikelihoodMIX"]) histLkMIX   = new TH1F( "MVA_LikelihoodMIX", "MVA_LikelihoodMIX", nbin,  0, 1 );
   if (Use["PDERS"])         histPD      = new TH1F( "MVA_PDERS",         "MVA_PDERS",         nbin,  0, 1 );
   if (Use["PDERSD"])        histPDD     = new TH1F( "MVA_PDERSD",        "MVA_PDERSD",        nbin,  0, 1 );
   if (Use["PDERSPCA"])      histPDPCA   = new TH1F( "MVA_PDERSPCA",      "MVA_PDERSPCA",      nbin,  0, 1 );
   if (Use["KNN"])           histKNN     = new TH1F( "MVA_KNN",           "MVA_KNN",           nbin,  0, 1 );
   if (Use["HMatrix"])       histHm      = new TH1F( "MVA_HMatrix",       "MVA_HMatrix",       nbin, -0.95, 1.55 );
   if (Use["Fisher"])        histFi      = new TH1F( "MVA_Fisher",        "MVA_Fisher",        nbin, -4, 4 );
   if (Use["FisherG"])       histFiG     = new TH1F( "MVA_FisherG",       "MVA_FisherG",       nbin, -1, 1 );
   if (Use["BoostedFisher"]) histFiB     = new TH1F( "MVA_BoostedFisher", "MVA_BoostedFisher", nbin, -2, 2 );
   if (Use["LD"])            histLD      = new TH1F( "MVA_LD",            "MVA_LD",            nbin, -2, 2 );
   if (Use["MLP"])           histNn      = new TH1F( "MVA_MLP",           "MVA_MLP",           nbin, -1.25, 1.5 );
   if (Use["MLPBFGS"])       histNnbfgs  = new TH1F( "MVA_MLPBFGS",       "MVA_MLPBFGS",       nbin, -1.25, 1.5 );
   if (Use["MLPBNN"])        histNnbnn   = new TH1F( "MVA_MLPBNN",        "MVA_MLPBNN",        nbin, -1.25, 1.5 );
   if (Use["CFMlpANN"])      histNnC     = new TH1F( "MVA_CFMlpANN",      "MVA_CFMlpANN",      nbin,  0, 1 );
   if (Use["TMlpANN"])       histNnT     = new TH1F( "MVA_TMlpANN",       "MVA_TMlpANN",       nbin, -1.3, 1.3 );
	if (Use["BDT"])     {
		histMattBdt     = new TH1F( "Matt_BDT",           "Matt_BDT",           15, -1.1, 0.35 );
		histTMVABdt     = new TH1F( "TMVA_BDT",           "TMVA_BDT",           36, -1.0, -0.1 );
	}
   if (Use["BDTD"])          histBdtD    = new TH1F( "MVA_BDTD",          "MVA_BDTD",          nbin, -0.8, 0.8 );
   if (Use["BDTG"])          histBdtG    = new TH1F( "MVA_BDTG",          "MVA_BDTG",          nbin, -1.0, 1.0 );
   if (Use["RuleFit"])       histRf      = new TH1F( "MVA_RuleFit",       "MVA_RuleFit",       nbin, -2.0, 2.0 );
   if (Use["SVM_Gauss"])     histSVMG    = new TH1F( "MVA_SVM_Gauss",     "MVA_SVM_Gauss",     nbin,  0.0, 1.0 );
   if (Use["SVM_Poly"])      histSVMP    = new TH1F( "MVA_SVM_Poly",      "MVA_SVM_Poly",      nbin,  0.0, 1.0 );
   if (Use["SVM_Lin"])       histSVML    = new TH1F( "MVA_SVM_Lin",       "MVA_SVM_Lin",       nbin,  0.0, 1.0 );
   if (Use["FDA_MT"])        histFDAMT   = new TH1F( "MVA_FDA_MT",        "MVA_FDA_MT",        nbin, -2.0, 3.0 );
   if (Use["FDA_GA"])        histFDAGA   = new TH1F( "MVA_FDA_GA",        "MVA_FDA_GA",        nbin, -2.0, 3.0 );
   if (Use["Category"])      histCat     = new TH1F( "MVA_Category",      "MVA_Category",      nbin, -2., 2. );
   if (Use["Plugin"])        histPBdt    = new TH1F( "MVA_PBDT",          "MVA_BDT",           nbin, -0.8, 0.8 );

   // PDEFoam also returns per-event error, fill in histogram, and also fill significance
   if (Use["PDEFoam"]) {
      histPDEFoam    = new TH1F( "MVA_PDEFoam",       "MVA_PDEFoam",              nbin,  0, 1 );
      histPDEFoamErr = new TH1F( "MVA_PDEFoamErr",    "MVA_PDEFoam error",        nbin,  0, 1 );
      histPDEFoamSig = new TH1F( "MVA_PDEFoamSig",    "MVA_PDEFoam significance", nbin,  0, 10 );
   }

   // Book example histogram for probability (the other methods are done similarly)
   TH1F *probHistFi(0), *rarityHistFi(0);
   if (Use["Fisher"]) {
      probHistFi   = new TH1F( "MVA_Fisher_Proba",  "MVA_Fisher_Proba",  nbin, 0, 1 );
      rarityHistFi = new TH1F( "MVA_Fisher_Rarity", "MVA_Fisher_Rarity", nbin, 0, 1 );
   }

   // Prepare input tree (this must be replaced by your data source)
   // in this example, there is a toy tree with signal and one with background events
   // we'll later on use only the "signal" events for the test in this example.
   //   
   TFile *input(0);
   TString fname = "/home/hep/wilken/taus/CMSSW_4_4_2_patch8/src/UserCode/wilken/V21/T_tW.root"; 
	double lumi = 4.457;
	Double_t  T_tW_weight = lumi/(lumiTtW/2.0); //T_tW_TuneZ2_7TeV_pythia6_tauola
  
   if (!gSystem->AccessPathName( fname )) 
      input = TFile::Open( fname ); // check if file in local directory exists
   else    
      input = TFile::Open( "http://root.cern.ch/files/tmva_class_example.root" ); // if not: download from ROOT server
   
   if (!input) {
      std::cout << "ERROR: could not open data file" << std::endl;
      exit(1);
   }
   std::cout << "--- TMVAClassificationApp    : Using input file: " << input->GetName() << std::endl;
   
   // --- Event loop

   // Prepare the event tree
   // - here the variable names have to corresponds to your tree
   // - you can use the same variables as above which is slightly faster,
   //   but of course you can use different ones and copy the values inside the event loop
   //
   std::cout << "--- Select signal sample" << std::endl;
   TTree* theTree = (TTree*)input->Get("BDT_tree");
   
    theTree->SetBranchAddress( "Hmass", &Hmass );
    theTree->SetBranchAddress( "Emumass", &Emumass );
    theTree->SetBranchAddress( "Hpt", &Hpt );
	theTree->SetBranchAddress( "Zpt", &Zpt );
	theTree->SetBranchAddress( "CSV0", &CSV0 );
	theTree->SetBranchAddress( "CSV1", &CSV1 );
	theTree->SetBranchAddress( "DeltaPhiHV", &DeltaPhiHV );
	theTree->SetBranchAddress( "DetaJJ", &DetaJJ );
	theTree->SetBranchAddress( "UnweightedEta", &UnweightedEta );
	theTree->SetBranchAddress( "lep0pt", &lep0pt );
	theTree->SetBranchAddress( "Ht", &Ht );
	theTree->SetBranchAddress( "EvntShpCircularity", &EvntShpCircularity );
	theTree->SetBranchAddress( "nJets", &nJets );
	theTree->SetBranchAddress( "naJets",&naJets);
	theTree->SetBranchAddress( "nSV", &nSV );
	theTree->SetBranchAddress( "lep1pt", &lep1pt );
	theTree->SetBranchAddress( "lep0pt", &lep0pt );
	theTree->SetBranchAddress( "EtaStandDev", &EtaStandDev );
	theTree->SetBranchAddress( "UnweightedEta", &UnweightedEta );
	theTree->SetBranchAddress( "jetCHF0", &jetCHF0 );
	theTree->SetBranchAddress( "jetCHF1", &jetCHF1 );
	theTree->SetBranchAddress( "lep_pfCombRelIso0", &lep_pfCombRelIso0 );
	theTree->SetBranchAddress( "lep_pfCombRelIso1", &lep_pfCombRelIso1 );
	theTree->SetBranchAddress( "DphiJJ", &DphiJJ );
	theTree->SetBranchAddress( "RMS_eta", &RMS_eta );
	theTree->SetBranchAddress( "PtbalZH", &PtbalZH );
	theTree->SetBranchAddress( "EventPt", &EventPt );
	theTree->SetBranchAddress( "AngleEMU", &AngleEMU );
	theTree->SetBranchAddress( "Centrality", &Centrality );
	theTree->SetBranchAddress( "EvntShpAplanarity", &EvntShpAplanarity );
	theTree->SetBranchAddress( "EvntShpSphericity", &EvntShpSphericity );
	theTree->SetBranchAddress( "EvntShpIsotropy", &EvntShpIsotropy );
	theTree->SetBranchAddress( "Trigweight", &Trigweight );
	theTree->SetBranchAddress( "B2011PUweight", &B2011PUweight );
	theTree->SetBranchAddress( "A2011PUweight", &A2011PUweight );
	theTree->SetBranchAddress( "btag2CSF", &btag2CSF );
	theTree->SetBranchAddress( "MET", &MET );
	theTree->SetBranchAddress( "Mte"      ,  &Mte     );
	theTree->SetBranchAddress( "Mtmu"      ,  &Mtmu     );
	theTree->SetBranchAddress( "dPhiHMET"      ,  &dPhiHMET     );
	theTree->SetBranchAddress( "DeltaPhijetMETmin"      ,  &DeltaPhijetMETmin     );
	theTree->SetBranchAddress( "delRjj"      ,  &delRjj     );
	theTree->SetBranchAddress( "delRemu"      ,  &delRemu     );
	theTree->SetBranchAddress( "DphiZMET"      ,  &DphiZMET     );
	theTree->SetBranchAddress( "Dphiemu"      ,  &Dphiemu     );
	theTree->SetBranchAddress( "DeltaPhijetMETmin"      ,  &DeltaPhijetMETmin     );
	theTree->SetBranchAddress("AngleHemu",&AngleHemu);
	theTree->SetBranchAddress("ProjVisT",&ProjVisT);
	theTree->SetBranchAddress("topMass",&topMass);
	theTree->SetBranchAddress("topPt",&topPt);
	theTree->SetBranchAddress("ZmassSVD",&ZmassSVD);
	theTree->SetBranchAddress("ZmassSVDnegSol",&ZmassSVDnegSol);
	theTree->SetBranchAddress("Zmass",&Zmass);
	theTree->SetBranchAddress("ZmassNegInclu",&ZmassNegInclu);
	theTree->SetBranchAddress("topMass",&topMass);
	theTree->SetBranchAddress("Mt",&Mt);
	
	theTree->SetBranchAddress( "eventFlavor", &eventFlavor );
	

   // Efficiency calculator for cut method
   Int_t    nSelCutsGA = 0;
   Double_t effS       = 0.7;

   std::vector<Float_t> vecVar(4); // vector for EvaluateMVA tests

   std::cout << "--- Processing: " << theTree->GetEntries() << " events" << std::endl;
   TStopwatch sw;
   sw.Start();
   int Nevents = 0, NpassBDT = 0;
   for (Long64_t ievt=0; ievt<theTree->GetEntries();ievt++) {
Nevents++;
      if (ievt%1000 == 0) std::cout << "--- ... Processing event: " << ievt << std::endl;

      theTree->GetEntry(ievt);

 //     var1 = userVar1 + userVar2;
   //   var2 = userVar1 - userVar2;

      // --- Return the MVA outputs and fill into histograms

      if (Use["CutsGA"]) {
         // Cuts is a special case: give the desired signal efficienciy
         Bool_t passed = reader->EvaluateMVA( "CutsGA method", effS );
         if (passed) nSelCutsGA++;
      }

      if (Use["Likelihood"   ])   histLk     ->Fill( reader->EvaluateMVA( "Likelihood method"    ) );
      if (Use["LikelihoodD"  ])   histLkD    ->Fill( reader->EvaluateMVA( "LikelihoodD method"   ) );
      if (Use["LikelihoodPCA"])   histLkPCA  ->Fill( reader->EvaluateMVA( "LikelihoodPCA method" ) );
      if (Use["LikelihoodKDE"])   histLkKDE  ->Fill( reader->EvaluateMVA( "LikelihoodKDE method" ) );
      if (Use["LikelihoodMIX"])   histLkMIX  ->Fill( reader->EvaluateMVA( "LikelihoodMIX method" ) );
      if (Use["PDERS"        ])   histPD     ->Fill( reader->EvaluateMVA( "PDERS method"         ) );
      if (Use["PDERSD"       ])   histPDD    ->Fill( reader->EvaluateMVA( "PDERSD method"        ) );
      if (Use["PDERSPCA"     ])   histPDPCA  ->Fill( reader->EvaluateMVA( "PDERSPCA method"      ) );
      if (Use["KNN"          ])   histKNN    ->Fill( reader->EvaluateMVA( "KNN method"           ) );
      if (Use["HMatrix"      ])   histHm     ->Fill( reader->EvaluateMVA( "HMatrix method"       ) );
      if (Use["Fisher"       ])   histFi     ->Fill( reader->EvaluateMVA( "Fisher method"        ) );
      if (Use["FisherG"      ])   histFiG    ->Fill( reader->EvaluateMVA( "FisherG method"       ) );
      if (Use["BoostedFisher"])   histFiB    ->Fill( reader->EvaluateMVA( "BoostedFisher method" ) );
      if (Use["LD"           ])   histLD     ->Fill( reader->EvaluateMVA( "LD method"            ) );
      if (Use["MLP"          ])   histNn     ->Fill( reader->EvaluateMVA( "MLP method"           ) );
      if (Use["MLPBFGS"      ])   histNnbfgs ->Fill( reader->EvaluateMVA( "MLPBFGS method"       ) );
      if (Use["MLPBNN"       ])   histNnbnn  ->Fill( reader->EvaluateMVA( "MLPBNN method"        ) );
      if (Use["CFMlpANN"     ])   histNnC    ->Fill( reader->EvaluateMVA( "CFMlpANN method"      ) );
      if (Use["TMlpANN"      ])   histNnT    ->Fill( reader->EvaluateMVA( "TMlpANN method"       ) );
	   if (Use["BDT"          ]) {
	   BDTvalue = reader->EvaluateMVA( "BDT method"           );
		   histMattBdt    ->Fill( BDTvalue,T_tW_weight*Trigweight*B2011PUweight );
		   histTMVABdt    ->Fill( BDTvalue,T_tW_weight*Trigweight*B2011PUweight );
	   }
      if (Use["BDTD"         ])   histBdtD   ->Fill( reader->EvaluateMVA( "BDTD method"          ) );
      if (Use["BDTG"         ])   histBdtG   ->Fill( reader->EvaluateMVA( "BDTG method"          ) );
      if (Use["RuleFit"      ])   histRf     ->Fill( reader->EvaluateMVA( "RuleFit method"       ) );
      if (Use["SVM_Gauss"    ])   histSVMG   ->Fill( reader->EvaluateMVA( "SVM_Gauss method"     ) );
      if (Use["SVM_Poly"     ])   histSVMP   ->Fill( reader->EvaluateMVA( "SVM_Poly method"      ) );
      if (Use["SVM_Lin"      ])   histSVML   ->Fill( reader->EvaluateMVA( "SVM_Lin method"       ) );
      if (Use["FDA_MT"       ])   histFDAMT  ->Fill( reader->EvaluateMVA( "FDA_MT method"        ) );
      if (Use["FDA_GA"       ])   histFDAGA  ->Fill( reader->EvaluateMVA( "FDA_GA method"        ) );
      if (Use["Category"     ])   histCat    ->Fill( reader->EvaluateMVA( "Category method"      ) );
      if (Use["Plugin"       ])   histPBdt   ->Fill( reader->EvaluateMVA( "P_BDT method"         ) );

      // Retrieve also per-event error
      if (Use["PDEFoam"]) {
         Double_t val = reader->EvaluateMVA( "PDEFoam method" );
         Double_t err = reader->GetMVAError();
         histPDEFoam   ->Fill( val );
         histPDEFoamErr->Fill( err );         
         if (err>1.e-50) histPDEFoamSig->Fill( val/err );
      }         

      // Retrieve probability instead of MVA output
      if (Use["Fisher"])   {
         probHistFi  ->Fill( reader->GetProba ( "Fisher method" ) );
         rarityHistFi->Fill( reader->GetRarity( "Fisher method" ) );
      }
   
	  // std::cout << "Ht is "<< Ht << endl;
	   if(BDTvalue>-1.50){
NpassBDT++;
hMjj_OpenSelection->Fill(Hmass,T_tW_weight*Trigweight*B2011PUweight );
hMmumu_OpenSelection->Fill(Emumass,T_tW_weight*Trigweight*B2011PUweight );
hPtjj_OpenSelection->Fill(Hpt,T_tW_weight*Trigweight*B2011PUweight );
hPtmumu_OpenSelection->Fill(Zpt,T_tW_weight*Trigweight*B2011PUweight );
hCSV0_OpenSelection->Fill(CSV0,T_tW_weight*Trigweight*B2011PUweight );
hCSV1_OpenSelection->Fill(CSV1,T_tW_weight*Trigweight*B2011PUweight );
hdphiVH_OpenSelection->Fill(DeltaPhiHV,T_tW_weight*Trigweight*B2011PUweight );
hdetaJJ_OpenSelection->Fill(DetaJJ,T_tW_weight*Trigweight*B2011PUweight );
hUnweightedEta_OpenSelection->Fill(UnweightedEta,T_tW_weight*Trigweight*B2011PUweight );
hPtmu0_OpenSelection->Fill(lep0pt,T_tW_weight*Trigweight*B2011PUweight );
	   hHt_OpenSelection->Fill(Ht,T_tW_weight*Trigweight*B2011PUweight );
	   hCircularity_OpenSelection->Fill(EvntShpCircularity,T_tW_weight*Trigweight*B2011PUweight );
	   hCHFb0_OpenSelection->Fill(jetCHF0, T_tW_weight*Trigweight*B2011PUweight );
	   hCHFb1_OpenSelection->Fill(jetCHF1, T_tW_weight*Trigweight*B2011PUweight );
	   hPtbalZH_OpenSelection->Fill(PtbalZH, T_tW_weight*Trigweight*B2011PUweight );
	   hPtmu1_OpenSelection->Fill(lep1pt, T_tW_weight*Trigweight*B2011PUweight );
	   hPFRelIsomu0_OpenSelection->Fill(lep_pfCombRelIso0, T_tW_weight*Trigweight*B2011PUweight );
	   hPFRelIsomu1_OpenSelection->Fill(lep_pfCombRelIso1, T_tW_weight*Trigweight*B2011PUweight );
	   hNjets_OpenSelection->Fill(nJets, T_tW_weight*Trigweight*B2011PUweight );
	   hRMSeta_OpenSelection->Fill(RMS_eta, T_tW_weight*Trigweight*B2011PUweight );
	   hStaDeveta_OpenSelection->Fill(EtaStandDev, T_tW_weight*Trigweight*B2011PUweight );
	   hdphiJJ_vect_OpenSelection->Fill(DphiJJ, T_tW_weight*Trigweight*B2011PUweight );
	   hCentrality_OpenSelection->Fill(Centrality, T_tW_weight*Trigweight*B2011PUweight );
	   hEventPt_OpenSelection->Fill(EventPt, T_tW_weight*Trigweight*B2011PUweight );
	   hAngleEMU_OpenSelection->Fill(AngleEMU, T_tW_weight*Trigweight*B2011PUweight );
	   hSphericity_OpenSelection->Fill(EvntShpSphericity, T_tW_weight*Trigweight*B2011PUweight );
	   hAplanarity_OpenSelection->Fill(EvntShpAplanarity, T_tW_weight*Trigweight*B2011PUweight );
	   hIsotropy_OpenSelection->Fill(EvntShpIsotropy, T_tW_weight*Trigweight*B2011PUweight );
		   hDphiDetajj_OpenSelection->Fill(DphiJJ, DetaJJ, T_tW_weight*Trigweight*B2011PUweight );
		   hMte_OpenSelection->Fill(Mte, T_tW_weight*Trigweight*B2011PUweight );
		   hMtmu_OpenSelection->Fill(Mtmu, T_tW_weight*Trigweight*B2011PUweight );
		   hdPhiHMET_OpenSelection->Fill(dPhiHMET, T_tW_weight*Trigweight*B2011PUweight );
		   hDphiemu_OpenSelection->Fill(DeltaPhijetMETmin, T_tW_weight*Trigweight*B2011PUweight );
		   hdelRjj_OpenSelection->Fill(delRjj, T_tW_weight*Trigweight*B2011PUweight );
		   hdelRemu_OpenSelection->Fill(delRemu, T_tW_weight*Trigweight*B2011PUweight );
		   hDphiZMET_OpenSelection->Fill(DphiZMET, T_tW_weight*Trigweight*B2011PUweight );
		   hDphiemu_OpenSelection->Fill(Dphiemu, T_tW_weight*Trigweight*B2011PUweight );
		   hDeltaPhijetMETmin_OpenSelection->Fill(DeltaPhijetMETmin, T_tW_weight*Trigweight*B2011PUweight );
		   hAngleHemu_OpenSelection->Fill(AngleHemu, T_tW_weight*Trigweight*B2011PUweight );
		   hProjVisT_OpenSelection->Fill(ProjVisT, T_tW_weight*Trigweight*B2011PUweight );
		   htopMass_OpenSelection->Fill(topMass, T_tW_weight*Trigweight*B2011PUweight );
		   htopPt_OpenSelection->Fill(topPt, T_tW_weight*Trigweight*B2011PUweight );
		   hVMt_OpenSelection->Fill(Mt, T_tW_weight*Trigweight*B2011PUweight );
		   hZmassSVD_OpenSelection->Fill(ZmassSVD, T_tW_weight*Trigweight*B2011PUweight );
		   hZmassSVDnegSol_OpenSelection->Fill(ZmassSVDnegSol, T_tW_weight*Trigweight*B2011PUweight );
		   hZmass_OpenSelection->Fill(Zmass, T_tW_weight*Trigweight*B2011PUweight );
		   hZmassNegInclu_OpenSelection->Fill(ZmassNegInclu, T_tW_weight*Trigweight*B2011PUweight );

	   }
	   treeWithBDT->Fill();
	   
   }//end event loop
   
   // Get elapsed time
   sw.Stop();
   std::cout << "--- End of event loop: "; sw.Print();
std::cout << "Number of Events: "<< Nevents << " Events passed BDT " << NpassBDT<< endl;
   // Get efficiency for cuts classifier
   
   if (Use["CutsGA"]) std::cout << "--- Efficiency for CutsGA method: " << double(nSelCutsGA)/theTree->GetEntries()
                                << " (for a required signal efficiency of " << effS << ")" << std::endl;

   if (Use["CutsGA"]) {

      // test: retrieve cuts for particular signal efficiency
      // CINT ignores dynamic_casts so we have to use a cuts-secific Reader function to acces the pointer  
      TMVA::MethodCuts* mcuts = reader->FindCutsMVA( "CutsGA method" ) ;

      if (mcuts) {      
         std::vector<Double_t> cutsMin;
         std::vector<Double_t> cutsMax;
         mcuts->GetCuts( 0.7, cutsMin, cutsMax );
         std::cout << "--- -------------------------------------------------------------" << std::endl;
         std::cout << "--- Retrieve cut values for signal efficiency of 0.7 from Reader" << std::endl;
         for (UInt_t ivar=0; ivar<cutsMin.size(); ivar++) {
            std::cout << "... Cut: " 
                      << cutsMin[ivar] 
                      << " < \"" 
                      << mcuts->GetInputVar(ivar)
                      << "\" <= " 
                      << cutsMax[ivar] << std::endl;
         }
         std::cout << "--- -------------------------------------------------------------" << std::endl;
      }
   }
	

   // --- Write histograms
   TFile *target  = new TFile( "BDTCut_T_tW.root","RECREATE" );
   if (Use["Likelihood"   ])   histLk     ->Write();
   if (Use["LikelihoodD"  ])   histLkD    ->Write();
   if (Use["LikelihoodPCA"])   histLkPCA  ->Write();
   if (Use["LikelihoodKDE"])   histLkKDE  ->Write();
   if (Use["LikelihoodMIX"])   histLkMIX  ->Write();
   if (Use["PDERS"        ])   histPD     ->Write();
   if (Use["PDERSD"       ])   histPDD    ->Write();
   if (Use["PDERSPCA"     ])   histPDPCA  ->Write();
   if (Use["KNN"          ])   histKNN    ->Write();
   if (Use["HMatrix"      ])   histHm     ->Write();
   if (Use["Fisher"       ])   histFi     ->Write();
   if (Use["FisherG"      ])   histFiG    ->Write();
   if (Use["BoostedFisher"])   histFiB    ->Write();
   if (Use["LD"           ])   histLD     ->Write();
   if (Use["MLP"          ])   histNn     ->Write();
   if (Use["MLPBFGS"      ])   histNnbfgs ->Write();
   if (Use["MLPBNN"       ])   histNnbnn  ->Write();
   if (Use["CFMlpANN"     ])   histNnC    ->Write();
   if (Use["TMlpANN"      ])   histNnT    ->Write();
	if (Use["BDT"          ])  {
		histMattBdt    ->Write();
		histTMVABdt    ->Write();
	}
   if (Use["BDTD"         ])   histBdtD   ->Write();
   if (Use["BDTG"         ])   histBdtG   ->Write(); 
   if (Use["RuleFit"      ])   histRf     ->Write();
   if (Use["SVM_Gauss"    ])   histSVMG   ->Write();
   if (Use["SVM_Poly"     ])   histSVMP   ->Write();
   if (Use["SVM_Lin"      ])   histSVML   ->Write();
   if (Use["FDA_MT"       ])   histFDAMT  ->Write();
   if (Use["FDA_GA"       ])   histFDAGA  ->Write();
   if (Use["Category"     ])   histCat    ->Write();
   if (Use["Plugin"       ])   histPBdt   ->Write();

   // Write also error and significance histos
   if (Use["PDEFoam"]) { histPDEFoam->Write(); histPDEFoamErr->Write(); histPDEFoamSig->Write(); }

   // Write also probability hists
   if (Use["Fisher"]) { if (probHistFi != 0) probHistFi->Write(); if (rarityHistFi != 0) rarityHistFi->Write(); }
   
   hCHFb0_OpenSelection->Write();
   hCHFb1_OpenSelection->Write();
   hPtbalZH_OpenSelection->Write();
   hPtmu1_OpenSelection->Write();
   hPFRelIsomu0_OpenSelection->Write();
   hPFRelIsomu1_OpenSelection->Write();
   hMjj_OpenSelection->Write();
   hMmumu_OpenSelection->Write();
   hPtjj_OpenSelection->Write();
   hPtmumu_OpenSelection->Write();
   hCSV0_OpenSelection->Write();
   hCSV1_OpenSelection->Write();
   hdphiVH_OpenSelection->Write();
   hdetaJJ_OpenSelection->Write();
   hUnweightedEta_OpenSelection->Write();
   hPtmu0_OpenSelection->Write();
	hCircularity_OpenSelection->Write();
	hHt_OpenSelection->Write();
	hNjets_OpenSelection->Write();
	hRMSeta_OpenSelection->Write();
hStaDeveta_OpenSelection->Write();
hdphiJJ_vect_OpenSelection->Write();
hCentrality_OpenSelection->Write();
hEventPt_OpenSelection->Write();
hAngleEMU_OpenSelection->Write();
hSphericity_OpenSelection->Write();
hAplanarity_OpenSelection->Write();
hIsotropy_OpenSelection->Write();
hDphiDetajj_OpenSelection->Write();
hMtmu_OpenSelection->Write();
hMte_OpenSelection->Write();
hdPhiHMET_OpenSelection->Write();
hDphiemu_OpenSelection->Write();
hdelRjj_OpenSelection->Write();
hdelRemu_OpenSelection->Write();
hDphiZMET_OpenSelection->Write();
hDeltaPhijetMETmin_OpenSelection->Write();
hAngleHemu_OpenSelection->Write();
hProjVisT_OpenSelection->Write();
htopMass_OpenSelection->Write();
htopPt_OpenSelection->Write();
hVMt_OpenSelection->Write();
hZmassSVD_OpenSelection->Write();
hZmassSVDnegSol_OpenSelection->Write();
hZmass_OpenSelection->Write();
hZmassNegInclu_OpenSelection->Write();

treeWithBDT->Write();

   
   target->Close();

  
   delete reader;

	hCHFb0_OpenSelection->Delete();
	hCHFb1_OpenSelection->Delete();
	hPtbalZH_OpenSelection->Delete();
	hPtmu1_OpenSelection->Delete();
	hPFRelIsomu0_OpenSelection->Delete();
	hPFRelIsomu1_OpenSelection->Delete();
	hMjj_OpenSelection->Delete();
	hMmumu_OpenSelection->Delete();
	hPtjj_OpenSelection->Delete();
	hPtmumu_OpenSelection->Delete();
	hCSV0_OpenSelection->Delete();
	hCSV1_OpenSelection->Delete();
	hdphiVH_OpenSelection->Delete();
	hdetaJJ_OpenSelection->Delete();
	hUnweightedEta_OpenSelection->Delete();
	hPtmu0_OpenSelection->Delete();
	hCircularity_OpenSelection->Delete();
	hHt_OpenSelection->Delete();
	hNjets_OpenSelection->Delete();
	hRMSeta_OpenSelection->Delete();
	hStaDeveta_OpenSelection->Delete();
	hdphiJJ_vect_OpenSelection->Delete();
	hCentrality_OpenSelection->Delete();
	hEventPt_OpenSelection->Delete();
	hAngleEMU_OpenSelection->Delete();
	hSphericity_OpenSelection->Delete();
	hAplanarity_OpenSelection->Delete();
	hIsotropy_OpenSelection->Delete();
	hDphiDetajj_OpenSelection->Delete();
	hMtmu_OpenSelection->Delete();
	hMte_OpenSelection->Delete();
	hdPhiHMET_OpenSelection->Delete();
	hDphiemu_OpenSelection->Delete();
	hdelRjj_OpenSelection->Delete();
	hdelRemu_OpenSelection->Delete();
	hDphiZMET_OpenSelection->Delete();
	hDeltaPhijetMETmin_OpenSelection->Delete();
	hAngleHemu_OpenSelection->Delete();
	hProjVisT_OpenSelection->Delete();
	htopMass_OpenSelection->Delete();
	htopPt_OpenSelection->Delete();
	hVMt_OpenSelection->Delete();
	hZmassSVD_OpenSelection->Delete();
	hZmassSVDnegSol_OpenSelection->Delete();
	hZmass_OpenSelection->Delete();
	hZmassNegInclu_OpenSelection->Delete();
	

	if (Use["BDT"          ])  {
		histMattBdt    ->Delete();
		histTMVABdt    ->Delete();
	}
	
	

treeWithBDT->Delete();

	
    
   std::cout << "==> T_tW_BDTCut is done!" << endl << std::endl;
   
   gROOT->ProcessLine(".q");
} 
