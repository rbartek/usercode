/**********************************************************************************
 * Project   : TMVA - a Root-integrated toolkit for multivariate data analysis    *
 * Package   : TMVA                                                               *
 * Exectuable: DYJetsToLL_PtZ_App_bflavor                                     *
 *                                                                                *
 * This macro provides a simple example on how to use the trained classifiers     *
 * within an analysis module                                                      *
 **********************************************************************************/

#include <cstdlib>
#include <vector>
#include <iostream>
#include <map>
#include <string>
#include "xsecnov10f.h"
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

void DYJetsToLL_PtZ_App_bflavor( TString myMethodList = "" ) 
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
   std::cout << "==> Start DYJetsToLL_PtZ_App_bflavor" << std::endl;

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
   Float_t Hmass, Zmass;
	Float_t Hpt, Zpt;
	Float_t CSV0, CSV1;
	Float_t DeltaPhiHV, DetaJJ;
	Int_t  nJets, eventFlavor; 
	Float_t Naj, nSV;
 	Float_t BDTvalue, Trigweight, B2011PUweight, btag2CSF, MET;
	Float_t alpha_j, qtb1, jetPhi0, jetPhi1, jetEta0, jetEta1, Zphi, Hphi;
	Float_t Ht, EvntShpCircularity, jetCHF0, jetCHF1;
	Float_t EtaStandDev, muonPFiso0, muonPFiso1, EvntShpIsotropy;
	Float_t mu0pt, mu1pt, UnweightedEta, DphiJJ, RMS_eta, EvntShpSphericity;
	Float_t PtbalZH, EventPt, Angle, Centrality, EvntShpAplanarity;
	

	reader->AddVariable( "Hmass", &Hmass );
	// reader->AddVariable( "Zmass", &Zmass );
	reader->AddVariable( "Hpt",                &Hpt );
   	reader->AddVariable( "CSV0",                &CSV0 );
	reader->AddVariable( "CSV1",                &CSV1 );
	reader->AddVariable( "Zpt",                &Zpt );
	reader->AddVariable( "DeltaPhiHV:= abs(deltaPhi(Hphi,Zphi))",                &DeltaPhiHV );
	reader->AddVariable( "DetaJJ:= abs(jetEta1-jetEta0)",                &DetaJJ );
	reader->AddVariable( "Naj", &Naj);
	reader->AddVariable( "UnweightedEta",                &UnweightedEta );
	reader->AddVariable( "EvntShpCircularity",                &EvntShpCircularity );
	reader->AddVariable( "alpha_j",                &alpha_j );
	reader->AddVariable( "qtb1",                &qtb1 );
	reader->AddVariable( "nSV",                &nSV );
	reader->AddVariable( "mu0pt",                &mu0pt );
	reader->AddVariable( "mu1pt",                &mu1pt );
	reader->AddVariable( "PtbalZH",                &PtbalZH );
	reader->AddVariable( "Angle",                &Angle );
	reader->AddVariable( "Centrality",                &Centrality );
	reader->AddVariable( "MET",                &MET );
	reader->AddVariable( "EvntShpAplanarity",                &EvntShpAplanarity );

   // Spectator variables declared in the training have to be added to the reader, too
  Float_t UnweightedEta, mu0pt;
//   reader->AddSpectator( "UnweightedEta",   &UnweightedEta );
  // reader->AddSpectator( "mu0pt",   &mu0pt );

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
    TH1F* hAngle_OpenSelection= new TH1F	("hAngle_OpenSelection", "Angle between H and Z", 45, 0, 4.5);
    TH1F* hSphericity_OpenSelection= new TH1F	("hSphericity_OpenSelection", "EventShapeVariables sphericity", 50, 0.0, 1);
    TH1F* hAplanarity_OpenSelection= new TH1F	("hAplanarity_OpenSelection", "EventShapeVariables Aplanarity", 50, -0.1, .4);
    TH1F* hIsotropy_OpenSelection= new TH1F	("hIsotropy_OpenSelection",  "EventShapeVariables isotropy", 30, 0.0, 1.3);
    TH2F* hDphiDetajj_OpenSelection= new TH2F	("hDphiDetajj_OpenSelection", "#Delta#phi vs #Delta#eta JJ", 25, -5, 5, 25, -5, 5);
	
	TTree *treeWithBDT = new TTree("treeWithBDT","Tree wiht BDT output");
	treeWithBDT->SetDirectory(0);
	treeWithBDT->Branch("nJets",&nJets, "nJets/I");
	treeWithBDT->Branch("Naj",&Naj, "Naj/F");
	treeWithBDT->Branch("eventFlavor",&eventFlavor, "eventFlavor/I");
	treeWithBDT->Branch("CSV0",&CSV0, "CSV0/F");
	treeWithBDT->Branch("CSV1",&CSV1, "CSV1/F");
	treeWithBDT->Branch("Zmass",&Zmass, "Zmass/F");
	treeWithBDT->Branch("Hmass",&Hmass, "Hmass/F");
	treeWithBDT->Branch("DeltaPhiHV",&DeltaPhiHV, "DeltaPhiHV/F");
	treeWithBDT->Branch("Hpt",&Hpt, "Hpt/F");
	treeWithBDT->Branch("Zpt",&Zpt, "Zpt/F");
	treeWithBDT->Branch("mu0pt",&mu0pt, "mu0pt/F");
	treeWithBDT->Branch("Ht",&Ht, "Ht/F");
	treeWithBDT->Branch("EtaStandDev",&EtaStandDev, "EtaStandDev/F");
	treeWithBDT->Branch("UnweightedEta",&UnweightedEta, "UnweightedEta/F");
	treeWithBDT->Branch("EvntShpCircularity",&EvntShpCircularity, "EvntShpCircularity/F");
	treeWithBDT->Branch("alpha_j",&alpha_j, "alpha_j/F");
	treeWithBDT->Branch("qtb1",&qtb1, "qtb1/F");
	treeWithBDT->Branch("nSV",&nSV, "nSV/F");
	treeWithBDT->Branch("Trigweight",&Trigweight, "Trigweight/F");
	treeWithBDT->Branch("B2011PUweight",&B2011PUweight, "B2011PUweight/F");
	treeWithBDT->Branch("btag2CSF",&btag2CSF, "btag2CSF/F");
	treeWithBDT->Branch("DetaJJ",&DetaJJ, "DetaJJ/F");
	treeWithBDT->Branch("jetCHF0",&jetCHF0, "jetCHF0/F");
	treeWithBDT->Branch("jetCHF1",&jetCHF1, "jetCHF1/F");
	treeWithBDT->Branch("jetPhi0",&jetPhi0, "jetPhi0/F");
	treeWithBDT->Branch("jetPhi1",&jetPhi1, "jetPhi1/F");
	treeWithBDT->Branch("jetEta0",&jetEta0, "jetEta0/F");
	treeWithBDT->Branch("jetEta1",&jetEta1, "jetEta1/F");
	treeWithBDT->Branch("mu1pt",&mu1pt, "mu1pt/F");
	treeWithBDT->Branch("muonPFiso0",&muonPFiso0, "muonPFiso0/F");
	treeWithBDT->Branch("muonPFiso1",&muonPFiso1, "muonPFiso1/F");
	treeWithBDT->Branch("DphiJJ",&DphiJJ, "DphiJJ/F");
	treeWithBDT->Branch("RMS_eta",&RMS_eta, "RMS_eta/F");
	treeWithBDT->Branch("PtbalZH",&PtbalZH, "PtbalZH/F");
	treeWithBDT->Branch("EventPt",&EventPt, "EventPt/F");
	treeWithBDT->Branch("Angle",&Angle, "Angle/F");
	treeWithBDT->Branch("Centrality",&Centrality, "Centrality/F");
	treeWithBDT->Branch("MET",&MET, "MET/F");
	treeWithBDT->Branch("EvntShpAplanarity",&EvntShpAplanarity, "EvntShpAplanarity/F");
	treeWithBDT->Branch("EvntShpSphericity",&EvntShpSphericity, "EvntShpSphericity/F");
	treeWithBDT->Branch("EvntShpIsotropy",&EvntShpIsotropy, "EvntShpIsotropy/F");
	treeWithBDT->Branch("Zphi",&Zphi, "Zphi/F");
	treeWithBDT->Branch("Hphi",&Hphi, "Hphi/F");
	treeWithBDT->Branch("BDTvalue",&BDTvalue, "BDTvalue/F");

	
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
	if (Use["BDT"]) {
		histMattBdt     = new TH1F( "Matt_BDT",           "Matt_BDT",           15, -1.1, 0.35 );
		histTMVABdt     = new TH1F( "TMVA_BDT",           "TMVA_BDT",           88, -1.0, 0.1 );
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
   TString fname = "/home/hep/wilken/Weights/CMSSW_4_2_8_patch3/src/UserCode/wilken/CSVShapeCorr/DYJetsToLL_PtZ.root"; 
	double lumi = 100.0;
	Double_t  DYJetsToLL_PtZ_weight = lumi/(lumiZJH/2.0) ; //DYJetsToLL_PtZ-100_TuneZ2_7TeV-madgraph-tauola
	//DYJetsToLL_PtZ_weight = DYJetsToLL_PtZ_weight* ( 4242 /46179.0) ;
  
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
   TTree* theTree = (TTree*)input->Get("BDT_btree");
   
    theTree->SetBranchAddress( "Hmass", &Hmass );
    theTree->SetBranchAddress( "Zmass", &Zmass );
    theTree->SetBranchAddress( "Hpt", &Hpt );
	theTree->SetBranchAddress( "Zpt", &Zpt );
	theTree->SetBranchAddress( "CSV0", &CSV0 );
	theTree->SetBranchAddress( "CSV1", &CSV1 );
	theTree->SetBranchAddress( "DeltaPhiHV", &DeltaPhiHV );
	theTree->SetBranchAddress( "DetaJJ", &DetaJJ );
	theTree->SetBranchAddress( "UnweightedEta", &UnweightedEta );
	theTree->SetBranchAddress( "mu0pt", &mu0pt );
	theTree->SetBranchAddress( "Ht", &Ht );
	theTree->SetBranchAddress( "EvntShpCircularity", &EvntShpCircularity );
	theTree->SetBranchAddress( "nJets", &nJets );
	theTree->SetBranchAddress( "Naj", &Naj );
	theTree->SetBranchAddress( "mu1pt", &mu1pt );
	theTree->SetBranchAddress( "mu0pt", &mu0pt );
	theTree->SetBranchAddress( "EtaStandDev", &EtaStandDev );
	theTree->SetBranchAddress( "UnweightedEta", &UnweightedEta );
	theTree->SetBranchAddress( "jetCHF0", &jetCHF0 );
	theTree->SetBranchAddress( "jetCHF1", &jetCHF1 );
	theTree->SetBranchAddress( "muonPFiso0", &muonPFiso0 );
	theTree->SetBranchAddress( "muonPFiso1", &muonPFiso1 );
	theTree->SetBranchAddress( "DphiJJ", &DphiJJ );
	theTree->SetBranchAddress( "RMS_eta", &RMS_eta );
	theTree->SetBranchAddress( "PtbalZH", &PtbalZH );
	theTree->SetBranchAddress( "EventPt", &EventPt );
	theTree->SetBranchAddress( "Angle", &Angle );
	theTree->SetBranchAddress( "Centrality", &Centrality );
	theTree->SetBranchAddress( "EvntShpAplanarity", &EvntShpAplanarity );
	theTree->SetBranchAddress( "EvntShpSphericity", &EvntShpSphericity );
	theTree->SetBranchAddress( "EvntShpIsotropy", &EvntShpIsotropy );
	theTree->SetBranchAddress( "Trigweight", &Trigweight );
	theTree->SetBranchAddress( "B2011PUweight", &B2011PUweight );
	theTree->SetBranchAddress( "btag2CSF", &btag2CSF );
	theTree->SetBranchAddress( "MET", &MET );
	theTree->SetBranchAddress( "eventFlavor", &eventFlavor );
	

   // Efficiency calculator for cut method
   Int_t    nSelCutsGA = 0;
   Double_t effS       = 0.7;

   std::vector<Float_t> vecVar(4); // vector for EvaluateMVA tests

	TTree* TMVATree = (TTree*)input->Get("TMVA_tree");
	int NumInTree = theTree->GetEntries();
	int NumberTMVAtree = TMVATree->GetEntries();
	DYJetsToLL_PtZ_weight = DYJetsToLL_PtZ_weight* ( NumInTree /float(NumberTMVAtree)) ;
	
	std::vector<Float_t> vecVar(4); // vector for EvaluateMVA tests
	
	std::cout << "--- Processing: " << theTree->GetEntries() << " events" << std::endl;
	TStopwatch sw;
	sw.Start();
	for (Long64_t ievt=0; ievt<NumInTree;ievt++) {
		
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
	   if (eventFlavor == 5 ) {
		   histMattBdt    ->Fill( BDTvalue,DYJetsToLL_PtZ_weight );
		   histTMVABdt    ->Fill( BDTvalue,DYJetsToLL_PtZ_weight );
hMjj_OpenSelection->Fill(Hmass,DYJetsToLL_PtZ_weight);
hMmumu_OpenSelection->Fill(Zmass,DYJetsToLL_PtZ_weight);
hPtjj_OpenSelection->Fill(Hpt,DYJetsToLL_PtZ_weight);
hPtmumu_OpenSelection->Fill(Zpt,DYJetsToLL_PtZ_weight);
hCSV0_OpenSelection->Fill(CSV0,DYJetsToLL_PtZ_weight);
hCSV1_OpenSelection->Fill(CSV1,DYJetsToLL_PtZ_weight);
hdphiVH_OpenSelection->Fill(DeltaPhiHV,DYJetsToLL_PtZ_weight);
hdetaJJ_OpenSelection->Fill(DetaJJ,DYJetsToLL_PtZ_weight);
hUnweightedEta_OpenSelection->Fill(UnweightedEta,DYJetsToLL_PtZ_weight);
hPtmu0_OpenSelection->Fill(mu0pt,DYJetsToLL_PtZ_weight);
	   hHt_OpenSelection->Fill(Ht,DYJetsToLL_PtZ_weight);
	   hCircularity_OpenSelection->Fill(EvntShpCircularity,DYJetsToLL_PtZ_weight);
	   hCHFb0_OpenSelection->Fill(jetCHF0, DYJetsToLL_PtZ_weight);
	   hCHFb1_OpenSelection->Fill(jetCHF1, DYJetsToLL_PtZ_weight);
	   hPtbalZH_OpenSelection->Fill(PtbalZH, DYJetsToLL_PtZ_weight);
	   hPtmu1_OpenSelection->Fill(mu1pt, DYJetsToLL_PtZ_weight);
	   hPFRelIsomu0_OpenSelection->Fill(muonPFiso0, DYJetsToLL_PtZ_weight);
	   hPFRelIsomu1_OpenSelection->Fill(muonPFiso1, DYJetsToLL_PtZ_weight);
	   hNjets_OpenSelection->Fill(nJets, DYJetsToLL_PtZ_weight);
	   hRMSeta_OpenSelection->Fill(RMS_eta, DYJetsToLL_PtZ_weight);
	   hStaDeveta_OpenSelection->Fill(EtaStandDev, DYJetsToLL_PtZ_weight);
	   hdphiJJ_vect_OpenSelection->Fill(DphiJJ, DYJetsToLL_PtZ_weight);
	   hCentrality_OpenSelection->Fill(Centrality, DYJetsToLL_PtZ_weight);
	   hEventPt_OpenSelection->Fill(EventPt, DYJetsToLL_PtZ_weight);
	   hAngle_OpenSelection->Fill(Angle, DYJetsToLL_PtZ_weight);
	   hSphericity_OpenSelection->Fill(EvntShpSphericity, DYJetsToLL_PtZ_weight);
	   hAplanarity_OpenSelection->Fill(EvntShpAplanarity, DYJetsToLL_PtZ_weight);
	   hIsotropy_OpenSelection->Fill(EvntShpIsotropy, DYJetsToLL_PtZ_weight);
	   hDphiDetajj_OpenSelection->Fill(DphiJJ, DetaJJ, DYJetsToLL_PtZ_weight);
	   }// is b quark

	   treeWithBDT->Fill();
   }//end event loop
   
   // Get elapsed time
   sw.Stop();
   std::cout << "--- End of event loop: "; sw.Print();

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

   TFile *target  = new TFile( "TMVApp_DYJetsToLL_PtZ_bflavor.root","RECREATE" );
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
hAngle_OpenSelection->Write();
hSphericity_OpenSelection->Write();
hAplanarity_OpenSelection->Write();
hIsotropy_OpenSelection->Write();
hDphiDetajj_OpenSelection->Write();

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
	hAngle_OpenSelection->Delete();
	hSphericity_OpenSelection->Delete();
	hAplanarity_OpenSelection->Delete();
	hIsotropy_OpenSelection->Delete();
	hDphiDetajj_OpenSelection->Delete();
	
treeWithBDT->Delete();
	
	if (Use["BDT"          ])  {
		histMattBdt    ->Delete();
		histTMVABdt    ->Delete();
	}
	
    
   std::cout << "==> DYJetsToLL_PtZ_App_bflavor is done!" << endl << std::endl;
   
   gROOT->ProcessLine(".q");
} 
