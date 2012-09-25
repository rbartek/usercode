/**********************************************************************************
 * Project   : TMVA - a Root-integrated toolkit for multivariate data analysis    *
 * Package   : TMVA                                                               *
 * Exectuable: TMVARegressionApplication                                          *
 *                                                                                *
 * This macro provides a simple example on how to use the trained regression MVAs *
 * within an analysis module                                                      *
 **********************************************************************************/

#include <cstdlib>
#include <vector>
#include <iostream>
#include <map>
#include <string>

#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TString.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TStopwatch.h"
#include "TLorentzVector.h"
#if not defined(__CINT__) || defined(__MAKECINT__)
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"



#endif

using namespace TMVA;

double evalEt( double pt, double eta, double phi, double e){
  TLorentzVector j;
  j.SetPtEtaPhiE(pt,eta,phi, e );
  return j.Et(); 

}


double evalMt( double pt, double eta, double phi, double e){
  TLorentzVector j;
  j.SetPtEtaPhiE(pt,eta,phi, e );
  return j.Mt(); 

}


double ptHiggs(double pt0, double eta0, double phi0, double e0, double pt1, double eta1, double phi1 , double e1, double oldpt0, double oldpt1) {
  TLorentzVector j0, j1 , H;

  j0.SetPtEtaPhiE(pt0,eta0,phi0, e0 * pt0/oldpt0);
  j1.SetPtEtaPhiE(pt1,eta1,phi1, e1  * pt1/oldpt1);
  H = j0 +j1;
  return H.Pt(); 


}


double massHiggs(double pt0, double eta0, double phi0, double e0, double pt1, double eta1, double phi1 , double e1, double oldpt0, double oldpt1) {
  TLorentzVector j0, j1 , H;

  j0.SetPtEtaPhiE(pt0,eta0,phi0, e0 * pt0/oldpt0);
  j1.SetPtEtaPhiE(pt1,eta1,phi1, e1  * pt1/oldpt1);
  H = j0 +j1;
  return fabs(H.M()); 


}




void Process( TString  *fname , bool doTheFirst, TString myMethodList = "" ) 
{
  // 0 or 1
  int isFirst =doTheFirst;
  int isSecond = !doTheFirst;
  if ( (isFirst * isSecond) ==1) { std::cout << "error, you should either run on the first or the second daugthers...."<<std::endl;}

   //---------------------------------------------------------------
   // This loads the library
   TMVA::Tools::Instance();

   // Default MVA methods to be trained + tested
   std::map<std::string,int> Use;

   // --- Mutidimensional likelihood and Nearest-Neighbour methods
   Use["PDERS"]           = 0;
   Use["PDEFoam"]         = 1; 
   Use["KNN"]             = 1;
   // 
   // --- Linear Discriminant Analysis
   Use["LD"]		        = 1;
   // 
   // --- Function Discriminant analysis
   Use["FDA_GA"]          = 1;
   Use["FDA_MC"]          = 0;
   Use["FDA_MT"]          = 0;
   Use["FDA_GAMT"]        = 0;
   // 
   // --- Neural Network
   Use["MLP"]             = 1; 
   // 
   // --- Support Vector Machine 
   Use["SVM"]             = 0;
   // 
   // --- Boosted Decision Trees
   Use["BDT"]             = 0;
   Use["BDTG"]            = 1;
   // ---------------------------------------------------------------

   std::cout << std::endl;
   std::cout << "==> Start TMVARegressionApplication" << std::endl;

   // Select methods (don't look at this code - not of interest)
   if (myMethodList != "") {
      for (std::map<std::string,int>::iterator it = Use.begin(); it != Use.end(); it++) it->second = 0;

      std::vector<TString> mlist = gTools().SplitString( myMethodList, ',' );
      for (UInt_t i=0; i<mlist.size(); i++) {
         std::string regMethod(mlist[i]);

         if (Use.find(regMethod) == Use.end()) {
            std::cout << "Method \"" << regMethod << "\" not known in TMVA under this name. Choose among the following:" << std::endl;
            for (std::map<std::string,int>::iterator it = Use.begin(); it != Use.end(); it++) std::cout << it->first << " ";
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
  


   float hjetpt[2],  hjetptraw[2],   hjetptleadtrack[2], hjetphi[2],  hjetgenpt[2], hjeteta[2], hjetet[2], hjete[2],  hjetmt[2],  hjetchf[2], hjetnhf[2], hjetcef[2], hjetnef[2], hjetnconstituents[2], hjetnch[2], hjetvtxmass[2], hjetvtxpt[2],  hjetvtx3dl[2], hjetvtx3del[2];  
   float hjete[2], hjetJECUnc[2];
   float rho25, minDeltaPhijetMET, hpfMET;
   if (isFirst){
   
   //................ adapt to the trainign variables....
 reader->AddVariable( "hJet_pt", &hjetpt[0] );
 //reader->AddVariable( "hJet_ptRaw", &hjetptraw[0] );
 //reader->AddVariable( "rho25", &rho25 );
 reader->AddVariable( "hJet_eta", &hjeteta[0] );
 //reader->AddVariable( "hJet_et:=evalEt(hJet_pt, hJet_eta, hJet_phi, hJet_e)", &hjetet[0]);
 //reader->AddVariable( "hJet_Mt:=evalMt(hJet_pt, hJet_eta, hJet_phi, hJet_e)", &hjetmt[0]);
 //reader->AddVariable( "hJet_ptLeadTrack", &hjetptleadtrack[0] );
 reader->AddVariable( "hJet_chf", &hjetchf[0] );
  //reader->AddVariable( "hJet_cef", &hjetcef[0] );
 reader->AddVariable( "hJet_nconstituents", &hjetnconstituents[0] );
 //reader->AddVariable( "hJet_nch", &hjetnch[0] );
 reader->AddVariable( "hJet_vtxPt", &hjetvtxpt[0] );
 reader->AddVariable( "hJet_vtx3dL",  &hjetvtx3dl[0] );
	   reader->AddVariable( "hJet_vtx3deL", &hjetvtx3del[0] );
	   reader->AddVariable( "hJet_JECUnc", &hjetJECUnc[0] );
	   reader->AddVariable( "hJet_e", &hjete[0] );
	   reader->AddVariable( "pfMET:=METtype1corr.et", &hpfMET );
	   reader->AddVariable( "minDeltaPhijetMET", &minDeltaPhijetMET );
   } else {
 reader->AddVariable( "hJet_pt", &hjetpt[1] );
 //reader->AddVariable( "hJet_ptRaw", &hjetptraw[1] );
 //reader->AddVariable( "rho25", &rho25 );

 reader->AddVariable( "hJet_eta", &hjeteta[1] );
// reader->AddVariable( "hJet_et:=evalEt(hJet_pt, hJet_eta, hJet_phi, hJet_e)", &hjetet[1]);
// reader->AddVariable( "hJet_Mt:=evalMt(hJet_pt, hJet_eta, hJet_phi, hJet_e)", &hjetmt[1]);
 //reader->AddVariable( "hJet_ptLeadTrack", &hjetptleadtrack[1] );
 reader->AddVariable( "hJet_chf", &hjetchf[1] );
  //reader->AddVariable( "hJet_cef", &hjetcef[1] );
 reader->AddVariable( "hJet_nconstituents", &hjetnconstituents[1] );
 //reader->AddVariable( "hJet_nch", &hjetnch[1] );
 reader->AddVariable( "hJet_vtxPt", &hjetvtxpt[1] );
 reader->AddVariable( "hJet_vtx3dL",  &hjetvtx3dl[1] );
 reader->AddVariable( "hJet_vtx3deL", &hjetvtx3del[1] );
 	   reader->AddVariable( "hJet_JECUnc", &hjetJECUnc[1] );
	   reader->AddVariable( "hJet_e", &hjete[1] );
	   reader->AddVariable( "pfMET:=METtype1corr.et", &hpfMET );
	   reader->AddVariable( "minDeltaPhijetMET", &minDeltaPhijetMET );
	   

   }

 


   // Spectator variables declared in the training have to be added to the reader, too
/*   Float_t spec1,spec2;
   reader->AddSpectator( "spec1:=var1*2",  &spec1 );
   reader->AddSpectator( "spec2:=var1*3",  &spec2 );
*/
   // --- Book the MVA methods

   TString dir    = "weights/";
   //  isFirst ? dir = "weights0/" : dir = "weights1/";

   TString prefix = "TMVARegression";

   // Book method(s)
   for (std::map<std::string,int>::iterator it = Use.begin(); it != Use.end(); it++) {
      if (it->second) {
         TString methodName = it->first + " method";
         TString weightfile = dir + prefix + "_" + TString(it->first) + ".weights.xml";
         reader->BookMVA( methodName, weightfile ); 
      }
   }
   
   // Book output histograms
   TH1* hists[100];
   Int_t nhists = -1;
   for (std::map<std::string,int>::iterator it = Use.begin(); it != Use.end(); it++) {
      TH1* h = new TH1F( it->first.c_str(), TString(it->first) + " method", 100, -100, 600 );
      if (it->second) hists[++nhists] = h;
   }
   nhists++;
   
   // Prepare input tree (this must be replaced by your data source)
   // in this example, there is a toy tree with signal and one with background events
   // we'll later on use only the "signal" events for the test in this example.
   //   



   TFile *input(0);// = TFile::Open("../mjj/TestZnunuHbb120PU32.root");


   //   TString fname = "../Updated_DiJetPt_ZH_ZToNuNu_HToBB_M-120_7TeV-powheg-herwigppProcV1_Fall11.root";
       //  TString fname = "TMVAReg.root";

   if (!gSystem->AccessPathName( *fname )) {
      input = TFile::Open( *fname ); // check if file in local directory exists
   } 
   else { 
      input = TFile::Open( "http://root.cern.ch/files/tmva_reg_example.root" ); // if not: download from ROOT server
   }
   
   if (!input) {
      std::cout << "ERROR: could not open data file" << std::endl;
      exit(1);
   }
   std::cout << "--- TMVARegressionApp        : Using input file: " << input->GetName() << std::endl;

   // --- Event loop

   // Prepare the tree
   // - here the variable names have to corresponds to your tree
   // - you can use the same variables as above which is slightly faster,
   //   but of course you can use different ones and copy the values inside the event loop
   //
   TTree*  theTree= (TTree*)input->Get("tautree");
   std::cout << "--- Select signal sample" << std::endl;

   float hJet_genPtReg0 , hJet_genPtReg1; 
   
   if (isSecond) theTree->SetBranchAddress( "hJet_genPtReg0", &hJet_genPtReg0 );	   
   theTree->SetBranchAddress( "hJet_pt", &hjetpt );	  
   theTree->SetBranchAddress( "rho25", &rho25 );	  
    theTree->SetBranchAddress( "hJet_ptRaw", &hjetptraw );	   
   //if (fname->Contains("MET"))   theTree->SetBranchAddress( "hJet_ptRaw", &hjetptraw );	   
   //       if (!fname->Contains("MET")) theTree->SetBranchAddress( "hJet_ptRaw", &hjetptraw );	   
   theTree->SetBranchAddress( "hJet_ptLeadTrack", &hjetptleadtrack );	   
   theTree->SetBranchAddress( "hJet_phi", &hjetphi );	   
   theTree->SetBranchAddress( "hJet_genPt", &hjetgenpt );	   
   theTree->SetBranchAddress( "hJet_eta", &hjeteta );	   
   theTree->SetBranchAddress( "hJet_e", &hjete);		   
   theTree->SetBranchAddress( "hJet_chf", &hjetchf );	   
   theTree->SetBranchAddress( "hJet_nch", &hjetnch );	   
   theTree->SetBranchAddress( "hJet_nconstituents", &hjetnconstituents );	   
   theTree->SetBranchAddress( "hJet_nhf", &hjetnhf );	   
   theTree->SetBranchAddress( "hJet_cef", &hjetcef );	   
   theTree->SetBranchAddress( "hJet_chf", &hjetchf );	   
   theTree->SetBranchAddress( "hJet_vtxMass", &hjetvtxmass );
   theTree->SetBranchAddress( "hJet_vtxPt", &hjetvtxpt );
   theTree->SetBranchAddress( "hJet_vtx3dL", &hjetvtx3dl );  
   theTree->SetBranchAddress( "hJet_vtx3deL", &hjetvtx3del );
   
   




   std::cout << "--- Processing: " << theTree->GetEntries() << " events" << std::endl;
   TStopwatch sw;
   sw.Start();

   float newHiggsPt, newHiggsMass;

   if (isFirst) fname->ReplaceAll("/home/hep/wilken/Regression/src/VHbbAnalysis/VHbbDataFormats/bin/AddedrootFilesWithRegression/", "");
   if (isFirst) fname->ReplaceAll("/home/hep/wilken/Regression/src/VHbbAnalysis/VHbbDataFormats/bin/AddedrootFilesWithRegression/", "");
   if (isFirst) fname->ReplaceAll("DiJetPt", "DiJetPtRegNew");
//   if (isFirst) fname->ReplaceAll("/", "");
   TFile *nfile2 = new TFile(fname->Data() ,"recreate");
   //TFile *nfile2 = new TFile(newname.c_str() ,"recreate");
    TTree *ntree = theTree->CopyTree("");
    
    TBranch * b_hJet_genPtReg0 ;
    TBranch * b_hJet_genPtReg1 ;
    TBranch * b_newHiggsPt;
    TBranch * b_newHiggsMass;
    if (isFirst) {b_hJet_genPtReg0 = ntree->Branch("hJet_genPtReg0"   ,  &hJet_genPtReg0            ,  "hJet_genPtReg0/F"       );        } else{
    b_hJet_genPtReg1 = ntree->Branch("hJet_genPtReg1"   ,  &hJet_genPtReg1            ,  "hJet_genPtReg1/F"       );          
    b_newHiggsPt = ntree->Branch("newHiggsPt"   ,  &newHiggsPt            ,  "newHiggsPt/F"       );          
    b_newHiggsMass = ntree->Branch("newHiggsMass"   ,  &newHiggsMass            ,  "newHiggsMass/F"       );          
    }
    //    TBranch * b_hjetpt = ntree->Branch("hjetpt" ,  &hjetpt            ,  "hjetpt[2]/F");        
    // TBranch * b_hjeteta = ntree->Branch("hjeteta" ,  &hjeteta            ,  "hjeteta[2]/F");            	


     for (Long64_t ievt=0; ievt<theTree->GetEntries();ievt++) {
       // for (Long64_t ievt=0; ievt<10000;ievt++) {


      //theTree->GetEntry(ievt);
      ntree->GetEntry(ievt);

      // Retrieve the MVA target values (regression outputs) and fill into histograms
      // NOTE: EvaluateRegression(..) returns a vector for multi-target regression

      for (Int_t ih=0; ih<nhists; ih++) {
         TString title = hists[ih]->GetTitle();
         hjetet[0] = evalEt(hjetpt[0], hjeteta[0], hjetphi[0], hjete[0]);
         hjetet[1] = evalEt(hjetpt[1], hjeteta[1], hjetphi[1], hjete[1]);
         hjetmt[0] = evalMt(hjetpt[0], hjeteta[0], hjetphi[0], hjete[0]);
         hjetmt[1] = evalMt(hjetpt[1], hjeteta[1], hjetphi[1], hjete[1]);

         Float_t val = (reader->EvaluateRegression( title ))[0];
         hists[ih]->Fill( val );
             
      }

      if (isFirst) {
	hJet_genPtReg0 = (reader->EvaluateRegression( "BDT method" )[0])     ; 
	b_hJet_genPtReg0->Fill();
      }else {
	hJet_genPtReg1 = (reader->EvaluateRegression( "BDT method" )[0])     ; 
	b_hJet_genPtReg1->Fill();
      newHiggsPt =     ptHiggs(hJet_genPtReg0, hjeteta[0], hjetphi[0], hjete[0], hJet_genPtReg1, hjeteta[1], hjetphi[1], hjete[1], hjetpt[0], hjetpt[1]);
      newHiggsMass = massHiggs(hJet_genPtReg0, hjeteta[0], hjetphi[0], hjete[0], hJet_genPtReg1, hjeteta[1], hjetphi[1], hjete[1], hjetpt[0], hjetpt[1]);
      
       b_newHiggsPt->Fill();
       b_newHiggsMass->Fill();
      }
     


      if (ievt%10000 == 0) {
         std::cout << "--- ... Processing event: " << ievt << std::endl;
      if (isFirst) 
	{std::cout << "regres output, hjet pt 0 , hJet gen pt 0  " << hJet_genPtReg0 << " , " << hjetpt[0] << " , " << hjetgenpt[0] << "\n";}
      if (isSecond) {
	std::cout << "regres output, hjet pt 0 , hJet gen pt 0  " << hJet_genPtReg0 << " , " << hjetpt[0] << " , " << hjetgenpt[0] << "\n";
      std::cout << "regres output, hjet pt 1 , hJet gen pt 1  " << hJet_genPtReg1 << " , " << hjetpt[1] << " , " << hjetgenpt[1] << "\n";
      std::cout << "new Higgs mass, pt  " << newHiggsMass << " , " << newHiggsPt  << "\n";
      }

      }
    

      
     }
     //  ntree->Fill(); 
  sw.Stop();
   std::cout << "--- End of event loop: "; sw.Print();
 
    
 
     
    ntree->Write("", TObject::kOverwrite);

    nfile2->Write();

   // --- Write histograms

   TFile *target  = new TFile( "TMVARegApp.root","RECREATE" );
   for (Int_t ih=0; ih<nhists; ih++) hists[ih]->Write();
   target->Close();

   std::cout << "--- Created root file: \"" << target->GetName() 
             << "\" containing the MVA output histograms" << std::endl;
  
   delete reader;
    
   std::cout << "==> TMVARegressionApplication is done!" << std::endl << std::endl;
}

void TMVARegressionApplicationProcessor(){

const int size =1;

string titles[size] = {

//  "../samples/07Jun2012Ntuples_Step1V32_Step2V2/07Jun2012Ntuples_Step1V32_Step2V2/DiJetPt_ZH_ZToNuNu_HToBB_M-115_8TeV-powheg-herwigpp.root",
  //"../samples/07Jun2012Ntuples_Step1V32_Step2V2/07Jun2012Ntuples_Step1V32_Step2V2/DiJetPt_ZH_ZToNuNu_HToBB_M-135_8TeV-powheg-herwigpp.root",

// "/home/hep/wilken/Regression/src/VHbbAnalysis/VHbbDataFormats/bin/AddedrootFilesWithRegression/RegTrig_DYPtZ.root",
 "/home/hep/wilken/Regression/src/VHbbAnalysis/VHbbDataFormats/bin/AddedrootFilesWithRegression/RegTrig_TTJets.root",
// "/home/hep/wilken/Regression/src/VHbbAnalysis/VHbbDataFormats/bin/AddedrootFilesWithRegression/RegTrig_T_tW.root",
// "/home/hep/wilken/Regression/src/VHbbAnalysis/VHbbDataFormats/bin/AddedrootFilesWithRegression/RegTrig_ZH125Fall11.root",

//	"/home/hep/wilken/Regression/src/VHbbAnalysis/VHbbDataFormats/bin/RegTrig_ZH125Fall11.root",
  /*
  "../../April25/TTbar/TTbar_1.root",
  "../../April25/TTbar/TTbar_2.root",
  "../../April25/TTbar/TTbar_3.root",
  "../../April25/TTbar/TTbar_4.root",
  "../../April25/TTbar/TTbar_5.root",
  "../../April25/TTbar/TTbar_6.root",*/
  /*
  "../../April25/SingleMu/SingleMu_1.root",
  "../../April25/SingleMu/SingleMu_2.root",
  "../../April25/SingleMu/SingleMu_3.root",
  "../../April25/SingleMu/SingleMu_4.root",
  "../../April25/SingleMu/SingleMu_5.root",
  "../../April25/SingleMu/SingleMu_6.root",
  "../../April25/SingleMu/SingleMu_7.root",
  "../../April25/SingleMu/SingleMu_8.root",
  "../../April25/SingleMu/SingleMu_9.root",
  "../../April25/SingleMu/SingleMu_10.root",
  "../../April25/SingleMu/SingleMu_11.root",
  "../../April25/SingleMu/SingleMu_12.root",
  "../../April25/SingleMu/SingleMu_13.root",
  "../../April25/SingleMu/SingleMu_14.root",
  "../../April25/SingleMu/SingleMu_15.root",
  "../../April25/SingleMu/SingleMu_16.root",
  
  
"../../April25/TTbar/TTbar_7.root",
  "../../April25/TTbar/TTbar_8.root",
  "../../April25/TTbar/TTbar_9.root",
  "../../April25/TTbar/TTbar_10.root",
  "../../April25/TTbar/TTbar_11.root",
  "../../April25/TTbar/TTbar_12.root",
  "../../April25/TTbar/TTbar_13.root",
  "../../April25/TTbar/TTbar_14.root",
  "../../April25/TTbar/TTbar_15.root",
  "../../April25/TTbar/TTbar_16.root",
  "../../April25/TTbar/TTbar_17.root",

  //"../samples/CorApril25WJetsPt100MadgraphWJets_Pt100Madgraph_Fall11.root",

  /* "../samples/CorApril25ZnunuH_115_summer11_1.root",
  "../samples/CorApril25ZnunuH_120_summer11_1.root",
  "../samples/CorApril25ZnunuH_125_summer11_1.root",
  */
  /*"../samples/CorApril25TTbarTTbar_1.root",
"../samples/CorApril25TTbarTTbar_10.root",
"../samples/CorApril25TTbarTTbar_11.root",
"../samples/CorApril25TTbarTTbar_12.root",
"../samples/CorApril25TTbarTTbar_13.root",
"../samples/CorApril25TTbarTTbar_14.root",
"../samples/CorApril25TTbarTTbar_15.root",
"../samples/CorApril25TTbarTTbar_16.root",
"../samples/CorApril25TTbarTTbar_17.root",
"../samples/CorApril25TTbarTTbar_2.root",
"../samples/CorApril25TTbarTTbar_3.root",
"../samples/CorApril25TTbarTTbar_4.root",
"../samples/CorApril25TTbarTTbar_5.root",
"../samples/CorApril25TTbarTTbar_6.root",
"../samples/CorApril25TTbarTTbar_7.root",
"../samples/CorApril25TTbarTTbar_8.root",
"../samples/CorApril25TTbarTTbar_9.root",
  */

  /*  
"../samples/CorApril25METBTagMETBtag_1.root",
"../samples/CorApril25METBTagMETBtag_2.root",
"../samples/CorApril25METMET_1.root",
"../samples/CorApril25METMET_2.root",
"../samples/CorApril25METMET_3.root",
"../samples/CorApril25METMET_4.root",
"../samples/CorApril25METMET_5.root",
"../samples/CorApril25METMET_6.root",
"../samples/CorApril25METMET_7.root",
"../samples/CorApril25METMET_8.root",
"../samples/CorApril25METMET_9.root",
  
"../samples/CorApril25SingleT_sChSingleT_sCh_Fall11.root",
"../samples/CorApril25SingleT_tChSingleT_tCh_Fall11.root",
"../samples/CorApril25SingleT_tWSingleT_tW_Fall11.root",
"../samples/CorApril25SingleTbar_sChSingleTbar_sCh_Fall11.root",
"../samples/CorApril25SingleTbar_tChSingleTbar_tCh_Fall11.root",
"../samples/CorApril25SingleTbar_tWSingleTbar_tW_Fall11.root",
"../samples/CorApril25WH_100.root",
"../samples/CorApril25WH_105.root",
"../samples/CorApril25WH_110.root",
"../samples/CorApril25WH_115.root",
"../samples/CorApril25WH_120.root",
"../samples/CorApril25WH_125.root",
"../samples/CorApril25WH_130.root",
"../samples/CorApril25WH_135.root",

"../samples/CorApril25WWWW_Fall11.root",
"../samples/CorApril25WZWZ_Fall11.root",
"../samples/CorApril25WbbMadgraphWbbMadgraph_Fall11.root",
"../samples/CorApril25ZH_100.root",
"../samples/CorApril25ZH_105.root",
"../samples/CorApril25ZH_110.root",

"../samples/CorApril25METBTagMETBtag_1.root",
  */
/*
"../samples/CorApril25ZH_115.root",
"../samples/CorApril25ZH_125.root",
"../samples/CorApril25ZZZZ_Fall11.root",
"../samples/CorApril25ZZZnunuJets_Pt100Madgraph_1.root",
"../samples/CorApril25ZnunuH_100.root",
"../samples/CorApril25ZnunuH_105.root",
*/
/*
"../samples/CorApril25ZnunuH_110.root",
"../samples/CorApril25ZnunuH_115.root",
"../samples/CorApril25ZnunuH_120.root",
"../samples/CorApril25ZnunuH_125.root",
"../samples/CorApril25ZnunuH_130.root",
"../samples/CorApril25ZnunuH_135.root",
"../samples/CorApril25ZnunuJetsPt100HerwigppZnunuJets_Pt100Herwigpp_1.root",
"../samples/CorApril25ZnunuJetsPt100HerwigppZnunuJets_Pt100Herwigpp_2.root",
"../samples/CorApril25ZnunuJetsPt100MadgraphZnunuJets_Pt100Madgraph_1.root",
"../samples/CorApril25ZnunuJets_HT100to200_MadgraphZnunuJets_HT_100to200_Madgraph_1.root",
"../samples/CorApril25ZnunuJets_HT100to200_MadgraphZnunuJets_HT_100to200_Madgraph_2.root",
"../samples/CorApril25ZnunuJets_HT200toInf_MadgraphZnunuJets_HT_200toInf_Madgraph_1.root",
"../samples/CorApril25ZnunuJets_HT200toInf_MadgraphZnunuJets_HT_200toInf_Madgraph_2.root",
"../samples/CorApril25ZnunuJets_HT50to100_MadgraphZnunuJets_HT_50to100_Madgraph_1.root",
"../samples/CorApril25ZnunuJets_HT50to100_MadgraphZnunuJets_HT_50to100_Madgraph_2.root",
"../samples/CorApril25ZnunuJets_HT50to100_MadgraphZnunuJets_HT_50to100_Madgraph_3.root",
"../samples/CorApril25ZnunuH_120_summer11_1.root"


  "../samples/CorApril25ZnunuH_125_summer11_1.root",
  "../samples/CorApril25ZnunuH_130_summer11_1.root",
  "../samples/CorApril25ZnunuH_135_summer11_1.root",
*/
};

for(int k=0;k<size;k++) {
  std::cout<< " file " << titles[k] << std::endl;
TString *fTT  = new TString(titles[k]);
 Process(fTT, 1, "BDT");
 Process(fTT, 0 , "BDT");

}
}
