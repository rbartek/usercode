#include <iostream>
#include <fstream>
#include <vector>


#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TMath.h"
#include "TLeaf.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TLegend.h"
#include "TProfile.h"
//#include "btagshape.h"

#include "/home/hep/wilken/taus/CMSSW_4_4_2_patch8/src/UserCode/wilken/interface/xsecV21.h" 
#include <cassert>

using std::cout;
using std::endl;
using std::vector;
using std::ofstream;

double lumi = 4.457; // fb-1
double lA = 2.219;
double lB = 2.238;

const float bjetsys=0.06;
const float cjetsys=0.12;
const float ljetsys=0.15;

TFile *inputS;
TFile *inputS2;
TFile* inputTT;
TFile* inputTtW;
TFile* inputTBARtW;
TFile *inputZZ;
TFile *inputWW;
TFile *inputWZ;
TFile *inputDY;
TFile *inputDY2;
TFile *inputWJ;
TFile *inputdata;

double bdtcut=-1.50;// //snippet not to cut on bdt -> put -1.5 
//directory where to save datacards

string dir="/home/hep/wilken/taus/CMSSW_4_4_2_patch8/src/UserCode/wilken/V21/";

enum {
	kl1pt, //0
	kl2pt, //1
	kpTjj, //2
	kpTZ, //3
	kcsv1, //4
	kcsv2, //5 
	kDphi, //6
	kNaj, //7 
	kmet, //8
	kDeta, //9
	kHt, //10
	kjchf1, //11
	kjchf2, //12
	kangleEMU, //13
	kueta, //14 
	kcentrality, //15
	kM, //16
	kDptHZ, //17
	hdrjj, //18
	kbdt, //19
	kdeltaPhiZMET,  //20
	kAplanarity, //21
	kSphericity, //22
	kCircularity,//23
	kIsotropy,//24
	kMte,   //25
	kMtmu   //26
};

enum {
	noboost,
	klowboosted,
	kboosted,
	khighboosted
};

bool PassBoost(double pth, double ptz, int boost, float Mjj)
{
	bool passes=false;
	
	if(boost==kboosted)
		if(ptz>100 && pth>100) passes=true;
	if(boost==klowboosted)
		if(ptz<=100 && pth<=100 && ptz>50 && pth>50 ) passes=true;
	if ( boost==khighboosted) 
		if( ptz>50 && pth>50 && Mjj>90 && Mjj<130)  passes=true;
	if(boost==noboost)passes=true;
	return passes;
}

const float Mjjcutl=90;
const float Mjjcuth=125;
const float csv1cut=0.9;
const float csv2cut=0.5;
const float Dphicut=3;
const float metcut=30;
const float hthcut=240; 
const float zhanglecut=2.4;
const float uetacut=3.5;

const float SFtt=1.10;
const float SFZl=0.95;
const float SFZh=1.;

//float mass2dcutl[4]={0,70,90,120};
//float mass2dcuth[4]={70,90,120,2000};

/*float mass2dcutl[N];
 float mass2dcuth[N];
 float bdt2dcutl[N];
 float bdt2dcuth[N];
 */

/*
 float mass2dcutl[12]={0,35,50,70,80,90,97,105,112,120,180,250};
 float mass2dcuth[12]={35,50,70,80,90,97,105,112,120,180,250,2000};
 float bdt2dcutl[12]={-1.5,-1.5,-1.5,-1.5,-1.5,-1.5,-1.5,-1.5,-1.5,-1.5,-1.5,-1.5};
 float bdt2dcuth[12]={1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5};
 */

//nuoivi 6 bins 
//float mass2dcutl[6]={0,  55,  80, 98,  132, 160};
//float mass2dcuth[6]={55, 80, 98, 132, 160, 250};
//float bdt2dcutl[6]={-1.5,-1.5,-1.5,-1.5,-1.5,-1.5};
//float bdt2dcuth[6]={1.5,1.5,1.5,1.5,1.5,1.5};

const int N=1; //total bins
const int Nm=1; //bins in mass 
//float mass2dcutl[4]={0,  70,  98,  132};
//float mass2dcuth[4]={70, 98,  132, 250};
//float bdt2dcutl[4]={-1.5,-1.5,-1.5,-1.5};
//float bdt2dcuth[4]={1.5,1.5,1.5,1.5};


float mass2dcutl[3]={0,   98,  132};
float mass2dcuth[3]={98,  132, 250};
float bdt2dcutl[3]={-1.5,-1.5,-1.5};
float bdt2dcuth[3]={1.5,1.5,1.5};

//float mass2dcutl[5]={0,  70,  98,  115, 132};
//float mass2dcuth[5]={70, 98,  115, 132, 250};
//float bdt2dcutl[5]={-1.5,-1.5,-1.5,-1.5,-1.5};
//float bdt2dcuth[5]={1.5,1.5,1.5,1.5,1.5};

//float mass2dcutl[6]={0,  70,  95,  105, 115, 132};
//float mass2dcuth[6]={70, 95,  105, 115, 132, 250};
//float bdt2dcutl[6]={-1.5,-1.5,-1.5,-1.5,-1.5,-1.5};
//float bdt2dcuth[6]={1.5,1.5,1.5,1.5,1.5,1.5};

//2D analysis, 12 bins
//float mass2dcutl[6]={0,  70,  95,  105, 115, 132};
//float mass2dcuth[6]={70, 95,  105, 115, 132, 250};
//float bdt2dcutl[2]={-1.5, -0.4};
//float bdt2dcuth[2]={-0.4, 1.5};

//float mass2dcutl[7]={0,  60,  85, 95,  105, 115, 132};
//float mass2dcuth[7]={60, 85,  95, 105, 115, 132, 250};
//float bdt2dcutl[7]={-1.5,-1.5,-1.5,-1.5,-1.5,-1.5,-1.5};
//float bdt2dcuth[7]={1.5,1.5,1.5,1.5,1.5,1.5,1.5};


//float mass2dcutl[8]={0,  60,  85, 95,  105, 115, 125, 135};
//float mass2dcuth[8]={60, 85,  95, 105, 115, 125, 135, 250};
//float bdt2dcutl[8]={-1.5,-1.5,-1.5,-1.5,-1.5,-1.5,-1.5,-1.5};
//float bdt2dcuth[8]={1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5};

//float mass2dcutl[9]={0,  50,  70, 80, 90,  100, 110, 120, 130};
//float mass2dcuth[9]={50, 70,  80, 90, 100, 110, 120, 130, 250};
//float bdt2dcutl[9]={-1.5,-1.5,-1.5,-1.5,-1.5,-1.5,-1.5,-1.5, -1.5};
//float bdt2dcuth[9]={1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5, 1.5};

//here 12 bins
//float mass2dcutl[4]={0,70,90,120};
//float mass2dcuth[4]={70,90,120,250};
//float bdt2dcutl[3]={-1.5,-0.9,-0.6};
//float bdt2dcuth[3]={-0.9,-0.6,1.1};

//here 8 bins, 2d
//float mass2dcutl[4]={0,70,90,120};
//float mass2dcuth[4]={70,90,120,2000};
//float bdt2dcutl[2]={-1.5,-0.6};
//float bdt2dcuth[2]={-0.6,1.1};

//here 12 bins
//float mass2dcutl[4]={0,70,90,120};
//float mass2dcuth[4]={70,90,120,2000};
//float bdt2dcutl[3]={-1.5,-0.9,-0.6};
//float bdt2dcuth[3]={-0.9,-0.6,1.1};

//here 16 bins
//float mass2dcutl[4]={0,70,90,120};
//float mass2dcuth[4]={70,90,130,2000};
//float bdt2dcutl[4]={-1.5,-0.9,-0.6,-0.3};
//float bdt2dcuth[4]={-0.9,-0.6,-0.3,1.1};

//here 20 bins
//float mass2dcutl[4]={0,70,90,120};
//float mass2dcuth[4]={70,90,130,2000};
//float bdt2dcutl[5]={-1.5,-0.9,-0.7,-0.5,-0.3};
//float bdt2dcuth[5]={-0.9,-0.7,-0.5,-0.3,1.1};

void makeDataCard(float eVH, float eWj, float eZjLF, float eZjHF,float eTT, float es_Top, float eVV,
				  float e_eVH, float e_eWj, float e_eZjLF, float e_eZjHF,float e_eTT, float e_es_Top,
				  float e_eVV, float bdt)
{
	
	ofstream myfile (TString::Format("%s/dataupdated%0.03f.txt",dir.c_str(),bdt).Data());
	
	myfile<<"imax 1  number of channels"<<endl;
	myfile<<"jmax 7  number of backgrounds"<<endl;;
	myfile<<"kmax *  number of nuisance parameters (sources of systematical uncertainties)"<<endl;;
	myfile<<"bin Zemu"<<endl;;
	myfile<<"observation 0"<<endl;
	myfile<<"bin \t Zemu \t Zemu  \t Zemu \t Zemu  \t Zemu \t Zemu \t Zemu  \t Zemu"<<endl;
	myfile<<"process \t VH \t Wj \t ZjLF \t ZjHF \t TT \t s_Top \t VV \t QCD"<<endl;
	myfile<<"process  \t 0 \t 1 \t 2  \t 3 \t 4  \t 5  \t 6  \t 7 "<<endl;
	//here put the expectations as computed from the macro
	myfile<<TString::Format("rate \t %0.3f  \t %0.3f \t %0.3f \t %0.3f \t %0.3f \t %0.3f \t %0.3f \t 0.000",
							eVH,eWj,eZjLF,eZjHF,eTT,es_Top,eVV).Data()<<endl;
	myfile<<"\n";
	myfile<<"lumi \t lnN \t 1.045 \t - \t - \t - \t - \t 1.045 \t 1.045 \t 1.045"<<endl; 
	myfile<<"pdf_qqbar \t lnN \t 1.01  \t - \t - \t - \t - \t - \t 1.01 \t -"<<endl; 
	myfile<<"pdf_gg \t lnN \t - \t - \t - \t - \t - \t 1.01 \t - \t 1.01"<<endl; 
	myfile<<"QCDscale_VH \t lnN \t 1.04  \t - \t - \t - \t - \t - \t - \t -"<<endl; 
	myfile<<"QCDscale_ttbar \t lnN \t -  \t - \t - \t - \t - \t 1.06 \t - \t -"<<endl;
	myfile<<"QCDscale_VV \t lnN \t -  \t - \t - \t - \t - \t - \t 1.035 \t -"<<endl; 
	myfile<<"QCDscale_QCD \t lnN \t -  \t - \t - \t - \t - \t - \t - \t 1.30"<<endl;
	myfile<<"CMS_vhbb_boost_EWK \t lnN  \t 1.05  \t - \t - \t - \t - \t - \t - \t - "<<endl; 
	myfile<<"CMS_vhbb_boost_QCD \t lnN \t 1.10 \t - \t - \t - \t - \t - \t - \t -"<<endl; 
	myfile<<"CMS_vhbb_ttbar \t lnN  \t - \t - \t - \t - \t - \t 1.29 \t - \t -"<<endl;       
	myfile<<"CMS_vhbb_VV \t lnN \t - \t - \t - \t - \t - \t - \t 1.295 \t -"<<endl;       
	myfile<<"CMS_trigger_m  \t lnN  \t 1.01  \t - \t - \t - \t - \t 1.01 \t 1.01 \t -"<<endl; 
	myfile<<"CMS_eff_m  \t lnN  1.04 \t -  \t - \t - \t - 1.04 1.04 1.04"<<endl; 
	myfile<<"CMS_eff_b  \t lnN  \t 1.10 \t 1.10 \t 1.10 \t 1.00 \t 1.10 \t 1.10 \t 1.10 \t 1.10"<<endl;  
	myfile<<"CMS_scale_j  \t lnN  \t 1.01 \t -  \t -  \t -  \t -  \t 1.01 \t 1.01 \t 1.01"<<endl;  
	myfile<<"CMS_res_j  \t lnN  \t 1.10 \t 1.10 \t 1.10 \t 1.10 \t 1.10 \t 1.10 \t 1.10 \t 1.10"<<endl;  
	
	//systematics associated to SF
	myfile<<"CMS_vhbb_ZjLF_SF \t lnN \t  -  \t - \t 1.049 \t - \t - \t - \t - \t -"<<endl; 
	myfile<<"CMS_vhbb_ZjHF_SF \t lnN  \t - \t - \t - \t 1.092 \t - \t - \t - \t -"<<endl; 
	myfile<<"CMS_vhbb_TT_SF  \t lnN  \t - \t - \t - \t - \t 1.053 \t - \t - \t -"<<endl; 
	
	//stat uncertainties on expectation (I think)
	myfile<<TString::Format("CMS_vhbb_stats_VH_Zemu \t lnN  \t %0.3f \t - \t - \t - \t - \t - \t - \t -", 
							1.+e_eVH/eVH).Data()<<endl;  
	myfile<<"CMS_vhbb_stats_Wj_Zemu \t lnN  \t - \t 1.00 \t - \t - \t - \t - \t - \t -"<<endl;
	myfile<<TString::Format("CMS_vhb_stats_ZjLF_Zemu \t lnN  \t - \t - \t %0.3f \t - \t - \t - \t - \t -",
							1.+e_eZjLF/eZjLF).Data()<<endl;  
	myfile<<TString::Format("CMS_vhbb_stats_ZjHF_Zemu \t lnN  \t - \t - \t - \t %0.3f \t - \t - \t - \t -",
							1.+e_eZjHF/eZjHF).Data()<<endl;    
	myfile<<TString::Format("CMS_vhbb_stats_TT_Zemu \t lnN   \t - \t - \t - \t - \t %0.3f \t - \t - \t -",
							1.+e_eTT/eTT).Data()<<endl;   
	myfile<<TString::Format("CMS_vhbb_stats_sT_Zemu \t lnN  \t - \t - \t - \t - \t - \t %0.3f \t - \t -",
							1+e_es_Top/es_Top).Data()<<endl;      
	// myfile<<"CMS_vhbb_stats_sT_Zemu \t lnN  \t - \t - \t - \t - \t - \t 2.00 \t - \t -"<<endl;      
	myfile<<TString::Format("CMS_vhbb_stats_VV_Zemu \t lnN  \t - \t - \t - \t - \t - \t - \t %0.3f \t -",
							1+e_eVV/eVV).Data()<<endl;       
	myfile.close();
}

void makeDataCardShape(float eVH, float eWj, float eZjLF, float eZjHF,float eTT, float es_Top, float eVV,
					   float e_eVH, float e_eWj, float e_eZjLF, float e_eZjHF,float e_eTT, float e_es_Top,
					   float e_eVV, float data, float bdt, string shapefile, string bs, int n2dbins=-1)
{
cout << "making shape datacard" << endl;
	ofstream myfile; 
	if(n2dbins==-1)  myfile.open(TString::Format("%sdataUpdatedShape%s%0.03f.txt", dir.c_str(),bs.c_str(),bdt).Data());
	else myfile.open(TString::Format("%sdataUpdatedShape%s.txt", 
									 dir.c_str(),bs.c_str()).Data());
	
	myfile<<"imax 1  number of channels"<<endl;
	myfile<<"jmax 7  number of backgrounds"<<endl;;
	myfile<<"kmax *  number of nuisance parameters (sources of systematical uncertainties)"<<endl;;
	myfile<<TString::Format("shapes * * %s $PROCESS $PROCESS_$SYSTEMATIC", shapefile.c_str()).Data()<<endl;
	myfile<<TString::Format("bin Zemu%s",bs.c_str()).Data()<<endl;;
	myfile<<TString::Format("observation %0.3f", data).Data()<<endl;
	cout<<"observation null " << endl;
	myfile<<TString::Format("bin \t Zemu%s \t Zemu%s \t Zemu%s  \t Zemu%s \t Zemu_noboosted \t Zemu_noboosted  \t Zemu_noboosted \t Zemu_noboosted",
							bs.c_str(),bs.c_str(),bs.c_str(),bs.c_str()).Data()<<endl;
	cout << "stupid line in data card" << endl;
	myfile<<"process \t VH \t Wj \t ZjLF \t ZjHF \t TT \t s_Top \t VV \t QCD"<<endl;
	myfile<<"process  \t 0 \t 1 \t 2  \t 3 \t 4  \t 5  \t 6  \t 7"<<endl;
	//here put the expectations as computed from the macro
	
	myfile<<TString::Format("rate \t %0.3f \t %0.3f \t %0.3f \t %0.3f \t %0.3f \t %0.3f \t %0.3f \t 0.000",
							eVH,eWj,eZjLF,eZjHF,eTT,es_Top,eVV).Data()<<endl;
	myfile<<"\n";
	myfile<<"lumi \t lnN \t 1.045 \t - \t - \t - \t - \t 1.045 \t 1.045 \t 1.045"<<endl; 
	
	myfile<<"pdf_qqbar \t lnN \t 1.01 \t - \t - \t - \t - \t - \t 1.01 \t -"<<endl; 
	myfile<<"pdf_gg \t lnN \t - \t - \t - \t - \t - \t 1.01 \t - \t 1.01"<<endl;
	
	myfile<<"QCDscale_VH \t lnN \t 1.04 \t - \t - \t - \t - \t - \t - \t -"<<endl; 
	myfile<<"QCDscale_ttbar \t lnN \t - \t - \t - \t - \t - \t 1.06 \t - \t -"<<endl;
	myfile<<"QCDscale_VV \t lnN \t - \t - \t - \t - \t - \t - \t 1.04 \t -"<<endl; 
	myfile<<"QCDscale_QCD \t lnN \t - \t - \t - \t - \t - \t - \t - \t 1.30"<<endl;
	
	myfile<<"CMS_vhbb_boost_EWK \t lnN  \t 1.05  \t - \t - \t - \t - \t - \t - \t - "<<endl; 
	myfile<<"CMS_vhbb_boost_QCD \t lnN \t 1.10 \t - \t - \t - \t - \t - \t - \t -"<<endl; 
	
	//myfile<<"CMS_vhbb_ttbar \t lnN \t - \t - \t - \t - \t - \t 1.30 \t - \t -"<<endl;       
	myfile<<"CMS_vhbb_ST \t lnN \t - \t - \t - \t - \t - \t 1.29 \t - \t -"<<endl;       
	myfile<<"CMS_vhbb_VV \t lnN \t - \t - \t - \t - \t - \t - \t 1.30 \t -"<<endl;       

	//systematics associated to SF
	myfile<<"CMS_vhbb_ZjLF_SF \t lnN \t  -  \t - \t 1.06 \t - \t - \t - \t - \t -"<<endl; 
	myfile<<"CMS_vhbb_ZjHF_SF \t lnN  \t -  \t - \t - \t 1.17 \t - \t - \t - \t -"<<endl; 
	myfile<<"CMS_vhbb_TT_SF  \t lnN  \t - \t - \t - \t - \t 1.21 \t - \t - \t -"<<endl; 
	
	myfile<<"CMS_trigger_m  \t lnN  \t 1.01  \t - \t - \t - \t - \t 1.01 \t 1.01 \t -"<<endl; 
	myfile<<"CMS_trigger_3  \t lnN  \t 1.02  \t - \t - \t - \t - \t 1.02 \t 1.02 \t -"<<endl; 
	myfile<<"CMS_eff_m  \t lnN  \t 1.04 \t - \t - \t - \t - \t 1.04 \t 1.04 \t 1.04"<<endl; 
	myfile<<"CMS_eff_e  \t lnN  \t 1.04 \t - \t - \t - \t - \t 1.04 \t 1.04 \t 1.04"<<endl; 
	myfile<<"CMS_eff_b  \t lnN  \t 1.11 \t 1.07 \t 1.07 \t 1 \t 1 \t 1.05 \t 1.11 \t -"<<endl; 
	myfile<<"CMS_fake_b  \t lnN  \t 1.05 \t 1.12 \t 1.12 \t 1 \t 1 \t 1.15 \t 1 \t -"<<endl; 
	myfile<<"CMS_scale_j  \t lnN  \t 1.03 \t -  \t -  \t -  \t -  \t 1.03 \t 1.03 \t -"<<endl;  
	myfile<<"CMS_res_j  \t lnN  \t 1.05 \t 1.04 \t 1.04 \t 1.04 \t 1.04 \t 1.04 \t 1.05 \t -"<<endl;  
	
	
	//stat uncertainties on expectation (I think)
	myfile<<TString::Format("CMS_vhbb_stats_VH_Zemu%s \t lnN  \t %0.3f \t - \t - \t - \t - \t - \t - \t -",bs.c_str(), 
							1.+e_eVH/eVH).Data()<<endl;  
	myfile<<TString::Format("CMS_vhbb_stats_Wj_Zemu%s \t lnN  \t - \t 1.00 \t - \t - \t - \t - \t - \t -",bs.c_str()).Data()<<endl;
	myfile<<TString::Format("CMS_vhb_stats_ZjLF_Zemu%s \t lnN  \t - \t - \t %0.3f \t - \t - \t - \t - \t -",bs.c_str(),
							1.+e_eZjLF/eZjLF).Data()<<endl;  
	myfile<<TString::Format("CMS_vhbb_stats_ZjHF_Zemu%s \t lnN  \t - \t - \t - \t %0.3f \t - \t - \t - \t -",bs.c_str(),
							1.+e_eZjHF/eZjHF).Data()<<endl;    
	myfile<<TString::Format("CMS_vhbb_stats_TT_Zemu%s \t lnN  \t -  \t - \t - \t - \t %0.3f \t - \t - \t -",bs.c_str(),
							1.+e_eTT/eTT).Data()<<endl;   
	myfile<<TString::Format("CMS_vhbb_stats_sT_Zemu%s \t lnN  \t - \t - \t - \t - \t - \t %0.3f \t - \t -",bs.c_str(),
							1+e_es_Top/es_Top).Data()<<endl;      
	//myfile<<TString::Format("CMS_vhbb_stats_sT_Zemu%s \t lnN  \t - \t - \t - \t - \t - \t 2.00 \t - \t -",bs.c_str()).Data()<<endl; 
	myfile<<TString::Format("CMS_vhbb_stats_VV_Zemu%s \t lnN  \t - \t - \t - \t - \t - \t - \t - \t %0.3f \t -",bs.c_str(),
							1+e_eVV/eVV).Data()<<endl;    
	//shape dependent systematics
	myfile<<"CMSeffb  \t shape \t 1. \t 1. \t 1. \t 0. \t 1. \t 1. \t 1. \t 1."<<endl;  
	// myfile<<TString::Format("CMS_vhbb_stats_Zemu%s \t shape  \t 1. \t - \t - \t - \t - \t - \t - \t -",bs.c_str()).Data()<<endl;  
	//myfile<<TString::Format("CMS_vhbb_stats_Zemu%s \t shape  \t - \t - \t 1. \t - \t - \t - \t - \t -",bs.c_str()).Data()<<endl;
	//myfile<<TString::Format("CMS_vhbb_stats_Zemu%s \t shape  \t -  \t - \t - \t 1. \t - \t - \t - \t -",bs.c_str()).Data()<<endl;
	//myfile<<TString::Format("CMS_vhbb_stats_Zemu%s \t shape  \t - \t - \t - \t - \t 1. \t - \t - \t -", bs.c_str()).Data()<<endl; 
	//myfile<<TString::Format("CMS_vhbb_stats_Zemu%s \t shape  \t - \t - \t - \t - \t - \t - \t 1. \t -",bs.c_str()).Data()<<endl;
	myfile.close();
	
}

void makeShapeHistos(string shapefile, TH1F* hVH, TH1F* hWj, TH1F* hZjLF, TH1F* hZjHF, 
					 TH1F* hTT, TH1F* hs_Top, TH1F* hVV, TH1F* hQCD, TH1F* hdata)
{
	TFile* file=new TFile(TString::Format("%s/%s",dir.c_str(),shapefile.c_str()).Data(),"RECREATE");
	// myfile<<"process \t VH \t Wj \t ZjLF \t ZjHF \t TT \t s_Top \t VV \t QCD"<<endl;
	//names have to match process names
	hVH->SetName("VH");   hVH->Write();
	hWj->SetName("Wj");  hWj->Write();
	hZjLF->SetName("ZjLF");  hZjLF->Write();
	hZjHF->SetName("ZjHF");  hZjHF->Write();
	hTT->SetName("TT");   hTT->Write();
	hs_Top->SetName("s_Top");   hs_Top->Write();
	hVV->SetName("VV");   hVV->Write();
	hQCD->SetName("QCD");   hQCD->Write();
	hdata->SetName("data_obs");   hdata->Write();
	file->Close();
}

void  AddSysHistos(string shapefile, TH1F* hVH_up, TH1F* hVH_down, 
				   TH1F* hWj_up, TH1F* hWj_down, TH1F* hZjLF_up,TH1F* hZjLF_down,
				   TH1F* hZjHF_up, TH1F* hZjHF_down, TH1F* hTT_up, 
				   TH1F* hTT_down, TH1F* hs_Top_up, TH1F* hs_Top_down,  
				   TH1F* hVV_up, TH1F* hVV_down, TH1F* hQCD_up,
				   TH1F* hQCD_down,  string sysUp="CMSeffbUp",
				   string sysDown="CMSeffbDown")
{
	cout << "made it into AddSysHistos" << endl;
	
	TFile* file=new TFile(TString::Format("%s/%s",dir.c_str(),
										  shapefile.c_str()).Data(),"UPDATE");
	// myfile<<"process \t VH \t Wj\t ZjLF \t ZjHF \t TT \t s_Top \t VV \t QCD"<<endl;
	//names have to match process names
	vector <string> process;
	process.push_back("VH");
	process.push_back("Wj");
	process.push_back("ZjLF");  
	process.push_back("ZjHF");  
	process.push_back("TT");   
	process.push_back("s_Top"); 
	process.push_back("VV");   
	process.push_back("QCD");  
	
	vector <TH1F*> hup;
	hup.push_back(hVH_up);
	hup.push_back(hWj_up);
	hup.push_back(hZjLF_up);
	hup.push_back(hZjHF_up); 
	hup.push_back(hTT_up); 
	hup.push_back(hs_Top_up); 
	hup.push_back(hVV_up); 
	hup.push_back(hQCD_up);
	vector <TH1F*> hdown;
	hdown.push_back(hVH_down);
	hdown.push_back(hWj_down);
	hdown.push_back(hZjLF_down);
	hdown.push_back(hZjHF_down); 
	hdown.push_back(hTT_down); 
	hdown.push_back(hs_Top_down); 
	hdown.push_back(hVV_down); 
	hdown.push_back(hQCD_down);
	cout << "array of histograms made" << endl;
	
	for(int i=0;i<7;i++){
		hup[i]->SetName(TString::Format("%s_%s",process[i].c_str(),sysUp.c_str())
						.Data()); 
		hup[i]->Write();
		hdown[i]->SetName(TString::Format("%s_%s",process[i].c_str(),sysDown.c_str())
						  .Data()); 
		hdown[i]->Write();
	}
	
	cout << "array of histograms written" << endl;
	
	file->Close();
}

bool Preselect(double csv1, double dphiZMET, double MZ, double CHFb0, double Mjj, double SVDZmass, double ptz)
{
	bool passes=false;
	
    if(csv1>0.244 && fabs(dphiZMET)<1.25  && MZ< 85 && MZ > 10 && ptz> 20//&& CHFb0 > 0.32  && Mjj<250. && Angleemu > 0.15
	   && (SVDZmass<0 ||(SVDZmass >30 && SVDZmass<100)) && Mjj>80 && Mjj<150 ) passes=true;
	
	return passes;
}

void FillBtagSys(float bdt, float w, int evf, TH1F* hup, TH1F* hdown)
{
	//b jets  
	if(abs(evf)==5){
		hup->Fill(bdt, w*(1.+bjetsys));  
		hdown->Fill(bdt, w*(1.-bjetsys));  
	}
	// c jets 
	if(abs(evf)==4){
		hup->Fill(bdt, w*(1.+cjetsys));  
		hdown->Fill(bdt, w*(1.-cjetsys));  
	}
	// light jets
	if(abs(evf)!=4 && abs(evf)!=5){
		hup->Fill(bdt, w*(1.+ljetsys));  
		hdown->Fill(bdt, w*(1.-ljetsys));  
	}
}

void FillBtagSys(float bdt, float w, int evf, 
				 TH1F* h0lup, TH1F* h0ldown, 
				 TH1F* h0hup, TH1F* h0hdown)
{
	//b jets  
	if(abs(evf)==5){
		h0hup->Fill(bdt, w*(1.+bjetsys));  
		h0hdown->Fill(bdt, w*(1.-bjetsys));  
	}
	
	// c jets 
	if(abs(evf)==4){
		h0lup->Fill(bdt, w*(1.+cjetsys));  
		h0ldown->Fill(bdt, w*(1.-cjetsys));  
	}
	// light jets
	if(abs(evf)!=4 && abs(evf)!=5){
		h0lup->Fill(bdt, w*(1.+ljetsys));  
		h0ldown->Fill(bdt, w*(1.-ljetsys));  
	}
}

void Make1dBinning(TTree* data, float scale, int boost, bool isMatt, int n2dbins=-1)
{
	
	TH1F* h=new TH1F("","",1000,0,250);
	h->Sumw2();
	
	float exp;
	float Hmass, Emumass;
	float Hpt, Zpt;
	float CSV0, CSV1;
	float jetCHF0, jetCHF1;
	float DphiJJ, DetaJJ, delRjj;
  	float Hphi, Zphi, DeltaPhiHV, Ht, MET, EventPt, PtbalZH, Mte, Mtmu, DphiZMET;
	float lep0pt, lep1pt, lep_pfCombRelIso0, lep_pfCombRelIso1, AngleEMU, delRemu;
	float EtaStandDev, RMS_eta, UnweightedEta;
	float Centrality, EvntShpAplanarity, EvntShpSphericity, EvntShpIsotropy, EvntShpCircularity;
	float Trigweight, B2011PUweight, A2011PUweight, btag2CSF, BDTvalue;
	float ZmassSVD;
	int nJets, eventFlavor, naJets, nSV;
	
	data->SetBranchAddress( "nJets", &nJets );
	data->SetBranchAddress( "naJets", &naJets );
	data->SetBranchAddress( "nSV", &nSV );	
	data->SetBranchAddress( "Hmass", &Hmass );
    data->SetBranchAddress( "Emumass", &Emumass );
    data->SetBranchAddress( "Hpt", &Hpt );
	data->SetBranchAddress( "Zpt", &Zpt );
	data->SetBranchAddress( "CSV0", &CSV0 );
	data->SetBranchAddress( "CSV1", &CSV1 );
	data->SetBranchAddress( "jetCHF0", &jetCHF0 );
	data->SetBranchAddress( "jetCHF1", &jetCHF1 );
	data->SetBranchAddress( "DetaJJ", &DetaJJ );
	data->SetBranchAddress( "DphiJJ", &DphiJJ );
	data->SetBranchAddress( "delRjj", &delRjj );
	data->SetBranchAddress( "Hphi", &Hphi );
	data->SetBranchAddress( "Zphi", &Zphi );
	data->SetBranchAddress( "DeltaPhiHV", &DeltaPhiHV );
	data->SetBranchAddress( "Ht", &Ht );
	data->SetBranchAddress( "MET", &MET );
	data->SetBranchAddress( "EventPt", &EventPt );
	data->SetBranchAddress( "PtbalZH", &PtbalZH );
	data->SetBranchAddress( "Mte", &Mte );
	data->SetBranchAddress( "Mtmu", &Mtmu );
	data->SetBranchAddress( "DphiZMET", &DphiZMET );
	data->SetBranchAddress( "Zphi", &Zphi );
	data->SetBranchAddress( "Hphi", &Hphi );
	data->SetBranchAddress( "lep_pfCombRelIso0", &lep_pfCombRelIso0 );
	data->SetBranchAddress( "lep_pfCombRelIso1", &lep_pfCombRelIso1 );
	data->SetBranchAddress( "lep0pt", &lep0pt );
	data->SetBranchAddress( "lep1pt", &lep1pt );
	data->SetBranchAddress( "AngleEMU", &AngleEMU );
	data->SetBranchAddress( "delRemu", &delRemu );
	data->SetBranchAddress( "RMS_eta", &RMS_eta );
	data->SetBranchAddress( "UnweightedEta", &UnweightedEta );
	data->SetBranchAddress( "EtaStandDev", &EtaStandDev );
	data->SetBranchAddress( "Centrality", &Centrality );
	data->SetBranchAddress( "EvntShpAplanarity", &EvntShpAplanarity );
	data->SetBranchAddress( "EvntShpSphericity", &EvntShpSphericity );
	data->SetBranchAddress( "EvntShpIsotropy", &EvntShpIsotropy );
	data->SetBranchAddress( "EvntShpCircularity", &EvntShpCircularity );
	data->SetBranchAddress( "eventFlavor", &eventFlavor );
	data->SetBranchAddress( "Trigweight", &Trigweight );
	data->SetBranchAddress( "B2011PUweight", &B2011PUweight );
	data->SetBranchAddress( "A2011PUweight", &A2011PUweight );
	data->SetBranchAddress( "btag2CSF", &btag2CSF );
	data->SetBranchAddress( "BDTvalue", &BDTvalue );
	data->SetBranchAddress( "ZmassSVD", &ZmassSVD );
	
	
	
	for(int i=0;i<data->GetEntries();i++){
		data->GetEntry(i);
		double bdt = BDTvalue;
		double wT_=Trigweight;
		//	  double wPU_=(wA2011PUweight*lA+B2011PUweight*lB)/(lA+lB);
		double wPU_= 1.0;
		wPU_= (lA*A2011PUweight + lB*B2011PUweight)/lumi;
		//double wb_=btag2CSF;
		double wb_=1.0;
		wT_ = 1.0;
		
		if(Preselect(CSV0, CSV1, Emumass, jetCHF0, Hmass, ZmassSVD , Hpt)){
			if(PassBoost(Hpt, Zpt, boost, Hmass)){
				h->Fill(Hmass,wT_*wPU_*wb_);
			}}//preselected and in boosted region
	} //end of loop
	cout<<"Make1dBinning myhisto before sclaing"<<myhisto->Integral(0,-1)<<endl;
	
	h->Scale(scale);
	h->SetLineWidth(2);
	
	exp=h->Integral(0,-1);
	//now that I have the integral of the stuff passing the selection
	float quantile= exp/(float)N;
	cout<<"total expected signal is "<<exp<<" quantile is  "<<quantile<<endl;
	
	int leftedge=0;
	int ncycle=0;
	for(int ibin=0;ibin<h->GetNbinsX();ibin++){
		//as soon as it gets more than quantile quantity, then accept it as quantile boundary 
		if( h->Integral(leftedge,ibin)>quantile){
			mass2dcutl[ncycle]=h->GetXaxis()->GetBinLowEdge(leftedge+1);
			leftedge=ibin;      
			mass2dcuth[ncycle]=h->GetXaxis()->GetBinLowEdge(leftedge+1);
			bdt2dcutl[ncycle]=-1.5;
			bdt2dcuth[ncycle]=1.5;      
			ncycle++;
		}
		if( ncycle==N-1){
			mass2dcutl[ncycle]=h->GetXaxis()->GetBinLowEdge(leftedge+1);
			mass2dcuth[ncycle]=h->GetXaxis()->GetBinLowEdge(h->GetNbinsX()+1);
			bdt2dcutl[ncycle]=-1.5;
			bdt2dcuth[ncycle]=1.5;      
			ibin=h->GetNbinsX();
		}
	}
	
	for(int i=0;i<N;i++) 
		cout<<"check if stuff has been set OK: " << i<< "   "<<mass2dcutl[i]<<"   "<< mass2dcuth[i]<<"   "<< 
		h->Integral(h->FindBin(mass2dcutl[i]), h->FindBin(mass2dcuth[i]))<<endl;
}

void HistoFiller(TTree* data, float scale, TH1F*& h, string var, int boost, 
				 float& exp, float& e_exp, int ivar, TH1F*& hup, TH1F*& hdown, TH1F*& hstatup, TH1F*& hstatdown,
				 bool isMatt, int n2dbins=-1)
{
	float histoFillValue = -99.99;
	h->Sumw2();
	
	if(ivar==kbdt){
		hup->Sumw2();
		hdown->Sumw2();
		hstatup->Sumw2();
		hstatdown->Sumw2();
	}
	
 	float Hmass, Emumass;
	float Hpt, Zpt;
	float CSV0, CSV1;
	float jetCHF0, jetCHF1;
	float DphiJJ, DetaJJ, delRjj;
  	float Hphi, Zphi, DeltaPhiHV, Ht, MET, EventPt, PtbalZH, Mte, Mtmu, DphiZMET;
	float lep0pt, lep1pt, lep_pfCombRelIso0, lep_pfCombRelIso1, AngleEMU, delRemu;
	float EtaStandDev, RMS_eta, UnweightedEta;
	float Centrality, EvntShpAplanarity, EvntShpSphericity, EvntShpIsotropy, EvntShpCircularity;
	float Trigweight, B2011PUweight, A2011PUweight, btag2CSF, BDTvalue;
	float ZmassSVD;
	int nJets, eventFlavor, naJets, nSV;
	
	data->SetBranchAddress( "nJets", &nJets );
	data->SetBranchAddress( "naJets", &naJets );
	data->SetBranchAddress( "nSV", &nSV );	
	data->SetBranchAddress( "Hmass", &Hmass );
    data->SetBranchAddress( "Emumass", &Emumass );
    data->SetBranchAddress( "Hpt", &Hpt );
	data->SetBranchAddress( "Zpt", &Zpt );
	data->SetBranchAddress( "CSV0", &CSV0 );
	data->SetBranchAddress( "CSV1", &CSV1 );
	data->SetBranchAddress( "jetCHF0", &jetCHF0 );
	data->SetBranchAddress( "jetCHF1", &jetCHF1 );
	data->SetBranchAddress( "DetaJJ", &DetaJJ );
	data->SetBranchAddress( "DphiJJ", &DphiJJ );
	data->SetBranchAddress( "delRjj", &delRjj );
	data->SetBranchAddress( "Hphi", &Hphi );
	data->SetBranchAddress( "Zphi", &Zphi );
	data->SetBranchAddress( "DeltaPhiHV", &DeltaPhiHV );
	data->SetBranchAddress( "Ht", &Ht );
	data->SetBranchAddress( "MET", &MET );
	data->SetBranchAddress( "EventPt", &EventPt );
	data->SetBranchAddress( "PtbalZH", &PtbalZH );
	data->SetBranchAddress( "Mte", &Mte );
	data->SetBranchAddress( "Mtmu", &Mtmu );
	data->SetBranchAddress( "DphiZMET", &DphiZMET );
	data->SetBranchAddress( "Zphi", &Zphi );
	data->SetBranchAddress( "Hphi", &Hphi );
	data->SetBranchAddress( "lep_pfCombRelIso0", &lep_pfCombRelIso0 );
	data->SetBranchAddress( "lep_pfCombRelIso1", &lep_pfCombRelIso1 );
	data->SetBranchAddress( "lep0pt", &lep0pt );
	data->SetBranchAddress( "lep1pt", &lep1pt );
	data->SetBranchAddress( "AngleEMU", &AngleEMU );
	data->SetBranchAddress( "delRemu", &delRemu );
	data->SetBranchAddress( "RMS_eta", &RMS_eta );
	data->SetBranchAddress( "UnweightedEta", &UnweightedEta );
	data->SetBranchAddress( "EtaStandDev", &EtaStandDev );
	data->SetBranchAddress( "Centrality", &Centrality );
	data->SetBranchAddress( "EvntShpAplanarity", &EvntShpAplanarity );
	data->SetBranchAddress( "EvntShpSphericity", &EvntShpSphericity );
	data->SetBranchAddress( "EvntShpIsotropy", &EvntShpIsotropy );
	data->SetBranchAddress( "EvntShpCircularity", &EvntShpCircularity );
	data->SetBranchAddress( "eventFlavor", &eventFlavor );
	data->SetBranchAddress( "Trigweight", &Trigweight );
	data->SetBranchAddress( "B2011PUweight", &B2011PUweight );
	data->SetBranchAddress( "A2011PUweight", &A2011PUweight );
	data->SetBranchAddress( "btag2CSF", &btag2CSF );
	data->SetBranchAddress( "BDTvalue", &BDTvalue );
	data->SetBranchAddress( "ZmassSVD", &ZmassSVD );
	
	
	for(int i=0;i<data->GetEntries();i++){
		data->GetEntry(i);
		histoFillValue = -99.99;
		
		if (ivar==kl1pt) histoFillValue = lep0pt;
		if (ivar==kl2pt) histoFillValue =  lep1pt; //1
		if (ivar==kpTjj) histoFillValue = Zpt;  //2
		if (ivar==kpTZ) histoFillValue =  Zpt; //3
		if (ivar==kcsv1) histoFillValue = CSV0;  //4
		if (ivar==kcsv2) histoFillValue = CSV1;  //5 
		if (ivar==kDphi) histoFillValue = DeltaPhiHV; //6
		if (ivar==kNaj) histoFillValue =  naJets; //7 
		if (ivar==kmet) histoFillValue = MET;  //8
		if (ivar==kDeta) histoFillValue = DetaJJ;   //9
		if (ivar==kHt) histoFillValue =  Ht; //10
		if (ivar==kjchf1) histoFillValue = jetCHF0;  //11
		if (ivar==kjchf2) histoFillValue = jetCHF1;  //12
		if (ivar==kangleEMU) histoFillValue =  AngleEMU; //13
		if (ivar==kueta) histoFillValue = UnweightedEta;  //14 
		if (ivar==kcentrality) histoFillValue = Centrality;  //15
		if (ivar==kM) histoFillValue =  Hmass;  //16
		if (ivar==kDptHZ) histoFillValue = PtbalZH;   //17
		if (ivar==hdrjj) histoFillValue =  delRjj;  //18
		if (ivar==kbdt) histoFillValue = BDTvalue;   //19
		if (ivar==kdeltaPhiZMET) histoFillValue = DphiZMET;   //20
		if (ivar==kAplanarity) histoFillValue = EvntShpAplanarity;   //21
		if (ivar==kSphericity) histoFillValue = EvntShpSphericity;  //22
		if (ivar==kCircularity) histoFillValue = EvntShpCircularity; //23
		if (ivar==kIsotropy) histoFillValue = EvntShpIsotropy;  //24
		if (ivar==kMte) histoFillValue = Mte;   //25
		if (ivar==kMtmu) histoFillValue = Mtmu;   //26
		
		
		double bdt = BDTvalue;
		double wT_=Trigweight;
		//	  double wPU_=(A2011PUweight*lA+B2011PUweight*lB)/(lA+lB);
		double wPU_= 1.0;
		wPU_= (lA*A2011PUweight + lB*B2011PUweight)/lumi;
		//double wb_=btag2CSF;
		double wb_=1.0;
		wT_ = 1.0;		
		if(Preselect(CSV0, CSV1, Emumass, jetCHF0, Hmass, ZmassSVD , Hpt)){
			if(PassBoost(Hpt, Zpt, boost, Hmass)){
				if( BDTvalue>bdtcut){
					/*  if( NewselectSignal(CSV0, CSV1, DeltaPhiHV,
					 Emumass, Hmass, 
					 MET, Ht, AngleEMU,
					 UnweightedEta))
					 */ 
					if(ivar!=kbdt) h->Fill(histoFillValue,wT_*wPU_*wb_);
					else h->Fill(BDTvalue,wT_*wPU_*wb_);
					
					if(ivar==kbdt) 
						FillBtagSys(BDTvalue, wT_*wPU_*wb_, eventFlavor, hup, hdown);
				}   
			}}//preselected and in boosted region 
	} //end of loop
	
	h->Scale(scale);
	
	h->SetLineWidth(2);
	if(ivar==kbdt){
		hup->Scale(scale);
		hdown->Scale(scale);
	}
	
	//stats error as an histogram, should work if scaling is done correctly
	if(ivar==kbdt) {
		for(int i=0;i<h->GetNbinsX()+1;i++){
			hstatup->SetBinContent(i, h->GetBinContent(i) + sqrt(h->GetBinContent(i)/scale)*scale);
			//h->GetBinError(i));
			hstatdown->SetBinContent(i, h->GetBinContent(i) - sqrt(h->GetBinContent(i)/scale)*scale);
			//-h->GetBinError(i));
		}
	}
	
	exp=h->Integral(0,-1);
	e_exp= sqrt(h->Integral(0,-1)/scale)*scale;
	
	/*
	 //try to remove 
	 if(ivar==kbdt){
	 if(exp<0.01 ){
	 exp=0;
	 e_exp=0;
	 for(int i=0;i<h->GetNbinsX()+1;i++){
	 hstatup->SetBinContent(i,0);
	 hstatdown->SetBinContent(i,0);
	 hup->SetBinContent(i,0);
	 hdown->SetBinContent(i,0);
	 }
	 }
	 }
	 */
	
}
//histo filler for data
void HistoFiller(TTree* data, TH1F*& h, string var, int boost, int ivar, bool isMatt, int n2dbins=-1) 
{
	float histoFillValue = -99.99;
	h->Sumw2();
	
	float Hmass, Emumass;
	float Hpt, Zpt;
	float CSV0, CSV1;
	float jetCHF0, jetCHF1;
	float DphiJJ, DetaJJ, delRjj;
  	float Hphi, Zphi, DeltaPhiHV, Ht, MET, EventPt, PtbalZH, Mte, Mtmu, DphiZMET;
	float lep0pt, lep1pt, lep_pfCombRelIso0, lep_pfCombRelIso1, AngleEMU, delRemu;
	float EtaStandDev, RMS_eta, UnweightedEta;
	float Centrality, EvntShpAplanarity, EvntShpSphericity, EvntShpIsotropy, EvntShpCircularity;
	float Trigweight, B2011PUweight, A2011PUweight, btag2CSF, BDTvalue;
	float ZmassSVD;
	int nJets, eventFlavor, naJets, nSV;
	
	data->SetBranchAddress( "nJets", &nJets );
	data->SetBranchAddress( "naJets", &naJets );
	data->SetBranchAddress( "nSV", &nSV );	
	data->SetBranchAddress( "Hmass", &Hmass );
    data->SetBranchAddress( "Emumass", &Emumass );
    data->SetBranchAddress( "Hpt", &Hpt );
	data->SetBranchAddress( "Zpt", &Zpt );
	data->SetBranchAddress( "CSV0", &CSV0 );
	data->SetBranchAddress( "CSV1", &CSV1 );
	data->SetBranchAddress( "jetCHF0", &jetCHF0 );
	data->SetBranchAddress( "jetCHF1", &jetCHF1 );
	data->SetBranchAddress( "DetaJJ", &DetaJJ );
	data->SetBranchAddress( "DphiJJ", &DphiJJ );
	data->SetBranchAddress( "delRjj", &delRjj );
	data->SetBranchAddress( "Hphi", &Hphi );
	data->SetBranchAddress( "Zphi", &Zphi );
	data->SetBranchAddress( "DeltaPhiHV", &DeltaPhiHV );
	data->SetBranchAddress( "Ht", &Ht );
	data->SetBranchAddress( "MET", &MET );
	data->SetBranchAddress( "EventPt", &EventPt );
	data->SetBranchAddress( "PtbalZH", &PtbalZH );
	data->SetBranchAddress( "Mte", &Mte );
	data->SetBranchAddress( "Mtmu", &Mtmu );
	data->SetBranchAddress( "DphiZMET", &DphiZMET );
	data->SetBranchAddress( "Zphi", &Zphi );
	data->SetBranchAddress( "Hphi", &Hphi );
	data->SetBranchAddress( "lep_pfCombRelIso0", &lep_pfCombRelIso0 );
	data->SetBranchAddress( "lep_pfCombRelIso1", &lep_pfCombRelIso1 );
	data->SetBranchAddress( "lep0pt", &lep0pt );
	data->SetBranchAddress( "lep1pt", &lep1pt );
	data->SetBranchAddress( "AngleEMU", &AngleEMU );
	data->SetBranchAddress( "delRemu", &delRemu );
	data->SetBranchAddress( "RMS_eta", &RMS_eta );
	data->SetBranchAddress( "UnweightedEta", &UnweightedEta );
	data->SetBranchAddress( "EtaStandDev", &EtaStandDev );
	data->SetBranchAddress( "Centrality", &Centrality );
	data->SetBranchAddress( "EvntShpAplanarity", &EvntShpAplanarity );
	data->SetBranchAddress( "EvntShpSphericity", &EvntShpSphericity );
	data->SetBranchAddress( "EvntShpIsotropy", &EvntShpIsotropy );
	data->SetBranchAddress( "EvntShpCircularity", &EvntShpCircularity );
	data->SetBranchAddress( "eventFlavor", &eventFlavor );
	data->SetBranchAddress( "Trigweight", &Trigweight );
	data->SetBranchAddress( "B2011PUweight", &B2011PUweight );
	data->SetBranchAddress( "A2011PUweight", &A2011PUweight );
	data->SetBranchAddress( "btag2CSF", &btag2CSF );
	data->SetBranchAddress( "BDTvalue", &BDTvalue );
	data->SetBranchAddress( "ZmassSVD", &ZmassSVD );
	
	
	for(int i=0;i<data->GetEntries();i++){
		data->GetEntry(i);
		
		histoFillValue = -99.99;
		
		if (ivar==kl1pt) histoFillValue = lep0pt;
		if (ivar==kl2pt) histoFillValue =  lep1pt; //1
		if (ivar==kpTjj) histoFillValue = Zpt;  //2
		if (ivar==kpTZ) histoFillValue =  Zpt; //3
		if (ivar==kcsv1) histoFillValue = CSV0;  //4
		if (ivar==kcsv2) histoFillValue = CSV1;  //5 
		if (ivar==kDphi) histoFillValue = DeltaPhiHV; //6
		if (ivar==kNaj) histoFillValue =  naJets; //7 
		if (ivar==kmet) histoFillValue = MET;  //8
		if (ivar==kDeta) histoFillValue = DetaJJ;   //9
		if (ivar==kHt) histoFillValue =  Ht; //10
		if (ivar==kjchf1) histoFillValue = jetCHF0;  //11
		if (ivar==kjchf2) histoFillValue = jetCHF1;  //12
		if (ivar==kangleEMU) histoFillValue =  AngleEMU; //13
		if (ivar==kueta) histoFillValue = UnweightedEta;  //14 
		if (ivar==kcentrality) histoFillValue = Centrality;  //15
		if (ivar==kM) histoFillValue =  Hmass;  //16
		if (ivar==kDptHZ) histoFillValue = PtbalZH;   //17
		if (ivar==hdrjj) histoFillValue =  delRjj;  //18
		if (ivar==kbdt) histoFillValue = BDTvalue;   //19
		if (ivar==kdeltaPhiZMET) histoFillValue = DphiZMET;   //20
		if (ivar==kAplanarity) histoFillValue = EvntShpAplanarity;   //21
		if (ivar==kSphericity) histoFillValue = EvntShpSphericity;  //22
		if (ivar==kCircularity) histoFillValue = EvntShpCircularity; //23
		if (ivar==kIsotropy) histoFillValue = EvntShpIsotropy;  //24
		if (ivar==kMte) histoFillValue = Mte;   //25
		if (ivar==kMtmu) histoFillValue = Mtmu;   //26
		
		
		if(Preselect(CSV0, CSV1, Emumass, jetCHF0, Hmass, ZmassSVD , Hpt)){
			if(PassBoost(Hpt, Zpt, boost, Hmass)){
				if( BDTvalue>bdtcut){
					/*  if( NewselectSignal(CSV0, CSV1, DeltaPhiHV,
					 Emumass, Hmass, 
					 MET, Ht, AngleEMU,
					 UnweightedEta))
					 */ 
					if(Hmass<95 || Hmass>125)
						if(ivar!=kbdt) h->Fill(histoFillValue);
						else  {
							if(BDTvalue<-0.25)
								h->Fill(BDTvalue);
						}
				}
			}}//preselected and in boosted region 
	} //end of loop
}


void HistoFillerLH(TTree* dataL, float wL, TTree* dataH,  float wH, TH1F*& hL, TH1F*& hH, string var, int boost,
				   float& expZL, float& expZH, float& e_expZL, float& e_expZH, 
				   int ivar, TH1F*& hL_up, TH1F*& hL_down, TH1F*& hH_up, 
				   TH1F*& hH_down, TH1F*& hstatLup, TH1F*& hstatLdown, TH1F*& hstatHup, TH1F*&  hstatHdown,
				   bool isMatt, int n2dbins=-1)
{
	float histoFillValue = -99.99;
	
	hL->Sumw2();
	hH->Sumw2();
	if(ivar==kbdt){
		hL_up->Sumw2();
		hH_up->Sumw2();
		hL_down->Sumw2();
		hH_down->Sumw2();
		hstatLup->Sumw2();	
		hstatLdown->Sumw2();	
		hstatHup->Sumw2();	
		hstatHdown->Sumw2();	
	}
	//n candidates 
	TH1F* hbdt_data0L=(TH1F*)hL->Clone("");
	TH1F* hbdt_data0H=(TH1F*)hL->Clone("");
	TH1F* hbdt_data4L=(TH1F*)hL->Clone("");
	TH1F* hbdt_data4H=(TH1F*)hL->Clone("");
	hbdt_data0L->Sumw2();
	hbdt_data0H->Sumw2();
	hbdt_data4L->Sumw2();
	hbdt_data4H->Sumw2();
	
	TH1F* hbdt_data0L_up=(TH1F*)hL_up->Clone("");
	TH1F* hbdt_data0H_up=(TH1F*)hL_up->Clone("");
	TH1F* hbdt_data0L_down=(TH1F*)hL_up->Clone("");
	TH1F* hbdt_data0H_down=(TH1F*)hL_up->Clone("");
	TH1F* hbdt_data4L_up=(TH1F*)hL_up->Clone("");
	TH1F* hbdt_data4H_up=(TH1F*)hL_up->Clone("");
	TH1F* hbdt_data4L_down=(TH1F*)hL_up->Clone("");
	TH1F* hbdt_data4H_down=(TH1F*)hL_up->Clone("");
	hbdt_data0L_up->Sumw2();
	hbdt_data0H_up->Sumw2();
	hbdt_data0L_down->Sumw2();
	hbdt_data0H_down->Sumw2();
	hbdt_data4L_up->Sumw2();
	hbdt_data4H_up->Sumw2();
	hbdt_data4L_down->Sumw2();
	hbdt_data4H_down->Sumw2();
	
	float Hmass, Emumass;
	float Hpt, Zpt;
	float CSV0, CSV1;
	float jetCHF0, jetCHF1;
	float DphiJJ, DetaJJ, delRjj;
  	float Hphi, Zphi, DeltaPhiHV, Ht, MET, EventPt, PtbalZH, Mte, Mtmu, DphiZMET;
	float lep0pt, lep1pt, lep_pfCombRelIso0, lep_pfCombRelIso1, AngleEMU, delRemu;
	float EtaStandDev, RMS_eta, UnweightedEta;
	float Centrality, EvntShpAplanarity, EvntShpSphericity, EvntShpIsotropy, EvntShpCircularity;
	float Trigweight, B2011PUweight, A2011PUweight, btag2CSF, BDTvalue;
	float ZmassSVD;
	int nJets, eventFlavor, naJets, nSV;
	
	dataL->SetBranchAddress( "nJets", &nJets );
	dataL->SetBranchAddress( "naJets", &naJets );
	dataL->SetBranchAddress( "nSV", &nSV );	
	dataL->SetBranchAddress( "Hmass", &Hmass );
    dataL->SetBranchAddress( "Emumass", &Emumass );
    dataL->SetBranchAddress( "Hpt", &Hpt );
	dataL->SetBranchAddress( "Zpt", &Zpt );
	dataL->SetBranchAddress( "CSV0", &CSV0 );
	dataL->SetBranchAddress( "CSV1", &CSV1 );
	dataL->SetBranchAddress( "jetCHF0", &jetCHF0 );
	dataL->SetBranchAddress( "jetCHF1", &jetCHF1 );
	dataL->SetBranchAddress( "DetaJJ", &DetaJJ );
	dataL->SetBranchAddress( "DphiJJ", &DphiJJ );
	dataL->SetBranchAddress( "delRjj", &delRjj );
	dataL->SetBranchAddress( "Hphi", &Hphi );
	dataL->SetBranchAddress( "Zphi", &Zphi );
	dataL->SetBranchAddress( "DeltaPhiHV", &DeltaPhiHV );
	dataL->SetBranchAddress( "Ht", &Ht );
	dataL->SetBranchAddress( "MET", &MET );
	dataL->SetBranchAddress( "EventPt", &EventPt );
	dataL->SetBranchAddress( "PtbalZH", &PtbalZH );
	dataL->SetBranchAddress( "Mte", &Mte );
	dataL->SetBranchAddress( "Mtmu", &Mtmu );
	dataL->SetBranchAddress( "DphiZMET", &DphiZMET );
	dataL->SetBranchAddress( "Zphi", &Zphi );
	dataL->SetBranchAddress( "Hphi", &Hphi );
	dataL->SetBranchAddress( "lep_pfCombRelIso0", &lep_pfCombRelIso0 );
	dataL->SetBranchAddress( "lep_pfCombRelIso1", &lep_pfCombRelIso1 );
	dataL->SetBranchAddress( "lep0pt", &lep0pt );
	dataL->SetBranchAddress( "lep1pt", &lep1pt );
	dataL->SetBranchAddress( "AngleEMU", &AngleEMU );
	dataL->SetBranchAddress( "delRemu", &delRemu );
	dataL->SetBranchAddress( "RMS_eta", &RMS_eta );
	dataL->SetBranchAddress( "UnweightedEta", &UnweightedEta );
	dataL->SetBranchAddress( "EtaStandDev", &EtaStandDev );
	dataL->SetBranchAddress( "Centrality", &Centrality );
	dataL->SetBranchAddress( "EvntShpAplanarity", &EvntShpAplanarity );
	dataL->SetBranchAddress( "EvntShpSphericity", &EvntShpSphericity );
	dataL->SetBranchAddress( "EvntShpIsotropy", &EvntShpIsotropy );
	dataL->SetBranchAddress( "EvntShpCircularity", &EvntShpCircularity );
	dataL->SetBranchAddress( "eventFlavor", &eventFlavor );
	dataL->SetBranchAddress( "Trigweight", &Trigweight );
	dataL->SetBranchAddress( "B2011PUweight", &B2011PUweight );
	dataL->SetBranchAddress( "A2011PUweight", &A2011PUweight );
	dataL->SetBranchAddress( "btag2CSF", &btag2CSF );
	dataL->SetBranchAddress( "BDTvalue", &BDTvalue );
	dataL->SetBranchAddress( "ZmassSVD", &ZmassSVD );
	
	for(int i=0;i<dataL->GetEntries();i++){
		dataL->GetEntry(i);    
		
		if (ivar==kl1pt) histoFillValue = lep0pt;
		if (ivar==kl2pt) histoFillValue =  lep1pt; //1
		if (ivar==kpTjj) histoFillValue = Zpt;  //2
		if (ivar==kpTZ) histoFillValue =  Zpt; //3
		if (ivar==kcsv1) histoFillValue = CSV0;  //4
		if (ivar==kcsv2) histoFillValue = CSV1;  //5 
		if (ivar==kDphi) histoFillValue = DeltaPhiHV; //6
		if (ivar==kNaj) histoFillValue =  naJets; //7 
		if (ivar==kmet) histoFillValue = MET;  //8
		if (ivar==kDeta) histoFillValue = DetaJJ;   //9
		if (ivar==kHt) histoFillValue =  Ht; //10
		if (ivar==kjchf1) histoFillValue = jetCHF0;  //11
		if (ivar==kjchf2) histoFillValue = jetCHF1;  //12
		if (ivar==kangleEMU) histoFillValue =  AngleEMU; //13
		if (ivar==kueta) histoFillValue = UnweightedEta;  //14 
		if (ivar==kcentrality) histoFillValue = Centrality;  //15
		if (ivar==kM) histoFillValue =  Hmass;  //16
		if (ivar==kDptHZ) histoFillValue = PtbalZH;   //17
		if (ivar==hdrjj) histoFillValue =  delRjj;  //18
		if (ivar==kbdt) histoFillValue = BDTvalue;   //19
		if (ivar==kdeltaPhiZMET) histoFillValue = DphiZMET;   //20
		if (ivar==kAplanarity) histoFillValue = EvntShpAplanarity;   //21
		if (ivar==kSphericity) histoFillValue = EvntShpSphericity;  //22
		if (ivar==kCircularity) histoFillValue = EvntShpCircularity; //23
		if (ivar==kIsotropy) histoFillValue = EvntShpIsotropy;  //24
		if (ivar==kMte) histoFillValue = Mte;   //25
		if (ivar==kMtmu) histoFillValue = Mtmu;   //26		
		
		
		double wT_=Trigweight;
		double wPU_= 1.0;
		wPU_= (lA*A2011PUweight + lB*B2011PUweight)/lumi;
		double wb_=1.0;
		wT_ = 1.0;
		int ef= eventFlavor;   
		
		if(Preselect(CSV0, CSV1, Emumass, jetCHF0, Hmass, ZmassSVD , Hpt)){
			if(PassBoost(Hpt, Zpt, boost, Hmass)) {
				if( BDTvalue>bdtcut){
					if (ef != 5){
						/*	if(NewselectSignal(CSV0, CSV1, DeltaPhiHV,
						 Emumass, Hmass, 
						 MET, Ht, AngleEMU,
						 UnweightedEta))
						 */
						if(ivar!=kbdt) hbdt_data0L->Fill(histoFillValue,wT_*wPU_*wb_);
						else hbdt_data0L->Fill(BDTvalue,wT_*wPU_*wb_);
						
					}
					//if(BDTvalue>bdtcut)
					if (ef == 5) {
						/* if( NewselectSignal(CSV0, CSV1, DeltaPhiHV,
						 Emumass, Hmass, 
						 MET, Ht, AngleEMU,
						 UnweightedEta))
						 */
						if(ivar!=kbdt) hbdt_data0H->Fill(histoFillValue,wT_*wPU_*wb_);
						else hbdt_data0H->Fill(BDTvalue,wT_*wPU_*wb_);
					}
					//fill btag systematics now
					if(ivar==kbdt)
						FillBtagSys(BDTvalue, wT_*wPU_*wb_, eventFlavor, 
									hbdt_data0L_up,  hbdt_data0L_down,  
									hbdt_data0H_up, hbdt_data0H_down);
				}
			}
		}
	}
	
	float Hmass, Emumass;
	float Hpt, Zpt;
	float CSV0, CSV1;
	float jetCHF0, jetCHF1;
	float DphiJJ, DetaJJ, delRjj;
  	float Hphi, Zphi, DeltaPhiHV, Ht, MET, EventPt, PtbalZH, Mte, Mtmu, DphiZMET;
	float lep0pt, lep1pt, lep_pfCombRelIso0, lep_pfCombRelIso1, AngleEMU, delRemu;
	float EtaStandDev, RMS_eta, UnweightedEta;
	float Centrality, EvntShpAplanarity, EvntShpSphericity, EvntShpIsotropy, EvntShpCircularity;
	float Trigweight, B2011PUweight, A2011PUweight, btag2CSF, BDTvalue;
	float ZmassSVD;
	int nJets, eventFlavor, naJets, nSV;
	
	dataH->SetBranchAddress( "nJets", &nJets );
	dataH->SetBranchAddress( "naJets", &naJets );
	dataH->SetBranchAddress( "nSV", &nSV );	
	dataH->SetBranchAddress( "Hmass", &Hmass );
    dataH->SetBranchAddress( "Emumass", &Emumass );
    dataH->SetBranchAddress( "Hpt", &Hpt );
	dataH->SetBranchAddress( "Zpt", &Zpt );
	dataH->SetBranchAddress( "CSV0", &CSV0 );
	dataH->SetBranchAddress( "CSV1", &CSV1 );
	dataH->SetBranchAddress( "jetCHF0", &jetCHF0 );
	dataH->SetBranchAddress( "jetCHF1", &jetCHF1 );
	dataH->SetBranchAddress( "DetaJJ", &DetaJJ );
	dataH->SetBranchAddress( "DphiJJ", &DphiJJ );
	dataH->SetBranchAddress( "delRjj", &delRjj );
	dataH->SetBranchAddress( "Hphi", &Hphi );
	dataH->SetBranchAddress( "Zphi", &Zphi );
	dataH->SetBranchAddress( "DeltaPhiHV", &DeltaPhiHV );
	dataH->SetBranchAddress( "Ht", &Ht );
	dataH->SetBranchAddress( "MET", &MET );
	dataH->SetBranchAddress( "EventPt", &EventPt );
	dataH->SetBranchAddress( "PtbalZH", &PtbalZH );
	dataH->SetBranchAddress( "Mte", &Mte );
	dataH->SetBranchAddress( "Mtmu", &Mtmu );
	dataH->SetBranchAddress( "DphiZMET", &DphiZMET );
	dataH->SetBranchAddress( "Zphi", &Zphi );
	dataH->SetBranchAddress( "Hphi", &Hphi );
	dataH->SetBranchAddress( "lep_pfCombRelIso0", &lep_pfCombRelIso0 );
	dataH->SetBranchAddress( "lep_pfCombRelIso1", &lep_pfCombRelIso1 );
	dataH->SetBranchAddress( "lep0pt", &lep0pt );
	dataH->SetBranchAddress( "lep1pt", &lep1pt );
	dataH->SetBranchAddress( "AngleEMU", &AngleEMU );
	dataH->SetBranchAddress( "delRemu", &delRemu );
	dataH->SetBranchAddress( "RMS_eta", &RMS_eta );
	dataH->SetBranchAddress( "UnweightedEta", &UnweightedEta );
	dataH->SetBranchAddress( "EtaStandDev", &EtaStandDev );
	dataH->SetBranchAddress( "Centrality", &Centrality );
	dataH->SetBranchAddress( "EvntShpAplanarity", &EvntShpAplanarity );
	dataH->SetBranchAddress( "EvntShpSphericity", &EvntShpSphericity );
	dataH->SetBranchAddress( "EvntShpIsotropy", &EvntShpIsotropy );
	dataH->SetBranchAddress( "EvntShpCircularity", &EvntShpCircularity );
	dataH->SetBranchAddress( "eventFlavor", &eventFlavor );
	dataH->SetBranchAddress( "Trigweight", &Trigweight );
	dataH->SetBranchAddress( "B2011PUweight", &B2011PUweight );
	dataH->SetBranchAddress( "A2011PUweight", &A2011PUweight );
	dataH->SetBranchAddress( "btag2CSF", &btag2CSF );
	dataH->SetBranchAddress( "BDTvalue", &BDTvalue );
	dataH->SetBranchAddress( "ZmassSVD", &ZmassSVD );
	
	
	for(int i=0;i<dataH->GetEntries();i++){
		dataH->GetEntry(i);
		
		if (ivar==kl1pt) histoFillValue = lep0pt;
		if (ivar==kl2pt) histoFillValue =  lep1pt; //1
		if (ivar==kpTjj) histoFillValue = Zpt;  //2
		if (ivar==kpTZ) histoFillValue =  Zpt; //3
		if (ivar==kcsv1) histoFillValue = CSV0;  //4
		if (ivar==kcsv2) histoFillValue = CSV1;  //5 
		if (ivar==kDphi) histoFillValue = DeltaPhiHV; //6
		if (ivar==kNaj) histoFillValue =  naJets; //7 
		if (ivar==kmet) histoFillValue = MET;  //8
		if (ivar==kDeta) histoFillValue = DetaJJ;   //9
		if (ivar==kHt) histoFillValue =  Ht; //10
		if (ivar==kjchf1) histoFillValue = jetCHF0;  //11
		if (ivar==kjchf2) histoFillValue = jetCHF1;  //12
		if (ivar==kangleEMU) histoFillValue =  AngleEMU; //13
		if (ivar==kueta) histoFillValue = UnweightedEta;  //14 
		if (ivar==kcentrality) histoFillValue = Centrality;  //15
		if (ivar==kM) histoFillValue =  Hmass;  //16
		if (ivar==kDptHZ) histoFillValue = PtbalZH;   //17
		if (ivar==hdrjj) histoFillValue =  delRjj;  //18
		if (ivar==kbdt) histoFillValue = BDTvalue;   //19
		if (ivar==kdeltaPhiZMET) histoFillValue = DphiZMET;   //20
		if (ivar==kAplanarity) histoFillValue = EvntShpAplanarity;   //21
		if (ivar==kSphericity) histoFillValue = EvntShpSphericity;  //22
		if (ivar==kCircularity) histoFillValue = EvntShpCircularity; //23
		if (ivar==kIsotropy) histoFillValue = EvntShpIsotropy;  //24
		if (ivar==kMte) histoFillValue = Mte;   //25
		if (ivar==kMtmu) histoFillValue = Mtmu;   //26		
		
		double wT_=Trigweight;
		double wPU_= 1.0;
		wPU_= (lA*A2011PUweight + lB*B2011PUweight)/lumi;
		double wb_=1.0;
		wT_ = 1.0;
		int ef= eventFlavor;   
		
		if(Preselect(CSV0, CSV1, Emumass, jetCHF0, Hmass, ZmassSVD , Hpt)){
			if(PassBoost(Hpt, Zpt, boost, Hmass)) {
				if( BDTvalue>bdtcut){
					if (ef != 5){
						/* if( NewselectSignal(CSV0, CSV1, DeltaPhiHV,
						 Emumass, Hmass, 
						 MET, Ht, AngleEMU,
						 UnweightedEta))
						 */ 
						if(ivar!=kbdt) hbdt_data4L->Fill(histoFillValue,wT_*wPU_*wb_);
						else hbdt_data4L->Fill(BDTvalue,wT_*wPU_*wb_);
					}
					if (ef == 5) {
						/*	if( NewselectSignal(CSV0, CSV1, DeltaPhiHV,
						 Emumass, Hmass, 
						 MET, Ht, AngleEMU,
						 UnweightedEta))
						 */ 
						if(ivar!=kbdt)  hbdt_data4H->Fill(histoFillValue,wT_*wPU_*wb_);
						else hbdt_data4H->Fill(BDTvalue,wT_*wPU_*wb_);
					}
					if(ivar==kbdt)
						FillBtagSys(BDTvalue, wT_*wPU_*wb_, eventFlavor, 
									hbdt_data4L_up,  hbdt_data4L_down,  
									hbdt_data4H_up, hbdt_data4H_down);	
				}
			}   
		}
	}//
	
	//now scale down properly
	hbdt_data0L->Scale(wL);
	hbdt_data0H->Scale(wL);
	hbdt_data4L->Scale(wH);
	hbdt_data4H->Scale(wH);
	
	hL->Add(hbdt_data0L);
	hL->Add(hbdt_data4L);
	hH->Add(hbdt_data0H);
	hH->Add(hbdt_data4H);
	
	hL->SetLineWidth(2);
	hH->SetLineWidth(2);
	
	if(ivar==kbdt){
		hbdt_data0L_up->Scale(wL);
		hbdt_data0H_up->Scale(wL);
		hbdt_data0L_down->Scale(wL);
		hbdt_data0H_down->Scale(wL);
		hbdt_data4L_up->Scale(wH);
		hbdt_data4H_up->Scale(wH);
		hbdt_data4L_down->Scale(wH);
		hbdt_data4H_down->Scale(wH);
		
		hL_up->Add(hbdt_data0L_up);
		hL_up->Add(hbdt_data4L_up);
		hL_down->Add(hbdt_data0L_down);
		hL_down->Add(hbdt_data4L_down);
		hH_up->Add(hbdt_data0H_up);
		hH_up->Add(hbdt_data4H_up);
		hH_down->Add(hbdt_data0H_down);
		hH_down->Add(hbdt_data4H_down);
		
		
		//stats error as an histogram, should work if scaling is done correctly
		for(int i=0;i<hbdt_data0L_up->GetNbinsX()+1;i++){
			hstatLup->SetBinContent(i, hbdt_data0L->GetBinContent(i)+ hbdt_data4L->GetBinContent(i)+
									sqrt(  sqrt(hbdt_data0L->GetBinContent(i)/wL)*wL * sqrt(hbdt_data0L->GetBinContent(i)/wL)*wL +
										 sqrt(hbdt_data4L->GetBinContent(i)/wH)*wH * sqrt(hbdt_data4L->GetBinContent(i)/wH)*wH ));
			hstatLdown->SetBinContent(i,  hbdt_data0L->GetBinContent(i)+ hbdt_data4L->GetBinContent(i) -
									  sqrt(  sqrt(hbdt_data0L->GetBinContent(i)/wL)*wL * sqrt(hbdt_data0L->GetBinContent(i)/wL)*wL +
										   sqrt(hbdt_data4L->GetBinContent(i)/wH)*wH * sqrt(hbdt_data4L->GetBinContent(i)/wH)*wH ));
			hstatHup->SetBinContent(i, hbdt_data0H->GetBinContent(i) + hbdt_data4H->GetBinContent(i) +
									sqrt( sqrt(hbdt_data0H->GetBinContent(i)/wL)*wL * sqrt(hbdt_data0H->GetBinContent(i)/wL)*wL +
										 + sqrt( hbdt_data4H->GetBinContent(i)/wH)*wH * sqrt( hbdt_data4H->GetBinContent(i)/wH)*wH ));
			hstatHdown->SetBinContent(i, hbdt_data0H->GetBinContent(i) + hbdt_data4H->GetBinContent(i) -
									  sqrt( sqrt(hbdt_data0H->GetBinContent(i)/wL)*wL * sqrt(hbdt_data0H->GetBinContent(i)/wL)*wL +
										   + sqrt( hbdt_data4H->GetBinContent(i)/wH)*wH * sqrt( hbdt_data4H->GetBinContent(i)/wH)*wH ));
		}
	} //ivar==kbdt
	
	expZL= hbdt_data0L->Integral(0,-1)+ hbdt_data4L->Integral(0,-1);
	expZH= hbdt_data0H->Integral(0,-1)+ hbdt_data4H->Integral(0,-1);
	e_expZL= sqrt ( sqrt(hbdt_data0L->Integral(0,-1)/wL)*wL* sqrt(hbdt_data0L->Integral(0,-1)/wL)*wL+ 
				   sqrt(hbdt_data4L->Integral(0,-1)/wH)*wH* sqrt(hbdt_data4L->Integral(0,-1)/wH)*wH);
	e_expZH= sqrt ( sqrt(hbdt_data0H->Integral(0,-1)/wL)*wL* sqrt(hbdt_data0H->Integral(0,-1)/wL)*wL+ 
				   sqrt(hbdt_data4H->Integral(0,-1)/wH)*wH* sqrt(hbdt_data4H->Integral(0,-1)/wH)*wH);
	
}


void plotvariables_updated(bool isMatt=true, int mass=115,  bool isFinal=false, int boost=noboost, bool EvalPunzi=false, bool twod=false) 
{
	gStyle->SetOptStat(false);
	
	inputS = TFile::Open("/home/hep/wilken/taus/CMSSW_4_4_2_patch8/src/UserCode/wilken/TMVA/BDTCut_Fall11_ZH115.root");
	inputTT = TFile::Open("/home/hep/wilken/taus/CMSSW_4_4_2_patch8/src/UserCode/wilken/TMVA/BDTCut_TTJets.root"); 
	inputTtW = TFile::Open("/home/hep/wilken/taus/CMSSW_4_4_2_patch8/src/UserCode/wilken/TMVA/BDTCut_T_tW.root");
	inputTBARtW = TFile::Open("/home/hep/wilken/taus/CMSSW_4_4_2_patch8/src/UserCode/wilken/TMVA/BDTCut_Tbar_tW.root");
	inputZZ  = TFile::Open("/home/hep/wilken/taus/CMSSW_4_4_2_patch8/src/UserCode/wilken/TMVA/BDTCut_ZZ.root"); 
	inputWW  = TFile::Open("/home/hep/wilken/taus/CMSSW_4_4_2_patch8/src/UserCode/wilken/TMVA/BDTCut_WW.root"); 
	inputWZ  = TFile::Open("/home/hep/wilken/taus/CMSSW_4_4_2_patch8/src/UserCode/wilken/TMVA/BDTCut_WZ.root"); 
	inputDY = TFile::Open("/home/hep/wilken/taus/CMSSW_4_4_2_patch8/src/UserCode/wilken/TMVA/BDTCut_DYM50_both.root");
	inputDY2 = TFile::Open("/home/hep/wilken/taus/CMSSW_4_4_2_patch8/src/UserCode/wilken/TMVA/BDTCut_DYPtZ_both.root");
	inputWJ  = TFile::Open("/home/hep/wilken/taus/CMSSW_4_4_2_patch8/src/UserCode/wilken/TMVA/BDTCut_WJets.root"); 
	
	//  inputdata  = TFile::Open("/scratch/tinti_g/BDT/data_all.root");
	
	//  TTree *data = (TTree*) inputdata->Get("treeWithBDT");
	TTree *signal = (TTree*) inputS->Get("treeWithBDT"); 
	TTree *bTT  = (TTree*) inputTT ->Get("treeWithBDT");   
	TTree* bTtW  = (TTree*) inputTtW ->Get("treeWithBDT");   
	TTree* bTBARtW  = (TTree*) inputTBARtW ->Get("treeWithBDT");   
	TTree* bZZ  = (TTree*) inputZZ ->Get("treeWithBDT");   
	TTree* bWW  = (TTree*) inputWW ->Get("treeWithBDT");  
	TTree* bWZ  = (TTree*) inputWZ ->Get("treeWithBDT");   
	TTree* bDY  = (TTree*) inputDY ->Get("treeWithBDT");   
	TTree* bDY2  = (TTree*) inputDY2 ->Get("treeWithBDT");   
	TTree* bWJ  = (TTree*) inputWJ ->Get("treeWithBDT");   
	
	vector <TTree*> t;
	t.push_back(bTT); //0
	t.push_back(bTtW); //1
	t.push_back(bTBARtW); //2
	t.push_back(bZZ); //3
	t.push_back(bWW); //4
	t.push_back(bWZ);//5 
	t.push_back(bDY);//6 
	t.push_back(bDY2);//7
	t.push_back(bWJ);//8 
	
	//bin0 in mass-> 0-70 GeV
	//bin1 in mass-> 70-90 GeV
	//bin2 in mass-> 90-120 GeV
	//bin3 in mass-> 1200-250 GeV
	
	//bin0 in bdt-> -1.-0.9
	//bin1 in bdt-> -0.9,-0.6
	//bin2 in bdt-> -0.6,0
	
	int n2dbins=-1;
	
	
    //bdtcut has leftmost cut
	//  const int bdtcutsteps=120;//80;     
	// for(int iscut=0; iscut<bdtcutsteps;iscut++){
	//if(iscut%2) continue; 	
    //want to scan bdtvalues from -0.5 to 
    //bdtcut=/*-0.4*/-1.+iscut*0.005;
	
	inputS2 =  TFile::Open("/home/hep/wilken/taus/CMSSW_4_4_2_patch8/src/UserCode/wilken/TMVA/BDTCut_Filteredemu.root");
	TTree *signal2 = (TTree*) inputS2->Get("treeWithBDT"); 
	cout<< signal->GetEntries()<<"   "<<signal2->GetEntries()<<"   "<< signal->GetEntries()/signal2->GetEntries()<<"   "<< signal->GetEntries()/(float)signal2->GetEntries()<< endl;
	
	//  assert(0);
	
	vector <double> sc;
	sc.push_back(lumi/(lumiTT/ 2.0));
	sc.push_back(lumi/( lumiTtW/ 2.0));
	sc.push_back(lumi/( lumiTtWb/2.0));
	sc.push_back(lumi/( lumiZZ/2.0));
	sc.push_back(lumi/( lumiWW/2.0));
	sc.push_back(lumi/( lumiWZ/2.0));
	sc.push_back(lumi/( lumiZJH/2.0));
	sc.push_back(lumi/( lumiZJL/2.0));
	sc.push_back(lumi/( lumiWJ/2.0));
	
	vector <string> ch;
	ch.push_back("TT"); //0
	ch.push_back("TtW"); //1
	ch.push_back("TtWb"); //2
	ch.push_back("ZZ"); //3
	ch.push_back("WW"); //4
	ch.push_back("WZ");  //5
	ch.push_back("ZJH"); //6 
	ch.push_back("ZJL"); //7
	ch.push_back("WJ");  //8
	
	//signal
	float scaleZH;
	if(mass==115)
//	float lumiFallBDT = 2.28079234375000000e+05/xsecbfZH115; //fall11 BDT_tree
//	scaleZH = lumi/(lumiFallBDT/2.0);
	scaleZH=lumi/(2.28079234375000000e+05/xsecbfZH115); //fall11 FOM_tree
	//scaleZH=lumi/(lumiZH115/2.0);
	if(mass==120)
		scaleZH=lumi/(numZH120/xsecbfZH120/(isFinal ? 1: 2));
	if(mass==125)
		scaleZH=lumi/(numZH125/xsecbfZH125/(isFinal ? 1: 2));
	if(mass==130)
		scaleZH=lumi/(numZH130/xsecbfZH130/(isFinal ? 1: 2));
	if(mass==135)
		scaleZH=lumi/(numZH135/xsecbfZH135/(isFinal ? 1: 2));
	if(mass==110)
		scaleZH=lumi/(numZH110/xsecbfZH110/(isFinal ? 1: 2));
	
	
	vector <string> var;
	vector <string> varname;
	var.push_back("lep0pt");  varname.push_back("lep1 p_{T} (GeV)"); //0
	var.push_back("lep1pt");  varname.push_back("lep2 p_{T} (GeV)"); //1
	var. push_back("Hpt");  varname.push_back("H p_{T} (GeV)"); //2
	var. push_back("Zpt");  varname.push_back("Z p_{T} (GeV)"); //3
	var. push_back("CSV0");  varname.push_back("csv1"); //4
	var. push_back("CSV1");  varname.push_back("csv2"); //5
	var. push_back("DeltaPhiHV");  varname.push_back("#Delta #phi (H,Z)"); //6
	var. push_back("naJets");  varname.push_back("# additional jets"); //7
	var. push_back("MET");  varname.push_back("M E_{T} (GeV)"); //8
	var. push_back("DetaJJ");  varname.push_back("#Delta #eta (j1,j2)"); //9
	var. push_back("Ht");  varname.push_back("#Sum component p_{T} (GeV)"); //10
	var. push_back("jetCHF0");  varname.push_back("charge fraction j1"); //11
	var. push_back("jetCHF1");  varname.push_back("charge fraction j2"); //12
	var. push_back("AngleEMU");  varname.push_back("angle(e,#mu)"); //13
	var. push_back("UnweightedEta");  varname.push_back("unweighted #eta"); //14
	var. push_back("Centrality");  varname.push_back("centrality"); //15
	var. push_back("Hmass");  varname.push_back("H mass (GeV)"); //16
	var. push_back("PtbalZH");  varname.push_back("H p_{T} - Z p_{T} (GeV)"); //17
	var. push_back("delRjj");  varname.push_back("dR (j,j)"); //18
	var. push_back("bdt");  //19
	varname.push_back("BDT output");
	var. push_back("DphiZMET");  varname.push_back("DphiZMET");
	var. push_back("EvntShpAplanarity");  varname.push_back("aplanarity");
	var. push_back("EvntShpSphericity");  varname.push_back("sphericity");
	var. push_back("EvntShpCircularity");  varname.push_back("circularity");
	var. push_back("EvntShpIsotropy");  varname.push_back("isotropy");
	var. push_back("Mte");  varname.push_back("Mte");
	var. push_back("Mtmu");  varname.push_back("Mtmu");
	
	const int NVar=27;
	if(var.size() != NVar) assert(0 && "mismatch between variables and plots: this should not happen \n");
	const int signalcolor=kBlack;
	const int Nbkg=10; //bkg
	
	int nbins[NVar]={12,12,13,13,10,10,30,10,12,50,50,10,10,20,10,20,25,40,50,/*14*/22,20,40,20,20,30,50,50};
	float bl[NVar]={10,10,10,10,0,0,0,1,10,0,10,0,0,0,0,0,0,-40,0,-1.1, -4,-0.1,0,0,0,0,0};
	float bh[NVar]={130,130,140,140,1,1,3.14,11,80,5,400,1,1,3.14,10,1,250,40,10,/*0.35*/1.1,4, 0.1,0.8,0.8,1,200,200};
	
	
	TH1F* hs[NVar];
	for(int itvar=0;itvar<NVar;itvar++)
		hs[itvar]=new TH1F(TString::Format("hs%d",itvar).Data(),
						   TString::Format("hs%d",itvar).Data(),
						   nbins[itvar],bl[itvar],bh[itvar]);
	
	TH1F* hdata[NVar];
	for(int itvar=0;itvar<NVar;itvar++)
		hdata[itvar]=new TH1F(TString::Format("hdata%d",itvar).Data(),
							  TString::Format("hdata%d",itvar).Data(),
							  nbins[itvar],bl[itvar],bh[itvar]);
	
	//only for variable BDT, create the 2histograms for systematics
	TH1F* hs_btag_up;
	TH1F* hs_btag_down;
	hs_btag_up=new TH1F("hs_btag_up", "hs_btag_up", nbins[kbdt],bl[kbdt],bh[kbdt]);
	hs_btag_down=new TH1F("hs_btag_down", "hs_btag_down", nbins[kbdt],bl[kbdt],bh[kbdt]);
	TH1F* hs_stat_up;
	TH1F* hs_stat_down;
	hs_stat_up=new TH1F("hs_stat_up", "hs_stat_up", nbins[kbdt],bl[kbdt],bh[kbdt]);
	hs_stat_down=new TH1F("hs_stat_down", "hs_stat_down", nbins[kbdt],bl[kbdt],bh[kbdt]);
	
	TH1F* h[NVar][Nbkg];
	for(int itvar=0;itvar<NVar;itvar++)
		for(int j=0;j<Nbkg;j++)
			h[itvar][j]=new TH1F(TString::Format("h%d_%d",itvar,j).Data(),
								 TString::Format("h%d_%d",itvar,j).Data(),
								 nbins[itvar],bl[itvar],bh[itvar]);
	
	TH1F* h_btag_up[Nbkg];
	TH1F* h_btag_down[Nbkg];
	TH1F* h_stat_up[Nbkg];
	TH1F* h_stat_down[Nbkg];
	for(int j=0;j<Nbkg;j++){
		h_btag_up[j]=new TH1F(TString::Format("h_b_tag_up%d",j).Data(),
							  TString::Format("h_b_tag_up%d",j).Data(),
							  nbins[kbdt],bl[kbdt],bh[kbdt]);
		h_btag_down[j]=new TH1F(TString::Format("h_b_tag_down%d",j).Data(),
								TString::Format("h_b_tag_down%d",j).Data(),
								nbins[kbdt],bl[kbdt],bh[kbdt]);
		h_stat_up[j]=new TH1F(TString::Format("h_stat_up%d",j).Data(),
							  TString::Format("h_stat_up%d",j).Data(),
							  nbins[kbdt],bl[kbdt],bh[kbdt]);
		h_stat_down[j]=new TH1F(TString::Format("h_stat_down%d",j).Data(),
								TString::Format("h_stat_down%d",j).Data(),
								nbins[kbdt],bl[kbdt],bh[kbdt]);
	}
	
	TH1F*  hL[NVar];
	TH1F*  hH[NVar];
	for(int itvar=0;itvar<NVar;itvar++){
		hL[itvar]=new TH1F(TString::Format("hL%d",itvar).Data(),
						   TString::Format("hL%d",itvar).Data(),
						   nbins[itvar],bl[itvar],bh[itvar]);
		hH[itvar]=new TH1F(TString::Format("hH%d",itvar).Data(),
						   TString::Format("hH%d",itvar).Data(),
						   nbins[itvar],bl[itvar],bh[itvar]);
	}
	
	TH1F*  hL_btag_up;
	TH1F*  hL_btag_down;
	TH1F*  hH_btag_up;
	TH1F*  hH_btag_down;
	TH1F*  hL_stat_up;
	TH1F*  hL_stat_down;
	TH1F*  hH_stat_up;
	TH1F*  hH_stat_down;
	hL_btag_up=new TH1F("hL_btag_up", "hL_btag_up", nbins[kbdt],bl[kbdt],bh[kbdt]);
	hL_btag_down=new TH1F("hL_btag_down", "hL_btag_down", nbins[kbdt],bl[kbdt],bh[kbdt]);
	hH_btag_up=new TH1F("hH_btag_up", "hH_btag_up", nbins[kbdt],bl[kbdt],bh[kbdt]);
	hH_btag_down=new TH1F("hH_btag_down", "hH_btag_down", nbins[kbdt],bl[kbdt],bh[kbdt]);
	hL_stat_up=new TH1F("hL_stat_up", "hL_stat_up", nbins[kbdt],bl[kbdt],bh[kbdt]);
	hL_stat_down=new TH1F("hL_stat_down", "hL_stat_down", nbins[kbdt],bl[kbdt],bh[kbdt]);
	hH_stat_up=new TH1F("hH_stat_up", "hH_stat_up", nbins[kbdt],bl[kbdt],bh[kbdt]);
	hH_stat_down=new TH1F("hH_stat_down", "hH_stat_down", nbins[kbdt],bl[kbdt],bh[kbdt]);
	
	//	int color[Nbkg]={kRed, kBlue, kGreen+2, kOrange, kCyan, kViolet, kGray, kRed-7, kYellow+2, kMagenta-5, kOrange+4, kGreen-2};
	int color[Nbkg]={kRed, kViolet, kGray, kRed-7, kYellow+2, kMagenta-5, kOrange+4, kGreen-2};
	
	TCanvas* c[NVar];
	for(int itvar=0;itvar<NVar;itvar++){
		c[itvar]=new TCanvas();
	}
	
	cout << "float rate and stat variables" << endl;
	
	float eVH, eWj, eZjLF, eZjHF,  eTT, es_Top, eVV, e_eVH, e_eWj, e_eZjLF, e_eZjHF, e_eTT, e_es_Top, e_eVV;
	vector <float> exp;
	vector <float> e_exp;
	cout << "what the hell is Nbkg-1 " << Nbkg-1 << endl;
	for(unsigned int jj=0;jj<Nbkg-1;jj++) {
		exp.push_back(-9999);
		e_exp.push_back(-9999);
	}
	cout << "rate and stat array initialized" << endl;
	
	//do the data first
	//  for(unsigned int ivar=0;ivar<var.size();ivar++)
	//  HistoFiller(data,hdata[ivar],var[ivar],boost, ivar, isMatt, n2dbins);
	
	//here do function to determine 2d/1d binning
	// Make1dBinning(signal,scaleZH, boost, isMatt, -1);
	
	for(unsigned int ivar=0;ivar<var.size();ivar++){
		HistoFiller(signal,scaleZH,hs[ivar],var[ivar], boost,eVH,e_eVH, ivar,
					hs_btag_up, hs_btag_down, hs_stat_up, hs_stat_down, 
					isMatt, n2dbins);
		cout << "signal histos filled using HistoFiller signalhisto "<< hs[0]->Integral(0,-1) << endl;
		
		for(int j=0;j<Nbkg-1;j++){
			if((j!=6) && (j!=7)) 
				HistoFiller(t[j],sc[j],h[ivar][j],var[ivar], boost, exp[j], e_exp[j],
							ivar, h_btag_up[j], h_btag_down[j],  h_stat_up[j], h_stat_down[j], 
							isMatt, n2dbins);
		}
		
		HistoFillerLH(t[7],sc[7],t[6],sc[6], hL[ivar], hH[ivar], var[ivar], 
					  boost,  eZjLF, eZjHF, e_eZjLF,  e_eZjHF, ivar, hL_btag_up,
					  hL_btag_down, hH_btag_up, hH_btag_down, hL_stat_up,
					  hL_stat_down, hH_stat_up, hH_stat_down, 
					  isMatt, n2dbins);
		
		
		//CONTROL REGION CORRECTIONS
		h[ivar][0]->Scale(SFtt);
		hL[ivar]->Scale(SFZl);
		hH[ivar]->Scale(SFZh);
		if(ivar==kbdt){
			h_btag_up[0]->Scale(SFtt);
			h_btag_down[0]->Scale(SFtt);    
			hL_btag_up->Scale(SFZl);
			hL_btag_down->Scale(SFZl);
			hH_btag_up->Scale(SFZh);
			hH_btag_down->Scale(SFZh);
			h_stat_up[0]->Scale(SFtt);
			h_stat_down[0]->Scale(SFtt);    
			hL_stat_up->Scale(SFZl);
			hL_stat_down->Scale(SFZl);
			hH_stat_up->Scale(SFZh);
			hH_stat_down->Scale(SFZh);
		}
		
		//now find expected bkgs and errors
		eZjLF*=SFZl;
		
		e_eZjLF*=SFZl;
		eZjHF*=SFZh;
		e_eZjHF*=SFZh;
		
		eTT=exp[0]*SFtt;
		e_eTT=e_exp[0]*SFtt;
		es_Top=exp[1]+exp[2];
		e_es_Top= sqrt( e_exp[1]*e_exp[1] + e_exp[2]*e_exp[2]);
		eVV= exp[3]+exp[4]+exp[5];
		e_eVV= sqrt(e_exp[3]*e_exp[3] + e_exp[4]*e_exp[4]+ e_exp[5]*e_exp[5]);
		
		
		eWj=exp[8];
		cout << "rate of Wjets " << eWj << endl;
		e_eWj=e_exp[8];
		
		c[ivar]->cd();
		hs[ivar]->GetXaxis()->SetTitle(varname[ivar].c_str());
		hs[ivar]->GetYaxis()->SetTitle("entries");
		hs[ivar]->SetLineColor(kBlack);
		hs[ivar]->SetFillColor(kBlack);
		hs[ivar]->Draw("hist e");
		float max=-9;
		
		if(ivar==kM) {
			double B=0;
			for(int j=0;j<Nbkg-1;j++){
				if((j!=6) && (j!=7))
					B+=h[ivar][j]->Integral(0,-1);
			}
			B+= hL[ivar]->Integral(0,-1);
			B+= hH[ivar]->Integral(0,-1);
			//     cout<<"punzi  "<<hs[ivar]->Integral(0,-1)/100./ (1.5 + TMath::Sqrt(B)+ 0.2*B)<<endl; 
		}
		
		if(EvalPunzi){
			if(ivar==kbdt) {
				//evaluate punzi
				
				const int bdtsteps=220;
				
				Double_t bdtlimit[bdtsteps];
				Double_t counter[bdtsteps];
				Double_t c_b[bdtsteps];  
				
				for(int is=0; is<bdtsteps;is++){
					bdtlimit[is]=-1+is*0.005;
					counter[is]=0.;
					c_b[is]=0.;
				}
				
				float Hmass, Emumass;
				float Hpt, Zpt;
				float CSV0, CSV1;
				float jetCHF0, jetCHF1;
				float DphiJJ, DetaJJ, delRjj;
				float Hphi, Zphi, DeltaPhiHV, Ht, MET, EventPt, PtbalZH, Mte, Mtmu, DphiZMET;
				float lep0pt, lep1pt, lep_pfCombRelIso0, lep_pfCombRelIso1, AngleEMU, delRemu;
				float EtaStandDev, RMS_eta, UnweightedEta;
				float Centrality, EvntShpAplanarity, EvntShpSphericity, EvntShpIsotropy, EvntShpCircularity;
				float Trigweight, B2011PUweight, A2011PUweight, btag2CSF, BDTvalue;
				float ZmassSVD;
				int nJets, eventFlavor, naJets, nSV;
				
				signal->SetBranchAddress( "nJets", &nJets );
				signal->SetBranchAddress( "naJets", &naJets );
				signal->SetBranchAddress( "nSV", &nSV );	
				signal->SetBranchAddress( "Hmass", &Hmass );
				signal->SetBranchAddress( "Emumass", &Emumass );
				signal->SetBranchAddress( "Hpt", &Hpt );
				signal->SetBranchAddress( "Zpt", &Zpt );
				signal->SetBranchAddress( "CSV0", &CSV0 );
				signal->SetBranchAddress( "CSV1", &CSV1 );
				signal->SetBranchAddress( "jetCHF0", &jetCHF0 );
				signal->SetBranchAddress( "jetCHF1", &jetCHF1 );
				signal->SetBranchAddress( "DetaJJ", &DetaJJ );
				signal->SetBranchAddress( "DphiJJ", &DphiJJ );
				signal->SetBranchAddress( "delRjj", &delRjj );
				signal->SetBranchAddress( "Hphi", &Hphi );
				signal->SetBranchAddress( "Zphi", &Zphi );
				signal->SetBranchAddress( "DeltaPhiHV", &DeltaPhiHV );
				signal->SetBranchAddress( "Ht", &Ht );
				signal->SetBranchAddress( "MET", &MET );
				signal->SetBranchAddress( "EventPt", &EventPt );
				signal->SetBranchAddress( "PtbalZH", &PtbalZH );
				signal->SetBranchAddress( "Mte", &Mte );
				signal->SetBranchAddress( "Mtmu", &Mtmu );
				signal->SetBranchAddress( "DphiZMET", &DphiZMET );
				signal->SetBranchAddress( "Zphi", &Zphi );
				signal->SetBranchAddress( "Hphi", &Hphi );
				signal->SetBranchAddress( "lep_pfCombRelIso0", &lep_pfCombRelIso0 );
				signal->SetBranchAddress( "lep_pfCombRelIso1", &lep_pfCombRelIso1 );
				signal->SetBranchAddress( "lep0pt", &lep0pt );
				signal->SetBranchAddress( "lep1pt", &lep1pt );
				signal->SetBranchAddress( "AngleEMU", &AngleEMU );
				signal->SetBranchAddress( "delRemu", &delRemu );
				signal->SetBranchAddress( "RMS_eta", &RMS_eta );
				signal->SetBranchAddress( "UnweightedEta", &UnweightedEta );
				signal->SetBranchAddress( "EtaStandDev", &EtaStandDev );
				signal->SetBranchAddress( "Centrality", &Centrality );
				signal->SetBranchAddress( "EvntShpAplanarity", &EvntShpAplanarity );
				signal->SetBranchAddress( "EvntShpSphericity", &EvntShpSphericity );
				signal->SetBranchAddress( "EvntShpIsotropy", &EvntShpIsotropy );
				signal->SetBranchAddress( "EvntShpCircularity", &EvntShpCircularity );
				signal->SetBranchAddress( "eventFlavor", &eventFlavor );
				signal->SetBranchAddress( "Trigweight", &Trigweight );
				signal->SetBranchAddress( "B2011PUweight", &B2011PUweight );
				signal->SetBranchAddress( "A2011PUweight", &A2011PUweight );
				signal->SetBranchAddress( "btag2CSF", &btag2CSF );
				signal->SetBranchAddress( "BDTvalue", &BDTvalue );
				signal->SetBranchAddress( "ZmassSVD", &ZmassSVD );
				
				
				for(int i=0;i<signal->GetEntries();i++){
					signal->GetEntry(i);
					
					double l = 	BDTvalue;					
					//double wPU_=(wPU->GetValue()*lA+B2011PUweight*lB)/(lA+lB);
					double wT_ = Trigweight;
					double wPU_ = B2011PUweight;
					wPU_= (lA*A2011PUweight + lB*B2011PUweight)/lumi;
					wT_ = 1.0;
					
					if(Preselect(CSV0, CSV1, Emumass, jetCHF0, Hmass, ZmassSVD , Hpt)){
						if(PassBoost(Hpt, Zpt, boost, Hmass))
							for(int is=0; is<bdtsteps;is++)
								if(l>bdtlimit[is]) 
									counter[is]+=scaleZH*wT_*wPU_*1.0;
					}
				} //loop on signal
				
				
				//start loop on background
				for(int j=0;j<Nbkg-1;j++){
					
					cout << "reading tree for bkg " << j << endl;
					
					float Hmass, Emumass;
					float Hpt, Zpt;
					float CSV0, CSV1;
					float jetCHF0, jetCHF1;
					float DphiJJ, DetaJJ, delRjj;
					float Hphi, Zphi, DeltaPhiHV, Ht, MET, EventPt, PtbalZH, Mte, Mtmu, DphiZMET;
					float lep0pt, lep1pt, lep_pfCombRelIso0, lep_pfCombRelIso1, AngleEMU, delRemu;
					float EtaStandDev, RMS_eta, UnweightedEta;
					float Centrality, EvntShpAplanarity, EvntShpSphericity, EvntShpIsotropy, EvntShpCircularity;
					float Trigweight, B2011PUweight, A2011PUweight, btag2CSF, BDTvalue;
					float ZmassSVD;
					int nJets, eventFlavor, naJets, nSV;
					
					t[j]->SetBranchAddress( "nJets", &nJets );
					t[j]->SetBranchAddress( "naJets", &naJets );
					t[j]->SetBranchAddress( "nSV", &nSV );	
					t[j]->SetBranchAddress( "Hmass", &Hmass );
					t[j]->SetBranchAddress( "Emumass", &Emumass );
					t[j]->SetBranchAddress( "Hpt", &Hpt );
					t[j]->SetBranchAddress( "Zpt", &Zpt );
					t[j]->SetBranchAddress( "CSV0", &CSV0 );
					t[j]->SetBranchAddress( "CSV1", &CSV1 );
					t[j]->SetBranchAddress( "jetCHF0", &jetCHF0 );
					t[j]->SetBranchAddress( "jetCHF1", &jetCHF1 );
					t[j]->SetBranchAddress( "DetaJJ", &DetaJJ );
					t[j]->SetBranchAddress( "DphiJJ", &DphiJJ );
					t[j]->SetBranchAddress( "delRjj", &delRjj );
					t[j]->SetBranchAddress( "Hphi", &Hphi );
					t[j]->SetBranchAddress( "Zphi", &Zphi );
					t[j]->SetBranchAddress( "DeltaPhiHV", &DeltaPhiHV );
					t[j]->SetBranchAddress( "Ht", &Ht );
					t[j]->SetBranchAddress( "MET", &MET );
					t[j]->SetBranchAddress( "EventPt", &EventPt );
					t[j]->SetBranchAddress( "PtbalZH", &PtbalZH );
					t[j]->SetBranchAddress( "Mte", &Mte );
					t[j]->SetBranchAddress( "Mtmu", &Mtmu );
					t[j]->SetBranchAddress( "DphiZMET", &DphiZMET );
					t[j]->SetBranchAddress( "Zphi", &Zphi );
					t[j]->SetBranchAddress( "Hphi", &Hphi );
					t[j]->SetBranchAddress( "lep_pfCombRelIso0", &lep_pfCombRelIso0 );
					t[j]->SetBranchAddress( "lep_pfCombRelIso1", &lep_pfCombRelIso1 );
					t[j]->SetBranchAddress( "lep0pt", &lep0pt );
					t[j]->SetBranchAddress( "lep1pt", &lep1pt );
					t[j]->SetBranchAddress( "AngleEMU", &AngleEMU );
					t[j]->SetBranchAddress( "delRemu", &delRemu );
					t[j]->SetBranchAddress( "RMS_eta", &RMS_eta );
					t[j]->SetBranchAddress( "UnweightedEta", &UnweightedEta );
					t[j]->SetBranchAddress( "EtaStandDev", &EtaStandDev );
					t[j]->SetBranchAddress( "Centrality", &Centrality );
					t[j]->SetBranchAddress( "EvntShpAplanarity", &EvntShpAplanarity );
					t[j]->SetBranchAddress( "EvntShpSphericity", &EvntShpSphericity );
					t[j]->SetBranchAddress( "EvntShpIsotropy", &EvntShpIsotropy );
					t[j]->SetBranchAddress( "EvntShpCircularity", &EvntShpCircularity );
					t[j]->SetBranchAddress( "eventFlavor", &eventFlavor );
					t[j]->SetBranchAddress( "Trigweight", &Trigweight );
					t[j]->SetBranchAddress( "B2011PUweight", &B2011PUweight );
					t[j]->SetBranchAddress( "A2011PUweight", &A2011PUweight );
					t[j]->SetBranchAddress( "btag2CSF", &btag2CSF );
					t[j]->SetBranchAddress( "BDTvalue", &BDTvalue );
					t[j]->SetBranchAddress( "ZmassSVD", &ZmassSVD );
					
					
					for(int i=0;i<t[j]->GetEntries();i++){
						t[j]->GetEntry(i);
						
						double l = BDTvalue;
						double wT_ = Trigweight;
						double wPU_= 1.0;
						wPU_= (lA*A2011PUweight + lB*B2011PUweight)/lumi;
						wT_ = 1.0;
						
						if(Preselect(CSV0, CSV1, Emumass, jetCHF0, Hmass, ZmassSVD , Hpt))
							if(PassBoost(Hpt, Zpt, boost, Hmass))
								for(int is=0; is<bdtsteps;is++){
									if(l>bdtlimit[is]){
										if(j!=0 && j!=2 && j!=3)
											c_b[is]+=sc[j]*wT_*wPU_*1.0;
										if(j==0)
											c_b[is]+=sc[j]*wT_*wPU_*1.0*SFtt;
										if(j==3){
											if(eventFlavor!=5)
												c_b[is]+=sc[j]*wT_*wPU_*1.0*SFZl;
											if(eventFlavor==5)
												c_b[is]+=sc[j]*wT_*wPU_*1.0*SFZh;
										}
										if(j==2){
											if(eventFlavor!=5)
												c_b[is]+=sc[j]*wT_*wPU_*1.0*SFZl;
											if(eventFlavor==5)
												c_b[is]+=sc[j]*wT_*wPU_*1.0*SFZh;
										}
									}
								}
					} //loop on signal
				} //loop on N
				
				vector <double> punzi;
				for(int is=0; is<bdtsteps;is++)
					punzi.push_back(counter[is]/(1.5+TMath::Sqrt(c_b[is])+0.2*c_b[is]));
				
				double maxpunzi=0;
				int maxis=-999;
				for(int is=0; is<bdtsteps;is++){
					if(punzi[is]>maxpunzi){
						maxpunzi=punzi[is];
						maxis=is;
					}
				}
				cout<<"**** maxpunzi is "<<maxpunzi<<"  at  "<<bdtlimit[maxis]<<"   signal "<<counter[maxis]<<"   and bkg "<<c_b[maxis] << endl;
			} //kbdt 
		} //evalPunzi
		
		TH1F* hTOP= new TH1F(TString::Format("hTOP%d",ivar).Data(),
							 TString::Format("hTOP%d",ivar).Data(),
							 nbins[ivar],bl[ivar],bh[ivar]);
		TH1F* hVV= new TH1F(TString::Format("hVV%d",ivar).Data(),
							TString::Format("hVV%d",ivar).Data(),
							nbins[ivar],bl[ivar],bh[ivar]);
		
		TH1F* htot= new TH1F(TString::Format("htot%d",ivar).Data(),
							 TString::Format("htot%d",ivar).Data(),
							 nbins[ivar],bl[ivar],bh[ivar]);
		hTOP->Sumw2();
		hVV->Sumw2();
		htot->Sumw2();
		
		hTOP->Add(h[ivar][0]);
		hTOP->Add(h[ivar][1]);
		hTOP->Add(h[ivar][2]);
		
		hVV->Add(h[ivar][3]);
		hVV->Add(h[ivar][4]);
		hVV->Add(h[ivar][5]);
		
		
		
		
		//    for(unsigned int j=0;j<t.size();j++){
		//   h[ivar][j]->SetLineColor(j+2);
		//   h[ivar][j]->Draw("hist e same");
		//   if(h[ivar][j]->GetMaximum()>max){
		//max=h[ivar][j]->GetMaximum();
		//  }
		// }
		max=hs[ivar]->GetMaximum();
		cout<<"signal "<<hs[ivar]->Integral(0,-1)<<endl;
		
		
		
		double B=0;
		for(int j=0;j<Nbkg-1;j++){
			if((j!=6) && (j!=7))
				B+=h[ivar][j]->Integral(0,-1);
			if (j==0) cout<<"TTJets "<<B<<endl;
			
		}
		B+= hL[ivar]->Integral(0,-1);
		B+= hH[ivar]->Integral(0,-1);
		cout<<"Background "<<B<<endl;
		
		hTOP->SetLineColor(kGreen+2);
		hTOP->Draw("hist e same");
		if(hTOP->GetMaximum()>max){
			max=hTOP->GetMaximum();
		}
		hVV->SetLineColor(kMagenta);
		hVV->Draw("hist e same");
		if(hVV->GetMaximum()>max){
			max=hVV->GetMaximum();
		}
		hL[ivar]->SetLineColor(kBlue);
		hL[ivar]->Draw("hist e same");
		cout<<"ZL "<<hL[ivar]->Integral(0,-1)<<endl;
		if(hL[ivar]->GetMaximum()>max){
			max=hL[ivar]->GetMaximum();
		}   
		hH[ivar]->SetLineColor(kRed);
		hH[ivar]->Draw("hist e same");
		if(hH[ivar]->GetMaximum()>max){
			max=hH[ivar]->GetMaximum();
		}
		cout<<"ZH "<<hH[ivar]->Integral(0,-1)<<endl;
		
		htot->Add(hTOP);
		htot->Add(hVV);
		htot->Add(hL[ivar]);
		htot->Add(hH[ivar]);
		htot->SetLineColor(kRed);
		htot->SetLineStyle(2);
		htot->Draw("hist e same");
		if(htot->GetMaximum()>max){
			max=htot->GetMaximum();
		}
		/*  hdata[ivar]->Draw("e same");
		 if(hdata[ivar]->GetMaximum()>max){
		 max=hdata[ivar]->GetMaximum();
		 }
		 */
		hs[ivar]->GetYaxis()->SetRangeUser(0, 1.1*max);
		TLegend* lg;
		//right legend  
		if(ivar == kueta || ivar==kDeta  || 
		   ivar==kmet || ivar== kNaj  ||
		   ivar==kcsv2 ||ivar== kpTjj 
		   || ivar==kl1pt || ivar==kl2pt
		   || ivar==kM || ivar==kDptHZ  ||
		   ivar==hdrjj  || ivar==kbdt || ivar==kdeltaPhiZMET  
		   || ivar==kAplanarity || ivar==kSphericity 
		   || ivar==kCircularity 
		   ||ivar==kMte   || ivar==kMtmu) 
			lg=new TLegend(.6,.6,.89,.89);
		//left lg  
		if(ivar== kDphi || ivar==kcentrality ||  ivar==kjchf1  || ivar==kjchf2
		   || ivar==kIsotropy || ivar== kHt  ||  ivar==kcsv1  || ivar==kpTZ ) 
			lg=new TLegend(.6-0.49,.6,.89-0.49,.89);
		//top central
		if(ivar==kangleEMU ) lg=new TLegend(.6-.3,.6,.89-.3,.89);
		
		
		//lg->AddEntry(hdata[ivar], "data", "PL");
		lg->AddEntry(hs[ivar],"signal ZH (m_{H}=115 GeV) ","F"); 
		lg->AddEntry(hTOP," TT+ sT","L");
		lg->AddEntry(hVV,"VV","L");
		lg->AddEntry(hL[ivar],"Z+light","L");
		lg->AddEntry(hH[ivar],"Z+bb","L");
		lg->AddEntry(htot,"All MC","L");
		lg->Draw();
		c[ivar]->Print(TString::Format((dir+"%s.gif").c_str(),var[ivar].c_str()).Data());
		
		
		
		
	} //loop on variables
	
	cout << "made it out of variable loop" << endl;
	
	ch.push_back("TT"); //0
	ch.push_back("TtW"); //5	//1
	ch.push_back("TtWb"); //6	//2
	ch.push_back("ZZ"); //7	//3
	ch.push_back("WW"); //8	//4
	ch.push_back("WZ");  //9	//5
	ch.push_back("ZJH"); //10		//6
	ch.push_back("ZJL"); //11	//7
	ch.push_back("WJ");  //12	//8
	
	//add all single tops together
	h[kbdt][1]->Add(h[kbdt][2]);
	// add all VV together
	h[kbdt][3]->Add(h[kbdt][4]);
	h[kbdt][3]->Add(h[kbdt][5]);
	
	//add all single tops together
	h_btag_up[1]->Add(h_btag_up[2]);
	h_btag_down[1]->Add(h_btag_down[2]);
	// add all VV together
	h_btag_up[3]->Add(h_btag_up[4]);
	h_btag_up[3]->Add(h_btag_up[5]);
	h_btag_down[3]->Add(h_btag_down[4]);
	h_btag_down[3]->Add(h_btag_down[5]);
	
	//add all single tops together
	h_stat_up[1]->Add(h_stat_up[2]);
	h_stat_down[1]->Add(h_stat_down[2]);
	// add all VV together
	h_stat_up[3]->Add(h_stat_up[4]);
	h_stat_up[3]->Add(h_stat_up[5]);
	h_stat_down[3]->Add(h_stat_down[4]);
	h_stat_down[3]->Add(h_stat_down[5]);
	
	
	//here do changes on VV
	//try to remove 
	/*
	 if(eVV<0.01){
	 eVV=0;
	 e_eVV=0;
	 for(int i=0;i<h[kbdt][7]->GetNbinsX()+1;i++){
	 h[kbdt][7]->SetBinContent(i, 0);
	 h_btag_up[7]->SetBinContent(i, 0);
	 h_btag_down[7]->SetBinContent(i, 0);
	 h_stat_up[7]->SetBinContent(i, 0);
	 h_stat_down[7]->SetBinContent(i, 0);
	 }
	 }
	 */
	 
	 cout << "histos combined" << endl;
	
	string bs;
	if(n2dbins==-1){
		if(boost==kboosted) bs="_boosted";
		if(boost==noboost) bs="_noboosted";
		if(boost==klowboosted) bs="_lowboosted";
		if(boost==khighboosted) bs="_highboosted";  
	} else{
		if(boost==kboosted) bs=TString::Format("_boosted%d",n2dbins);
		if(boost==klowboosted) bs=TString::Format("_lowboosted%d",n2dbins);
		if(boost==khighboosted) bs=TString::Format("_highboosted%d",n2dbins);  
	}
	
	cout << "bs decided upon" << endl;
	
	string shapefile;
	shapefile = TString::Format("dataUpdatedShape%s.root",bs.c_str()).Data();
	
	
    if(!EvalPunzi){
		//fudge, use empty histogram for data for now  
		makeShapeHistos(shapefile, hs[kbdt], h[kbdt][8], hL[kbdt], hH[kbdt], h[kbdt][0], h[kbdt][1], h[kbdt][3], h[kbdt][9], hdata[kbdt]);
		//		makeShapeHistos(shapefile, hs[kbdt], h[kbdt][12], h[kbdt][12], hL[kbdt], hH[kbdt], h[kbdt][0], h[kbdt][1], h[kbdt][7], h[kbdt][13], hdata[kbdt]);
		cout << "shape histos made" << endl;
		
		//addSystematics with shape <- only 
		AddSysHistos(shapefile, hs_btag_up, hs_btag_down, h_btag_up[8], 
					 h_btag_down[8], hL_btag_up, hL_btag_down, hH_btag_up, 
					 hH_btag_down, h_btag_up[0], h_btag_down[0], h_btag_up[1],
					 h_btag_down[1], h_btag_up[3], h_btag_down[3], h_btag_up[9],
					 h_btag_down[9], "CMSeffbUp","CMSeffbDown");
		/*    AddSysHistos(shapefile, hs_stat_up, hs_stat_down, h_stat_up[8], 
		 h_stat_down[8], hL_stat_up, hL_stat_down, hH_stat_up, 
		 hH_stat_down, h_stat_up[0], h_stat_down[0], h_stat_up[1],
		 h_stat_down[1], h_stat_up[7], h_stat_down[7], h_stat_up[9],
		 h_stat_down[9], TString::Format("CMS_vhbb_stats_Zemu%sUp",bs.c_str()).Data(),
		 TString::Format("CMS_vhbb_stats_Zemu%sDown",bs.c_str()).Data());
		 */
		cout << "systematics added to histos" << endl;
		
		float e_data = 0.0;
		makeDataCardShape(eVH, eWj, eZjLF, eZjHF, eTT,  es_Top, eVV,
						  e_eVH, e_eWj, e_eZjLF, e_eZjHF, e_eTT, e_es_Top, e_eVV, e_data,
						  //					  e_eVH, e_eWjLF, e_eWjHF, e_eZjLF, e_eZjHF, e_eTT, e_es_Top, e_eVV, hdata[kbdt]->Integral(0,-1),
						  bdtcut, shapefile, bs, n2dbins);
	}
	
	cout<<"last line plot variables" << endl;
	
}
