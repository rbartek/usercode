//Form: k factor (when needed) * production cross section * BF1 * BF2 (BFs when needed). All results in fb 

//ZH to llbb BF = 0.101 
const double xsecbfZH110 = 0.5869*0.745*0.101*1000;
const double xsecbfZH115 = 0.5117*0.704*0.101*1000;
const double xsecbfZH120 = 0.4483*0.648*0.101*1000;
const double xsecbfZH125 = 0.3943*0.577*0.101*1000;            
const double xsecbfZH130 = 0.3473*0.493*0.101*1000;           
const double xsecbfZH135 = 0.3074*0.403*0.101*1000;            

//Z to Leptons
const double xsecbfZJL  = 3503.71*1000; //Inclusive
const double xsecbfZJH  = 40.5005*1000; //ZPt>100
const double xsecbfZJ50_70 = 111.406*1000; //Z PT 50-70
const double xsecbfZJ70_100 = 62.128*1000; //Z PT 70-100
const double xsecbfZJ200_400HT = 23.4332878559322*1000; //HT 200-400
const double xsecbfZJ400HT = 3.35643541016949*1000; //HT > 400
const double xsecbfZoneJet  = 666.298749152542*1000; //DYJetsToLL1
const double xsecbfZtwoJet  = 214.973393220339*1000; //DYJetsToLL2
const double xsecbfZthreeJet  = 60.6913833898305*1000; //DYJetsToLL3
const double xsecbfZfourJet  = 27.3645689491525*1000; //DYJetsToLL4
const double xsecbfWJ   = 36257.2*1000; //WJetsToLNu_TuneZ2Star_8TeV-madgraph-tarball/Summer12_DR53X-PU_S10_START53_V7A-v*/AODSIM
const double xsecbfWZ   = 32.32*1000;
const double xsecbfWW   = 57.11*1000;
const double xsecbfZZ   = 8.297*1000; 
const double xsecbfTT   = 234*1000; 
const double xsecbfTs   = 3.79*1000; 
const double xsecbfTt   = 56.4*1000; 
const double xsecbfTtW  = 11.1*1000; 
const double xsecbfTsb  = 1.76*1000; 
const double xsecbfTtb  = 30.7*1000; 
const double xsecbfTtWb = 11.1*1000; 
const double xsecbfQCD120 = 115613.74*1000; 
const double xsecbfQCD170 = 24297.5*1000; 
const double xsecbfQCD300 = 1169.5762*1000; 

const double numWJ = 8.13240960000000000e+07;
const double numZJH = 1.13720900000000000e+06;
const double numZJL = 3.61874080000000000e+07;
const double numTT = 5.84488960000000000e+07;
//TTbar_M60
//TTbar 3803462.25;
const double numTs = 1111111111;
const double numTt = 1111111111111;
const double numTtW = 8.14224000000000000e+05;
const double numTsb = 11111111111;
const double numTtb = 1992339.75;
const double numTtWb = 8.09875000000000000e+05;
const double numWW = 4.22534900000000000e+06;
const double numWZ = 4.26467100000000000e+06;
const double numZH100 = 1111111111;
const double numZH105 = 219999.0;
const double numZH110 = 218777.0;
const double numZH115 = 1.09994200000000000e+06;
const double numZH120 = 219999.0;
const double numZH125 = 1.09993800000000000e+06;
const double numZH125Filtered = 2.26620000000000000e+04;
const double numZH125Fall = 2.16744000000000000e+05;
const double numZH130 = 111111111111;
const double numZH135 = 1090594.0;
const double numZZ = 9.79892300000000000e+06;
const double numQCD120 = 1111111111;
const double numQCD170 = 11111111111;
const double numQCD300 = 111111111111;

const double lumiZH115  = numZH115/  xsecbfZH115  ;     
const double lumiZH125  = numZH125/  xsecbfZH125  ;
//const double lumiZH125Filtered  = numZH125Filtered/  xsecbfZH125  ;     
const double lumiZH125Filtered  = numZH125/  xsecbfZH125  ;     
const double lumiZH125Fall  = numZH125Fall/  xsecbfZH125  ;          
const double lumiZJL = numZJL  /  xsecbfZJL  ;    
const double lumiZJH = numZJH  /  xsecbfZJH  ;    
const double lumiWJ  = numWJ   /  xsecbfWJ  ;     
const double lumiWZ  = numWZ   /  xsecbfWZ  ;     
const double lumiWW  = numWW   /  xsecbfWW  ;     
const double lumiZZ  = numZZ   /  xsecbfZZ  ;     
const double lumiTT  = numTT   /  xsecbfTT  ;     
const double lumiTs  = numTs   /  xsecbfTs  ;     
const double lumiTt  = numTt   /  xsecbfTt  ;     
const double lumiTtW = numTtW  /  xsecbfTtW ;    
const double lumiTsb = numTsb  /  xsecbfTsb ;     
const double lumiTtb = numTtb  /  xsecbfTtb ;      
const double lumiTtWb= numTtWb /  xsecbfTtWb;    
const double lumiQCD120= numQCD120 /  xsecbfQCD120;    
const double lumiQCD170= numQCD170 /  xsecbfQCD170;    
const double lumiQCD300= numQCD300 /  xsecbfQCD300;    

