//Form: k factor (when needed) * production cross section * BF1 * BF2 (BFs when needed). All results in fb 

//ZH to llbb BF = 0.101 
const double xsecbfZH115 = 0.4107*0.704*0.101*1000;
const double xsecbfZH120 = 0.3598*0.648*0.101*1000;
const double xsecbfZH125 = 0.3158*0.577*0.101*1000;            
const double xsecbfZH130 = 0.2778*0.493*0.101*1000;           
const double xsecbfZH135 = 0.2453*0.403*0.101*1000;            

//Z to Leptons
const double xsecbfZJL  = 3048*1000; 
const double xsecbfZJH  = 1.3*25.8*1000; 
const double xsecbfWJ   = 31314*1000; 
const double xsecbfWZ   = 18.2*1000;
const double xsecbfWW   = 43*1000;
const double xsecbfZZ   = 5.9*1000; 
const double xsecbfTT   = 165*1000; 
const double xsecbfTs   = 3.19*1000; 
const double xsecbfTt   = 41.92*1000; 
const double xsecbfTtW  = 7.87*1000; 
const double xsecbfTsb  = 1.44*1000; 
const double xsecbfTtb  = 22.65*1000; 
const double xsecbfTtWb = 7.87*1000; 

const double numWJ = 9.82263680000000000e+07;
const double numZJH = 1.18055337500000000e+06;
const double numZJL = 3.62730200000000000e+07;
const double numTT = 3.24947580000000000e+07; 
//TTbar_M60
//TTbar 3803462.25;
const double numTs = 268272.21875;
const double numTt = 3.56846950000000000e+06;
const double numTtW = 8.41851250000000000e+05;
const double numTsb = 139683.328125;
const double numTtb = 1992339.75;
const double numTtWb = 3.38248250000000000e+05;
const double numWW = 3.04573050000000000e+06;
const double numWZ = 4.41649050000000000e+06;
const double numZH100 = 219999.0;
const double numZH105 = 219999.0;
const double numZH110 = 218777.0;
const double numZH115 = 1.0397629E+06;
const double numZH120 = 219999.0;
const double numZH125 = 216759.0;
const double numZH130 = 225991.546875;
const double numZH135 = 189999.0;
const double numZZ = 4.19944500000000000e+06;


const double lumiZH115  = numZH115/  xsecbfZH115  ;     
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

