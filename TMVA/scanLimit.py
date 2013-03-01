import sys
import os
import commands
import string


TTT = 100
EEE = 20
DDD = 2

while DDD < 5:
 while EEE < 120:
  while TTT < 700:
    os.system("sed -i 's/NTrees=TTT/NTrees=%i/g' /home/hep/wilken/BDTSelection/src/UserCode/wilken/TMVA/OptiBDTparm.C" % TTT)
    os.system("sed -i 's/nEventsMin=EEE/nEventsMin=%i/g' /home/hep/wilken/BDTSelection/src/UserCode/wilken/TMVA/OptiBDTparm.C" % EEE)
    os.system("sed -i 's/MaxDepth=DDD/MaxDepth=%i/g' /home/hep/wilken/BDTSelection/src/UserCode/wilken/TMVA/OptiBDTparm.C" % DDD)
    os.chdir("/home/hep/wilken/BDTSelection/src/UserCode/wilken/TMVA/")
    os.system("pwd")
    os.system("eval `scramv1 runtime -sh`")
    os.system("root -l OptiBDTparm.C")
    os.system("cp weights/TMVAClassification_BDT.class.C /home/hep/wilken/CMSSW_5_2_5/src/.")
    os.chdir("/home/hep/wilken/CMSSW_5_2_5/src/")
    os.system("pwd")
#    os.system("cmsenv")
    os.system("eval `scramv1 runtime -sh`")
    os.system(" root -b -q runSeparation.C | tee logoutput%i_%i_%i.out"% (TTT,EEE,DDD))
#    os.system("combineCards.py vhbb_ZmmHighPt_8TeV.txt vhbb_ZeeHighPt_8TeV.txt > high.txt")
#    os.system("grep \"The best limit is\" logoutput.out | awk '{print $5}' > cresult")
#    fres = open("cresult", "r")
#    tres = (fres.read().rstrip('\n'))
#    res = float(tres)
    os.chdir("/home/hep/wilken/")
    os.system("pwd")
#    os.system("eval `scramv1 runtime -sh`")
    os.system("grep 'BoostType=AdaBoost:AdaBoostBeta=0.3:SeparationType' /home/hep/wilken/BDTSelection/src/UserCode/wilken/TMVA/OptiBDTparm.C")
##    os.system("echo 'NTrees=%i:nEventsMin=%i:MaxDepth=%i      %f' >> FinalResults.txt" % (TTT,EEE,DDD,res))
    os.system("echo 'NTrees=%i:nEventsMin=%i:MaxDepth=%i' >> FinalResults.txt" % (TTT,EEE,DDD))
    os.system("sed -i 's/NTrees=%i/NTrees=TTT/g' /home/hep/wilken/BDTSelection/src/UserCode/wilken/TMVA/OptiBDTparm.C" % TTT)
    os.system("sed -i 's/nEventsMin=%i/nEventsMin=EEE/g' /home/hep/wilken/BDTSelection/src/UserCode/wilken/TMVA/OptiBDTparm.C" % EEE)
    os.system("sed -i 's/MaxDepth=%i/MaxDepth=DDD/g' /home/hep/wilken/BDTSelection/src/UserCode/wilken/TMVA/OptiBDTparm.C" % DDD)
    TTT = TTT + 100
  EEE = EEE + 20
  TTT = 100
 DDD = DDD + 1
 EEE = 20
 TTT = 100
