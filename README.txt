cmsrel CMSSW_7_1_3
cd CMSSW_7_1_3/src
cmsenv
git clone https://github.com/jgran/HLTBabyMaker HLTStudy/HLTBabyMaker/
scram b
cd HLTStudy/HLTBabyMaker
cmsRun python/hlt_tk1b.py
