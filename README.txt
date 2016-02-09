cmsrel CMSSW_8_0_0_pre5
cd CMSSW_8_0_0_pre5/src
cmsenv
git clone git@github.com:olivito/HLTBabyMaker.git HLTStudy/HLTBabyMaker/
scram b
cd HLTStudy/HLTBabyMaker
cmsRun python/hlt_tk1b.py
