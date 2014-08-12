{

  gROOT->ProcessLine(".L libTools.so");
  gROOT->ProcessLine(".L ScanChain.C+");

  TChain *ch = new TChain("Events"); 
  ch->Add("/hadoop/cms/store/user/jgran/HLT_Study_jgran/QCD_Pt-80to120_Tune4C_13TeV_pythia8/merged/*.root");
  ch->Add("/hadoop/cms/store/user/jgran/HLT_Study_jgran/QCD_Pt-120to170_Tune4C_13TeV_pythia8/merged/*.root");
  ch->Add("/hadoop/cms/store/user/jgran/HLT_Study_jgran/QCD_Pt-170to300_Tune4C_13TeV_pythia8/merged/*.root");
  ch->Add("/hadoop/cms/store/user/jgran/HLT_Study_jgran/QCD_Pt-300to470_Tune4C_13TeV_pythia8/merged/*.root");
  ch->Add("/hadoop/cms/store/user/jgran/HLT_Study_jgran/QCD_Pt-470to600_Tune4C_13TeV_pythia8/merged/*.root");
  ch->Add("/hadoop/cms/store/user/jgran/HLT_Study_jgran/QCD_Pt-600to800_Tune4C_13TeV_pythia8/merged/*.root");
  ch->Add("/hadoop/cms/store/user/jgran/HLT_Study_jgran/QCD_Pt-800to1000_Tune4C_13TeV_pythia8/merged/*.root");
  ch->Add("/home/users/jgran/HLT/CMSSW_7_1_3/src/HLTStudy/HLTBabyMaker/merging/merged/*.root");
  ch->Add("/hadoop/cms/store/user/jgran/HLT_Study_jgran/QCD_Pt-1400to1800_Tune4C_13TeV_pythia8/merged/*.root");
  ScanChain(ch); 
}
