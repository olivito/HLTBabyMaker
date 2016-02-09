#include <iostream>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/Common/interface/TriggerNames.h"
#include "FWCore/Common/interface/TriggerResultsByName.h"

#include "HLTStudy/HLTBabyMaker/interface/TrigBoolMaker.h"

//____________________________________________________________________________
TrigBoolMaker::TrigBoolMaker(const edm::ParameterSet& iConfig) {

    //-------------- produces statements  -----------------
    produces<int> ("HLTL1HTT125").setBranchAlias("HLT_L1HTT125");
    produces<int> ("HLTL1HTT150").setBranchAlias("HLT_L1HTT150");
    produces<int> ("HLTL1HTT175").setBranchAlias("HLT_L1HTT175");
    produces<int> ("HLTL1HTT200").setBranchAlias("HLT_L1HTT200");

    produces<int> ("HLTL1ETM50").setBranchAlias("HLT_L1ETM50");
    produces<int> ("HLTL1ETM60").setBranchAlias("HLT_L1ETM60");
    produces<int> ("HLTL1ETM70").setBranchAlias("HLT_L1ETM70");

    produces<int> ("HLTPFHT350PFMET100").setBranchAlias("HLT_PFHT350_PFMET100");
    produces<int> ("HLTPFHT800").setBranchAlias("HLT_PFHT800");
    produces<int> ("HLTPFMETNoMu90PFMHTNoMu90IDTight").setBranchAlias("HLT_PFMETNoMu90_PFMHTNoMu90_IDTight");
    produces<int> ("HLTPFMET90PFMHT90IDTight").setBranchAlias("HLT_PFMET90_PFMHT90_IDTight");
    
    //-------------- consumes statements  -----------------

    triggerResultsToken = consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("triggerResultsInputTag_"));

    processName = iConfig.getParameter<std::string>("processName_");
}

//____________________________________________________________________________
TrigBoolMaker::~TrigBoolMaker() {}

//____________________________________________________________________________
void TrigBoolMaker::beginJob() {}

//____________________________________________________________________________
void TrigBoolMaker::beginRun(edm::Run const & iRun, edm::EventSetup const& iSetup) {
  
  using namespace std;
  using namespace edm;

  bool changed(true);
  hltConfig.init(iRun,iSetup,processName,changed);
}

//____________________________________________________________________________
void TrigBoolMaker::endJob() {}

//____________________________________________________________________________
// ------------ method called to produce the data  ------------
void TrigBoolMaker::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  
    //-------------- branch auto_ptr defs  -----------------
  
    std::auto_ptr<int> HLT_L1HTT125     (new int);
    std::auto_ptr<int> HLT_L1HTT150     (new int);
    std::auto_ptr<int> HLT_L1HTT175     (new int);
    std::auto_ptr<int> HLT_L1HTT200     (new int);
    
    std::auto_ptr<int> HLT_L1ETM50     (new int);
    std::auto_ptr<int> HLT_L1ETM60     (new int);
    std::auto_ptr<int> HLT_L1ETM70     (new int);
    
    std::auto_ptr<int> HLT_PFHT350_PFMET100     (new int);
    std::auto_ptr<int> HLT_PFHT800     (new int);
    std::auto_ptr<int> HLT_PFMETNoMu90_PFMHTNoMu90_IDTight     (new int);
    std::auto_ptr<int> HLT_PFMET90_PFMHT90_IDTight     (new int);
    
    //-------------- getByToken  -----------------

    iEvent.getByToken(triggerResultsToken, triggerResultsHandle);

    //-------------- retrieving decisions  -----------------
    
    *HLT_L1HTT125 = analyzeTrigger(iEvent,iSetup,"HLT_L1HTT125");
    *HLT_L1HTT150 = analyzeTrigger(iEvent,iSetup,"HLT_L1HTT150");
    *HLT_L1HTT175 = analyzeTrigger(iEvent,iSetup,"HLT_L1HTT175");
    *HLT_L1HTT200 = analyzeTrigger(iEvent,iSetup,"HLT_L1HTT200");
    
    *HLT_L1ETM50 = analyzeTrigger(iEvent,iSetup,"HLT_L1ETM50");
    *HLT_L1ETM60 = analyzeTrigger(iEvent,iSetup,"HLT_L1ETM60");
    *HLT_L1ETM70 = analyzeTrigger(iEvent,iSetup,"HLT_L1ETM70");
    
    *HLT_PFHT350_PFMET100 = analyzeTrigger(iEvent,iSetup,"HLT_PFHT350_PFMET100_v1");
    *HLT_PFHT800 = analyzeTrigger(iEvent,iSetup,"HLT_PFHT800_v2");
    *HLT_PFMETNoMu90_PFMHTNoMu90_IDTight = analyzeTrigger(iEvent,iSetup,"HLT_PFMETNoMu90_PFMHTNoMu90_IDTight_v1");
    *HLT_PFMET90_PFMHT90_IDTight = analyzeTrigger(iEvent,iSetup,"HLT_PFMET90_PFMHT90_IDTight_v2");

    //-------------- put statements -----------------

    iEvent.put(HLT_L1HTT125,   "HLTL1HTT125" );
    iEvent.put(HLT_L1HTT150,   "HLTL1HTT150" );
    iEvent.put(HLT_L1HTT175,   "HLTL1HTT175" );
    iEvent.put(HLT_L1HTT200,   "HLTL1HTT200" );
    
    iEvent.put(HLT_L1ETM50,   "HLTL1ETM50" );
    iEvent.put(HLT_L1ETM60,   "HLTL1ETM60" );
    iEvent.put(HLT_L1ETM70,   "HLTL1ETM70" );

    iEvent.put(HLT_PFHT350_PFMET100,   "HLTPFHT350PFMET100" );
    iEvent.put(HLT_PFHT800,   "HLTPFHT800" );
    iEvent.put(HLT_PFMETNoMu90_PFMHTNoMu90_IDTight,   "HLTPFMETNoMu90PFMHTNoMu90IDTight" );
    iEvent.put(HLT_PFMET90_PFMHT90_IDTight,   "HLTPFMET90PFMHT90IDTight" );

}

//____________________________________________________________________________
bool TrigBoolMaker::analyzeTrigger(const edm::Event& iEvent, const edm::EventSetup& iSetup, const std::string& triggerName) {
  
  using namespace std;
  using namespace edm;
  using namespace trigger;

  const bool verbose = false;
  
  if (verbose) cout << endl;

  const unsigned int ntrigs(hltConfig.size());
  const unsigned int triggerIndex(hltConfig.triggerIndex(triggerName));
  assert(triggerIndex==iEvent.triggerNames(*triggerResultsHandle).triggerIndex(triggerName));

  // abort on invalid trigger name
  if (triggerIndex>=ntrigs) {
    cout << "OverlapAnalyzer::analyzeTrigger: path "
	 << triggerName << " - not found!" << endl;
    return false;
  }

  if (verbose) {
    cout << "OverlapAnalyzer::analyzeTrigger: path "
	 << triggerName << " [" << triggerIndex << "]" << endl;
  }
  // modules on this trigger path
  const unsigned int m(hltConfig.size(triggerIndex));
  const vector<string>& moduleLabels(hltConfig.moduleLabels(triggerIndex));

  bool wasRun = triggerResultsHandle->wasrun(triggerIndex);
  bool accept = triggerResultsHandle->accept(triggerIndex);
  bool error = triggerResultsHandle->error(triggerIndex);
  const unsigned int moduleIndex(triggerResultsHandle->index(triggerIndex));
  // Results from TriggerResults product
  if (verbose) {
    cout << " Trigger path status:"
	 << " WasRun=" << wasRun
	 << " Accept=" << accept
	 << " Error =" << error
	 << endl;
    cout << " Last active module - label/type: "
	 << moduleLabels[moduleIndex] << "/" << hltConfig.moduleType(moduleLabels[moduleIndex])
	 << " [" << moduleIndex << " out of 0-" << (m-1) << " on this path]"
	 << endl;
  }
  assert (moduleIndex<m);

  if (!wasRun || !accept || error) return false;
  return true;
}

 
//define this as a plug-in
DEFINE_FWK_MODULE(TrigBoolMaker);
