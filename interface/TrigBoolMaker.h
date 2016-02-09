#ifndef TrigBoolMaker_H
#define TrigBoolMaker_H

// system include files
#include <memory>
#include <string>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/Common/interface/TriggerResults.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"

class TrigBoolMaker : public edm::EDProducer {
public:
    explicit TrigBoolMaker (const edm::ParameterSet&);
    ~TrigBoolMaker();

private:
    virtual void beginJob() ;
    virtual void beginRun(edm::Run const &, edm::EventSetup const&);
    virtual void produce(edm::Event&, const edm::EventSetup&);
    virtual void endJob() ;
    bool analyzeTrigger(const edm::Event& iEvent, const edm::EventSetup& iSetup, const std::string& triggerName);

    // ----------member data ---------------------------
    edm::EDGetTokenT<edm::TriggerResults> triggerResultsToken;
    std::string processName;
  
    edm::Handle<edm::TriggerResults>           triggerResultsHandle;
    HLTConfigProvider hltConfig;

};


#endif
