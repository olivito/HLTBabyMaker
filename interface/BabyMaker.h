#ifndef BabyMaker_H
#define BabyMaker_H

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

class BabyMaker : public edm::EDProducer {
public:
    explicit BabyMaker (const edm::ParameterSet&);
    ~BabyMaker();

private:
    virtual void beginJob() ;
    virtual void produce(edm::Event&, const edm::EventSetup&);
    virtual void endJob() ;

    // ----------member data ---------------------------
    edm::InputTag pfJetsInputTag;
    edm::InputTag pfMetInputTag;
    edm::InputTag pfHTInputTag;
    edm::InputTag caloHTInputTag;
    edm::InputTag hemInputTag;
    edm::InputTag genJetsInputTag;
    edm::InputTag caloJetsInputTag;
    
};


#endif
