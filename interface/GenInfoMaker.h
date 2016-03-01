#ifndef GenInfoMaker_H
#define GenInfoMaker_H

// system include files
#include <memory>
#include <string>
#include <vector>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"


class GenInfoMaker : public edm::EDProducer {
public:
    explicit GenInfoMaker (const edm::ParameterSet&);
    ~GenInfoMaker();

private:
    virtual void beginJob() ;
    virtual void produce(edm::Event&, const edm::EventSetup&);
    virtual void endJob() ;

    // ----------member data ---------------------------
  //    edm::EDGetTokenT<edm::View<PileupSummaryInfo> > pileupSummaryToken;
    edm::EDGetTokenT<std::vector<PileupSummaryInfo> > pileupSummaryToken;
    edm::EDGetTokenT<GenEventInfoProduct> genEvtInfoToken;  
};


#endif
