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

#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/JetReco/interface/CaloJet.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/METReco/interface/GenMET.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

class BabyMaker : public edm::EDProducer {
public:
    explicit BabyMaker (const edm::ParameterSet&);
    ~BabyMaker();

private:
    virtual void beginJob() ;
    virtual void produce(edm::Event&, const edm::EventSetup&);
    virtual void endJob() ;

    // ----------member data ---------------------------
    edm::EDGetTokenT<edm::View<reco::PFJet> > pfJetsToken;
    edm::EDGetTokenT<edm::View<reco::MET> > pfMetToken;
    edm::EDGetTokenT<edm::View<reco::MET> > pfHTToken;
    edm::EDGetTokenT<edm::View<reco::MET> > pfHTTightIDToken;
    edm::EDGetTokenT<edm::View<reco::MET> > caloMetToken;
    edm::EDGetTokenT<edm::View<reco::MET> > caloHTToken;
    edm::EDGetTokenT<edm::View<reco::GenJet> > genJetsToken;
    edm::EDGetTokenT<edm::View<reco::CaloJet> > caloJetsToken;
    edm::EDGetTokenT<edm::View<reco::GenMET> > genMETToken;
    edm::EDGetTokenT<edm::View<PileupSummaryInfo> > pileupSummaryToken;
    
};


#endif
