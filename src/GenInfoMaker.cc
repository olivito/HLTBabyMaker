#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "TMath.h"

#include "HLTStudy/HLTBabyMaker/interface/GenInfoMaker.h"

//____________________________________________________________________________
GenInfoMaker::GenInfoMaker(const edm::ParameterSet& iConfig) {

    //-------------- produces statements  -----------------
  
    produces<float> ("mcweight").setBranchAlias("mc_weight");
    produces<float> ("pthat").setBranchAlias("pt_hat");
    produces<float> ("pthatpumax").setBranchAlias("pt_hat_pumax");

    produces<std::vector<float> > ("TrueNumInteractions").setBranchAlias("TrueNumInteractions");
    produces<std::vector<int> >   ("numPUvertices").setBranchAlias("num_PU_vertices");


    //-------------- consumes statements  -----------------

    //    pileupSummaryToken = consumes<edm::View<PileupSummaryInfo> >(iConfig.getParameter<edm::InputTag>("PileupSummaryInputTag_"));
    pileupSummaryToken = consumes<std::vector<PileupSummaryInfo> >(iConfig.getParameter<edm::InputTag>("PileupSummaryInputTag_"));
    genEvtInfoToken = consumes<GenEventInfoProduct>(iConfig.getParameter<edm::InputTag>("genEvtInfoInputTag_"));
}

//____________________________________________________________________________
GenInfoMaker::~GenInfoMaker() {}

//____________________________________________________________________________
void GenInfoMaker::beginJob() {}

//____________________________________________________________________________
void GenInfoMaker::endJob() {}

//____________________________________________________________________________
// ------------ method called to produce the data  ------------
void GenInfoMaker::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  
    //-------------- branch auto_ptr defs  -----------------

    std::auto_ptr<float> mc_weight     (new float);
    std::auto_ptr<float> pt_hat     (new float);
    std::auto_ptr<float> pt_hat_pumax     (new float);

    std::auto_ptr<std::vector<float> > TrueNumInteractions (new std::vector<float>);
    std::auto_ptr<std::vector<int> > num_PU_vertices       (new std::vector<int>);

    //-------------- getByToken  -----------------

    //    edm::Handle<edm::View<PileupSummaryInfo> > pusummary_h;
    edm::Handle<std::vector<PileupSummaryInfo> > pusummary_h;
    iEvent.getByToken(pileupSummaryToken, pusummary_h);

    edm::Handle<GenEventInfoProduct> genEvtInfo_h;
    iEvent.getByToken(genEvtInfoToken, genEvtInfo_h);

    //-------------- retrieving objects  -----------------
    
    *pt_hat_pumax = 0.;
    if (pusummary_h.isValid()) {
      //      for(edm::View<PileupSummaryInfo>::const_iterator pusummary_it = pusummary_h->begin(); pusummary_it != pusummary_h->end(); pusummary_it++){
      for(std::vector<PileupSummaryInfo>::const_iterator pusummary_it = pusummary_h->begin(); pusummary_it != pusummary_h->end(); pusummary_it++){
	TrueNumInteractions->push_back(pusummary_it->getTrueNumInteractions());
	num_PU_vertices    ->push_back(pusummary_it->getPU_NumInteractions());
	// get largest ptHat for a PU event
	//  based on https://github.com/cms-steam/RemovePileUpDominatedEvents
	//  but not as sophisticated: not using genJets from PU collisions..
	for(const auto& ptHat : pusummary_it->getPU_pT_hats()) {
	  if (ptHat > *pt_hat_pumax) *pt_hat_pumax = ptHat;
	}
      }
    }

    if (genEvtInfo_h.isValid()) {
      *mc_weight = genEvtInfo_h->weight();
      *pt_hat = genEvtInfo_h->qScale();
    }

    //-------------- put statements -----------------


    if (genEvtInfo_h.isValid()) {
      iEvent.put(mc_weight,  "mcweight" );
      iEvent.put(pt_hat,  "pthat" );
    }

    if (pusummary_h.isValid()) {
      iEvent.put(TrueNumInteractions, "TrueNumInteractions");
      iEvent.put(num_PU_vertices, "numPUvertices");
      iEvent.put(pt_hat_pumax, "pthatpumax");
    }
}

 
//define this as a plug-in
DEFINE_FWK_MODULE(GenInfoMaker);
