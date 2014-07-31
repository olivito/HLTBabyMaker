#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/PatCandidates/interface/MET.h"

#include "HLTStudy/HLTBabyMaker/interface/BabyMaker.h"
#include "HLTStudy/HLTBabyMaker/interface/tools.h"
#include "HLTStudy/HLTBabyMaker/interface/MT2.h"
#include "HLTStudy/HLTBabyMaker/interface/hemJet.h"

typedef math::XYZTLorentzVector LorentzVector_;


BabyMaker::BabyMaker(const edm::ParameterSet& iConfig) {

    produces<std::vector<float> > ("pfjetspt").setBranchAlias("pfjets_pt");
    produces<std::vector<float> > ("pfjetseta").setBranchAlias("pfjets_eta");
    produces<std::vector<float> > ("pfjetsphi").setBranchAlias("pfjets_phi");

    produces<std::vector<float> > ("genjetspt").setBranchAlias("genjets_pt");
    produces<std::vector<float> > ("genjetseta").setBranchAlias("genjets_eta");
    produces<std::vector<float> > ("genjetsphi").setBranchAlias("genjets_phi");

    produces<float> ("metpt").setBranchAlias("met_pt");
    produces<float> ("meteta").setBranchAlias("met_eta");
    produces<float> ("metphi").setBranchAlias("met_phi");

    produces<float> ("htpt").setBranchAlias("ht_pt");
    produces<float> ("hteta").setBranchAlias("ht_eta");
    produces<float> ("htphi").setBranchAlias("ht_phi");

    produces<LorentzVector> ("hem1p4hlt").setBranchAlias("hem1_p4_hlt");
    produces<LorentzVector> ("hem2p4hlt").setBranchAlias("hem2_p4_hlt");
    produces<LorentzVector> ("hem1p4snt").setBranchAlias("hem1_p4_snt");
    produces<LorentzVector> ("hem2p4snt").setBranchAlias("hem2_p4_snt");

    produces<float> ("pseudoMT2hlt").setBranchAlias("pseudoMT2_hlt");
    produces<float> ("pseudoMT2snt").setBranchAlias("pseudoMT2_snt");
    produces<float> ("mt2hlt").setBranchAlias("mt2_hlt");
    produces<float> ("mt2snt").setBranchAlias("mt2_snt");

    produces<float> ("dPhiHemhlt").setBranchAlias("dPhiHem_hlt");
    produces<float> ("dPhiHemsnt").setBranchAlias("dPhiHem_snt");
    produces<float> ("dPhijj").setBranchAlias("dPhijj");

    pfJetsInputTag = iConfig.getParameter<edm::InputTag>("pfJetsInputTag_");
    pfMetInputTag = iConfig.getParameter<edm::InputTag>("pfMetInputTag_");
    pfHTInputTag = iConfig.getParameter<edm::InputTag>("pfHTInputTag_");
    hemInputTag = iConfig.getParameter<edm::InputTag>("hemInputTag_");
    genJetsInputTag = iConfig.getParameter<edm::InputTag>("genJetsInputTag_");
}


BabyMaker::~BabyMaker() {}

void  BabyMaker::beginJob() {
}

void BabyMaker::endJob() {
}


// ------------ method called to produce the data  ------------
void BabyMaker::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  
    std::auto_ptr<std::vector<float> > pfjets_pt   (new std::vector<float>);
    std::auto_ptr<std::vector<float> > pfjets_eta  (new std::vector<float>);
    std::auto_ptr<std::vector<float> > pfjets_phi  (new std::vector<float>);

    std::auto_ptr<std::vector<float> > genjets_pt   (new std::vector<float>);
    std::auto_ptr<std::vector<float> > genjets_eta  (new std::vector<float>);
    std::auto_ptr<std::vector<float> > genjets_phi  (new std::vector<float>);

    std::auto_ptr<float> met_pt   (new float);
    std::auto_ptr<float> met_eta  (new float);
    std::auto_ptr<float> met_phi  (new float);

    std::auto_ptr<float> ht_pt   (new float);
    std::auto_ptr<float> ht_eta  (new float);
    std::auto_ptr<float> ht_phi  (new float);

    std::auto_ptr<LorentzVector> hem1_p4_hlt  (new LorentzVector);
    std::auto_ptr<LorentzVector> hem2_p4_hlt  (new LorentzVector);
    std::auto_ptr<LorentzVector> hem1_p4_snt  (new LorentzVector);
    std::auto_ptr<LorentzVector> hem2_p4_snt  (new LorentzVector);

    std::auto_ptr<float> pseudoMT2_hlt  (new float);
    std::auto_ptr<float> pseudoMT2_snt  (new float);
    std::auto_ptr<float> mt2_hlt  (new float);
    std::auto_ptr<float> mt2_snt  (new float);

    std::auto_ptr<float> dPhiHem_hlt  (new float);
    std::auto_ptr<float> dPhiHem_snt  (new float);
    std::auto_ptr<float> dPhijj       (new float);

    edm::Handle<edm::View<reco::PFJet> > jet_h;
    iEvent.getByLabel(pfJetsInputTag, jet_h);

    edm::Handle<edm::View<reco::MET> > ht_h;
    iEvent.getByLabel(pfHTInputTag, ht_h);

    edm::Handle<edm::View<reco::MET> > met_h;
    iEvent.getByLabel(pfMetInputTag, met_h);

    edm::Handle<std::vector<LorentzVector_> > hem_h;
    iEvent.getByLabel(hemInputTag, hem_h);

    edm::Handle<edm::View<reco::GenJet> > genjet_h;
    iEvent.getByLabel(genJetsInputTag, genjet_h);

    
    std::vector<LorentzVector> goodJets;

    for(edm::View<reco::PFJet>::const_iterator jet_it = jet_h->begin(); jet_it != jet_h->end(); jet_it++){

      if(jet_it->pt() < 20.0) continue;

      pfjets_pt  ->push_back(jet_it->pt());
      pfjets_eta ->push_back(jet_it->eta());
      pfjets_phi ->push_back(jet_it->phi());

      if(jet_it->pt() < 40.0) continue;
      if(fabs(jet_it->eta()) > 3.0) continue;
      if(goodJets.size() == 7) continue;
      goodJets.push_back(LorentzVector(jet_it->p4()));
  
    } 

    for(edm::View<reco::GenJet>::const_iterator genjet_it = genjet_h->begin(); genjet_it != genjet_h->end(); genjet_it++){

      if(genjet_it->pt() < 20.0) continue;

      genjets_pt  ->push_back(genjet_it->pt());
      genjets_eta ->push_back(genjet_it->eta());
      genjets_phi ->push_back(genjet_it->phi());

    } 

    *met_pt   = (met_h->front()).pt();
    *met_eta  = (met_h->front()).eta();
    *met_phi  = (met_h->front()).phi();

    *ht_pt   = (ht_h->front()).pt();
    *ht_eta  = (ht_h->front()).eta();
    *ht_phi  = (ht_h->front()).phi();

    if(hem_h->size() > 1){
      *hem1_p4_hlt = hem_h->at(0);
      *hem2_p4_hlt = hem_h->at(1);
    }
    else{
      *hem1_p4_hlt = LorentzVector(0,0,0,0); 
      *hem2_p4_hlt = LorentzVector(0,0,0,0); 
    }
    
    float dPhiHem = DeltaPhi( (*hem1_p4_hlt).phi(), (*hem2_p4_hlt).phi() );
    *pseudoMT2_hlt = (*hem1_p4_hlt).pt() * (*hem2_p4_hlt).pt() * (1 + std::cos(dPhiHem));
    *mt2_hlt = HemMT2(*met_pt, *met_phi, *hem1_p4_hlt, *hem2_p4_hlt);
    *dPhiHem_hlt = dPhiHem;


    if(goodJets.size() > 1){
      
      std::vector<LorentzVector> hemJets = getHemJets(goodJets);

      *hem1_p4_snt = hemJets.at(0);
      *hem2_p4_snt = hemJets.at(1);

      dPhiHem = DeltaPhi( (*hem1_p4_snt).phi(), (*hem2_p4_snt).phi() );
      *pseudoMT2_snt = (*hem1_p4_snt).pt() * (*hem2_p4_snt).pt() * (1 + std::cos(dPhiHem));
      *mt2_snt = HemMT2(*met_pt, *met_phi, *hem1_p4_snt, *hem2_p4_snt);

      *dPhiHem_snt = dPhiHem;
      *dPhijj = DeltaPhi(goodJets.at(0).phi(), goodJets.at(1).phi());
    }
    else{
      *hem1_p4_snt = LorentzVector(0,0,0,0); 
      *hem2_p4_snt = LorentzVector(0,0,0,0); 
      *pseudoMT2_snt = -999.9;
      *mt2_snt = -999.9;
      *dPhijj = -999.9;
    }
    

    iEvent.put(pfjets_pt,   "pfjetspt" );
    iEvent.put(pfjets_eta,  "pfjetseta" );
    iEvent.put(pfjets_phi,  "pfjetsphi" );

    iEvent.put(genjets_pt,   "genjetspt" );
    iEvent.put(genjets_eta,  "genjetseta" );
    iEvent.put(genjets_phi,  "genjetsphi" );

    iEvent.put(met_pt,   "metpt" );
    iEvent.put(met_eta,  "meteta" );
    iEvent.put(met_phi,  "metphi" );

    iEvent.put(ht_pt,   "htpt" );
    iEvent.put(ht_eta,  "hteta" );
    iEvent.put(ht_phi,  "htphi" );

    iEvent.put(hem1_p4_hlt,  "hem1p4hlt" );
    iEvent.put(hem2_p4_hlt,  "hem2p4hlt" );
    iEvent.put(hem1_p4_snt,  "hem1p4snt" );
    iEvent.put(hem2_p4_snt,  "hem2p4snt" );

    iEvent.put(pseudoMT2_hlt,  "pseudoMT2hlt" );
    iEvent.put(pseudoMT2_snt,  "pseudoMT2snt" );
    iEvent.put(mt2_hlt,  "mt2hlt" );
    iEvent.put(mt2_snt,  "mt2snt" );

    iEvent.put(dPhiHem_hlt,  "dPhiHemhlt" );
    iEvent.put(dPhiHem_snt,  "dPhiHemsnt" );
    iEvent.put(dPhijj,  "dPhijj" );
}

//define this as a plug-in
DEFINE_FWK_MODULE(BabyMaker);
