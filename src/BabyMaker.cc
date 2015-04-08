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

    produces<std::vector<float> > ("calojetspt").setBranchAlias("calojets_pt");
    produces<std::vector<float> > ("calojetseta").setBranchAlias("calojets_eta");
    produces<std::vector<float> > ("calojetsphi").setBranchAlias("calojets_phi");

    produces<float> ("metpt").setBranchAlias("met_pt");
    produces<float> ("metphi").setBranchAlias("met_phi");

    produces<float> ("calometpt").setBranchAlias("calomet_pt");
    produces<float> ("calometphi").setBranchAlias("calomet_phi");

    produces<float> ("genmetpt").setBranchAlias("genmet_pt");
    produces<float> ("genmetphi").setBranchAlias("genmet_phi");

    produces<float> ("pfht").setBranchAlias("pf_ht");
    produces<float> ("caloht").setBranchAlias("calo_ht");

    produces<float> ("pfmhtpt").setBranchAlias("pf_mht_pt");
    produces<float> ("pfmhteta").setBranchAlias("pf_mht_eta");
    produces<float> ("pfmhtphi").setBranchAlias("pf_mht_phi");

    produces<float> ("calomhtpt").setBranchAlias("calo_mht_pt");
    produces<float> ("calomhteta").setBranchAlias("calo_mht_eta");
    produces<float> ("calomhtphi").setBranchAlias("calo_mht_phi");

/*
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
*/

    produces<float> ("scale1fb").setBranchAlias("scale1fb");

    pfJetsInputTag = iConfig.getParameter<edm::InputTag>("pfJetsInputTag_");
    pfMetInputTag = iConfig.getParameter<edm::InputTag>("pfMetInputTag_");
    pfHTInputTag = iConfig.getParameter<edm::InputTag>("pfHTInputTag_");
    caloMetInputTag = iConfig.getParameter<edm::InputTag>("caloMetInputTag_");
    caloHTInputTag = iConfig.getParameter<edm::InputTag>("caloHTInputTag_");
    //hemInputTag = iConfig.getParameter<edm::InputTag>("hemInputTag_");
    genJetsInputTag = iConfig.getParameter<edm::InputTag>("genJetsInputTag_");
    caloJetsInputTag = iConfig.getParameter<edm::InputTag>("caloJetsInputTag_");
    genMETInputTag = iConfig.getParameter<edm::InputTag>("genMETInputTag_");
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

    std::auto_ptr<std::vector<float> > calojets_pt   (new std::vector<float>);
    std::auto_ptr<std::vector<float> > calojets_eta  (new std::vector<float>);
    std::auto_ptr<std::vector<float> > calojets_phi  (new std::vector<float>);

    std::auto_ptr<float> met_pt   (new float);
    std::auto_ptr<float> met_phi  (new float);

    std::auto_ptr<float> calomet_pt   (new float);
    std::auto_ptr<float> calomet_phi  (new float);

    std::auto_ptr<float> genmet_pt   (new float);
    std::auto_ptr<float> genmet_phi  (new float);

    std::auto_ptr<float> pf_ht   (new float);
    std::auto_ptr<float> calo_ht   (new float);

    std::auto_ptr<float> pf_mht_pt   (new float);
    std::auto_ptr<float> pf_mht_eta   (new float);
    std::auto_ptr<float> pf_mht_phi   (new float);

    std::auto_ptr<float> calo_mht_pt   (new float);
    std::auto_ptr<float> calo_mht_eta   (new float);
    std::auto_ptr<float> calo_mht_phi   (new float);

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

    std::auto_ptr<float> scale1fb     (new float);

    edm::Handle<edm::View<reco::PFJet> > jet_h;
    iEvent.getByLabel(pfJetsInputTag, jet_h);

    edm::Handle<edm::View<reco::MET> > pfht_h;
    iEvent.getByLabel(pfHTInputTag, pfht_h);

    edm::Handle<edm::View<reco::MET> > caloht_h;
    iEvent.getByLabel(caloHTInputTag, caloht_h);

    edm::Handle<edm::View<reco::MET> > met_h;
    iEvent.getByLabel(pfMetInputTag, met_h);

    edm::Handle<edm::View<reco::MET> > calomet_h;
    iEvent.getByLabel(caloMetInputTag, calomet_h);

    edm::Handle<edm::View<reco::GenMET> > genmet_h;
    iEvent.getByLabel(genMETInputTag, genmet_h);

    //edm::Handle<std::vector<LorentzVector_> > hem_h;
    //iEvent.getByLabel(hemInputTag, hem_h);

    edm::Handle<edm::View<reco::GenJet> > genjet_h;
    iEvent.getByLabel(genJetsInputTag, genjet_h);

    edm::Handle<edm::View<reco::CaloJet> > calojet_h;
    iEvent.getByLabel(caloJetsInputTag, calojet_h);


    *scale1fb = 1.0;
    
    std::vector<LorentzVector> goodJets;

    for(edm::View<reco::PFJet>::const_iterator jet_it = jet_h->begin(); jet_it != jet_h->end(); jet_it++){

      if(jet_it->pt() < 35.0) continue;

      pfjets_pt  ->push_back(jet_it->pt());
      pfjets_eta ->push_back(jet_it->eta());
      pfjets_phi ->push_back(jet_it->phi());

      if(jet_it->pt() < 35.0) continue;
      if(fabs(jet_it->eta()) > 3.0) continue;
      if(goodJets.size() == 7) continue;
      goodJets.push_back(LorentzVector(jet_it->p4()));
  
    } 

    for(edm::View<reco::GenJet>::const_iterator genjet_it = genjet_h->begin(); genjet_it != genjet_h->end(); genjet_it++){

      if(genjet_it->pt() < 35.0) continue;

      genjets_pt  ->push_back(genjet_it->pt());
      genjets_eta ->push_back(genjet_it->eta());
      genjets_phi ->push_back(genjet_it->phi());

    } 

    for(edm::View<reco::CaloJet>::const_iterator calojet_it = calojet_h->begin(); calojet_it != calojet_h->end(); calojet_it++){

      if(calojet_it->pt() < 35.0) continue;

      calojets_pt  ->push_back(calojet_it->pt());
      calojets_eta ->push_back(calojet_it->eta());
      calojets_phi ->push_back(calojet_it->phi());

    } 

    *met_pt   = (met_h->front()).pt();
    *met_phi  = (met_h->front()).phi();

    *calomet_pt   = (calomet_h->front()).pt();
    *calomet_phi  = (calomet_h->front()).phi();

    *genmet_pt   = (genmet_h->front()).pt();
    *genmet_phi  = (genmet_h->front()).phi();

    *pf_ht    = (pfht_h->front()).sumEt();
    *calo_ht  = (caloht_h->front()).sumEt();

    *pf_mht_pt    = (pfht_h->front()).pt();
    *pf_mht_eta   = (pfht_h->front()).eta();
    *pf_mht_phi   = (pfht_h->front()).phi();

    *calo_mht_pt    = (caloht_h->front()).pt();
    *calo_mht_eta   = (caloht_h->front()).eta();
    *calo_mht_phi   = (caloht_h->front()).phi();

/*
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
*/
    

    iEvent.put(pfjets_pt,   "pfjetspt" );
    iEvent.put(pfjets_eta,  "pfjetseta" );
    iEvent.put(pfjets_phi,  "pfjetsphi" );

    iEvent.put(genjets_pt,   "genjetspt" );
    iEvent.put(genjets_eta,  "genjetseta" );
    iEvent.put(genjets_phi,  "genjetsphi" );

    iEvent.put(calojets_pt,   "calojetspt" );
    iEvent.put(calojets_eta,  "calojetseta" );
    iEvent.put(calojets_phi,  "calojetsphi" );

    iEvent.put(met_pt,   "metpt" );
    iEvent.put(met_phi,  "metphi" );

    iEvent.put(calomet_pt,   "calometpt" );
    iEvent.put(calomet_phi,  "calometphi" );

    iEvent.put(genmet_pt,   "genmetpt" );
    iEvent.put(genmet_phi,  "genmetphi" );

    iEvent.put(pf_ht,   "pfht" );
    iEvent.put(calo_ht,   "caloht" );

    iEvent.put(pf_mht_pt,   "pfmhtpt" );
    iEvent.put(pf_mht_eta,  "pfmhteta" );
    iEvent.put(pf_mht_phi,  "pfmhtphi" );

    iEvent.put(calo_mht_pt,   "calomhtpt" );
    iEvent.put(calo_mht_eta,  "calomhteta" );
    iEvent.put(calo_mht_phi,  "calomhtphi" );

/*
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
*/

    iEvent.put(scale1fb,  "scale1fb" );
}

//define this as a plug-in
DEFINE_FWK_MODULE(BabyMaker);
