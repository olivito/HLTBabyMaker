#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "TMath.h"

#include "HLTStudy/HLTBabyMaker/interface/BabyMaker.h"

typedef math::XYZTLorentzVector LorentzVector;

//____________________________________________________________________________
BabyMaker::BabyMaker(const edm::ParameterSet& iConfig) {

    //-------------- produces statements  -----------------
  
    produces<std::vector<float> > ("pfjetspt").setBranchAlias("pfjets_pt");
    produces<std::vector<float> > ("pfjetseta").setBranchAlias("pfjets_eta");
    produces<std::vector<float> > ("pfjetsphi").setBranchAlias("pfjets_phi");

    produces<std::vector<float> > ("genjetspt").setBranchAlias("genjets_pt");
    produces<std::vector<float> > ("genjetseta").setBranchAlias("genjets_eta");
    produces<std::vector<float> > ("genjetsphi").setBranchAlias("genjets_phi");

    produces<std::vector<float> > ("calojetspt").setBranchAlias("calojets_pt");
    produces<std::vector<float> > ("calojetseta").setBranchAlias("calojets_eta");
    produces<std::vector<float> > ("calojetsphi").setBranchAlias("calojets_phi");

    produces<float> ("pfmetpt").setBranchAlias("pfmet_pt");
    produces<float> ("pfmetphi").setBranchAlias("pfmet_phi");

    produces<float> ("pfmetnomupt").setBranchAlias("pfmet_nomu_pt");
    produces<float> ("pfmetnomuphi").setBranchAlias("pfmet_nomu_phi");

    produces<float> ("calometpt").setBranchAlias("calomet_pt");
    produces<float> ("calometphi").setBranchAlias("calomet_phi");

    produces<float> ("genmetpt").setBranchAlias("genmet_pt");
    produces<float> ("genmetphi").setBranchAlias("genmet_phi");

    produces<float> ("genht").setBranchAlias("gen_ht");
    
    produces<float> ("pfht").setBranchAlias("pf_ht");
    produces<float> ("caloht").setBranchAlias("calo_ht");

    produces<float> ("pfmhtpt").setBranchAlias("pf_mht_pt");
    produces<float> ("pfmhtphi").setBranchAlias("pf_mht_phi");

    produces<float> ("pfmhtnomutightidpt").setBranchAlias("pf_mht_nomu_tightid_pt");
    produces<float> ("pfmhtnomutightidphi").setBranchAlias("pf_mht_nomu_tightid_phi");

    produces<float> ("pfmhttightidpt").setBranchAlias("pf_mht_tightid_pt");
    produces<float> ("pfmhttightidphi").setBranchAlias("pf_mht_tightid_phi");

    produces<float> ("calomhtpt").setBranchAlias("calo_mht_pt");
    produces<float> ("calomhtphi").setBranchAlias("calo_mht_phi");

    produces<int> ("npixvtx").setBranchAlias("n_pix_vtx");

    produces<float> ("scale1fb").setBranchAlias("scale1fb");

    produces<std::vector<float> > ("TrueNumInteractions").setBranchAlias("TrueNumInteractions");
    produces<std::vector<int> >   ("numPUvertices").setBranchAlias("num_PU_vertices");

    produces<std::vector<float> > ("pfjetsofflinept").setBranchAlias("pfjets_offline_pt");
    produces<std::vector<float> > ("pfjetsofflineeta").setBranchAlias("pfjets_offline_eta");
    produces<std::vector<float> > ("pfjetsofflinephi").setBranchAlias("pfjets_offline_phi");

    produces<float> ("pfofflineht").setBranchAlias("pf_offline_ht");
    produces<float> ("pfmetofflinept").setBranchAlias("pfmet_offline_pt");
    produces<float> ("pfmetofflinephi").setBranchAlias("pfmet_offline_phi");

    //-------------- consumes statements  -----------------

    pfJetsToken = consumes<edm::View<reco::PFJet> >(iConfig.getParameter<edm::InputTag>("pfJetsInputTag_"));
    pfMetToken = consumes<edm::View<reco::MET> >(iConfig.getParameter<edm::InputTag>("pfMetInputTag_"));
    pfMetNoMuToken = consumes<edm::View<reco::MET> >(iConfig.getParameter<edm::InputTag>("pfMetNoMuInputTag_"));
    pfHTToken = consumes<edm::View<reco::MET> >(iConfig.getParameter<edm::InputTag>("pfHTInputTag_"));
    pfHTNoMuTightIDToken = consumes<edm::View<reco::MET> >(iConfig.getParameter<edm::InputTag>("pfHTNoMuTightIDInputTag_"));
    pfHTTightIDToken = consumes<edm::View<reco::MET> >(iConfig.getParameter<edm::InputTag>("pfHTTightIDInputTag_"));
    caloMetToken = consumes<edm::View<reco::MET> >(iConfig.getParameter<edm::InputTag>("caloMetInputTag_"));
    caloHTToken = consumes<edm::View<reco::MET> >(iConfig.getParameter<edm::InputTag>("caloHTInputTag_"));
    genJetsToken = consumes<edm::View<reco::GenJet> >(iConfig.getParameter<edm::InputTag>("genJetsInputTag_"));
    caloJetsToken = consumes<edm::View<reco::CaloJet> >(iConfig.getParameter<edm::InputTag>("caloJetsInputTag_"));
    genMETToken = consumes<edm::View<reco::GenMET> >(iConfig.getParameter<edm::InputTag>("genMETInputTag_"));
    pileupSummaryToken = consumes<edm::View<PileupSummaryInfo> >(iConfig.getParameter<edm::InputTag>("PileupSummaryInputTag_"));
    pixelVerticesToken = consumes<edm::View<reco::Vertex> >(iConfig.getParameter<edm::InputTag>("pixelVerticesInputTag_"));
    
    pfJetsOfflineToken = consumes<edm::View<pat::Jet> >(iConfig.getParameter<edm::InputTag>("pfJetsOfflineInputTag_"));
    pfMetOfflineToken = consumes<edm::View<pat::MET> >(iConfig.getParameter<edm::InputTag>("pfMetOfflineInputTag_"));
}

//____________________________________________________________________________
BabyMaker::~BabyMaker() {}

//____________________________________________________________________________
void BabyMaker::beginJob() {}

//____________________________________________________________________________
void BabyMaker::endJob() {}

//____________________________________________________________________________
// ------------ method called to produce the data  ------------
void BabyMaker::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  
    //-------------- branch auto_ptr defs  -----------------

    std::auto_ptr<std::vector<float> > pfjets_pt   (new std::vector<float>);
    std::auto_ptr<std::vector<float> > pfjets_eta  (new std::vector<float>);
    std::auto_ptr<std::vector<float> > pfjets_phi  (new std::vector<float>);

    std::auto_ptr<std::vector<float> > genjets_pt   (new std::vector<float>);
    std::auto_ptr<std::vector<float> > genjets_eta  (new std::vector<float>);
    std::auto_ptr<std::vector<float> > genjets_phi  (new std::vector<float>);

    std::auto_ptr<std::vector<float> > calojets_pt   (new std::vector<float>);
    std::auto_ptr<std::vector<float> > calojets_eta  (new std::vector<float>);
    std::auto_ptr<std::vector<float> > calojets_phi  (new std::vector<float>);

    std::auto_ptr<float> pfmet_pt   (new float);
    std::auto_ptr<float> pfmet_phi  (new float);

    std::auto_ptr<float> pfmet_nomu_pt   (new float);
    std::auto_ptr<float> pfmet_nomu_phi  (new float);

    std::auto_ptr<float> calomet_pt   (new float);
    std::auto_ptr<float> calomet_phi  (new float);

    std::auto_ptr<float> genmet_pt   (new float);
    std::auto_ptr<float> genmet_phi  (new float);

    std::auto_ptr<float> gen_ht   (new float);
    
    std::auto_ptr<float> pf_ht   (new float);
    std::auto_ptr<float> calo_ht   (new float);

    std::auto_ptr<float> pf_mht_pt   (new float);
    std::auto_ptr<float> pf_mht_phi   (new float);

    std::auto_ptr<float> pf_mht_nomu_tightid_pt   (new float);
    std::auto_ptr<float> pf_mht_nomu_tightid_phi   (new float);

    std::auto_ptr<float> pf_mht_tightid_pt   (new float);
    std::auto_ptr<float> pf_mht_tightid_phi   (new float);

    std::auto_ptr<float> calo_mht_pt   (new float);
    std::auto_ptr<float> calo_mht_phi   (new float);

    std::auto_ptr<int> n_pix_vtx     (new int);

    std::auto_ptr<float> scale1fb     (new float);

    std::auto_ptr<std::vector<float> > TrueNumInteractions (new std::vector<float>);
    std::auto_ptr<std::vector<int> > num_PU_vertices       (new std::vector<int>);
    
    std::auto_ptr<std::vector<float> > pfjets_offline_pt   (new std::vector<float>);
    std::auto_ptr<std::vector<float> > pfjets_offline_eta  (new std::vector<float>);
    std::auto_ptr<std::vector<float> > pfjets_offline_phi  (new std::vector<float>);

    std::auto_ptr<float> pf_offline_ht   (new float);
    std::auto_ptr<float> pfmet_offline_pt   (new float);
    std::auto_ptr<float> pfmet_offline_phi  (new float);

    //-------------- getByToken  -----------------

    edm::Handle<edm::View<reco::PFJet> > jet_h;
    iEvent.getByToken(pfJetsToken, jet_h);

    edm::Handle<edm::View<reco::MET> > pfht_h;
    iEvent.getByToken(pfHTToken, pfht_h);

    edm::Handle<edm::View<reco::MET> > pfhtnomutightid_h;
    iEvent.getByToken(pfHTNoMuTightIDToken, pfhtnomutightid_h);

    edm::Handle<edm::View<reco::MET> > pfhttightid_h;
    iEvent.getByToken(pfHTTightIDToken, pfhttightid_h);

    edm::Handle<edm::View<reco::MET> > caloht_h;
    iEvent.getByToken(caloHTToken, caloht_h);

    edm::Handle<edm::View<reco::MET> > pfmet_h;
    iEvent.getByToken(pfMetToken, pfmet_h);

    edm::Handle<edm::View<reco::MET> > pfmetnomu_h;
    iEvent.getByToken(pfMetNoMuToken, pfmetnomu_h);

    edm::Handle<edm::View<reco::MET> > calomet_h;
    iEvent.getByToken(caloMetToken, calomet_h);

    edm::Handle<edm::View<reco::GenMET> > genmet_h;
    iEvent.getByToken(genMETToken, genmet_h);

    edm::Handle<edm::View<reco::GenJet> > genjet_h;
    iEvent.getByToken(genJetsToken, genjet_h);

    edm::Handle<edm::View<reco::CaloJet> > calojet_h;
    iEvent.getByToken(caloJetsToken, calojet_h);

    edm::Handle<edm::View<PileupSummaryInfo> > pusummary_h;
    iEvent.getByToken(pileupSummaryToken, pusummary_h);

    edm::Handle<edm::View<reco::Vertex> > pixelvertices_h;
    iEvent.getByToken(pixelVerticesToken, pixelvertices_h);

    edm::Handle<edm::View<pat::Jet> > jet_offline_h;
    iEvent.getByToken(pfJetsOfflineToken, jet_offline_h);

    edm::Handle<edm::View<pat::MET> > pfmet_offline_h;
    iEvent.getByToken(pfMetOfflineToken, pfmet_offline_h);

    //-------------- retrieving objects  -----------------

    *scale1fb = 1.0;
    
    std::vector<LorentzVector> goodJets;

    for(edm::View<reco::PFJet>::const_iterator jet_it = jet_h->begin(); jet_it != jet_h->end(); jet_it++){

      if(jet_it->pt() < 40.0) continue;

      pfjets_pt  ->push_back(jet_it->pt());
      pfjets_eta ->push_back(jet_it->eta());
      pfjets_phi ->push_back(jet_it->phi());

      if(jet_it->pt() < 40.0) continue;
      if(fabs(jet_it->eta()) > 3.0) continue;
      if(goodJets.size() == 7) continue;
      goodJets.push_back(LorentzVector(jet_it->p4()));
  
    } 

    *gen_ht = 0;
    if (genjet_h.isValid()) {
      for(edm::View<reco::GenJet>::const_iterator genjet_it = genjet_h->begin(); genjet_it != genjet_h->end(); genjet_it++){

	if(genjet_it->pt() < 40.0) continue;

	genjets_pt  ->push_back(genjet_it->pt());
	genjets_eta ->push_back(genjet_it->eta());
	genjets_phi ->push_back(genjet_it->phi());

	if (fabs(genjet_it->eta()) < 3.0) *gen_ht += genjet_it->pt();
      }
    }

    for(edm::View<reco::CaloJet>::const_iterator calojet_it = calojet_h->begin(); calojet_it != calojet_h->end(); calojet_it++){

      if(calojet_it->pt() < 40.0) continue;

      calojets_pt  ->push_back(calojet_it->pt());
      calojets_eta ->push_back(calojet_it->eta());
      calojets_phi ->push_back(calojet_it->phi());

    }

    *pfmet_pt   = (pfmet_h->front()).pt();
    *pfmet_phi  = (pfmet_h->front()).phi();

    *pfmet_nomu_pt   = (pfmetnomu_h->front()).pt();
    *pfmet_nomu_phi  = (pfmetnomu_h->front()).phi();

    *calomet_pt   = (calomet_h->front()).pt();
    *calomet_phi  = (calomet_h->front()).phi();

    if (genmet_h.isValid()) {
      *genmet_pt   = (genmet_h->front()).pt();
      *genmet_phi  = (genmet_h->front()).phi();
    }

    *pf_ht    = (pfht_h->front()).sumEt();
    *calo_ht  = (caloht_h->front()).sumEt();

    *pf_mht_pt    = (pfht_h->front()).pt();
    *pf_mht_phi   = (pfht_h->front()).phi();

    *pf_mht_nomu_tightid_pt    = (pfhtnomutightid_h->front()).pt();
    *pf_mht_nomu_tightid_phi   = (pfhtnomutightid_h->front()).phi();

    *pf_mht_tightid_pt    = (pfhttightid_h->front()).pt();
    *pf_mht_tightid_phi   = (pfhttightid_h->front()).phi();

    *calo_mht_pt    = (caloht_h->front()).pt();
    *calo_mht_phi   = (caloht_h->front()).phi();

    *n_pix_vtx = pixelvertices_h->size();
    
    if (pusummary_h.isValid()) {
      for(edm::View<PileupSummaryInfo>::const_iterator pusummary_it = pusummary_h->begin(); pusummary_it != pusummary_h->end(); pusummary_it++){
	TrueNumInteractions->push_back(pusummary_it->getTrueNumInteractions());
	num_PU_vertices    ->push_back(pusummary_it->getPU_NumInteractions());
      }
    }

    *pf_offline_ht = 0.;
    if (jet_offline_h.isValid()) {
      for(edm::View<pat::Jet>::const_iterator jet_offline_it = jet_offline_h->begin(); jet_offline_it != jet_offline_h->end(); jet_offline_it++){

	if(jet_offline_it->pt() < 30.0) continue;

	pfjets_offline_pt  ->push_back(jet_offline_it->pt());
	pfjets_offline_eta ->push_back(jet_offline_it->eta());
	pfjets_offline_phi ->push_back(jet_offline_it->phi());

	if(fabs(jet_offline_it->eta()) > 3.0) continue;
	*pf_offline_ht += jet_offline_it->pt();
      }
    }

    if (pfmet_offline_h.isValid()) {
      *pfmet_offline_pt   = (pfmet_offline_h->front()).pt();
      *pfmet_offline_phi  = (pfmet_offline_h->front()).phi();
    }

    //-------------- put statements -----------------

    iEvent.put(pfjets_pt,   "pfjetspt" );
    iEvent.put(pfjets_eta,  "pfjetseta" );
    iEvent.put(pfjets_phi,  "pfjetsphi" );

    if (genjet_h.isValid()) {
      iEvent.put(genjets_pt,   "genjetspt" );
      iEvent.put(genjets_eta,  "genjetseta" );
      iEvent.put(genjets_phi,  "genjetsphi" );
    }

    iEvent.put(calojets_pt,   "calojetspt" );
    iEvent.put(calojets_eta,  "calojetseta" );
    iEvent.put(calojets_phi,  "calojetsphi" );

    iEvent.put(pfmet_pt,   "pfmetpt" );
    iEvent.put(pfmet_phi,  "pfmetphi" );

    iEvent.put(pfmet_nomu_pt,   "pfmetnomupt" );
    iEvent.put(pfmet_nomu_phi,  "pfmetnomuphi" );

    iEvent.put(calomet_pt,   "calometpt" );
    iEvent.put(calomet_phi,  "calometphi" );

    if (genmet_h.isValid()) {
      iEvent.put(genmet_pt,   "genmetpt" );
      iEvent.put(genmet_phi,  "genmetphi" );
    }

    if (genjet_h.isValid()) iEvent.put(gen_ht, "genht");

    iEvent.put(pf_ht,   "pfht" );
    iEvent.put(calo_ht,   "caloht" );

    iEvent.put(pf_mht_pt,   "pfmhtpt" );
    iEvent.put(pf_mht_phi,  "pfmhtphi" );

    iEvent.put(pf_mht_nomu_tightid_pt,   "pfmhtnomutightidpt" );
    iEvent.put(pf_mht_nomu_tightid_phi,  "pfmhtnomutightidphi" );

    iEvent.put(pf_mht_tightid_pt,   "pfmhttightidpt" );
    iEvent.put(pf_mht_tightid_phi,  "pfmhttightidphi" );

    iEvent.put(calo_mht_pt,   "calomhtpt" );
    iEvent.put(calo_mht_phi,  "calomhtphi" );

    iEvent.put(n_pix_vtx,  "npixvtx" );
    
    iEvent.put(scale1fb,  "scale1fb" );

    if (pusummary_h.isValid()) {
      iEvent.put(TrueNumInteractions, "TrueNumInteractions");
      iEvent.put(num_PU_vertices, "numPUvertices");
    }

    if (jet_offline_h.isValid()) {
      iEvent.put(pfjets_offline_pt,   "pfjetsofflinept" );
      iEvent.put(pfjets_offline_eta,  "pfjetsofflineeta" );
      iEvent.put(pfjets_offline_phi,  "pfjetsofflinephi" );

      iEvent.put(pf_offline_ht,   "pfofflineht" );
    }

    if (pfmet_offline_h.isValid()) {
      iEvent.put(pfmet_offline_pt,   "pfmetofflinept" );
      iEvent.put(pfmet_offline_phi,  "pfmetofflinephi" );
    }
}

 
//define this as a plug-in
DEFINE_FWK_MODULE(BabyMaker);
