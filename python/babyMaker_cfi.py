import FWCore.ParameterSet.Config as cms

babyMaker = cms.EDProducer("BabyMaker",
                            pfJetsInputTag_   = cms.InputTag("hltAK4PFJetsCorrected"),
                            pfMetInputTag_    = cms.InputTag("hltPFMETProducer"),
                            pfMetNoMuInputTag_    = cms.InputTag("hltPFMETNoMuProducer"),
                            pfHTInputTag_     = cms.InputTag("hltPFHT"),
                            pfHTNoMuTightIDInputTag_  = cms.InputTag("hltPFMHTNoMuTightID"),
                            pfHTTightIDInputTag_  = cms.InputTag("hltPFMHTTightID"),
                            caloMetInputTag_  = cms.InputTag("hltMet"),
                            caloHTInputTag_   = cms.InputTag("hltHtMht"),
                            genJetsInputTag_  = cms.InputTag("ak4GenJetsNoNu"),
                            caloJetsInputTag_ = cms.InputTag("hltAK4CaloJetsCorrected"),
                            genMETInputTag_   = cms.InputTag("genMetTrue"),
                            PileupSummaryInputTag_   = cms.InputTag("addPileupInfo"),
)


