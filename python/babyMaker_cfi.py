import FWCore.ParameterSet.Config as cms

babyMaker = cms.EDProducer("BabyMaker",
                            pfJetsInputTag_   = cms.InputTag("hltAK4PFJetsCorrected"),
                            pfMetInputTag_    = cms.InputTag("hltPFMETProducer"),
                            pfHTInputTag_     = cms.InputTag("hltPFHT"),
                            caloMetInputTag_  = cms.InputTag("hltMet"),
                            caloHTInputTag_   = cms.InputTag("hltHtMht"),
                            genJetsInputTag_  = cms.InputTag("ak4GenJetsNoNu"),
                            caloJetsInputTag_ = cms.InputTag("hltAK4CaloJetsCorrected"),
                            genMETInputTag_   = cms.InputTag("genMetTrue"),
)


