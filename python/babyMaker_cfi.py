import FWCore.ParameterSet.Config as cms

babyMaker = cms.EDProducer("BabyMaker",
                            pfJetsInputTag_   = cms.InputTag("hltAK4PFJetL1FastL2L3Corrected"),
                            pfMetInputTag_    = cms.InputTag("hltPFMETProducer"),
                            pfHTInputTag_     = cms.InputTag("hltPFHT"),
                            caloHTInputTag_   = cms.InputTag("hltHtMht"),
                            hemInputTag_      = cms.InputTag("hltRHemisphere"),
                            genJetsInputTag_  = cms.InputTag("ak4GenJetsNoNu"),
                            caloJetsInputTag_ = cms.InputTag("hltCaloJetL1FastJetCorrected"),
                            genMETInputTag_  = cms.InputTag("genMetTrue"),
)


