import FWCore.ParameterSet.Config as cms

babyMaker = cms.EDProducer("BabyMaker",
                            #aliasPrefix = cms.untracked.string("evt"),
                            pfJetsInputTag_ = cms.InputTag("hltAK4PFJetL1FastL2L3Corrected"),
                            pfMetInputTag_  = cms.InputTag("hltPFMETProducer"),
                            pfHTInputTag_   = cms.InputTag("hltPFHT"),
                            hemInputTag_    = cms.InputTag("hltRHemisphere")
)


