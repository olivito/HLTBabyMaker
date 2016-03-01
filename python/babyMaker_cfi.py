import FWCore.ParameterSet.Config as cms

babyMaker = cms.EDProducer("BabyMaker",
                            pfJetsInputTag_   = cms.InputTag("hltAK4PFJetsCorrected"),
                            pfMetInputTag_    = cms.InputTag("hltPFMETProducer"),
                            pfMetNoMuInputTag_    = cms.InputTag("hltPFMETNoMuProducer"),
                            pfMetNoHFInputTag_    = cms.InputTag("hltPFMETNoHFProducer"),
                            pfHTInputTag_     = cms.InputTag("hltPFHT"),
                            pfHTNoMuTightIDInputTag_  = cms.InputTag("hltPFMHTNoMuTightID"),
                            pfHTTightIDInputTag_  = cms.InputTag("hltPFMHTTightID"),
                            pfHTTightIDCentralInputTag_  = cms.InputTag("hltPFMHTTightIDCentral"),
                            pfHTNoHFTightIDCentralInputTag_  = cms.InputTag("hltPFMHTNoHFTightIDCentral"),
                            caloMetInputTag_  = cms.InputTag("hltMet"),
                            caloMetHBHECleanedInputTag_  = cms.InputTag("hltMetClean"),
                            caloHTInputTag_   = cms.InputTag("hltHtMht"),
                            genJetsInputTag_  = cms.InputTag("ak4GenJetsNoNu"),
                            caloJetsInputTag_ = cms.InputTag("hltAK4CaloJetsCorrected"),
                            genMETInputTag_   = cms.InputTag("genMetTrue"),
                            pixelVerticesInputTag_   = cms.InputTag("hltPixelVertices"),
                            pfJetsOfflineInputTag_   = cms.InputTag("slimmedJets","","PAT"),
                            pfMetOfflineInputTag_    = cms.InputTag("slimmedMETs","","PAT"),
)


