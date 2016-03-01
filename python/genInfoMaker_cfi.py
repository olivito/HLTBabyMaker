import FWCore.ParameterSet.Config as cms

genInfoMaker = cms.EDProducer("GenInfoMaker",
                            PileupSummaryInputTag_   = cms.InputTag("addPileupInfo"),
                            genEvtInfoInputTag_   = cms.InputTag("generator"),
)


