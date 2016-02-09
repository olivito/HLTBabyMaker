import FWCore.ParameterSet.Config as cms

trigBoolMaker = cms.EDProducer("TrigBoolMaker",
                            triggerResultsInputTag_   = cms.InputTag("TriggerResults","","reHLT"),
                            processName_   = cms.string("reHLT"),
)


