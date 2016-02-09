import FWCore.ParameterSet.Config as cms

process = cms.Process('HLTTRIGBOOLS')

# load event level configurations
process.load('Configuration/EventContent/EventContent_cff')
process.load("Configuration.StandardSequences.Services_cff")
process.load('Configuration.Geometry.GeometryRecoDB_cff')
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff")
process.load("Configuration.StandardSequences.GeometryRecoDB_cff")

# services
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 100000

# global tag
from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag as customiseGlobalTag
process.GlobalTag = customiseGlobalTag(process.GlobalTag, globaltag = 'auto:run2_mc_GRun')
process.GlobalTag.connect   = 'frontier://FrontierProd/CMS_CONDITIONS'
process.GlobalTag.pfnPrefix = cms.untracked.string('frontier://FrontierProd/')
for pset in process.GlobalTag.toGet.value():
    pset.connect = pset.connect.value().replace('frontier://FrontierProd/', 'frontier://FrontierProd/')
#process.GlobalTag.globaltag = "auto:hltonline"

process.load('HLTStudy.HLTBabyMaker.trigBoolMaker_cfi')


process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)


# Input source
process.source = cms.Source("PoolSource",
   fileNames = cms.untracked.vstring('file:ntuple_hlt_mc_v5.root'),
)


process.out = cms.OutputModule("PoolOutputModule",
  fileName     = cms.untracked.string('ntuple_trigbools.root'),
  dropMetaData = cms.untracked.string("NONE")
)

process.out.outputCommands = cms.untracked.vstring( 'drop *' )
process.out.outputCommands.extend(cms.untracked.vstring('keep *_*babyMaker*_*_*'))
process.out.outputCommands.extend(cms.untracked.vstring('keep *_*trigBoolMaker*_*_*'))
process.out.outputCommands.extend(cms.untracked.vstring('keep *_TriggerResults_*_HLT'))
process.out.outputCommands.extend(cms.untracked.vstring('keep *_TriggerResults_*_reHLT'))


# Path and EndPath definitions
process.maker = cms.Path(process.trigBoolMaker)
process.outpath = cms.EndPath(process.out)
