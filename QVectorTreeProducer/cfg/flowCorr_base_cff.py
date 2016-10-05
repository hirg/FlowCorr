import FWCore.ParameterSet.Config as cms

process = cms.Process("FlowCorr")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1

process.maxEvents = cms.untracked.PSet( 
   input = cms.untracked.int32(100)
)

process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(False),
    SkipEvent = cms.untracked.vstring('ProductNotFound')
)

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        'root://cms-xrd-global.cern.ch///store/user/qwang/HIMinimumBias3/HIMinBias_v2/160128_201904/0000/HIMinBias_1.root'
    )
)


process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '75X_dataRun2_v12', '')

process.load("HLTrigger.HLTfilters.hltHighLevel_cfi")
process.hltMB = process.hltHighLevel.clone()
process.hltMB.HLTPaths = [
               "HLT_HIL1MinimumBiasHF1AND*",
               "HLT_HIL1MinimumBiasHF2AND*"
              ]
process.hltMB.andOr = cms.bool(True)
process.hltMB.throw = cms.bool(False)

### centrality
process.load("RecoHI.HiCentralityAlgos.CentralityBin_cfi")
process.newCentralityBin = process.centralityBin.clone()

process.TFileService = cms.Service("TFileService",
                                  fileName=cms.string("flowCorr_data.root")
                                  )

#Analyzer
process.load("FlowCorr.QVectorTreeProducer.flowCorr_cfi")
process.testcumu = process.flowCorr.clone()

process.p = cms.Path(process.newCentralityBin * process.hltMB * process.testcumu)

