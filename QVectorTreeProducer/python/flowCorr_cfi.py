import FWCore.ParameterSet.Config as cms

flowCorr = cms.EDAnalyzer('QVectorTreeProducer',
     #tracks
     tracks       = cms.InputTag('hiGeneralTracks'),
     trackQuality = cms.untracked.string('highpurity'),
     noffmin      = cms.untracked.int32(0),
     noffmax      = cms.untracked.int32(10000),
      #--- offline tracks
     minEtaOff    = cms.untracked.double(-2.4),
     maxEtaOff    = cms.untracked.double(2.4),
     minPtOff     = cms.untracked.double(0.4),
     maxPtOff     = cms.untracked.double(99999.),
     chargeOff    = cms.untracked.vint32(-1,1),
     isPixTrkOff  = cms.untracked.bool(False),
     dzdzErrorOff = cms.untracked.double(3.0),
     d0d0ErrorOff = cms.untracked.double(3.0),
     ptErrorPtOff = cms.untracked.double(0.1),
     chi2nOff     = cms.untracked.double(0.15),
     nHitsOff     = cms.untracked.int32(11),
     trkAlgoOff   = cms.untracked.vint32(4,5,6,7),
      #--- reference tracks
     minEtaRef    = cms.untracked.double(-2.4),
     maxEtaRef    = cms.untracked.double(2.4),
     minPtRef     = cms.untracked.double(0.3),
     maxPtRef     = cms.untracked.double(3.0),
     chargeRef    = cms.untracked.vint32(-1,1),
     isPixTrkRef  = cms.untracked.bool(False),
     dzdzErrorRef = cms.untracked.double(3.0),
     d0d0ErrorRef = cms.untracked.double(3.0),
     ptErrorPtRef = cms.untracked.double(0.1),
     chi2nRef     = cms.untracked.double(0.15),
     nHitsRef     = cms.untracked.int32(11),
     trkAlgoRef   = cms.untracked.vint32(4,5,6,7),

     #vertex
     vertex  = cms.InputTag('hiSelectedVertex'),
     minVz   = cms.untracked.double(-15.0),
     maxVz   = cms.untracked.double(15.0),
     maxRho  = cms.untracked.double(0.2),
     nVtxMax = cms.untracked.int32(9999),

     #calotower
     caloTower = cms.InputTag('towerMaker'),

     #centrality
     centralitySrc    = cms.InputTag('hiCentrality'),
     centralityBinSrc = cms.untracked.InputTag('centralityBin','HFtowers'),
     minCent          = cms.untracked.int32(-1),
     maxCent          = cms.untracked.int32(-1),

     #event plane
     evtPlane     = cms.InputTag('hiEvtPlane'),
     evtPlaneFlat = cms.InputTag('hiEvtPlaneFlat',''),  
     epLvl        = cms.untracked.int32(0),

     #global
)
