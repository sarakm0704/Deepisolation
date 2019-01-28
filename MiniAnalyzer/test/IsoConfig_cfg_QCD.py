import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.MessageLogger.cerr.FwkReport.reportEvery = 50000

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
'/store/mc/RunIISummer16MiniAODv2/QCD_Pt-30to50_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/120000/062CCE36-3DB7-E611-9FE8-44A842CFD64D.root',
             )    
)

process.demo = cms.EDAnalyzer('IsoAnalyzer',
    muonLabel =  cms.InputTag('slimmedMuons'),
    pfSrc = cms.InputTag("packedPFCandidates"), 
    VertexTag = cms.InputTag("offlineSlimmedPrimaryVertices"),
    genLabel = cms.InputTag("genParticles"),
    dRch = cms.untracked.double(0.4),
    dRnh = cms.untracked.double(0.4),
    dRph = cms.untracked.double(0.4),
    vetoCHThreshold =   cms.untracked.double(0.0),
    vetoNHThreshold =   cms.untracked.double(0.5),
    vetoPHThreshold =   cms.untracked.double(0.5),
    vetoCHdR =   cms.untracked.double(0.0),
    vetoNHdR =   cms.untracked.double(0.0),
    vetoPHdR =   cms.untracked.double(0.0),
    vetoCHre = cms.untracked.vdouble(0.0,0.0,0.0,0.0),
    vetoNHre = cms.untracked.vdouble(0.0,0.0,0.0,0.0),
    vetoPHre = cms.untracked.vdouble(-0.1,0.1,-0.2,0.2),
)

process.TFileService = cms.Service("TFileService",
     fileName = cms.string('IsoNtuple_QCD.root')
 )

process.p = cms.Path(process.demo)
