import FWCore.ParameterSet.Config as cms

from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing ('analysis')
options.register('Labels', 0, VarParsing.multiplicity.singleton, VarParsing.varType.int, "1: DY, 0:QCD")
options.parseArguments()

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.MessageLogger.cerr.FwkReport.reportEvery = 50000

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
'/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext2-v1/110000/005ED0EB-79F1-E611-B6DA-02163E011C2B.root',
    )
)

process.demo = cms.EDAnalyzer('IsoAnalyzer',
    muonLabel =  cms.InputTag('slimmedMuons'),
    pfSrc = cms.InputTag("packedPFCandidates"), 
    VertexTag = cms.InputTag("offlineSlimmedPrimaryVertices"),
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
    label = cms.untracked.int32(options.Labels),
)

process.TFileService = cms.Service("TFileService",
     fileName = cms.string('22_DY.root')
 )

process.p = cms.Path(process.demo)
