import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
        #'file:/afs/cern.ch/cms/Tutorials/TWIKI_DATA/TTJets_8TeV_53X.root'
        #'/store/data/Run2016H/DoubleMuon/MINIAOD/PromptReco-v3/000/284/068/00000/08D18E09-6B9F-E611-ABE6-02163E0144F3.root'
        #'/store/mc/RunIIFall17MiniAOD/ZToMuMu_NNPDF31_13TeV-powheg_M_50_120/MINIAODSIM/94X_mc2017_realistic_v10-v1/20000/4CFC5DA8-F904-E811-A134-FA163E3BAEF5.root'
    #  QCD   'root://cms-xrd-global.cern.ch//store/mc/RunIISummer16MiniAODv2/QCD_Pt-30to50_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/120000/062CCE36-3DB7-E611-9FE8-44A842CFD64D.root'
        'root://cms-xrd-global.cern.ch//store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext2-v1/110000/005ED0EB-79F1-E611-B6DA-02163E011C2B.root'
    )
)

process.demo = cms.EDAnalyzer('IsoAnalyzer',
    muonLabel =  cms.InputTag('slimmedMuons'),
    #pfCandidateLabel = cms.InputTag('particleFlow'),
    pfSrc = cms.InputTag("packedPFCandidates"), 
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
     fileName = cms.string('IsoNtuple_Z.root')
 )

process.p = cms.Path(process.demo)
