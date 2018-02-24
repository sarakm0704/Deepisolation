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
#IsoNtuple_QCD.root
        #'root://cms-xrd-global.cern.ch//store/mc/RunIISummer16MiniAODv2/QCD_Pt-30to50_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/120000/062CCE36-3DB7-E611-9FE8-44A842CFD64D.root',
#'root://cms-xrd-global.cern.ch//store/mc/RunIISummer16MiniAODv2/QCD_Pt-30to50_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/120000/08C906C8-C5B6-E611-B825-B499BAAC055E.root',
#'root://cms-xrd-global.cern.ch//store/mc/RunIISummer16MiniAODv2/QCD_Pt-30to50_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/120000/5A97FFBE-CEB6-E611-AB9B-484D7E8DF0D3.root',
#'root://cms-xrd-global.cern.ch//store/mc/RunIISummer16MiniAODv2/QCD_Pt-30to50_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/120000/640A347C-C9B6-E611-9C29-6C3BE5B564A8.root',
#'root://cms-xrd-global.cern.ch//store/mc/RunIISummer16MiniAODv2/QCD_Pt-30to50_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/120000/6416D62E-C3B6-E611-8CEA-484D7E8DF0D3.root',

#IsoNtuple_QCD_1.root
#'root://cms-xrd-global.cern.ch//store/mc/RunIISummer16MiniAODv2/QCD_Pt-30to50_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/120000/76EA4211-3DB7-E611-AFA6-0CC47A7C345E.root',
#'root://cms-xrd-global.cern.ch//store/mc/RunIISummer16MiniAODv2/QCD_Pt-30to50_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/120000/7C7CBDDC-C3B6-E611-A981-484D7E8DF0C6.root',

#IsoNtuple_QCD_2.root_pass
#'root://cms-xrd-global.cern.ch//store/mc/RunIISummer16MiniAODv2/QCD_Pt-30to50_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/120000/96C689B9-C5B6-E611-93E9-484D7E8DF0D3.root',

#IsoNtuple_QCD_3.root_pass
#'root://cms-xrd-global.cern.ch//store/mc/RunIISummer16MiniAODv2/QCD_Pt-30to50_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/120000/BA2B4893-C7B6-E611-91EA-44A842CFC9CC.root',

#IsoNtuple_QCD_4.root_pass
#'root://cms-xrd-global.cern.ch//store/mc/RunIISummer16MiniAODv2/QCD_Pt-30to50_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/120000/D4DD3606-CCB6-E611-8C45-484D7E8DF0ED.root'
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
     fileName = cms.string('IsoNtuple_QCD_2.root')
 )

process.p = cms.Path(process.demo)
