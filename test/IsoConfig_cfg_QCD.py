import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.MessageLogger.cerr.FwkReport.reportEvery = 50000

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
'/store/mc/RunIISummer16MiniAODv2/QCD_Pt-30to50_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/120000/062CCE36-3DB7-E611-9FE8-44A842CFD64D.root',
'/store/mc/RunIISummer16MiniAODv2/QCD_Pt-30to50_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/120000/08C906C8-C5B6-E611-B825-B499BAAC055E.root',
'/store/mc/RunIISummer16MiniAODv2/QCD_Pt-30to50_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/120000/5A97FFBE-CEB6-E611-AB9B-484D7E8DF0D3.root',
'/store/mc/RunIISummer16MiniAODv2/QCD_Pt-30to50_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/120000/640A347C-C9B6-E611-9C29-6C3BE5B564A8.root',
'/store/mc/RunIISummer16MiniAODv2/QCD_Pt-30to50_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/120000/6416D62E-C3B6-E611-8CEA-484D7E8DF0D3.root',
'/store/mc/RunIISummer16MiniAODv2/QCD_Pt-30to50_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/120000/76EA4211-3DB7-E611-AFA6-0CC47A7C345E.root',
'/store/mc/RunIISummer16MiniAODv2/QCD_Pt-30to50_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/120000/7C7CBDDC-C3B6-E611-A981-484D7E8DF0C6.root',
'/store/mc/RunIISummer16MiniAODv2/QCD_Pt-30to50_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/120000/96C689B9-C5B6-E611-93E9-484D7E8DF0D3.root',
'/store/mc/RunIISummer16MiniAODv2/QCD_Pt-30to50_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/120000/BA2B4893-C7B6-E611-91EA-44A842CFC9CC.root',
'/store/mc/RunIISummer16MiniAODv2/QCD_Pt-30to50_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/120000/D4DD3606-CCB6-E611-8C45-484D7E8DF0ED.root',
'/store/mc/RunIISummer16MiniAODv2/QCD_Pt-30to50_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/120000/DA5421DC-CAB6-E611-BE2E-484D7E8DF0D3.root',
'/store/mc/RunIISummer16MiniAODv2/QCD_Pt-30to50_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/120000/DA617CC7-CCB6-E611-8D92-B499BAAC039C.root',
'/store/mc/RunIISummer16MiniAODv2/QCD_Pt-30to50_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/04E77FA9-1FB6-E611-82A4-24BE05C44BB1.root',
'/store/mc/RunIISummer16MiniAODv2/QCD_Pt-30to50_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/0E7EC73D-0CB6-E611-99CC-24BE05C616E1.root',
'/store/mc/RunIISummer16MiniAODv2/QCD_Pt-30to50_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/0EEC81C5-1AB6-E611-A7AD-5065F381F272.root',
'/store/mc/RunIISummer16MiniAODv2/QCD_Pt-30to50_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/100B1AD8-1CB6-E611-A692-A0000420FE80.root',
'/store/mc/RunIISummer16MiniAODv2/QCD_Pt-30to50_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/162ED122-15B6-E611-B9C3-A0000420FE80.root',
'/store/mc/RunIISummer16MiniAODv2/QCD_Pt-30to50_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/18F932E1-22B6-E611-9B9D-A0000420FE80.root',
'/store/mc/RunIISummer16MiniAODv2/QCD_Pt-30to50_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/1C0A4BF6-21B6-E611-A25E-A0000420FE80.root',
'/store/mc/RunIISummer16MiniAODv2/QCD_Pt-30to50_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/1C53DCC6-21B6-E611-AA9D-A0000420FE80.root',
'/store/mc/RunIISummer16MiniAODv2/QCD_Pt-30to50_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/1E462460-14B6-E611-BDC2-A0000420FE80.root',
'/store/mc/RunIISummer16MiniAODv2/QCD_Pt-30to50_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/227ABB34-25B6-E611-A6B0-A0000420FE80.root',
'/store/mc/RunIISummer16MiniAODv2/QCD_Pt-30to50_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/2A07B59D-19B6-E611-B8EA-B8CA3A70A5E8.root',
'/store/mc/RunIISummer16MiniAODv2/QCD_Pt-30to50_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/2AD55091-20B6-E611-9E31-5065F3820351.root',
'/store/mc/RunIISummer16MiniAODv2/QCD_Pt-30to50_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/2AF7BA1B-19B6-E611-AD3F-24BE05CECBD1.root',
'/store/mc/RunIISummer16MiniAODv2/QCD_Pt-30to50_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/30266E6C-C1B6-E611-8323-0CC47A4C8EB6.root',
'/store/mc/RunIISummer16MiniAODv2/QCD_Pt-30to50_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/309E5F65-0FB6-E611-9297-A0000420FE80.root',
'/store/mc/RunIISummer16MiniAODv2/QCD_Pt-30to50_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/368F4C84-20B6-E611-88F6-A0000420FE80.root',
'/store/mc/RunIISummer16MiniAODv2/QCD_Pt-30to50_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/3C662563-09B6-E611-AC87-24BE05C6C7E1.root',
'/store/mc/RunIISummer16MiniAODv2/QCD_Pt-30to50_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/40101E86-0DB6-E611-AD0B-0CC47A4D76D0.root',
'/store/mc/RunIISummer16MiniAODv2/QCD_Pt-30to50_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/42094C23-15B6-E611-9E80-A0000420FE80.root',
'/store/mc/RunIISummer16MiniAODv2/QCD_Pt-30to50_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/4625B831-25B6-E611-AF68-24BE05C4D8C1.root',
'/store/mc/RunIISummer16MiniAODv2/QCD_Pt-30to50_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/464C7048-11B6-E611-BBEE-A0000420FE80.root',
'/store/mc/RunIISummer16MiniAODv2/QCD_Pt-30to50_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/5037E11B-19B6-E611-8570-24BE05CE2EE1.root',
'/store/mc/RunIISummer16MiniAODv2/QCD_Pt-30to50_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/549BB963-18B6-E611-94EA-B8CA3A70BAC8.root',
'/store/mc/RunIISummer16MiniAODv2/QCD_Pt-30to50_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/5C6803DA-1EB6-E611-A8C0-A0000420FE80.root',
'/store/mc/RunIISummer16MiniAODv2/QCD_Pt-30to50_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/5EC896F5-12B6-E611-B481-A0000420FE80.root',
'/store/mc/RunIISummer16MiniAODv2/QCD_Pt-30to50_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/6220120E-23B6-E611-A548-A0000420FE80.root',
'/store/mc/RunIISummer16MiniAODv2/QCD_Pt-30to50_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/64A6514D-C1B6-E611-BB73-5065F38122D1.root',
'/store/mc/RunIISummer16MiniAODv2/QCD_Pt-30to50_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/664020F8-15B6-E611-BDF5-A0000420FE80.root',
'/store/mc/RunIISummer16MiniAODv2/QCD_Pt-30to50_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/68CD698B-27B6-E611-BAC6-A0369F3016EC.root',
'/store/mc/RunIISummer16MiniAODv2/QCD_Pt-30to50_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/6AE1B31B-19B6-E611-BB15-A0000420FE80.root',
'/store/mc/RunIISummer16MiniAODv2/QCD_Pt-30to50_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/6C4B834A-0EB6-E611-88B8-0025905A60DA.root',
'/store/mc/RunIISummer16MiniAODv2/QCD_Pt-30to50_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/6E9DFDFE-1CB6-E611-A498-A0000420FE80.root',
'/store/mc/RunIISummer16MiniAODv2/QCD_Pt-30to50_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/74043C6C-18B6-E611-A14F-A0000420FE80.root',
'/store/mc/RunIISummer16MiniAODv2/QCD_Pt-30to50_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/78468EF8-16B6-E611-A45F-A0000420FE80.root',
'/store/mc/RunIISummer16MiniAODv2/QCD_Pt-30to50_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/7A6B2841-1CB6-E611-8CFE-A0000420FE80.root',
'/store/mc/RunIISummer16MiniAODv2/QCD_Pt-30to50_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/8487829F-1AB6-E611-A621-A0000420FE80.root',
'/store/mc/RunIISummer16MiniAODv2/QCD_Pt-30to50_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/84F518A7-13B6-E611-B9C4-24BE05CE4D91.root',
'/store/mc/RunIISummer16MiniAODv2/QCD_Pt-30to50_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/8AB96A1C-19B6-E611-823E-A0000420FE80.root',
'/store/mc/RunIISummer16MiniAODv2/QCD_Pt-30to50_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/8CE84119-23B6-E611-8C4D-A0369F3016EC.root',
'/store/mc/RunIISummer16MiniAODv2/QCD_Pt-30to50_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/8E5F9325-15B6-E611-AF66-B8CA3A70BAC8.root',
'/store/mc/RunIISummer16MiniAODv2/QCD_Pt-30to50_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/92551110-12B6-E611-9B15-24BE05C46B01.root',
'/store/mc/RunIISummer16MiniAODv2/QCD_Pt-30to50_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/96AA35D8-0DB6-E611-B5C5-A0000420FE80.root',
'/store/mc/RunIISummer16MiniAODv2/QCD_Pt-30to50_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/A65AB32D-11B6-E611-B820-A0000420FE80.root',
'/store/mc/RunIISummer16MiniAODv2/QCD_Pt-30to50_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/A8AC686E-C1B6-E611-9A14-0025905B85A0.root',
'/store/mc/RunIISummer16MiniAODv2/QCD_Pt-30to50_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/B602DCF6-16B6-E611-A21D-A0000420FE80.root',
'/store/mc/RunIISummer16MiniAODv2/QCD_Pt-30to50_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/B6582EED-1DB6-E611-B51B-A0000420FE80.root',
'/store/mc/RunIISummer16MiniAODv2/QCD_Pt-30to50_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/B8EA81CD-1EB6-E611-A269-A0000420FE80.root',
'/store/mc/RunIISummer16MiniAODv2/QCD_Pt-30to50_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/C691A17E-25B6-E611-8DDD-A0000420FE80.root',
'/store/mc/RunIISummer16MiniAODv2/QCD_Pt-30to50_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/C847FDEC-1BB6-E611-A6DC-24BE05CE2EE1.root',
'/store/mc/RunIISummer16MiniAODv2/QCD_Pt-30to50_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/CC8A2539-24B6-E611-8E99-A0000420FE80.root',
'/store/mc/RunIISummer16MiniAODv2/QCD_Pt-30to50_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/CE336F38-21B6-E611-9B91-B8CA3A70A5E8.root',
'/store/mc/RunIISummer16MiniAODv2/QCD_Pt-30to50_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/D2975A39-18B6-E611-B237-A0000420FE80.root',
'/store/mc/RunIISummer16MiniAODv2/QCD_Pt-30to50_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/DA8A4E96-26B6-E611-9BA0-A0000420FE80.root',
'/store/mc/RunIISummer16MiniAODv2/QCD_Pt-30to50_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/DC9C5437-25B6-E611-A623-A0000420FE80.root',
'/store/mc/RunIISummer16MiniAODv2/QCD_Pt-30to50_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/E44FF8BC-1DB6-E611-BB44-B8CA3A70A5E8.root',
'/store/mc/RunIISummer16MiniAODv2/QCD_Pt-30to50_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/F653362B-24B6-E611-BD09-A0000420FE80.root',
'/store/mc/RunIISummer16MiniAODv2/QCD_Pt-30to50_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/FCA026D3-1AB6-E611-B9BA-A0000420FE80.root',
'/store/mc/RunIISummer16MiniAODv2/QCD_Pt-30to50_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/80000/00087D12-10B6-E611-A34C-B499BAABD28C.root'
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
     fileName = cms.string('IsoNtuple_QCD.root')
 )

process.p = cms.Path(process.demo)
