// -*- C++ -*-
//
// Package:    Iso/IsoAnalyzer
// Class:      IsoAnalyzer
// 
/**\class IsoAnalyzer IsoAnalyzer.cc Iso/IsoAnalyzer/plugins/IsoAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Ji Eun Choi
//         Created:  Mon, 05 Feb 2018 09:35:36 GMT
//
//


// system include files
#include <memory>
#include <iostream>
// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/PatCandidates/interface/Muon.h"


#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "CommonTools/Utils/interface/StringObjectFunction.h"
#include "DataFormats/PatCandidates/interface/LookupTableRecord.h"

#include "DataFormats/RecoCandidate/interface/IsoDepositDirection.h"
#include "DataFormats/RecoCandidate/interface/IsoDeposit.h"
#include "DataFormats/RecoCandidate/interface/IsoDepositVetos.h"
#include "DataFormats/PatCandidates/interface/Isolation.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"

#include "TH1.h"
#include "TTree.h"

//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<> and also remove the line from
// constructor "usesResource("TFileService");"
// This will improve performance in multithreaded jobs.

using namespace std;

class IsoAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit IsoAnalyzer(const edm::ParameterSet&);
      ~IsoAnalyzer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      // ----------member data ---------------------------
      //edm::InputTag minTracks_; //used to select what tracks to read from configuration file
      edm::EDGetTokenT<edm::View<pat::Muon> > muonLabel_;
      //edm::EDGetTokenT<reco::PFCandidateCollection> pfCandidateLabel_;
      edm::EDGetTokenT<pat::PackedCandidateCollection> pfSrc_;

      TTree *tree;

      double dRch_, dRnh_, dRph_;
      double vetoCHThreshold_, vetoNHThreshold_, vetoPHThreshold_;
      double vetoCHdR_, vetoNHdR_, vetoPHdR_;
      vector<double> vetoCHre_;
      vector<double> vetoNHre_;
      vector<double> vetoPHre_;

      unsigned int EVENT, RUN, LUMI;
      double pt, eta, phi;

      double chIso0_005, chIso005_01, chIso01_015, chIso015_02, chIso02_025, chIso025_03;
      double chIso, schIso;
      double nhIso0_005, nhIso005_01, nhIso01_015, nhIso015_02, nhIso02_025, nhIso025_03;
      double nhIso;
      double phIso0_005, phIso005_01, phIso01_015, phIso015_02, phIso02_025, phIso025_03;
      double phIso;
      double Pileup, Pileup0_005, Pileup005_01, Pileup01_015, Pileup015_02, Pileup02_025, Pileup025_03;
      double absIso, relIso;

      double chIso_org, puChIso_org, nhIso_org, phIso_org, absIso_org, relIso_org;

      int chN, nhN, phN;
      double ch_dR, ch_dphi, ch_deta;
      double nh_dR, nh_dphi, nh_deta;
      double ph_dR, ph_dphi, ph_deta;

      int chn, nhn, phn;
      double chR, nhR, phR;
      double chPhi, nhPhi, phPhi;
      double chEta, nhEta, phEta;
 
      double chiso0_005, chiso005_01, chiso01_015, chiso015_02, chiso02_025, chiso025_03;
      double chiso, schiso;
      double nhiso0_005, nhiso005_01, nhiso01_015, nhiso015_02, nhiso02_025, nhiso025_03;
      double nhiso;
      double phiso0_005, phiso005_01, phiso01_015, phiso015_02, phiso02_025, phiso025_03;
      double phiso;
      double pileup0_005, pileup005_01, pileup01_015, pileup015_02, pileup02_025, pileup025_03;
      double pileup;
      double absiso, reliso;
   
      double chiso_org, puChiso_org, nhiso_org, phiso_org, absiso_org, reliso_org, neutral_org;
      double neutral;

      TH1D *IsoNtuple; 
};


IsoAnalyzer::IsoAnalyzer(const edm::ParameterSet& iConfig)
 :
  muonLabel_(consumes<edm::View<pat::Muon> >(iConfig.getParameter<edm::InputTag>( "muonLabel" ))),
  //pfCandidateLabel_(consumes<pat::PFCandidateCollection>(iConfig.getParameter<edm::InputTag>("pfCandidateLabel"))),
  pfSrc_(consumes<pat::PackedCandidateCollection>(iConfig.getParameter<edm::InputTag>("pfSrc"))),
  dRch_(iConfig.getUntrackedParameter<double>("dRch",0.4)),
  dRnh_(iConfig.getUntrackedParameter<double>("dRnh",0.4)),
  dRph_(iConfig.getUntrackedParameter<double>("dRph",0.4)),
  vetoCHThreshold_(iConfig.getUntrackedParameter<double>("vetoCHThreshold",0.0)),
  vetoNHThreshold_(iConfig.getUntrackedParameter<double>("vetoNHThreshold",0.0)),
  vetoPHThreshold_(iConfig.getUntrackedParameter<double>("vetoPHThreshold",0.0)),
  vetoCHdR_(iConfig.getUntrackedParameter<double>("vetoCHdR",0.0)),
  vetoNHdR_(iConfig.getUntrackedParameter<double>("vetoNHdR",0.0)),
  vetoPHdR_(iConfig.getUntrackedParameter<double>("vetoPHdR",0.0)),
  vetoCHre_(iConfig.getUntrackedParameter< std::vector<double> >("vetoCHre")),
  vetoNHre_(iConfig.getUntrackedParameter< std::vector<double> >("vetoNHre")),
  vetoPHre_(iConfig.getUntrackedParameter< std::vector<double> >("vetoPHre"))

{
  usesResource("TFileService");
  edm::Service<TFileService> fs;
  tree = fs->make<TTree>("tree", "Tree for isolation study");

  IsoNtuple = fs->make<TH1D>("tracks" , "Tracks" , 100 , 0 , 5000 );
}

IsoAnalyzer::~IsoAnalyzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


// member functions
// ------------ method called for each event  ------------
void IsoAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

  using namespace reco;
  using namespace isodeposit;
  using namespace edm;

  pt = eta = phi =0;

  chIso0_005 = 0;
  chIso005_01 = 0;
  chIso01_015 = 0;
  chIso015_02 = 0;
  chIso02_025 = 0;
  chIso025_03 = 0;
  chIso = 0;
  schIso = 0;

  nhIso0_005 = 0;
  nhIso005_01 = 0;
  nhIso01_015 = 0;
  nhIso015_02 = 0;
  nhIso02_025 = 0;
  nhIso025_03 = 0;
  nhIso = 0;

  phIso0_005 = 0;
  phIso005_01 = 0;
  phIso01_015 = 0;
  phIso015_02 = 0;
  phIso02_025 = 0;
  phIso025_03 = 0;
  phIso = 0;

  absIso = 0;
  relIso = 0;

  Pileup = 0;
  Pileup0_005 = 0;
  Pileup005_01 = 0;
  Pileup01_015 = 0;
  Pileup015_02 = 0;
  Pileup02_025 = 0;
  Pileup025_03 = 0;

  chIso_org = 0;
  puChIso_org = 0;
  nhIso_org = 0;
  phIso_org = 0;
  absIso_org = 0;
  relIso_org = 0;

  chN = 0;
  nhN = 0;
  phN = 0;

  ch_dR = 0;
  ch_deta = 0;
  ch_dphi = 0;

  nh_dR = 0;
  nh_deta = 0;
  nh_dphi = 0;

  ph_dR = 0;
  ph_deta = 0;
  ph_dphi = 0;

  edm::Handle<edm::View<pat::Muon> > muons;
  iEvent.getByToken(muonLabel_, muons);

  //edm::Handle<reco::PFCandidateCollection> pfCandidates_; 
  //iEvent.getByToken(pfCandidateLabel_, pfCandidates_);

  edm::Handle<pat::PackedCandidateCollection> pfcands;
  iEvent.getByToken(pfSrc_, pfcands);

  bool GoodMuon = false;

  for (unsigned i = 0; i != muons->size(); i++) {
    pat::Muon muon = muons->at(i);

    if( muon.pt() < 20) continue;
    if( muon.eta() < -2.1) continue;
    if( muon.eta() > 2.1) continue;

    if( muon.pt() > 20 && std::abs(muon.eta()) < 2.1) GoodMuon = true;

    EVENT  = iEvent.id().event();
    RUN    = iEvent.id().run();
    LUMI   = iEvent.id().luminosityBlock();

    pt = muon.pt();
    eta = muon.eta();
    phi = muon.phi();
 
    IsoDeposit::Direction Dir = Direction(muon.eta(),muon.phi());

    IsoDeposit::AbsVetos vetos_ch;
    vetos_ch.push_back(new ConeVeto( Dir, vetoCHdR_ ));
    vetos_ch.push_back(new ThresholdVeto( vetoCHThreshold_ ));
    vetos_ch.push_back(new RectangularEtaPhiVeto( Dir, vetoCHre_[0], vetoCHre_[1], vetoCHre_[2], vetoCHre_[3]));

    IsoDeposit::AbsVetos vetos_nh;
    vetos_nh.push_back(new ConeVeto( Dir, vetoNHdR_ ));
    vetos_nh.push_back(new ThresholdVeto( vetoNHThreshold_ ));
    vetos_nh.push_back(new RectangularEtaPhiVeto( Dir, vetoNHre_[0], vetoNHre_[1], vetoNHre_[2], vetoNHre_[3]));

    IsoDeposit::AbsVetos vetos_ph;
    vetos_ph.push_back(new ConeVeto( Dir, vetoPHdR_ ));
    vetos_ph.push_back(new ThresholdVeto( vetoPHThreshold_ ));
    vetos_ph.push_back(new RectangularEtaPhiVeto( Dir, vetoPHre_[0], vetoPHre_[1], vetoPHre_[2], vetoPHre_[3]));


    chn = 0;
    nhn = 0;
    phn = 0;

    chR = 0;
    nhR = 0;
    phR = 0;
  
    chPhi = 0;
    nhPhi = 0;
    phPhi = 0;

    chEta = 0;
    nhEta = 0;
    phEta = 0;
 
    chiso0_005 = 0;
    chiso005_01 = 0;
    chiso01_015 = 0;
    chiso015_02 = 0;
    chiso02_025 = 0;
    chiso025_03 = 0;
    chiso = 0;
    schiso = 0;

    nhiso0_005 = 0;
    nhiso005_01 = 0;
    nhiso01_015 = 0;
    nhiso015_02 = 0;
    nhiso02_025 = 0;
    nhiso025_03 = 0;
    nhiso = 0;

    phiso0_005 = 0;
    phiso005_01 = 0;
    phiso01_015 = 0;
    phiso015_02 = 0;
    phiso02_025 = 0;
    phiso025_03 = 0;
    phiso = 0;

    absiso = 0;
    reliso = 0;
   
    chiso_org = 0;
    puChiso_org = 0;
    nhiso_org = 0;
    phiso_org = 0;
    absiso_org = 0;
    reliso_org = 0;
 
    pileup = 0; 
    pileup0_005 = 0; 
    pileup005_01 = 0; 
    pileup01_015 = 0; 
    pileup015_02 = 0; 
    pileup02_025 = 0; 
    pileup025_03 = 0; 
    neutral = 0;
    neutral_org = 0;

    std::vector<reco::CandidatePtr> footprint; 
    for (unsigned int i = 0, n=muon.numberOfSourceCandidatePtrs(); i < n; ++i) {
        footprint.push_back(muon.sourceCandidatePtr(i));  
    }

    for (unsigned int i = 0, n = pfcands->size(); i < n; ++i) {
      const pat::PackedCandidate &pfc = (*pfcands)[i];
      const unsigned int absId = std::abs(pfc.pdgId());
      double dR = deltaR( muon.p4(), pfc.p4() );
      double dPhi = deltaPhi( muon.phi(), pfc.phi() );
      double dEta = abs( muon.eta() - pfc.eta() );
      double iso = pfc.pt();

      if( absId ==  211 ) {  // charged hadron
        if( dR < 0.3 ) {
          //testing footprint
          if(std::find(footprint.begin(), footprint.end(),  reco::CandidatePtr(pfcands,i)) != footprint.end()) {
             continue;
          } 
          //end footprint 
          if( pfc.fromPV() >= 2 ){ 
            chn++;
            chR = chR + dR; 
            chPhi = chPhi + dPhi; 
            chEta = chEta + dEta;
            chiso = chiso + iso;
            if( dR < 0.05 ){
              chiso0_005 = chiso0_005 + iso;
            }
            else if( dR < 0.1 ){
              chiso005_01 = chiso005_01 + iso;
            }
            else if( dR < 0.15 ){
              chiso01_015 = chiso01_015 + iso;
            }
            else if( dR < 0.2 ){
              chiso015_02 = chiso015_02 + iso;
            }
            else if( dR < 0.25 ){
              chiso02_025 = chiso02_025 + iso;
            }
            else{
              chiso025_03 = chiso025_03 + iso;
            }
            schiso = chiso0_005 + chiso005_01 + chiso01_015 + chiso015_02 + chiso02_025 + chiso025_03;
          
          }
          else{
            pileup = pileup + iso;
            if( dR < 0.05 ){
              pileup0_005 = pileup0_005 + iso;
            }
            else if( dR < 0.1 ){
              pileup005_01 = pileup005_01 + iso;
            }
            else if( dR < 0.15 ){
              pileup01_015 = pileup01_015 + iso;
            }
            else if( dR < 0.2 ){
              pileup015_02 = pileup015_02 + iso;
            }
            else if( dR < 0.25 ){
              pileup02_025 = pileup02_025 + iso;
            }
            else{
              pileup025_03 = pileup025_03 + iso;
            } 
          
          }
        }
      } 

      if( absId ==  130 ) {  // neutral hadron
        if( dR < 0.3 ) {
          nhn++;
          nhR = nhR + dR;
          nhPhi = nhPhi + dPhi; 
          nhEta = nhEta + dEta;
          nhiso = nhiso + iso;
            if( dR < 0.05 ){
              nhiso0_005 = nhiso0_005 + iso;
            }
            else if( dR < 0.1 ){
              nhiso005_01 = nhiso005_01 + iso;
            }
            else if( dR < 0.15 ){
              nhiso01_015 = nhiso01_015 + iso;
            }
            else if( dR < 0.2 ){
              nhiso015_02 = nhiso015_02 + iso;
            }
            else if( dR < 0.25 ){
              nhiso02_025 = nhiso02_025 + iso;
            }
            else{
              nhiso025_03 = nhiso025_03 + iso;
            }
        }
      }

      if( absId == 22 ) {  // photon 
        if( dR < 0.3 ) {
          phn++;
          phR = phR + dR;
          //cout << phR << " " << dR << endl;
          phPhi = phPhi + dPhi; 
          phEta = phEta + dEta;
          phiso = phiso + iso;
            if( dR < 0.05 ){
              phiso0_005 = phiso0_005 + iso;
            }
            else if( dR < 0.1 ){
              phiso005_01 = phiso005_01 + iso;
            }
            else if( dR < 0.15 ){
              phiso01_015 = phiso01_015 + iso;
            }
            else if( dR < 0.2 ){
              phiso015_02 = phiso015_02 + iso;
            }
            else if( dR < 0.25 ){
              phiso02_025 = phiso02_025 + iso;
            }
            else{
              phiso025_03 = phiso025_03 + iso;
            }
        }
      }
     
     neutral = nhiso + phiso;
     absiso = chiso + std::max(0.0, neutral-0.5*pileup);
     reliso = absiso / muon.pt();

    } //end of pf candidates loop

    chiso_org = muon.chargedHadronIso();
    puChiso_org = muon.puChargedHadronIso();
    nhiso_org = muon.neutralHadronIso();
    phiso_org = muon.photonIso();
    neutral_org = nhiso_org + phiso_org;
    absiso_org = chiso_org + std::max(0.0, neutral_org-0.5*puChiso_org);
    reliso_org = absiso_org / muon.pt();

    if ( chn != 0 ){
      chR = chR/chn;
      chPhi = chPhi/chn;
      chEta = chEta/chn;
    }

    if(nhn != 0){
      nhR = nhR/nhn;
      nhPhi = nhPhi/nhn;
      nhEta = nhEta/nhn;
    }

    if(phn != 0){
      phR = phR/phn;
      phPhi = phPhi/phn;
      phEta = phEta/phn;
    }


/*
    cout << "here 1" << endl;
    double chiso = muon.isoDeposit(pat::PfChargedHadronIso)->depositAndCountWithin(dRch_, vetos_ch).first;
    cout << "here 2" << endl;
    double nhiso = muon.isoDeposit(pat::PfNeutralHadronIso)->depositAndCountWithin(dRnh_, vetos_nh).first;
    cout << "here 3 " << endl;
    double phiso = muon.isoDeposit(pat::PfGammaIso)->depositAndCountWithin(dRph_, vetos_ph).first;
    cout << "here 4 " << endl;
    double chn = muon.isoDeposit(pat::PfChargedHadronIso)->depositAndCountWithin(dRch_, vetos_ch).second;
    cout << "here 5 " << endl;
    double nhn = muon.isoDeposit(pat::PfNeutralHadronIso)->depositAndCountWithin(dRnh_, vetos_nh).second;
    cout << "here 6 " << endl;
    double phn = muon.isoDeposit(pat::PfGammaIso)->depositAndCountWithin(dRph_, vetos_ph).second;
    cout << "here 7 " << endl;

*/
    //to study surrounding particles
    typedef reco::IsoDeposit::const_iterator IM;

//    cout << "pt = " << muon.pt() << endl;
  } 


  if( GoodMuon ) tree->Fill();

}


// ------------ method called once each job just before starting event loop  ------------
void IsoAnalyzer::beginJob(){

   tree->Branch("EVENT",&EVENT,"EVENT/i");
   tree->Branch("RUN"  ,&RUN  ,"RUN/i");
   tree->Branch("LUMI" ,&LUMI ,"LUMI/i");

   tree->Branch("pt" ,&pt ,"pt/d");
   tree->Branch("eta",&eta,"eta/d");
   tree->Branch("phi",&phi,"phi/d");

   tree->Branch("chIso0_005",&chiso0_005,"chiso0_005/d");
   tree->Branch("chIso005_01",&chiso005_01,"chiso005_01/d");
   tree->Branch("chIso01_015",&chiso01_015,"chiso01_015/d");
   tree->Branch("chIso015_02",&chiso015_02,"chiso015_02/d");
   tree->Branch("chIso02_025",&chiso02_025,"chiso02_025/d");
   tree->Branch("chIso025_03",&chiso025_03,"chiso025_03/d");
   tree->Branch("chIso",&chiso,"chiso/d");
   tree->Branch("schIso",&schiso,"schiso/d");

   tree->Branch("nhIso0_005",&nhiso0_005,"nhiso0_005/d");
   tree->Branch("nhIso005_01",&nhiso005_01,"nhiso005_01/d");
   tree->Branch("nhIso01_015",&nhiso01_015,"nhiso01_015/d");
   tree->Branch("nhIso015_02",&nhiso015_02,"nhiso015_02/d");
   tree->Branch("nhIso02_025",&nhiso02_025,"nhiso02_025/d");
   tree->Branch("nhIso025_03",&nhiso025_03,"nhiso025_03/d");
   tree->Branch("nhIso",&nhiso,"nhiso/d");

   tree->Branch("phIso0_005",&phiso0_005,"phiso0_005/d");
   tree->Branch("phIso005_01",&phiso005_01,"phiso005_01/d");
   tree->Branch("phIso01_015",&phiso01_015,"phiso01_015/d");
   tree->Branch("phIso015_02",&phiso015_02,"phiso015_02/d");
   tree->Branch("phIso02_025",&phiso02_025,"phiso02_025/d");
   tree->Branch("phIso025_03",&phiso025_03,"phiso025_03/d");
   tree->Branch("phIso",&phiso,"phiso/d");

   tree->Branch("Pileup0_005",&pileup0_005,"pileup0_005/d");
   tree->Branch("Pileup005_01",&pileup005_01,"pileup005_01/d");
   tree->Branch("Pileup01_015",&pileup01_015,"pileup01_015/d");
   tree->Branch("Pileup015_02",&pileup015_02,"pileup015_02/d");
   tree->Branch("Pileup02_025",&pileup02_025,"pileup02_025/d");
   tree->Branch("Pileup025_03",&pileup025_03,"pileup025_03/d");
   tree->Branch("Pileup",&pileup,"pileup/d");

   tree->Branch("abIso",&absiso,"absiso/d");
   tree->Branch("relIso",&reliso,"reliso/d");

   tree->Branch("chIso_org",&chiso_org,"chiso_org/d");
   tree->Branch("puChIso_org",&puChiso_org,"puChiso_org/d");
   tree->Branch("nhIso_org",&nhiso_org,"nhiso_org/d");
   tree->Branch("phIso_org",&phiso_org,"phiso_org/d");
   tree->Branch("absIso_org",&absiso_org,"absiso_org/d");
   tree->Branch("relIso_org",&reliso_org,"reliso_org/d");

   tree->Branch("chN",&chn,"chn/I");
   tree->Branch("nhN",&nhn,"nhn/I");
   tree->Branch("phN",&phn,"phn/I");

   tree->Branch("ch_dR",&chR,"chR/d");
   tree->Branch("nh_dR",&nhR,"nhR/d");
   tree->Branch("ph_dR",&phR,"phR/d");

   tree->Branch("ch_dPhi",&ch_dphi,"ch_dphi/d");
   tree->Branch("nh_dPhi",&nh_dphi,"nh_dphi/d");
   tree->Branch("ph_dPhi",&ph_dphi,"ph_dphi/d");

   tree->Branch("ch_dEta",&ch_deta,"ch_deta/d");
   tree->Branch("nh_dEta",&nh_deta,"nh_deta/d");
   tree->Branch("ph_dEta",&ph_deta,"ph_deta/d");


}

// ------------ method called once each job just after ending the event loop  ------------
void IsoAnalyzer::endJob() 
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void IsoAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);

  //Specify that only 'tracks' is allowed
  //To use, remove the default given above and uncomment below
  //ParameterSetDescription desc;
  //desc.addUntracked<edm::InputTag>("tracks","ctfWithMaterialTracks");
  //descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(IsoAnalyzer);
