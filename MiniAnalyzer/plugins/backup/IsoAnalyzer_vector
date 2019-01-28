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

       double dRch_;
       double dRnh_;
       double dRph_;
       double vetoCHThreshold_; 
       double vetoNHThreshold_;
       double vetoPHThreshold_;
       double vetoCHdR_;
       double vetoNHdR_;
       double vetoPHdR_;
       vector<double> vetoCHre_;
       vector<double> vetoNHre_;
       vector<double> vetoPHre_;

       TTree *tree;

       unsigned int EVENT;
       unsigned int RUN;
       unsigned int LUMI;

       std::vector<double>* eta;
       std::vector<double>* phi;
       std::vector<double>* pt;

       std::vector<double>* chIso0_005;
       std::vector<double>* chIso005_01;
       std::vector<double>* chIso01_015;
       std::vector<double>* chIso015_02;
       std::vector<double>* chIso02_025;
       std::vector<double>* chIso025_03;
       std::vector<double>* chIso;
       std::vector<double>* schIso;

       std::vector<double>* nhIso0_005;
       std::vector<double>* nhIso005_01;
       std::vector<double>* nhIso01_015;
       std::vector<double>* nhIso015_02;
       std::vector<double>* nhIso02_025;
       std::vector<double>* nhIso025_03;
       std::vector<double>* nhIso;

       std::vector<double>* phIso0_005;
       std::vector<double>* phIso005_01;
       std::vector<double>* phIso01_015;
       std::vector<double>* phIso015_02;
       std::vector<double>* phIso02_025;
       std::vector<double>* phIso025_03;
       std::vector<double>* phIso;

       std::vector<double>* absIso;
       std::vector<double>* relIso;

       std::vector<double>* Pileup;
       std::vector<double>* Pileup0_005;
       std::vector<double>* Pileup005_01;
       std::vector<double>* Pileup01_015;
       std::vector<double>* Pileup015_02;
       std::vector<double>* Pileup02_025;
       std::vector<double>* Pileup025_03;

       std::vector<double>* chIso_org;
       std::vector<double>* puChIso_org;
       std::vector<double>* nhIso_org;
       std::vector<double>* phIso_org;
       std::vector<double>* absIso_org;
       std::vector<double>* relIso_org;

       std::vector<double>* chN;
       std::vector<double>* nhN;
       std::vector<double>* phN;

       std::vector<double>* ch_dR;
       std::vector<double>* ch_th;
       std::vector<double>* ch_dphi;
       std::vector<double>* ch_deta;

       std::vector<double>* nh_dR;
       std::vector<double>* nh_th;
       std::vector<double>* nh_dphi;
       std::vector<double>* nh_deta;

       std::vector<double>* ph_dR;
       std::vector<double>* ph_th;
       std::vector<double>* ph_dphi;
       std::vector<double>* ph_deta;

//       std::vector<double>* chargedhadron_pt;
//       std::vector<double>* chargedhadron_eta;
//       std::vector<double>* chargedhadron_phi;
 
       TH1D * IsoNtuple; 
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
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

   pt = new std::vector<double>();
   eta = new std::vector<double>();
   phi = new std::vector<double>();

   chIso0_005 = new std::vector<double>();
   chIso005_01 = new std::vector<double>();
   chIso01_015 = new std::vector<double>();
   chIso015_02 = new std::vector<double>();
   chIso02_025 = new std::vector<double>();
   chIso025_03 = new std::vector<double>();
   chIso = new std::vector<double>();
   schIso = new std::vector<double>();

   nhIso0_005 = new std::vector<double>();
   nhIso005_01 = new std::vector<double>();
   nhIso01_015 = new std::vector<double>();
   nhIso015_02 = new std::vector<double>();
   nhIso02_025 = new std::vector<double>();
   nhIso025_03 = new std::vector<double>();
   nhIso = new std::vector<double>();

   phIso0_005 = new std::vector<double>();
   phIso005_01 = new std::vector<double>();
   phIso01_015 = new std::vector<double>();
   phIso015_02 = new std::vector<double>();
   phIso02_025 = new std::vector<double>();
   phIso025_03 = new std::vector<double>();
   phIso = new std::vector<double>();

   absIso = new std::vector<double>();
   relIso = new std::vector<double>();

   Pileup = new std::vector<double>();
   Pileup0_005 = new std::vector<double>();
   Pileup005_01 = new std::vector<double>();
   Pileup01_015 = new std::vector<double>();
   Pileup015_02 = new std::vector<double>();
   Pileup02_025 = new std::vector<double>();
   Pileup025_03 = new std::vector<double>();

   chIso_org = new std::vector<double>();
   puChIso_org = new std::vector<double>();
   nhIso_org = new std::vector<double>();
   phIso_org = new std::vector<double>();
   absIso_org = new std::vector<double>();
   relIso_org = new std::vector<double>();

   chN = new std::vector<double>();
   nhN = new std::vector<double>();
   phN = new std::vector<double>();

   ch_dR = new std::vector<double>();
   ch_th = new std::vector<double>();
   ch_deta = new std::vector<double>();
   ch_dphi = new std::vector<double>();

   nh_dR = new std::vector<double>();
   nh_th = new std::vector<double>();
   nh_deta = new std::vector<double>();
   nh_dphi = new std::vector<double>();

   ph_dR = new std::vector<double>();
   ph_th = new std::vector<double>();
   ph_deta = new std::vector<double>();
   ph_dphi = new std::vector<double>();

//   chargedhadron_pt = new std::vector<double>();
//   chargedhadron_eta = new std::vector<double>();
//   chargedhadron_phi = new std::vector<double>();

   IsoNtuple = fs->make<TH1D>("tracks" , "Tracks" , 100 , 0 , 5000 );

}


IsoAnalyzer::~IsoAnalyzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
IsoAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  using namespace reco;
  using namespace isodeposit;
  using namespace edm;
 
  pt->clear();
  eta->clear();
  phi->clear();

  chIso0_005->clear();
  chIso005_01->clear();
  chIso01_015->clear();
  chIso015_02->clear();
  chIso02_025->clear();
  chIso025_03->clear();
  chIso->clear();
  schIso->clear();

  nhIso0_005->clear();
  nhIso005_01->clear();
  nhIso01_015->clear();
  nhIso015_02->clear();
  nhIso02_025->clear();
  nhIso025_03->clear();
  nhIso->clear();

  phIso0_005->clear();
  phIso005_01->clear();
  phIso01_015->clear();
  phIso015_02->clear();
  phIso02_025->clear();
  phIso025_03->clear();
  phIso->clear();

  absIso->clear();
  relIso->clear();

  Pileup->clear();
  Pileup0_005->clear();
  Pileup005_01->clear();
  Pileup01_015->clear();
  Pileup015_02->clear();
  Pileup02_025->clear();
  Pileup025_03->clear();

  chIso_org->clear();
  puChIso_org->clear();
  nhIso_org->clear();
  phIso_org->clear();
  absIso_org->clear();
  relIso_org->clear();

  chN->clear();
  nhN->clear();
  phN->clear();

  ch_dR->clear();
  ch_th->clear();
  ch_deta->clear();
  ch_dphi->clear();

  nh_dR->clear();
  nh_th->clear();
  nh_deta->clear();
  nh_dphi->clear();

  ph_dR->clear();
  ph_th->clear();
  ph_deta->clear();
  ph_dphi->clear();

//  chargedhadron_pt->clear();
//  chargedhadron_eta->clear();
//  chargedhadron_phi->clear();

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

    pt->push_back(muon.pt());
    eta->push_back(muon.eta());
    phi->push_back(muon.phi());
 
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

    float chn = 0;
    float nhn = 0;
    float phn = 0;

    float chR = 0;
    float nhR = 0;
    float phR = 0;
  
    float chPhi = 0;
    float nhPhi = 0;
    float phPhi = 0;

    float chEta = 0;
    float nhEta = 0;
    float phEta = 0;
 
    float chiso0_005 = 0;
    float chiso005_01 = 0;
    float chiso01_015 = 0;
    float chiso015_02 = 0;
    float chiso02_025 = 0;
    float chiso025_03 = 0;
    float chiso = 0;
    float schiso = 0;

    float nhiso0_005 = 0;
    float nhiso005_01 = 0;
    float nhiso01_015 = 0;
    float nhiso015_02 = 0;
    float nhiso02_025 = 0;
    float nhiso025_03 = 0;
    float nhiso = 0;

    float phiso0_005 = 0;
    float phiso005_01 = 0;
    float phiso01_015 = 0;
    float phiso015_02 = 0;
    float phiso02_025 = 0;
    float phiso025_03 = 0;
    float phiso = 0;

    float absiso = 0;
    float reliso = 0;
   
    float chiso_org = 0;
    float puChiso_org = 0;
    float nhiso_org = 0;
    float phiso_org = 0;
    float absiso_org = 0;
    float reliso_org = 0;
 
    float pileup = 0; 
    float pileup0_005 = 0; 
    float pileup005_01 = 0; 
    float pileup01_015 = 0; 
    float pileup015_02 = 0; 
    float pileup02_025 = 0; 
    float pileup025_03 = 0; 
    float neutral = 0;
    float neutral_org = 0;

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
      if( absId == 22) {  // photon 
        if( dR < 0.3 ) {
          phn++;
          phR = phR + dR;
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

    chN->push_back(chn);
    nhN->push_back(nhn);
    phN->push_back(phn);
  
    chIso0_005->push_back(chiso0_005);   
    chIso005_01->push_back(chiso005_01);   
    chIso01_015->push_back(chiso01_015);   
    chIso015_02->push_back(chiso015_02);   
    chIso02_025->push_back(chiso02_025);   
    chIso025_03->push_back(chiso025_03);   
    chIso->push_back(chiso);  
    schIso->push_back(schiso);  
 
    nhIso0_005->push_back(nhiso0_005);   
    nhIso005_01->push_back(nhiso005_01);   
    nhIso01_015->push_back(nhiso01_015);   
    nhIso015_02->push_back(nhiso015_02);   
    nhIso02_025->push_back(nhiso02_025);   
    nhIso025_03->push_back(nhiso025_03);   
    nhIso->push_back(nhiso);   

    phIso0_005->push_back(phiso0_005);   
    phIso005_01->push_back(phiso005_01);   
    phIso01_015->push_back(phiso01_015);   
    phIso015_02->push_back(phiso015_02);   
    phIso02_025->push_back(phiso02_025);   
    phIso025_03->push_back(phiso025_03);   
    phIso->push_back(phiso); 
 
    Pileup->push_back(pileup); 
    Pileup0_005->push_back(pileup0_005); 
    Pileup005_01->push_back(pileup005_01); 
    Pileup01_015->push_back(pileup01_015); 
    Pileup015_02->push_back(pileup015_02); 
    Pileup02_025->push_back(pileup02_025); 
    Pileup025_03->push_back(pileup025_03); 
    absIso->push_back(absiso);   
    relIso->push_back(reliso);    

    chiso_org = muon.chargedHadronIso();
    puChiso_org = muon.puChargedHadronIso();
    nhiso_org = muon.neutralHadronIso();
    phiso_org = muon.photonIso();
    neutral_org = nhiso_org + phiso_org;
    absiso_org = chiso_org + std::max(0.0, neutral_org-0.5*puChiso_org);
    reliso_org = absiso_org / muon.pt();

    chIso_org->push_back(chiso_org);
    puChIso_org->push_back(puChiso_org);
    nhIso_org->push_back(nhiso_org);
    phIso_org->push_back(phiso_org);
    relIso_org->push_back(reliso_org);
    absIso_org->push_back(absiso_org);

    //dR
    float chR_avg = 0;
    if ( chn != 0 ) chR_avg = chR/chn;
    ch_dR->push_back(chR_avg);
    nh_dR->push_back(nhR/nhn);
    ph_dR->push_back(phR/phn);

    //dphi
    float chPhi_avg = 0;
    if ( chn != 0 ) chPhi_avg = chPhi/chn;
    ch_dphi->push_back(chPhi_avg);
    nh_dphi->push_back(nhPhi/nhn);
    ph_dphi->push_back(phPhi/phn);

    //deta
    float chEta_avg = 0;
    if ( chn != 0 ) chEta_avg = chEta/chn;
    ch_deta->push_back(chEta_avg);
    nh_deta->push_back(nhEta/nhn);
    ph_deta->push_back(phEta/phn);

//    typedef pat::PackedCandidateCollection::const_iterator CI;
//    for(CI ci = pfcands->begin(); ci!=pfcands_->end(); ++ci) {
//      const pat::PackedCandidate& pfc = *ci;
//      if( pfc.particleId() == 1 ) {
//        chargedhadron_pt->push_back( pfc.pt());
//        chargedhadron_eta->push_back( pfc.eta());
//        chargedhadron_phi->push_back( pfc.phi());
//      }
//    }

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


    chN->push_back(chn);
    nhN->push_back(nhn);
    phN->push_back(phn);    
    chIso->push_back(chiso);
    nhIso->push_back(nhiso);
    phIso->push_back(phiso);

*/
    //to study surrounding particles
    typedef reco::IsoDeposit::const_iterator IM;

//    cout << "pt = " << muon.pt() << endl;
  } 


  if( GoodMuon ) tree->Fill();

}


// ------------ method called once each job just before starting event loop  ------------
void 
IsoAnalyzer::beginJob()
{

   tree->Branch("EVENT",&EVENT,"EVENT/i");
   tree->Branch("RUN",&RUN,"RUN/i");
   tree->Branch("LUMI",&LUMI,"LUMI/i");

   tree->Branch("pt","std::vector<double>",&pt);
   tree->Branch("eta","std::vector<double>",&eta);
   tree->Branch("phi","std::vector<double>",&phi);

   tree->Branch("chIso0_005","std::vector<double>",&chIso0_005);
   tree->Branch("chIso005_01","std::vector<double>",&chIso005_01);
   tree->Branch("chIso01_015","std::vector<double>",&chIso01_015);
   tree->Branch("chIso015_02","std::vector<double>",&chIso015_02);
   tree->Branch("chIso02_025","std::vector<double>",&chIso02_025);
   tree->Branch("chIso025_03","std::vector<double>",&chIso025_03);
   tree->Branch("chIso","std::vector<double>",&chIso);
   tree->Branch("schIso","std::vector<double>",&schIso);

   tree->Branch("nhIso0_005","std::vector<double>",&nhIso0_005);
   tree->Branch("nhIso005_01","std::vector<double>",&nhIso005_01);
   tree->Branch("nhIso01_015","std::vector<double>",&nhIso01_015);
   tree->Branch("nhIso015_02","std::vector<double>",&nhIso015_02);
   tree->Branch("nhIso02_025","std::vector<double>",&nhIso02_025);
   tree->Branch("nhIso025_03","std::vector<double>",&nhIso025_03);
   tree->Branch("nhIso","std::vector<double>",&nhIso);

   tree->Branch("phIso0_005","std::vector<double>",&phIso0_005);
   tree->Branch("phIso005_01","std::vector<double>",&phIso005_01);
   tree->Branch("phIso01_015","std::vector<double>",&phIso01_015);
   tree->Branch("phIso015_02","std::vector<double>",&phIso015_02);
   tree->Branch("phIso02_025","std::vector<double>",&phIso02_025);
   tree->Branch("phIso025_03","std::vector<double>",&phIso025_03);
   tree->Branch("phIso","std::vector<double>",&phIso);

   tree->Branch("Pileup","std::vector<double>",&Pileup);
   tree->Branch("Pileup0_005","std::vector<double>",&Pileup0_005);
   tree->Branch("Pileup005_01","std::vector<double>",&Pileup005_01);
   tree->Branch("Pileup01_015","std::vector<double>",&Pileup01_015);
   tree->Branch("Pileup015_02","std::vector<double>",&Pileup015_02);
   tree->Branch("Pileup02_025","std::vector<double>",&Pileup02_025);
   tree->Branch("Pileup025_03","std::vector<double>",&Pileup025_03);
   tree->Branch("absIso","std::vector<double>",&absIso);
   tree->Branch("relIso","std::vector<double>",&relIso);

   tree->Branch("chIso_org","std::vector<double>",&chIso_org);
   tree->Branch("puChIso_org","std::vector<double>",&puChIso_org);
   tree->Branch("nhIso_org","std::vector<double>",&nhIso_org);
   tree->Branch("phIso_org","std::vector<double>",&phIso_org);
   tree->Branch("absIso_org","std::vector<double>",&absIso_org);
   tree->Branch("relIso_org","std::vector<double>",&relIso_org);

   tree->Branch("chN","std::vector<double>",&chN);
   tree->Branch("nhN","std::vector<double>",&nhN);
   tree->Branch("phN","std::vector<double>",&phN);

   tree->Branch("ch_dR","std::vector<double>",&ch_dR);
   tree->Branch("ch_th","std::vector<double>",&ch_th);
   tree->Branch("ch_deta","std::vector<double>",&ch_deta);
   tree->Branch("ch_dphi","std::vector<double>",&ch_dphi);

   tree->Branch("nh_dR","std::vector<double>",&nh_dR);
   tree->Branch("nh_th","std::vector<double>",&nh_th);
   tree->Branch("nh_deta","std::vector<double>",&nh_deta);
   tree->Branch("nh_dphi","std::vector<double>",&nh_dphi);

   tree->Branch("ph_dR","std::vector<double>",&ph_dR);
   tree->Branch("ph_th","std::vector<double>",&ph_th);
   tree->Branch("ph_deta","std::vector<double>",&ph_deta);
   tree->Branch("ph_dphi","std::vector<double>",&ph_dphi);

//   tree->Branch("chargedhadron_pt","std::vector<double>",&chargedhadron_pt);
//   tree->Branch("chargedhadron_phi","std::vector<double>",&chargedhadron_phi);
//   tree->Branch("chargedhadron_eta","std::vector<double>",&chargedhadron_eta);


}

// ------------ method called once each job just after ending the event loop  ------------
void 
IsoAnalyzer::endJob() 
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
IsoAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
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
