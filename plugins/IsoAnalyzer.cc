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

       std::vector<double>* chIso0_1;
       std::vector<double>* chIso1_2;
       std::vector<double>* chIso2_3;
       std::vector<double>* chIso;
       std::vector<double>* schIso;

       std::vector<double>* nhIso0_1;
       std::vector<double>* nhIso1_2;
       std::vector<double>* nhIso2_3;
       std::vector<double>* nhIso;

       std::vector<double>* phIso0_1;
       std::vector<double>* phIso1_2;
       std::vector<double>* phIso2_3;
       std::vector<double>* phIso;

       std::vector<double>* absIso;
       std::vector<double>* relIso;
       std::vector<double>* Pileup;

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

       std::vector<double>* chargedhadron_pt;
       std::vector<double>* chargedhadron_eta;
       std::vector<double>* chargedhadron_phi;
 
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

   chIso0_1 = new std::vector<double>();
   chIso1_2 = new std::vector<double>();
   chIso2_3 = new std::vector<double>();
   chIso = new std::vector<double>();
   schIso = new std::vector<double>();

   nhIso0_1 = new std::vector<double>();
   nhIso1_2 = new std::vector<double>();
   nhIso2_3 = new std::vector<double>();
   nhIso = new std::vector<double>();

   phIso0_1 = new std::vector<double>();
   phIso1_2 = new std::vector<double>();
   phIso2_3 = new std::vector<double>();
   phIso = new std::vector<double>();

   absIso = new std::vector<double>();
   relIso = new std::vector<double>();
   Pileup = new std::vector<double>();

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

   chargedhadron_pt = new std::vector<double>();
   chargedhadron_eta = new std::vector<double>();
   chargedhadron_phi = new std::vector<double>();

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

  chIso0_1->clear();
  chIso1_2->clear();
  chIso2_3->clear();
  chIso->clear();
  schIso->clear();

  nhIso0_1->clear();
  nhIso1_2->clear();
  nhIso2_3->clear();
  nhIso->clear();

  phIso0_1->clear();
  phIso1_2->clear();
  phIso2_3->clear();
  phIso->clear();
  absIso->clear();
  relIso->clear();
  Pileup->clear();

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

  chargedhadron_pt->clear();
  chargedhadron_eta->clear();
  chargedhadron_phi->clear();

  edm::Handle<edm::View<pat::Muon> > muons;
  iEvent.getByToken(muonLabel_, muons);

  //edm::Handle<reco::PFCandidateCollection> pfCandidates_; 
  //iEvent.getByToken(pfCandidateLabel_, pfCandidates_);

  edm::Handle<pat::PackedCandidateCollection> pfcands;
  iEvent.getByToken(pfSrc_, pfcands);

  EVENT  = iEvent.id().event();
  RUN    = iEvent.id().run();
  LUMI   = iEvent.id().luminosityBlock();

  for (unsigned i = 0; i != muons->size(); i++) {
    pat::Muon muon = muons->at(i);

    if( muon.pt() < 20) continue;
    if( muon.eta() < -2.1) continue;
    if( muon.eta() > 2.1) continue;

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
 
    float chiso0_1 = 0;
    float chiso1_2 = 0;
    float chiso2_3 = 0;
    float chiso = 0;
    float schiso = 0;

    float nhiso0_1 = 0;
    float nhiso1_2 = 0;
    float nhiso2_3 = 0;
    float nhiso = 0;

    float phiso0_1 = 0;
    float phiso1_2 = 0;
    float phiso2_3 = 0;
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
            if( dR < 0.1 ){
              chiso0_1 = chiso0_1 + iso;
            }
            else if( dR < 0.2 ){
              chiso1_2 = chiso1_2 + iso;
            }
            else{
              chiso2_3 = chiso2_3 + iso;
            }
            schiso = chiso0_1 + chiso1_2 + chiso2_3;

          }
          else{
            pileup = pileup + iso;
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
            if( dR < 0.1 ){
              nhiso0_1 = nhiso0_1 + iso;
            }
            else if( dR < 0.2 ){
              nhiso1_2 = nhiso1_2 + iso;
            }
            else{
              nhiso2_3 = nhiso2_3 + iso;
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
            if( dR < 0.1 ){
              phiso0_1 = phiso0_1 + iso;
            }
            else if( dR < 0.2 ){
              phiso1_2 = phiso1_2 + iso;
            }
            else{
              phiso2_3 = phiso2_3 + iso;
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
  
    chIso0_1->push_back(chiso0_1);   
    chIso1_2->push_back(chiso1_2);   
    chIso2_3->push_back(chiso2_3);   
    chIso->push_back(chiso);  
    schIso->push_back(schiso);  
 
    nhIso0_1->push_back(nhiso0_1);   
    nhIso1_2->push_back(nhiso1_2);   
    nhIso2_3->push_back(nhiso2_3);   
    nhIso->push_back(nhiso);   

    phIso0_1->push_back(phiso0_1);   
    phIso1_2->push_back(phiso1_2);   
    phIso2_3->push_back(phiso2_3);   
    phIso->push_back(phiso); 
 
    Pileup->push_back(pileup); 
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


  tree->Fill();

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

   tree->Branch("chIso0_1","std::vector<double>",&chIso0_1);
   tree->Branch("chIso1_2","std::vector<double>",&chIso1_2);
   tree->Branch("chIso2_3","std::vector<double>",&chIso2_3);
   tree->Branch("chIso","std::vector<double>",&chIso);
   tree->Branch("schIso","std::vector<double>",&schIso);

   tree->Branch("nhIso0_1","std::vector<double>",&nhIso0_1);
   tree->Branch("nhIso1_2","std::vector<double>",&nhIso1_2);
   tree->Branch("nhIso2_3","std::vector<double>",&nhIso2_3);
   tree->Branch("nhIso","std::vector<double>",&nhIso);

   tree->Branch("phIso0_1","std::vector<double>",&phIso0_1);
   tree->Branch("phIso1_2","std::vector<double>",&phIso1_2);
   tree->Branch("phIso2_3","std::vector<double>",&phIso2_3);
   tree->Branch("phIso","std::vector<double>",&phIso);

   tree->Branch("Pileup","std::vector<double>",&Pileup);
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

   tree->Branch("chargedhadron_pt","std::vector<double>",&chargedhadron_pt);
   tree->Branch("chargedhadron_phi","std::vector<double>",&chargedhadron_phi);
   tree->Branch("chargedhadron_eta","std::vector<double>",&chargedhadron_eta);


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
