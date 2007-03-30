#include "RecoParticleFlow/PFTracking/test/TestPreId.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "SimDataFormats/Track/interface/SimTrackContainer.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/ParticleFlowReco/interface/PFRecTrack.h"
#include "DataFormats/ParticleFlowReco/interface/PFRecTrackFwd.h"
#include "SimTracker/TrackAssociation/interface/TrackAssociatorBase.h"
#include "SimTracker/Records/interface/TrackAssociatorRecord.h"
#include "FWCore/Framework/interface/ESHandle.h"
//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//

using namespace edm;
using namespace std;
using namespace reco;
TestPreId::TestPreId(const edm::ParameterSet& iConfig):
  conf_(iConfig)
{
   //now do what ever initialization is needed
  refitLabel_= 
    iConfig.getParameter<InputTag>("RefitModuleLabel");
  mcLabel_ = iConfig.getParameter<InputTag>("SimHits");
  
  ttLabel_ = iConfig.getParameter<InputTag>("TrTruth");

  pfTrackLabel_= 
    iConfig.getParameter<InputTag>("PFRecTrackLabel");

  outputfile_ = conf_.getParameter<std::string>("OutputFile");
}


TestPreId::~TestPreId()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to for each event  ------------
void
TestPreId::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{


   Handle<SimTrackContainer> simTracks;
   iEvent.getByLabel(mcLabel_,simTracks);
   const SimTrackContainer simTC = *(simTracks.product());

   Handle<TrackingParticleCollection>  TPCollectionH ;
   iEvent.getByLabel(ttLabel_,TPCollectionH);


   Handle<TrackCollection> tkCollection;
   iEvent.getByLabel(refitLabel_, tkCollection);
   const  reco::TrackCollection  tC = *(tkCollection.product());


   Handle<PFRecTrackCollection> thePfRecTrackCollection;
   iEvent.getByLabel(pfTrackLabel_,thePfRecTrackCollection);
   const PFRecTrackCollection PfRTkColl = *(thePfRecTrackCollection.product());

   ESHandle<TrackAssociatorBase> theHitsAssociator;
   iSetup.get<TrackAssociatorRecord>().get("TrackAssociatorByHits",theHitsAssociator);
   TrackAssociatorBase * associatorByHits = (TrackAssociatorBase *) theHitsAssociator.product();

   reco::RecoToSimCollection p = 
     associatorByHits->associateRecoToSim (tkCollection,TPCollectionH,&iEvent );


   PFRecTrackCollection::const_iterator ipftrak;
   for (ipftrak=PfRTkColl.begin(); ipftrak!=PfRTkColl.end(); ipftrak++){
     int itk=(*ipftrak).recTrackId();

     TrackRef track(tkCollection, itk);
     try{ 
       vector<pair<TrackingParticleRef, double> > tp = p[track];
       
       if (tp.size()>0){
	 TrackingParticleRef tpr =tp[0].first;
	 int code= (*tpr).pdgId();
	 if (abs(code)==11){
	   hpt_elec->Fill(tC[itk].pt());
	   heta_elec->Fill(fabs(tC[itk].eta()));
	   if((*ipftrak).algoType()==PFRecTrack::KF_ELCAND){
	     hpt_elec_sel->Fill(tC[itk].pt());
	     heta_elec_sel->Fill(fabs(tC[itk].eta()));
	   }
	 }
	 if (abs(code)==211){
	   hpt_pion->Fill(tC[itk].pt());
	   heta_pion->Fill(fabs(tC[itk].eta()));
	   if((*ipftrak).algoType()==PFRecTrack::KF_ELCAND){
	     hpt_pion_sel->Fill(tC[itk].pt());
	     heta_pion_sel->Fill(fabs(tC[itk].eta()));
	   }
	 }

       }
     }
     catch (Exception event) {
     }
   }
   

   
}


// ------------ method called once each job just before starting event loop  ------------
void 
TestPreId::beginJob(const edm::EventSetup&)
{
  tf1 = new TFile(outputfile_.c_str(), "RECREATE");
  hpt_pion = new TH1F("hpt_pion","hpt_pion",50,0.,50.);
  hpt_pion_sel = new TH1F("hpt_pion_sel","hpt_pion_sel",50,0.,50.);
  hpt_elec = new TH1F("hpt_elec","hpt_elec",50,0.,50.);
  hpt_elec_sel = new TH1F("hpt_elec_sel","hpt_elec_sel",50,0.,50.);
  heta_pion = new TH1F("heta_pion","heta_pion",30,0.,2.4);
  heta_pion_sel = new TH1F("heta_pion_sel","heta_pion_sel",30,0.,2.4);
  heta_elec = new TH1F("heta_elec","heta_elec",30,0.,2.4);
  heta_elec_sel = new TH1F("heta_elec_sel","heta_elec_sel",30,0.,2.4);
}

// ------------ method called once each job just after ending the event loop  ------------
void 
TestPreId::endJob() {
  tf1->Write();
  tf1->Close();
}

