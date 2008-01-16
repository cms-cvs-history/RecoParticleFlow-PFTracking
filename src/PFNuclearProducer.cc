#include <memory>
#include "RecoParticleFlow/PFTracking/interface/PFNuclearProducer.h"
#include "RecoParticleFlow/PFTracking/interface/PFTrackTransformer.h"
#include "DataFormats/ParticleFlowReco/interface/PFNuclearInteraction.h"
#include "TrackingTools/PatternTools/interface/Trajectory.h"
using namespace std;
using namespace edm;
PFNuclearProducer::PFNuclearProducer(const ParameterSet& iConfig):
  pfTransformer_(0)
{
  produces<reco::PFRecTrackCollection>();
  produces<reco::PFNuclearInteractionCollection>();

  nuclearContainers_ = 
    iConfig.getParameter< vector < InputTag > >("nuclearColList");
}

PFNuclearProducer::~PFNuclearProducer()
{
  delete pfTransformer_;
}

void
PFNuclearProducer::produce(Event& iEvent, const EventSetup& iSetup)
{
  typedef reco::NuclearInteraction::trackRef_iterator trackRef_iterator;

  //create the empty collections 
  auto_ptr< reco::PFNuclearInteractionCollection > 
    pfNuclearColl (new reco::PFNuclearInteractionCollection);
  auto_ptr< reco::PFRecTrackCollection > 
    pfNuclearRecTrackColl (new reco::PFRecTrackCollection);
  
  reco::PFRecTrackRefProd pfTrackRefProd = iEvent.getRefBeforePut<reco::PFRecTrackCollection>();
  int hid=0;

  // loop on the nuclear interaction collections
  for (uint istr=0; istr<nuclearContainers_.size();istr++){
    
    Handle<reco::NuclearInteractionCollection> nuclCollH;
    iEvent.getByLabel(nuclearContainers_[istr], nuclCollH);
    const reco::NuclearInteractionCollection& nuclColl = *(nuclCollH.product());

    // loop on all NuclearInteraction 
    for( uint icoll=0; icoll < nuclColl.size(); icoll++) {
      reco::PFRecTrackRefVector pfRecTkcoll;

      // convert the secondary tracks
      for(trackRef_iterator it = nuclColl[icoll].secondaryTracks_begin(); it!=nuclColl[icoll].secondaryTracks_end(); it++){
        reco::PFRecTrack pftrack( (*it)->charge(), 
       				reco::PFRecTrack::KF, 
       				it->key(), (reco::TrackRef)((*it).castTo<reco::TrackRef>()) );
        Trajectory FakeTraj;
        bool valid = pfTransformer_->addPoints( pftrack, **it, FakeTraj);
        if(valid) {
	  pfRecTkcoll.push_back(reco::PFRecTrackRef( pfTrackRefProd, hid++ ));	
          pfNuclearRecTrackColl->push_back(pftrack);
        }
      }
      reco::NuclearInteractionRef niRef(nuclCollH, icoll);
      pfNuclearColl->push_back( reco::PFNuclearInteraction( niRef, pfRecTkcoll ));
    }
  }
  iEvent.put(pfNuclearRecTrackColl);
  iEvent.put(pfNuclearColl);
}

// ------------ method called once each job just before starting event loop  ------------
void 
PFNuclearProducer::beginJob(const EventSetup& iSetup)
{
  pfTransformer_= new PFTrackTransformer();
  pfTransformer_->OnlyProp();
}

// ------------ method called once each job just after ending the event loop  ------------
void 
PFNuclearProducer::endJob() {
}