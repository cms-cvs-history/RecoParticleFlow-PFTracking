#include "RecoParticleFlow/PFTracking/interface/VertexFilter.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackExtra.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "FWCore/ParameterSet/interface/InputTag.h"
#include "TrackingTools/TransientTrackingRecHit/interface/TransientTrackingRecHit.h"
#include "TrackingTools/PatternTools/interface/Trajectory.h"
#include "TrackingTools/PatternTools/interface/TrajTrackAssociation.h"
//
// class decleration
//




VertexFilter::VertexFilter(const edm::ParameterSet& iConfig):
  conf_(iConfig)
{
  dist = conf_.getParameter<double>("DistFromVertex");
  produces<reco::TrackCollection>();
  produces<TrackingRecHitCollection>();
  produces<reco::TrackExtraCollection>();
  produces<std::vector<Trajectory> >();
  produces<TrajTrackAssociationCollection>();

}


VertexFilter::~VertexFilter()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to produce the data  ------------
void
VertexFilter::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;
  using namespace reco;
  using namespace std;

  auto_ptr<TrackCollection> selTracks(new TrackCollection);
  auto_ptr<TrackingRecHitCollection> selHits(new TrackingRecHitCollection);
  auto_ptr<TrackExtraCollection> selTrackExtras(new TrackExtraCollection);
  auto_ptr<vector<Trajectory> > outputTJ(new vector<Trajectory> );
  auto_ptr<TrajTrackAssociationCollection> trajTrackMap( new TrajTrackAssociationCollection() );

  Ref<TrackExtraCollection>::key_type idx = 0;
  Ref<TrackCollection>::key_type iTkRef = 0;

  TrackExtraRefProd rTrackExtras = iEvent.getRefBeforePut<TrackExtraCollection>();
  TrackingRecHitRefProd rHits = iEvent.getRefBeforePut<TrackingRecHitCollection>();
  Handle<std::vector<Trajectory> > TrajectoryCollection;
  Handle<TrajTrackAssociationCollection> assoMap;

  InputTag tkTag             = conf_.getParameter<InputTag>("recTracks");
  InputTag vtxTag            = conf_.getParameter<InputTag>("recVertices");
  
  Handle<TrackCollection> tkCollection;
  Handle<VertexCollection> vtxCollection;
  
  iEvent.getByLabel( tkTag             , tkCollection);
  iEvent.getByLabel( vtxTag            , vtxCollection);
  TrackCollection  tC = *(tkCollection.product());
  const VertexCollection  vxC = *(vtxCollection.product());
  iEvent.getByLabel(tkTag,TrajectoryCollection);
  iEvent.getByLabel(tkTag,assoMap);
  // TrackCollection::iterator itc;
  VertexCollection::const_iterator ivxc;
  for(TrajTrackAssociationCollection::const_iterator it = assoMap->begin();
      it != assoMap->end(); ++it){
    const Ref<vector<Trajectory> > traj = it->key;
    const reco::TrackRef itc = it->val;
    
    math::XYZPoint gp1=(*itc).vertex();
    bool hasAVertex=false;
    
    for (ivxc= vxC.begin();ivxc!=vxC.end();ivxc++){
      math::XYZPoint gp2=(*ivxc).position();
	if ((gp1-gp2).Mag2()<dist) hasAVertex=true;
    }
    
    if (hasAVertex){
      Track track =(*itc);
      selTracks->push_back( track );
      outputTJ->push_back(*traj);
  

      iTkRef++;
      TrackExtraRef teref=TrackExtraRef ( rTrackExtras, idx ++ );
      track.setExtra( teref );

      TrackExtra & tx =const_cast<TrackExtra &>((*(*itc).extra()));

      
      selTrackExtras->push_back(tx);
      trackingRecHit_iterator irhit=(*itc).recHitsBegin();
      size_t i = 0;
      for (; irhit!=(*itc).recHitsEnd(); irhit++){
	TrackingRecHit *hit=(*(*irhit)).clone();
	track.setHitPattern( * hit, i ++ );
	selHits->push_back(hit );
      }
    }
  }
//   for (int i=0; i<10;i++){
//     trajTrackMap->insert(edm::Ref<std::vector<Trajectory> > (outputTJ, i),
// 			 edm::Ref<reco::TrackCollection>    (selTracks, i)); 
//   }
  iEvent.put(selTracks);
  iEvent.put(selTrackExtras);    
  iEvent.put(selHits);
  iEvent.put(outputTJ);
  iEvent.put( trajTrackMap );
}

// ------------ method called once each job just before starting event loop  ------------
void 
VertexFilter::beginJob(const edm::EventSetup&)
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
VertexFilter::endJob() {
}


