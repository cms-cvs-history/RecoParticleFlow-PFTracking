#ifndef ConvBremGsfTrackMapProducer_H
#define ConvBremGsfTrackMapProducer_H

#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "TrackingTools/GsfTools/interface/MultiTrajectoryStateTransform.h"

/// \brief Abstract
/*!
\author Daniele Benedetti
\date April 2013
  ConvBremGsfTrackMapProducer produces a map between KF tracks from converted brems and GSF track. 
*/


class PFTrackTransformer;
class ConvBremTrackFinder;

class ConvBremGsfTrackMapProducer : public edm::EDProducer {
public:
  
  ///Constructor
  explicit ConvBremGsfTrackMapProducer(const edm::ParameterSet&);
  
  ///Destructor
  ~ConvBremGsfTrackMapProducer();
  
private:
  virtual void beginRun(edm::Run&,const edm::EventSetup&) ;
  virtual void endRun() ;
  
  ///Produce the Track collection
  virtual void produce(edm::Event&, const edm::EventSetup&);
  

  /// InputTags
  edm::InputTag trackCollection_;
  edm::InputTag gsfTrackLabel_;  
  edm::InputTag vtx_h;
  edm::InputTag pfEcalClusters_;

  /// track quality
  bool useQuality_;
  reco::TrackBase::TrackQuality trackQuality_;

  /// Conv Brem Finder
  bool useConvBremFinder_;
  double mvaConvBremFinderID_;
  std::string path_mvaWeightFileConvBrem_;

  /// Transformers  
  PFTrackTransformer *pfTransformer_; 
  MultiTrajectoryStateTransform  mtsTransform_;
  ConvBremTrackFinder *convBremTrackFinder_;

};
#endif
