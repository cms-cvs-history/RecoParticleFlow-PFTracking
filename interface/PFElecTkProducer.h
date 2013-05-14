#ifndef PFElecTkProducer_H
#define PFElecTkProducer_H

#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "DataFormats/ParticleFlowReco/interface/GsfPFRecTrackFwd.h"
#include "DataFormats/ParticleFlowReco/interface/GsfPFRecTrack.h"
#include "DataFormats/ParticleFlowReco/interface/PFRecTrack.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "TrackingTools/GsfTools/interface/MultiTrajectoryStateMode.h"
#include "TrackingTools/GsfTools/interface/MultiTrajectoryStateTransform.h"
#include "DataFormats/EgammaReco/interface/ElectronSeed.h"
#include "DataFormats/ParticleFlowReco/interface/PFCluster.h"
#include "DataFormats/ParticleFlowReco/interface/PFClusterFwd.h"

// the new map of gsftracks to brem cconv tracks stored in the event (AA)
#include "DataFormats/EgammaTrackReco/interface/GsfTrackToBremConvTracksMapFwd.h"

class PFTrackTransformer;
class GsfTrack;
class MagneticField;
class TrackerGeometry;
class CaloGeometry;
// not needed when we save gsf- conv map (AA)
//class ConvBremPFTrackFinder;

/// \brief Abstract
/*!
\author Michele Pioppi, Daniele Benedetti
\date January 2007

 PFElecTkProducer reads the merged GsfTracks collection
 built with the TrackerDriven and EcalDriven seeds 
 and transform them in PFGsfRecTracks.
*/

class PFElecTkProducer : public edm::EDProducer {
 public:
  
     ///Constructor
     explicit PFElecTkProducer(const edm::ParameterSet&);

     ///Destructor
     ~PFElecTkProducer();

	  void setCaloGeometry( const CaloGeometry* geo) { caloGeometry_ = geo; }

   private:
      virtual void beginRun(const edm::Run&,const edm::EventSetup&) override;
      virtual void endRun(const edm::Run&,const edm::EventSetup&) override;

      ///Produce the PFRecTrack collection
      virtual void produce(edm::Event&, const edm::EventSetup&) override;

    
//      int FindPfRef(const reco::PFRecTrackCollection & PfRTkColl, 
//        reco::GsfTrack, bool);

      /// get the index of the PFRecTrack corresponding to tkRef (AA)  
      int FindPfTkIndex(const reco::PFRecTrackCollection & PfRTkColl, 
        reco::TrackRef & tkRef);
        
        
        
      
//      bool isFifthStep(reco::PFRecTrackRef pfKfTrack);
      bool isFifthStep(reco::TrackRef kfTrack);  // for modified gsfEctra (AA)

      bool applySelection(reco::GsfTrack);
      
      bool resolveGsfTracks(const std::vector<reco::GsfPFRecTrack> &GsfPFVec,
			    unsigned int ngsf,
			    std::vector<unsigned int> &secondaries,
			    const reco::PFClusterCollection & theEClus);

      float minTangDist(const reco::GsfPFRecTrack& primGsf,
			const reco::GsfPFRecTrack& secGsf); 
      
      bool isSameEgSC(const reco::ElectronSeedRef& nSeedRef,
		      const reco::ElectronSeedRef& iSeedRef,
		      bool& bothGsfEcalDriven,
		      float& SCEnergy);

      bool isSharingEcalEnergyWithEgSC(const reco::GsfPFRecTrack& nGsfPFRecTrack,
				       const reco::GsfPFRecTrack& iGsfPFRecTrack,
				       const reco::ElectronSeedRef& nSeedRef,
				       const reco::ElectronSeedRef& iSeedRef,
				       const reco::PFClusterCollection& theEClus,
				       bool& bothGsfTrackerDriven,
				       bool& nEcalDriven,
				       bool& iEcalDriven,
				       float& nEnergy,
				       float& iEnergy);
      
      bool isInnerMost(const reco::GsfTrackRef& nGsfTrack,
		       const reco::GsfTrackRef& iGsfTrack,
		       bool& sameLayer);
      
      bool isInnerMostWithLostHits(const reco::GsfTrackRef& nGsfTrack,
				   const reco::GsfTrackRef& iGsfTrack,
				   bool& sameLayer);
      
      void createGsfPFRecTrackRef(const edm::OrphanHandle<reco::GsfPFRecTrackCollection>& gsfPfHandle,
				  std::vector<reco::GsfPFRecTrack>& gsfPFRecTrackPrimary,
				  const std::map<unsigned int, std::vector<reco::GsfPFRecTrack> >& MapPrimSec);
      
      // helper: convert a vector of TrackRefs to a vector PFRecTrackRefs (AA)
 //     const std::vector<reco::PFRecTrackRef> trkRefVecToPFTrkRefVec(const std::vector<reco::TrackRef>& convTrkRefVector);
        const std::vector<reco::PFRecTrackRef> trkRefVecToPFTrkRefVec(const edm::RefVector<std::vector<reco::Track> >& convTrkRefVector,
                                                                      edm::Handle<reco::PFRecTrackCollection>& thePfRecTrackCollection);
      
          
      // ----------member data ---------------------------
      reco::GsfPFRecTrack pftrack_;
      reco::GsfPFRecTrack secpftrack_;
      edm::ParameterSet conf_;
      edm::InputTag gsfTrackLabel_;
      edm::InputTag pfTrackLabel_;
      edm::InputTag primVtxLabel_;
      edm::InputTag pfEcalClusters_;
      edm::InputTag pfNuclear_;
      edm::InputTag pfConv_;
      edm::InputTag pfV0_;
// for the new map (AA)
      edm::InputTag gsfTrackToBremConvTracksMapTag_;
      
      bool useNuclear_;
      bool useConversions_;
      bool useV0_;
      bool applyAngularGsfClean_;
      double detaCutGsfClean_;
      double dphiCutGsfClean_;

      ///PFTrackTransformer
      PFTrackTransformer *pfTransformer_;     
      const MultiTrajectoryStateMode *mtsMode_;
      MultiTrajectoryStateTransform  mtsTransform_;
// not needed when we save gsf- conv map (AA)
//      ConvBremPFTrackFinder *convBremFinder_;


      ///Trajectory of GSfTracks in the event?
      bool trajinev_;
      bool modemomentum_;
      bool applySel_;
      bool applyGsfClean_;
      bool useFifthStep_;
      bool useFifthStepForEcalDriven_;
      bool useFifthStepForTrackDriven_;
      //   bool useFifthStepSec_;
      bool debugGsfCleaning_;
      double SCEne_;
      double detaGsfSC_;
      double dphiGsfSC_;
      double maxPtConvReco_; // does not seem to be used anywhere (AA)

// not needed when we save gsf- conv map (AA)
/*
      /// Conv Brem Finder
      bool useConvBremFinder_;
      double mvaConvBremFinderID_;
      std::string path_mvaWeightFileConvBrem_;
*/
		// needed for linking by DetID when running on AOD files (AA)
		const CaloGeometry* caloGeometry_;

		// switch for AOD inputs (AA)
		bool runOnAOD_;

};
#endif
