#include "RecoParticleFlow/PFTracking/interface/ConvBremTrackFinder.h"
#include "FWCore/Framework/interface/ESHandle.h"


#include "RecoParticleFlow/PFProducer/interface/Utils.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "DataFormats/EgammaReco/interface/ElectronSeed.h"
#include "DataFormats/EgammaReco/interface/ElectronSeedFwd.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h" 
#include "DataFormats/GsfTrackReco/interface/GsfTrackFwd.h" 
#include "TrackingTools/IPTools/interface/IPTools.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "RecoParticleFlow/PFClusterTools/interface/PFEnergyCalibration.h"
#include "TMath.h"
#include "RecoParticleFlow/PFClusterTools/interface/LinkByRecHit.h"
#include "DataFormats/ParticleFlowReco/interface/PFRecTrack.h"
#include "TrackingTools/PatternTools/interface/Trajectory.h"


using namespace edm;
using namespace std;
using namespace reco;

ConvBremTrackFinder::ConvBremTrackFinder(const TransientTrackBuilder& builder,
					 double mvaBremConvCut,
					 string mvaWeightFileConvBrem,
					 const PFTrackTransformer& pfTkTransformer,
					 const MultiTrajectoryStateTransform& mtjstate,
					 const TrackBase::TrackQuality& trackQuality):
  builder_(builder),
  mvaBremConvCut_(mvaBremConvCut),
  mvaWeightFileConvBrem_(mvaWeightFileConvBrem),
  pfTkTransformer_(pfTkTransformer),
  mtjstate_(mtjstate),
  trackQuality_(trackQuality)
{
  tmvaReader_ = new TMVA::Reader("!Color:Silent");
  tmvaReader_->AddVariable("secR",&secR);
  tmvaReader_->AddVariable("sTIP",&sTIP);
  tmvaReader_->AddVariable("nHITS1",&nHITS1);
  tmvaReader_->AddVariable("secPin",&secPin);
  tmvaReader_->AddVariable("Epout",&Epout);
  tmvaReader_->AddVariable("detaBremKF",&detaBremKF);
  tmvaReader_->AddVariable("ptRatioGsfKF",&ptRatioGsfKF);
  tmvaReader_->BookMVA("BDT",mvaWeightFileConvBrem.c_str());

  pfcalib_ = new PFEnergyCalibration();

}
ConvBremTrackFinder::~ConvBremTrackFinder(){delete tmvaReader_; delete pfcalib_; }

void
ConvBremTrackFinder::runConvBremFinder(const Handle<TrackCollection>& theTrackColl,
				       const Handle<VertexCollection>& primaryVertex,
				       const reco::PFClusterCollection & theEClus,
				       reco::GsfTrack gsftrack)
{
  

  found_ = false;
  bool debug = false;
  trackRef_vec_.clear();
  
  if(debug)
    cout << "runConvBremFinder:: Entering " << endl;
  


  // Access Track collections 
  const TrackCollection& TkColl = *(theTrackColl.product());
  reco::TrackCollection::const_iterator tkit = TkColl.begin();
  reco::TrackCollection::const_iterator tkitend= TkColl.end();



  // Access Gsf Track inner/outer momentums and tangents 
  const MultiTrajectoryStateMode *mtsMode_(0);
  TrajectoryStateOnSurface inTSOS = mtjstate_.innerStateOnSurface(gsftrack);
  TrajectoryStateOnSurface outTSOS = mtjstate_.outerStateOnSurface(gsftrack);

  GlobalVector InMom;
  if(inTSOS.isValid()) {
    mtsMode_->momentumFromModeCartesian(inTSOS,InMom);
  }
  else {
    InMom = GlobalVector(gsftrack.pxMode(),gsftrack.pyMode(),gsftrack.pzMode());
  }
  
  GlobalVector OutMom;
  bool isOutFailed = false;
  if(outTSOS.isValid()) {
    mtsMode_->momentumFromModeCartesian(outTSOS,OutMom);
  }
  else {
    isOutFailed = true;
  }
    


   // Compute the gsf radius of the inner hit
  float gsfR = sqrt(gsftrack.innerPosition().x()*gsftrack.innerPosition().x() + 
		    gsftrack.innerPosition().y()*gsftrack.innerPosition().y() );	

  // Access the gsf track tangents
  vector<GsfTangent> gsftang = gsftrack.gsfExtra()->tangents();
  ElectronSeedRef eleSeedRef= gsftrack.extra()->seedRef().castTo<ElectronSeedRef>();


  // Pass the track quality 
  // reco::TrackBase::TrackQuality trackQuality = reco::TrackBase::qualityByName("highPurity");


  // loop over the KF tracks 
  unsigned int i = 0;
  for(;tkit!=tkitend;++tkit,i++){
    
    if(!((TkColl[i]).quality(trackQuality_))) 
      continue;


    reco::TrackRef trackRef(theTrackColl, i);


    // NOTA BENE: add the GSF track ref filter implemented by Anton to identify  
    // the track from the ECAL driven seeded electrons. 

    if (eleSeedRef->ctfTrack().isNonnull()){
      if(trackRef == eleSeedRef->ctfTrack()) {
	if (debug) 
	  cout << " ** SAME REF SEED ** " << endl;
	continue;
      }
    }


    double dphi_presel = fabs(tkit->phi()-gsftrack.phi()); 
    if (dphi_presel > TMath::Pi()) dphi_presel -= TMath::TwoPi();
    double deta_presel = fabs(tkit->eta()-gsftrack.eta());
   
    // run a very loose preselection

    if( fabs(dphi_presel)> 1.0  || fabs(deta_presel) > 0.4) continue;
  
    double minDEtaBremKF = 1000.;
    double minDPhiBremKF = 1000.;
    double minDRBremKF = 1000.;
    double minDEtaBremKFPos = 1000.;
    double minDPhiBremKFPos = 1000.;
    double minDRBremKFPos = 1000.;
    
    double secEta = tkit->innerMomentum().eta();
    double secPhi = tkit->innerMomentum().phi();
    double secEtaPos = tkit->innerPosition().eta();
    double secPhiPos = tkit->innerPosition().phi();

    // The first Brem tangent is taken from the direction of the innermomentum
    double bremEta = InMom.eta();
    double bremPhi = InMom.phi();
    double deta = fabs(bremEta - secEta);
    double dphi = fabs(bremPhi - secPhi);
    if (dphi>TMath::Pi()) dphi-= TMath::TwoPi();
    double DR = sqrt(deta*deta + dphi*dphi);
    
    
    double detaPos = fabs(bremEta - secEtaPos);
    double dphiPos = fabs(bremPhi - secPhiPos);
    if (dphiPos>TMath::Pi()) dphiPos-= TMath::TwoPi();
    double DRPos = sqrt(detaPos*detaPos + dphiPos*dphiPos);
    
    if(DR < minDRBremKF) {
      
      minDRBremKF = DR;
      minDEtaBremKF = deta;
      minDPhiBremKF = fabs(dphi);
    }
    
    if(DRPos < minDRBremKFPos) {
      minDRBremKFPos = DR;
      minDEtaBremKFPos = detaPos;
      minDPhiBremKFPos = fabs(dphiPos);
    }

    for(unsigned int iTang = 0; iTang < gsftrack.gsfExtra()->tangentsSize(); iTang++) {
      
      if ((sqrt(gsftang[iTang].position().x()*gsftang[iTang].position().x() 
		+ gsftang[iTang].position().y()*gsftang[iTang].position().y())>110) 
	  ||(fabs(gsftang[iTang].position().z())>280)) continue;    
    
      bremEta = gsftang[iTang].momentum().eta();
      bremPhi = gsftang[iTang].momentum().phi();
      deta = fabs(bremEta - secEta);
      dphi = fabs(bremPhi - secPhi);
      if (dphi>TMath::Pi()) dphi-= TMath::TwoPi();
      DR = sqrt(deta*deta + dphi*dphi);
      
      
      detaPos = fabs(bremEta - secEtaPos);
      dphiPos = fabs(bremPhi - secPhiPos);
      if (dphiPos>TMath::Pi()) dphiPos-= TMath::TwoPi();
      DRPos = sqrt(detaPos*detaPos + dphiPos*dphiPos);
      
      if(DR < minDRBremKF) {
	
	minDRBremKF = DR;
	minDEtaBremKF = deta;
	minDPhiBremKF = fabs(dphi);
      }
      
      if(DRPos < minDRBremKFPos) {
	minDRBremKFPos = DR;
	minDEtaBremKFPos = detaPos;
	minDPhiBremKFPos = fabs(dphiPos);
      }
    }  // end tangent loop
    // The last Brem tangent is taken from the direction of the outermomentum
    if(isOutFailed == false) {
      bremEta = OutMom.eta();
      bremPhi = OutMom.phi();
      deta = fabs(bremEta - secEta);
      dphi = fabs(bremPhi - secPhi);
      if (dphi>TMath::Pi()) dphi-= TMath::TwoPi();
      DR = sqrt(deta*deta + dphi*dphi);
      
      
      detaPos = fabs(bremEta - secEtaPos);
      dphiPos = fabs(bremPhi - secPhiPos);
      if (dphiPos>TMath::Pi()) dphiPos-= TMath::TwoPi();
      DRPos = sqrt(detaPos*detaPos + dphiPos*dphiPos);
      
      if(DR < minDRBremKF) {
	
	minDRBremKF = DR;
	minDEtaBremKF = deta;
	minDPhiBremKF = fabs(dphi);
      }
      
      if(DRPos < minDRBremKFPos) {
	minDRBremKFPos = DR;
	minDEtaBremKFPos = detaPos;
	minDPhiBremKFPos = fabs(dphiPos);
      }
    } // end of last brem tangent
    
    // secR
    secR = sqrt(tkit->innerPosition().x()*tkit->innerPosition().x() + 
		tkit->innerPosition().y()*tkit->innerPosition().y() );   
    

    // apply loose selection (to be parallel) between the secondary track and brem-tangents.
    // Moreover if the secR is internal with respect to the GSF track by two pixel layers discard it.
    if( (minDPhiBremKF < 0.1 || minDPhiBremKFPos < 0.1) &&
	(minDEtaBremKF < 0.02 ||  minDEtaBremKFPos < 0.02)&&
	secR > (gsfR-8)) {

      if(debug)
	cout << "runConvBremFinder:: OK Find track and BREM close " 
	     << " MinDphi " << minDPhiBremKF << " MinDeta " << minDEtaBremKF  << endl;
     
      reco::PFRecTrack pfrectrack( trackRef->charge(), 
				   reco::PFRecTrack::KF, 
				   i, trackRef );
      bool valid = false;
      Trajectory FakeTraj;
      valid = pfTkTransformer_.addPoints( pfrectrack, *trackRef, FakeTraj);
      if(valid) {
	

	float MinDist = 100000.;
	float EE_calib = 0.; 
	pfrectrack.calculatePositionREP();
	for (PFClusterCollection::const_iterator clus = theEClus.begin();
	     clus != theEClus.end();
	     clus++ ) {
	  // Removed unusd variable, left this in case it has side effects
	  clus->position();
	  double dist = -1.;
	  PFCluster clust = *clus;
	  clust.calculatePositionREP();
	  dist =  pfrectrack.extrapolatedPoint( reco::PFTrajectoryPoint::ECALShowerMax ).isValid() ?
	    LinkByRecHit::testTrackAndClusterByRecHit(pfrectrack , clust ) : -1.;
	  
	  if(dist > 0.) {
	    bool applyCrackCorrections = false;
	    vector<double> ps1Ene(0);
	    vector<double> ps2Ene(0);
	    double ps1,ps2;
	    ps1=ps2=0.;
	    if(dist < MinDist) {
	      MinDist = dist;
	      EE_calib = pfcalib_->energyEm(*clus,ps1Ene,ps2Ene,ps1,ps2,applyCrackCorrections);
	    }
	  }
	}
	if(MinDist > 0. && MinDist < 100000.) {
	  
	  // compute all the input variables for conv brem selection
	  
	  secPout = sqrt(tkit->outerMomentum().x()*tkit->outerMomentum().x() +
			 tkit->outerMomentum().y()*tkit->outerMomentum().y() +
			 tkit->outerMomentum().z()*tkit->outerMomentum().z());
	  
	  secPin = sqrt(tkit->innerMomentum().x()*tkit->innerMomentum().x() +
			tkit->innerMomentum().y()*tkit->innerMomentum().y() +
			tkit->innerMomentum().z()*tkit->innerMomentum().z());
	  
	  
	  // maybe put innter momentum pt? 
	  ptRatioGsfKF = tkit->pt()/(gsftrack.ptMode());
	  
	  Vertex dummy;
	  const Vertex *pv = &dummy;
	  edm::Ref<VertexCollection> pvRef;
	  if (primaryVertex->size() != 0) {
	    pv = &*primaryVertex->begin();
	    // we always use the first vertex (at the moment)
	    pvRef = edm::Ref<VertexCollection>(primaryVertex, 0);
	  } else { // create a dummy PV
	    Vertex::Error e;
	    e(0, 0) = 0.0015 * 0.0015;
	    e(1, 1) = 0.0015 * 0.0015;
	    e(2, 2) = 15. * 15.;
	    Vertex::Point p(0, 0, 0);
	    dummy = Vertex(p, e, 0, 0, 0);
	  }
	  
	  
	  // direction of the Gsf track
	  GlobalVector direction(gsftrack.innerMomentum().x(), 
				 gsftrack.innerMomentum().y(), 
				 gsftrack.innerMomentum().z());
	  
	  TransientTrack transientTrack = builder_.build(*trackRef);	 
	  sTIP = IPTools::signedTransverseImpactParameter(transientTrack, direction, *pv).second.significance();
	  
	  
	  Epout = EE_calib/secPout;
	  
	  // eta distance brem-secondary kf track
	  detaBremKF = minDEtaBremKF;
	  
	  // Number of commont hits
	  trackingRecHit_iterator  nhit=gsftrack.recHitsBegin();
	  trackingRecHit_iterator  nhit_end=gsftrack.recHitsEnd();
	  unsigned int tmp_sh = 0;
	  //uint ish=0;
	  
	  for (;nhit!=nhit_end;++nhit){
	    if ((*nhit)->isValid()){
	      trackingRecHit_iterator  ihit=tkit->recHitsBegin();
	      trackingRecHit_iterator  ihit_end=tkit->recHitsEnd();
	      for (;ihit!=ihit_end;++ihit){
		if ((*ihit)->isValid()) {
		  // method 1
		  if((*nhit)->sharesInput(&*(*ihit),TrackingRecHit::all))  tmp_sh++;
		  
		  // method 2 to switch in case of problem with rechit collections
		  //  if(((*ihit)->geographicalId()==(*nhit)->geographicalId())&&
		  //  (((*nhit)->localPosition()-(*ihit)->localPosition()).mag()<0.01)) ish++;
		  
		  
		}
	      }
	    }
	  }
	  
	  nHITS1 = tmp_sh;
	  
	  double mvaValue = tmvaReader_->EvaluateMVA("BDT");
	  
	  if(debug) 
	    cout << " The imput variables for conv brem tracks identification method 2" << endl
		 << " secR          " << secR << " gsfR " << gsfR  << endl
		 << " N shared hits " << nHITS1 << endl
		 << " sTIP          " << sTIP << endl
		 << " detaBremKF    " << detaBremKF << endl
		 << " E/pout        " << Epout << endl
		 << " pin           " << secPin << endl
		 << " ptRatioKFGsf  " << ptRatioGsfKF << endl
		 << " ***** MVA ***** " << mvaValue << endl;
	  
	  if(mvaValue > mvaBremConvCut_) {
	    found_ = true;
	    trackRef_vec_.push_back(trackRef);
	  }
	  
	}
      } // end if valid
    }
  }
    
}


