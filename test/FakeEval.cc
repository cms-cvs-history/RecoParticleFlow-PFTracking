#include "RecoParticleFlow/PFTracking/test/FakeEval.h"
#include "SimTracker/TrackAssociation/interface/TrackAssociatorBase.h"
#include "SimTracker/TrackerHitAssociation/interface/TrackerHitAssociator.h"
#include "SimTracker/Records/interface/TrackAssociatorRecord.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "SimDataFormats/HepMCProduct/interface/HepMCProduct.h"





//
// class decleration
//
using namespace edm;
using namespace std;
using namespace reco;


//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
FakeEval::FakeEval(const edm::ParameterSet& iConfig): 
  conf_(iConfig)
{
   //now do what ever initialization is needed

}


FakeEval::~FakeEval()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to for each event  ------------
void
FakeEval::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  evt++;
  bool debug = conf_.getParameter<bool>("DEBUG");
  int Iter = conf_.getParameter<int32_t>("Iter");

  vector< int > newtracks;
  vector< int > oldtracks;
  //COLLECTIONS TO BE ANALYZED
  //INPUT_TAG  
  InputTag simTag            = conf_.getParameter<InputTag>("simTracks");
  InputTag tpTag             = conf_.getParameter<InputTag>("trackingParticles");
  double ptcut  = conf_.getParameter<double>("PT_CUT");
  //HANDLE
  Handle<SimTrackContainer> simTracks;
  Handle<TrackingParticleCollection>  tpCollection ;

  //GET COLLECTIONS
  iEvent.getByLabel( simTag            , simTracks);
  iEvent.getByLabel( tpTag             , tpCollection);
  //
  const TrackingParticleCollection tPC   = *(tpCollection.product());
  const SimTrackContainer simTC = *(simTracks.product());



  ////ANALISI HIT ////////////////////
 
  
  Handle<PSimHitContainer> PixelBarrelHitsLowTof;
  Handle<PSimHitContainer> PixelBarrelHitsHighTof;
  Handle<PSimHitContainer> PixelEndcapHitsLowTof;
  Handle<PSimHitContainer> PixelEndcapHitsHighTof;
  Handle<PSimHitContainer> TIBHitsLowTof;
  Handle<PSimHitContainer> TIBHitsHighTof;
  Handle<PSimHitContainer> TIDHitsLowTof;
  Handle<PSimHitContainer> TIDHitsHighTof;
  Handle<PSimHitContainer> TOBHitsLowTof;
  Handle<PSimHitContainer> TOBHitsHighTof;
  Handle<PSimHitContainer> TECHitsLowTof;
  Handle<PSimHitContainer> TECHitsHighTof;
  
  iEvent.getByLabel("g4SimHits","TrackerHitsPixelBarrelLowTof", PixelBarrelHitsLowTof);
  iEvent.getByLabel("g4SimHits","TrackerHitsPixelBarrelHighTof", PixelBarrelHitsHighTof);
  iEvent.getByLabel("g4SimHits","TrackerHitsPixelEndcapLowTof", PixelEndcapHitsLowTof);
  iEvent.getByLabel("g4SimHits","TrackerHitsPixelEndcapHighTof", PixelEndcapHitsHighTof);
  iEvent.getByLabel("g4SimHits","TrackerHitsTIBLowTof", TIBHitsLowTof);
  iEvent.getByLabel("g4SimHits","TrackerHitsTIBHighTof", TIBHitsHighTof);
  iEvent.getByLabel("g4SimHits","TrackerHitsTIDLowTof", TIDHitsLowTof);
  iEvent.getByLabel("g4SimHits","TrackerHitsTIDHighTof", TIDHitsHighTof);
  iEvent.getByLabel("g4SimHits","TrackerHitsTOBLowTof", TOBHitsLowTof);
  iEvent.getByLabel("g4SimHits","TrackerHitsTOBHighTof", TOBHitsHighTof);
  iEvent.getByLabel("g4SimHits","TrackerHitsTECLowTof", TECHitsLowTof);
  iEvent.getByLabel("g4SimHits","TrackerHitsTECHighTof", TECHitsHighTof);
  
  vector<PSimHit> theTrackerHits; 
  
  theTrackerHits.insert(theTrackerHits.end(), PixelBarrelHitsLowTof->begin(), PixelBarrelHitsLowTof->end()); 
  theTrackerHits.insert(theTrackerHits.end(), PixelBarrelHitsHighTof->begin(), PixelBarrelHitsHighTof->end());
  theTrackerHits.insert(theTrackerHits.end(), PixelEndcapHitsLowTof->begin(), PixelEndcapHitsLowTof->end()); 
  theTrackerHits.insert(theTrackerHits.end(), PixelEndcapHitsHighTof->begin(), PixelEndcapHitsHighTof->end());
  theTrackerHits.insert(theTrackerHits.end(), TIBHitsLowTof->begin(), TIBHitsLowTof->end()); 
  theTrackerHits.insert(theTrackerHits.end(), TIBHitsHighTof->begin(), TIBHitsHighTof->end());
  theTrackerHits.insert(theTrackerHits.end(), TIDHitsLowTof->begin(), TIDHitsLowTof->end()); 
  theTrackerHits.insert(theTrackerHits.end(), TIDHitsHighTof->begin(), TIDHitsHighTof->end());
  theTrackerHits.insert(theTrackerHits.end(), TOBHitsLowTof->begin(), TOBHitsLowTof->end()); 
  theTrackerHits.insert(theTrackerHits.end(), TOBHitsHighTof->begin(), TOBHitsHighTof->end());
  theTrackerHits.insert(theTrackerHits.end(), TECHitsLowTof->begin(), TECHitsLowTof->end()); 
  theTrackerHits.insert(theTrackerHits.end(), TECHitsHighTof->begin(), TECHitsHighTof->end());

  //SIMTRACK COLLECTION
  float ptsim=0;
  float ptrec=0;
  float ptrec_gen=0;
  //  float ptrec_orig=0;
  //  float ptrec_orig_gen=0;

  uint Minhit = conf_.getParameter<uint32_t>("MIN_HIT");




  ESHandle<TrackAssociatorBase> theHitsAssociator;
  iSetup.get<TrackAssociatorRecord>().get("TrackAssociatorByHits",theHitsAssociator);
  TrackAssociatorBase * associatorByHits = (TrackAssociatorBase *) theHitsAssociator.product();
  
  InputTag tkTagFirst    = conf_.getParameter<InputTag>("FirstTrackCollection");
  Handle<TrackCollection> tkFirstCollection;
  iEvent.getByLabel( tkTagFirst             , tkFirstCollection);

  InputTag tkTagSecond    = conf_.getParameter<InputTag>("SecondTrackCollection");
  Handle<TrackCollection> tkSecondCollection;
  iEvent.getByLabel( tkTagSecond             , tkSecondCollection);

  InputTag tkTagThird    = conf_.getParameter<InputTag>("ThirdTrackCollection");
  Handle<TrackCollection> tkThirdCollection;
  iEvent.getByLabel( tkTagThird             , tkThirdCollection);


  float ptuntrack=0;
  vector < pair<uint,pair<uint,float> > > simsim;
 //  for(SimTrackContainer::const_iterator it=simTC.begin();it!=simTC.end();it++){
  for(SimTrackContainer::size_type i=0; i<simTC.size(); ++i){
    uint isimhit=0;
    //   bool tracked=false;
    //   i++;
    TrackingParticleRef tp (tpCollection,i); 

    for (vector<PSimHit>::iterator isim = theTrackerHits.begin();
	 isim != theTrackerHits.end(); ++isim){
    
      if (isim->trackId()-1==tp.index()) isimhit++;
    }


    if ((!simTC[i].noGenpart())&&
	(simTC[i].charge()!=0)&&
	(simTC[i].momentum().Pt()>ptcut)&&
	(isimhit>Minhit)){ 
  
      ptsim+=simTC[i].momentum().Pt(); 
      hpt_sim->Fill(simTC[i].momentum().Pt());
      hpt_sim_red->Fill(simTC[i].momentum().Pt());
      simsim.push_back(make_pair(tp.index(),make_pair(0,simTC[i].momentum().Pt())));

    }
  }
  //SIM-REC TRACK ASSOCIATION

 

  float ptrec1 =0;
  float ptrec2 =0;
  float ptrec3 =0;
  float ptrec4 =0;
  //CONTRIBUTION TO CHARGED COMPONENT
  float ptgood=0;
  float ptsec=0;
  float ptfake=0;
 
  float ptwrong=0;
  float ptdouble=0;
  vector<uint> alreadyasso;
  /////////////////////////////
  //FIRST COLLECTION
  if (Iter>0){

    const TrackCollection  tCFirst = *(tkFirstCollection.product()); 

    RecoToSimCollection pFirst = 
      associatorByHits->associateRecoToSim (tkFirstCollection,tpCollection,&iEvent );
    
    tot1+=tCFirst.size(); 
    
    for(TrackCollection::size_type i=0; i<tCFirst.size(); ++i) {
     hpt_rec1->Fill(tCFirst[i].pt());
     hpt_rec1_red->Fill(tCFirst[i].pt());
     hpt_rec2->Fill(tCFirst[i].pt());
     hpt_rec2_red->Fill(tCFirst[i].pt());
     hpt_rec3->Fill(tCFirst[i].pt());
     hpt_rec3_red->Fill(tCFirst[i].pt());
     hpt_rec4->Fill(tCFirst[i].pt());
     hpt_rec4_red->Fill(tCFirst[i].pt());
     ptrec+=tCFirst[i].pt();
     ptrec1+=tCFirst[i].pt();
     TrackRef track(tkFirstCollection, i);
      
      
      try{ 
	std::vector<std::pair<TrackingParticleRef, double> > tp = pFirst[track];
	if(debug)    cout << "Reco Track " << setw(2) << track.index() 
			  << " pT: "  << setw(6) << track->pt() 
			  <<  " matched to " << tp.size() 
			  << " MC Tracks" << std::endl;
	if (tp.size()>0) ptrec_gen+=track->pt();
	uint izo=0;
	for (std::vector<std::pair<TrackingParticleRef, double> >::const_iterator it = tp.begin(); 
	     it != tp.end(); ++it) {
	  TrackingParticleRef tpr = it->first;
	  newtracks.push_back(tpr.index());
	  //	cout<<"ff "<<tpr->genParticle().size()<<endl;
	  if (tpr->genParticle().size()>0){
	    gen1++;
	    bool firsttime=true;
	 
	    for (uint iass=0;iass<alreadyasso.size();iass++){
	      //	      if (tpr.index()==alreadyasso[iass]) firsttime=false;
	      if (tpr.index()==alreadyasso[iass]) firsttime=false;
	 
	    }

	    for (uint iass=0;iass<simsim.size();iass++){
	      if (tpr.index()==simsim[iass].first) simsim[iass].second.first++;
	    }
	    double assocChi2 = it->second;
	    if(debug)  cout << "\t\tMCTrack " << setw(2) << tpr.index() 
			    << " pT: " << setw(6) << tpr->pt()  
			    << " NShared: " << assocChi2 << endl;
	    float df=tpr->pt()-track->pt();
	    if (firsttime){
	      alreadyasso.push_back(tpr.index());
	      
	      if (izo==0) {

		float res = 0.055+(tpr->pt()*0.017);
		if (fabs(df/tpr->pt())<res) ptgood+=df;
		else ptwrong+=df;
		//	      cout<<"DIFF "<<df<<" "<<res<<endl;

		hpt_gen->Fill(tpr->pt());
		hid_gen->Fill(IdTranslator(tpr->pdgId()));
		//	    if (track->found()>5) gen5++;
	    



	      }
	    }else ptdouble+=track->pt();
	    izo++;
	  }
	  if (tpr->genParticle().size()==0){
	    ptsec+=track->pt();
	    sec1++;
	    hpt_sec->Fill(tpr->pt());
	    hid_sec->Fill(IdTranslator(tpr->pdgId()));
	  }
	}
      } catch (Exception event) {

	ptfake+=track->pt();
	fak1++;
	if(debug) cout << "->   Track " << setw(2) << track.index() << " pT: " 
		       << setprecision(2) << setw(6) << track->pt() 
		       <<  " matched to 0  MC Tracks" << endl;
	
	
      }
    }
  }
  

  //SECOND COLLECTION
  if (Iter>1){
    const TrackCollection  tCSecond = *(tkSecondCollection.product()); 

    RecoToSimCollection pSecond = 
      associatorByHits->associateRecoToSim (tkSecondCollection,tpCollection,&iEvent );
    
    tot2+=tCSecond.size(); 
    
    for(TrackCollection::size_type i=0; i<tCSecond.size(); ++i) {
     hpt_rec2->Fill(tCSecond[i].pt());
     hpt_rec2_red->Fill(tCSecond[i].pt());
     hpt_rec3->Fill(tCSecond[i].pt());
     hpt_rec3_red->Fill(tCSecond[i].pt());
     hpt_rec4->Fill(tCSecond[i].pt());
     hpt_rec4_red->Fill(tCSecond[i].pt());
     ptrec+=tCSecond[i].pt();
     ptrec2+=tCSecond[i].pt();
     TrackRef track(tkSecondCollection, i);
      
      
      try{ 
	std::vector<std::pair<TrackingParticleRef, double> > tp = pSecond[track];
	if(debug)    cout << "Reco Track " << setw(2) << track.index() 
			  << " pT: "  << setw(6) << track->pt() 
			  <<  " matched to " << tp.size() 
			  << " MC Tracks" << std::endl;
	if (tp.size()>0) ptrec_gen+=track->pt();
	uint izo=0;
	for (std::vector<std::pair<TrackingParticleRef, double> >::const_iterator it = tp.begin(); 
	     it != tp.end(); ++it) {
	  TrackingParticleRef tpr = it->first;
	  newtracks.push_back(tpr.index());
	  //	cout<<"ff "<<tpr->genParticle().size()<<endl;
	  if (tpr->genParticle().size()>0){
	    gen2++;
	    bool firsttime=true;
	    for (uint iass=0;iass<alreadyasso.size();iass++){
	      if (tpr.index()==alreadyasso[iass]) firsttime=false;
	    }
	    for (uint iass=0;iass<simsim.size();iass++){
	      if (tpr.index()==simsim[iass].first) simsim[iass].second.first++;
	    }
	    double assocChi2 = it->second;
	    if(debug)  cout << "\t\tMCTrack " << setw(2) << tpr.index() 
			    << " pT: " << setw(6) << tpr->pt()  
			    << " NShared: " << assocChi2 << endl;
	    float df=tpr->pt()-track->pt();
	    if (firsttime){
	      alreadyasso.push_back(tpr.index());
	      if (izo==0) {

		float res = 0.055+(tpr->pt()*0.017);
		if (fabs(df/tpr->pt())<res) ptgood+=df;
		else ptwrong+=df;
		//	      cout<<"DIFF "<<df<<" "<<res<<endl;

		hpt_gen->Fill(tpr->pt());
		hid_gen->Fill(IdTranslator(tpr->pdgId()));
		//	    if (track->found()>5) gen5++;
	    
	      }
	    }else ptdouble+=track->pt();
	    izo++;
	  }
	  if (tpr->genParticle().size()==0){
	    ptsec+=track->pt();
	    sec2++;
	    hpt_sec->Fill(tpr->pt());
	    hid_sec->Fill(IdTranslator(tpr->pdgId()));
	  }

	}
      } catch (Exception event) {
	//     cout<<"hu" <<endl;
	ptfake+=track->pt();
	fak2++;
	if(debug) cout << "->   Track " << setw(2) << track.index() << " pT: " 
		       << setprecision(2) << setw(6) << track->pt() 
		       <<  " matched to 0  MC Tracks" << endl;
	
	
      }
    }
  }

  //THIRD COLLECTION
  if (Iter>2){

    const TrackCollection  tCThird = *(tkThirdCollection.product()); 

    RecoToSimCollection pThird = 
      associatorByHits->associateRecoToSim (tkThirdCollection,tpCollection,&iEvent );
    
    tot3+=tCThird.size(); 
    
    for(TrackCollection::size_type i=0; i<tCThird.size(); ++i) {
     hpt_rec3->Fill(tCThird[i].pt());
     hpt_rec3_red->Fill(tCThird[i].pt());
     hpt_rec4->Fill(tCThird[i].pt());
     hpt_rec4_red->Fill(tCThird[i].pt());
     ptrec+=tCThird[i].pt();
     ptrec3+=tCThird[i].pt();
     TrackRef track(tkThirdCollection, i);
      
      
      try{ 
	std::vector<std::pair<TrackingParticleRef, double> > tp = pThird[track];
	if(debug)    cout << "Reco Track " << setw(2) << track.index() 
			  << " pT: "  << setw(6) << track->pt() 
			  <<  " matched to " << tp.size() 
			  << " MC Tracks" << std::endl;
	if (tp.size()>0) ptrec_gen+=track->pt();
	uint izo=0;
	for (std::vector<std::pair<TrackingParticleRef, double> >::const_iterator it = tp.begin(); 
	     it != tp.end(); ++it) {
	  TrackingParticleRef tpr = it->first;
	  newtracks.push_back(tpr.index());
	  //	cout<<"ff "<<tpr->genParticle().size()<<endl;
	  if (tpr->genParticle().size()>0){
	    gen3++;
	    bool firsttime=true;
	    for (uint iass=0;iass<alreadyasso.size();iass++){
	      if (tpr.index()==alreadyasso[iass]) firsttime=false;
	    }
	    for (uint iass=0;iass<simsim.size();iass++){
	      if (tpr.index()==simsim[iass].first) simsim[iass].second.first++;
	    }
	    double assocChi2 = it->second;
	    if(debug)  cout << "\t\tMCTrack " << setw(2) << tpr.index() 
			    << " pT: " << setw(6) << tpr->pt()  
			    << " NShared: " << assocChi2 << endl;
	    float df=tpr->pt()-track->pt();
	    if (firsttime){
	      alreadyasso.push_back(tpr.index());
	      if (izo==0) {

		float res = 0.055+(tpr->pt()*0.017);
		if (fabs(df/tpr->pt())<res) ptgood+=df;
		else ptwrong+=df;
		//	      cout<<"DIFF "<<df<<" "<<res<<endl;

		hpt_gen->Fill(tpr->pt());
		hid_gen->Fill(IdTranslator(tpr->pdgId()));
		//	    if (track->found()>5) gen5++;
	    



	      }
	    }else ptdouble+=track->pt();
	    izo++;
	  }
	  if (tpr->genParticle().size()==0){
	    ptsec+=track->pt();
	    sec3++;
	    hpt_sec->Fill(tpr->pt());
	    hid_sec->Fill(IdTranslator(tpr->pdgId()));
	  }
	}
      } catch (Exception event) {

	ptfake+=track->pt();
	fak3++;
	if(debug) cout << "->   Track " << setw(2) << track.index() << " pT: " 
		       << setprecision(2) << setw(6) << track->pt() 
		       <<  " matched to 0  MC Tracks" << endl;
	
	
      }
    }
  }

  //FOURTH COLLECTION
  if (Iter>3){
    InputTag tkTagFourth    = conf_.getParameter<InputTag>("FourthTrackCollection");
    Handle<TrackCollection> tkFourthCollection;
    iEvent.getByLabel( tkTagFourth             , tkFourthCollection);
    const TrackCollection  tCFourth = *(tkFourthCollection.product()); 

    RecoToSimCollection pFourth = 
      associatorByHits->associateRecoToSim (tkFourthCollection,tpCollection,&iEvent );
    
    tot4+=tCFourth.size(); 
    
    for(TrackCollection::size_type i=0; i<tCFourth.size(); ++i) {
     hpt_rec4->Fill(tCFourth[i].pt());
     hpt_rec4_red->Fill(tCFourth[i].pt());
     ptrec+=tCFourth[i].pt();
     ptrec4+=tCFourth[i].pt();
     TrackRef track(tkFourthCollection, i);
      
      
      try{ 
	std::vector<std::pair<TrackingParticleRef, double> > tp = pFourth[track];
	if(debug)    cout << "Reco Track " << setw(2) << track.index() 
			  << " pT: "  << setw(6) << track->pt() 
			  <<  " matched to " << tp.size() 
			  << " MC Tracks" << std::endl;
	if (tp.size()>0) ptrec_gen+=track->pt();
	for (std::vector<std::pair<TrackingParticleRef, double> >::const_iterator it = tp.begin(); 
	     it != tp.end(); ++it) {
	  TrackingParticleRef tpr = it->first;
	  newtracks.push_back(tpr.index());

	  if (tpr->genParticle().size()>0){
	    gen4++;
	    
	    hpt_gen->Fill(tpr->pt());
	    hid_gen->Fill(IdTranslator(tpr->pdgId()));

	  }
	  if (tpr->genParticle().size()==0){
	    sec4++;
	    hpt_sec->Fill(tpr->pt());
	    hid_sec->Fill(IdTranslator(tpr->pdgId()));

	  }
	  double assocChi2 = it->second;
	  if(debug)  cout << "\t\tMCTrack " << setw(2) << tpr.index() 
			  << " pT: " << setw(6) << tpr->pt()  
			  << " NShared: " << assocChi2 << endl;
	}
      } catch (Exception event) {
	fak4++;
	if(debug) cout << "->   Track " << setw(2) << track.index() << " pT: " 
		       << setprecision(2) << setw(6) << track->pt() 
		       <<  " matched to 0  MC Tracks" << endl;
	
	
      }
    }
  }



  for (uint isims=0;isims<simsim.size();isims++){
 
    if (simsim[isims].second.first==0) ptuntrack+=simsim[isims].second.second;
  }


  hpt_gen_res1->Fill(ptsim-ptrec1);
  hpt_gen_res2->Fill(ptsim-(ptrec1+ptrec2));
  hpt_gen_res3->Fill(ptsim-(ptrec1+ptrec2+ptrec3));
  hpt_gen_res4->Fill(ptsim-(ptrec1+ptrec2+ptrec3+ptrec4));
  if (debug){
    cout<<"TOTAL  = "<<ptsim-(ptrec1+ptrec2+ptrec3+ptrec4)<<" "<<(ptgood+ptuntrack+ptwrong)-(ptfake+ptsec+ptdouble)<<endl;
    cout<<"GOOD    = "<<ptgood<<endl;
    cout<<"SEC     = "<<ptsec<<endl;
    cout<<"FAKE    = "<<ptfake<<endl;
    cout<<"DOUBLE  = "<<ptdouble<<endl;
    cout<<"WRONG   = "<<ptwrong<<endl;
    cout<<"UNTRACK = "<<ptuntrack<<endl;
  }

  good_contr->Fill(ptgood);
  fake_contr->Fill(ptgood-ptfake);
  untrack_contr->Fill(ptgood+ptuntrack);
  sec_contr->Fill(ptgood-ptsec);
  wrong_contr->Fill(ptgood+ptwrong);
  double_contr->Fill(ptgood-ptdouble);
  tot_contr->Fill((ptgood+ptuntrack+ptwrong)-(ptfake+ptsec+ptdouble));
}


// ------------ method called once each job just before starting event loop  ------------
void 
FakeEval::beginJob(const edm::EventSetup&)
{
  string OutputFile = conf_.getParameter<string>("OutputFile");
  tf1 = new TFile(OutputFile.c_str(), "RECREATE");
  hpt_sec= new TH1F("hpt_sec","hpt_sec",200,0.0,20);
  hid_sec= new TH1F("hid_sec","hid_sec",13,-6.5,6.5);
  hpt_gen= new TH1F("hpt_gen","hpt_gen",200,0.0,20);
  hid_gen= new TH1F("hid_gen","hid_gen",13,-6.5,6.5);
  hpt_res= new TH1F("hpt_res","hpt_res",100,-50,50);

  hpt_gen_res1= new TH1F("hpt_gen_res1","hpt_gen_res1",60,-180.,180.);
  hpt_gen_res2= new TH1F("hpt_gen_res2","hpt_gen_res2",60,-180.,180.);
  hpt_gen_res3= new TH1F("hpt_gen_res3","hpt_gen_res3",60,-180.,180.);
  hpt_gen_res4= new TH1F("hpt_gen_res4","hpt_gen_res4",60,-180.,180.);

  tot_contr  = new TH1F("tot_contr","tot_contr",60,-60.,60.);
  good_contr = new TH1F("good_contr","good_contr",60,-60.,60.);
  fake_contr = new TH1F("fake_contr","fake_contr",60,-60.,60.);
  untrack_contr = new TH1F("untrack_contr","untrack_contr",60,-60.,60.);
  sec_contr = new TH1F("sec_contr","sec_contr",60,-60.,60.);
  wrong_contr = new TH1F("wrong_contr","wrong_contr",60,-60.,60.);
  double_contr = new TH1F("double_contr","double_contr",60,-60.,60.); 

  hpt_sim = new TH1F("hpt_sim","hpt_sim",30,0,30);
  hpt_sim_red = new TH1F("hpt_sim_red","hpt_sim_red",30,0,3);
  hpt_rec1 = new TH1F("hpt_rec1","hpt_rec1",30,0,30);
  hpt_rec1_red = new TH1F("hpt_rec1_red","hpt_rec1_red",30,0,3);
  hpt_rec2 = new TH1F("hpt_rec2","hpt_rec2",30,0,30);
  hpt_rec2_red = new TH1F("hpt_rec2_red","hpt_rec2_red",30,0,3);
  hpt_rec3 = new TH1F("hpt_rec3","hpt_rec3",30,0,30);
  hpt_rec3_red = new TH1F("hpt_rec3_red","hpt_rec3_red",30,0,3);
  hpt_rec4 = new TH1F("hpt_rec4","hpt_rec4",30,0,30);
  hpt_rec4_red = new TH1F("hpt_rec4_red","hpt_rec4_red",30,0,3);
  tot1=0;  tot2=0; tot3=0; tot4=0;
  fak1=0;  fak2=0; fak3=0; fak4=0;

  sec1=0;  sec2=0;   sec3=0;  sec4=0;
  gen1=0;  gen2=0;   gen3=0;  gen4=0;
  altk=0;
  origfak=0;origfak5=0;
  evt=0; origtot=0; origtot5=0;
}

// ------------ method called once each job just after ending the event loop  ------------
void 
FakeEval::endJob() {
  tf1->Write();
  tf1->Close();
  cout<<"EVENTI "<<evt<<endl;
  float ev= float(evt);
  cout<<"RECO "<<tot1/ev<<" "<<tot2/ev<<" "<<tot3/ev<<" "<<tot4/ev <<endl;
  cout<<"FAKE "<<fak1/ev<<" "<<fak2/ev<<" "<<fak3/ev<<" "<<fak4/ev <<endl;
  cout<<"GEN  "<<gen1/ev<<" "<<gen2/ev<<" "<<gen3/ev<<" "<<gen4/ev <<endl;
  cout<<"SEC  "<<sec1/ev<<" "<<sec2/ev<<" "<<sec3/ev<<" "<<sec4/ev <<endl;



//   cout<<endl<<"ALREADY TRACKED IN THE FIRST LEVEL "<<100*float(altk)/float(tot)<<"%"<<endl;
//   cout<<"FAKE "<<endl
//       <<"ORIG "<<origfak<<" more than 5 "<<origfak5<<endl
//       <<"NEW  "<<fak<<" more than 5 "<<fak5<<endl;
    
}
int FakeEval::IdTranslator(int id){
  int newid=0;
  if (id==11) newid=-1;
  if (id==-11) newid=1;
  if (id==13) newid=-2;
  if (id==-13) newid=2;
  if (id==211) newid=3;
  if (id==-211) newid=-3;
  if (id==321) newid=4;
  if (id==-321) newid=-4;
  if (id==2212) newid=5;
  if (id==-2212) newid=-5;
 
  if (id==22) newid=6;
  return newid;
}

