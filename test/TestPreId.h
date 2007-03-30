// -*- C++ -*-
//
// Package:    TestPreId
// Class:      TestPreId
// 
/**\class TestPreId TestPreId.cc fg/TestPreId/src/TestPreId.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Michele Pioppi
//         Created:  Wed Mar 28 19:40:20 CEST 2007
// $Id$
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include <TROOT.h>
#include <TFile.h>
#include <TH1F.h>


//
// class decleration
//

class TestPreId : public edm::EDAnalyzer {
 public:
  explicit TestPreId(const edm::ParameterSet&);
  ~TestPreId();
  
  
 private:
  virtual void beginJob(const edm::EventSetup&) ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
  
  // ----------member data ---------------------------
 private:
  edm::ParameterSet conf_;
  edm::InputTag refitLabel_;
  edm::InputTag mcLabel_;
  edm::InputTag ttLabel_;
  edm::InputTag pfTrackLabel_;

  std::string outputfile_;
  TFile *tf1;
  TH1F  *hpt_pion, *hpt_elec , *heta_pion, *heta_elec;
  TH1F  *hpt_pion_sel, *hpt_elec_sel , *heta_pion_sel, *heta_elec_sel;
};
