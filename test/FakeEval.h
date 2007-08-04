#include <memory>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include <TROOT.h>
#include <TFile.h>
#include <TH1F.h>
class FakeEval : public edm::EDAnalyzer {
   public:
      explicit FakeEval(const edm::ParameterSet&);
      ~FakeEval();


   private:
  virtual void beginJob(const edm::EventSetup&) ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
  int IdTranslator(int);
 private:
  edm::ParameterSet conf_;
  TFile *tf1;
  TH1F *hpt_sec,*hid_sec;
  TH1F *hpt_gen,*hid_gen;
  TH1F *hpt_res,*hpt_gen_res,*hpt_orig_res,*hpt_orig_gen_res;
  TH1F *hpt_sim,*hpt_sim_red,*hpt_rec1,*hpt_rec1_red,*hpt_rec2,*hpt_rec2_red,
    *hpt_rec3,*hpt_rec3_red,*hpt_rec4,*hpt_rec4_red;
  TH1F *hpt_gen_res1,*hpt_gen_res2,*hpt_gen_res3,*hpt_gen_res4;
  TH1F *good_contr,*fake_contr,*untrack_contr,*sec_contr,*wrong_contr,*double_contr,*tot_contr;
  int tot, tot5;
  int tot1,tot2,tot3,tot4;
  int fak1,fak2,fak3,fak4;
  int sec1,sec2,sec3,sec4;
  int gen1,gen2,gen3,gen4;
  int fak, fak5;
//   int sec, sec5;
//   int gen, gen5;
  int altk;
  int origfak,origfak5;
  int evt,origtot,origtot5;

      // ----------member data ---------------------------
};
