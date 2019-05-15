// -*- C++ -*-
//
// Package:    Analyzer/JHanalyzer
// Class:      JHanalyzer
// 
/**\class JHanalyzer JHanalyzer.cc Analyzer/JHanalyzer/plugins/JHanalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  JunHo Choi
//         Created:  Fri, 23 Mar 2018 03:24:18 GMT
//
//

//vector<reco::GenParticle>             "prunedGenParticles"        ""                "PAT"     

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"



using namespace edm;
using namespace reco;
using namespace std;

#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenRunInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHERunInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Run.h"//to use edm::Run                                                                                 


#include <TTree.h>
#include <TFile.h>
#include <TLorentzVector.h>


//
// Class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<> and also remove the line from
// constructor "usesResource("TFileService");"
// This will improve performance in multithreaded jobs.

class JHanalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit JHanalyzer(const edm::ParameterSet&);
      ~JHanalyzer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

  //GenEventInfoProduct                   "generator"                 ""                "SIM"   

  edm::EDGetTokenT<GenParticleCollection> genParticles_Token;
  edm::EDGetTokenT<GenEventInfoProduct> genInfo_Token;
  edm::EDGetTokenT<LHEEventProduct> LHEInfo_Token;




  TTree *tree1;
  double Z_pt;
      // ----------member data ---------------------------
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
JHanalyzer::JHanalyzer(const edm::ParameterSet& iConfig)

{
   //now do what ever initialization is needed
 
  //vector<reco::GenParticle>             "genParticles"              ""                "SIM"     

  usesResource("TFileService");
   genParticles_Token = consumes<GenParticleCollection>(edm::InputTag("genParticles"));
   genInfo_Token = consumes<GenEventInfoProduct>(edm::InputTag("generator"));
   LHEInfo_Token = consumes<LHEEventProduct>(edm::InputTag("externalLHEProducer"));

}


JHanalyzer::~JHanalyzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
JHanalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  ////////////initialize/////////////


   using namespace edm;
   edm::Handle<LHEEventProduct> LHEInfo;
   iEvent.getByToken(LHEInfo_Token, LHEInfo);


   //---LHE info and reweight---//
   int lheinfoweightsize= LHEInfo->weights().size();
   int lheinfocommentssize = LHEInfo->comments_size();
   double w0=LHEInfo->originalXWGTUP();
   for (int i_lhe =0; i_lhe < lheinfoweightsize; i_lhe++){
     //weight id in lhe file//
     //cout<<"weight_id="<<LHEInfo->weights()[i_lhe].id<<endl;

     //event weight//
     //cout<<LHEInfo->weights()[i_lhe].wgt/w0<<endl;                                                                                      
     
     
   }

   for (int i =0; i < lheinfocommentssize; i++){                                                                                          
     //cout<<"comment i ="<<i<<"=" << LHEInfo->getComment(i)<<endl;                                                                       
   }                                                                                                                                        
   

   const lhef::HEPEUP& lheEvent = LHEInfo->hepeup();
   std::vector<lhef::HEPEUP::FiveVector> lheParticles = lheEvent.PUP;
   Int_t nLHEParticle = lheParticles.size();
   for( Int_t idxParticle = 0; idxParticle < nLHEParticle; ++idxParticle ){

     Int_t id = lheEvent.IDUP[idxParticle];
     Int_t status = lheEvent.ISTUP[idxParticle];
     double px = lheParticles[idxParticle][0];
     double py = lheParticles[idxParticle][1];
     double pz = lheParticles[idxParticle][2];
     double ee = lheParticles[idxParticle][3];
     double mm = lheParticles[idxParticle][4];
     cout<<"idxParticle="<<idxParticle<<" id="<<id<<" status="<<status<<" px="<<px<<" py="<<py<<" pz="<<pz<<" ee="<<ee<<" mm="<<mm<<endl;
     
   }




   tree1->Fill();


#ifdef THIS_IS_AN_EVENT_EXAMPLE
   //Handle<ExampleData> pIn;
   // iEvent.getByLabel("example",pIn);
#endif
   
#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
   // ESHandle<SetupData> pSetup;
   //iSetup.get<SetupRecord>().get(pSetup);
#endif
}



/*
 fsr_t->Fill();
   fsr_p->Fill();
   fsr_pl->Fill();
   photon_t->Fill();
   isr_t->Fill();

   fsr01_1->Fill();
   fsr01_2->Fill();

   fsr02_1->Fill();
   fsr02_2->Fill();

   fsr03_1->Fill();
   fsr03_2->Fill();

   fsr04_1->Fill();
   fsr04_2->Fill();
 */

// ------------ method called once each job just before starting event loop  ------------
void 
JHanalyzer::beginJob()
{
  cout<<"begin job"<<endl;

  edm::Service<TFileService> fs;
  tree1=fs->make<TTree>("tree1","tree1");
  tree1->Branch("Z_pt",&Z_pt,"Z_pt/D");
  cout<<"end of beginjob"<<endl;


}

// ------------ method called once each job just after ending the event loop  ------------
void 
JHanalyzer::endJob() 
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
JHanalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(JHanalyzer);
