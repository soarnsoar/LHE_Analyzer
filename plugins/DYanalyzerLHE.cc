// -*- C++ -*-
//
// Package:    Analyzer/DYanalyzerLHE
// Class:      DYanalyzerLHE
// 
/**\class DYanalyzerLHE DYanalyzerLHE.cc Analyzer/DYanalyzerLHE/plugins/DYanalyzerLHE.cc

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

class DYanalyzerLHE : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit DYanalyzerLHE(const edm::ParameterSet&);
      ~DYanalyzerLHE();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;
  TH1F * h_dimuon_pt, * h_dimuon_eta, * h_dimuon_phi, *h_dimuon_mass;
  TH1F * h_muonp_pt, * h_muonp_eta, * h_muonp_phi, *h_muonp_mass;
  TH1F * h_muonm_pt, * h_muonm_eta, * h_muonm_phi, *h_muonm_mass;

  TH1F * h_dielectron_pt, * h_dielectron_eta, * h_dielectron_phi, *h_dielectron_mass;
  TH1F * h_electronp_pt, * h_electronp_eta, * h_electronp_phi, *h_electronp_mass;
  TH1F * h_electronm_pt, * h_electronm_eta, * h_electronm_phi, *h_electronm_mass;

  //GenEventInfoProduct                   "generator"                 ""                "SIM"   

  edm::EDGetTokenT<GenParticleCollection> genParticles_Token;
  edm::EDGetTokenT<GenEventInfoProduct> genInfo_Token;
  edm::EDGetTokenT<LHEEventProduct> LHEInfo_Token;





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
DYanalyzerLHE::DYanalyzerLHE(const edm::ParameterSet& iConfig)

{
   //now do what ever initialization is needed
 
  //vector<reco::GenParticle>             "genParticles"              ""                "SIM"     

  usesResource("TFileService");
   genParticles_Token = consumes<GenParticleCollection>(edm::InputTag("genParticles"));
   genInfo_Token = consumes<GenEventInfoProduct>(edm::InputTag("generator"));
   LHEInfo_Token = consumes<LHEEventProduct>(edm::InputTag("externalLHEProducer"));

}


DYanalyzerLHE::~DYanalyzerLHE()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
DYanalyzerLHE::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
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
   


   int i_muonp=-1, i_muonm=-1;
   int i_electronp=-1, i_electronm=-1;


   const lhef::HEPEUP& lheEvent = LHEInfo->hepeup();
   std::vector<lhef::HEPEUP::FiveVector> lheParticles = lheEvent.PUP;
   Int_t nLHEParticle = lheParticles.size();
   for( Int_t idxParticle = 0; idxParticle < nLHEParticle; ++idxParticle ){

     Int_t id = lheEvent.IDUP[idxParticle];
     Int_t status = lheEvent.ISTUP[idxParticle];
     //double px = lheParticles[idxParticle][0];
     //double py = lheParticles[idxParticle][1];
     //double pz = lheParticles[idxParticle][2];
     //double ee = lheParticles[idxParticle][3];
     //double mm = lheParticles[idxParticle][4];
     //cout<<"idxParticle="<<idxParticle<<" id="<<id<<" status="<<status<<" px="<<px<<" py="<<py<<" pz="<<pz<<" ee="<<ee<<" mm="<<mm<<endl;
     if(status==1){
       if(id==13){
         i_muonp=idxParticle;
       }
       else if(id==-13){
         i_muonm=idxParticle;
       }
       else if(id==11){
         i_electronp=idxParticle;
       }
       else if(id==-11){
         i_electronm=idxParticle;
       }

     }
   }


   TLorentzVector vp,vm,vdilep;

   if(i_muonp!=-1 && i_muonm!=-1){ //if dimuon events                                                                                                                             


     vp.SetPxPyPzE(lheParticles[i_muonp][0], lheParticles[i_muonp][1],lheParticles[i_muonp][2],lheParticles[i_muonp][3]  );
     vm.SetPxPyPzE(lheParticles[i_muonm][0], lheParticles[i_muonm][1],lheParticles[i_muonm][2],lheParticles[i_muonm][3]  );

     vdilep=vp+vm;


     h_dimuon_pt->Fill(  vdilep.Perp(),w0 );
     h_dimuon_eta->Fill(  vdilep.Eta() ,w0);
     h_dimuon_phi->Fill(  vdilep.Phi() ,w0);
     h_dimuon_mass->Fill(  vdilep.M() ,w0);

     h_muonp_pt->Fill(vp.Perp(),w0);
     h_muonp_eta->Fill(vp.Eta(),w0);
     h_muonp_phi->Fill(vp.Phi(),w0);

     h_muonm_pt->Fill(vm.Perp(),w0);
     h_muonm_eta->Fill(vm.Eta(),w0);
     h_muonm_phi->Fill(vm.Phi(),w0);

   }
   else if(i_electronp!=-1 && i_electronm != -1){ //if dielectron events                                                                                                          

     vp.SetPxPyPzE(lheParticles[i_electronp][0], lheParticles[i_electronp][1],lheParticles[i_electronp][2],lheParticles[i_electronp][3]  );
     vm.SetPxPyPzE(lheParticles[i_electronm][0], lheParticles[i_electronm][1],lheParticles[i_electronm][2],lheParticles[i_electronm][3]  );


     vdilep=vp+vm;


     h_dielectron_pt->Fill(  vdilep.Perp() ,w0);
     h_dielectron_eta->Fill(  vdilep.Eta() ,w0);
     h_dielectron_phi->Fill( vdilep.Phi() ,w0);
     h_dielectron_mass->Fill( vdilep.M() ,w0);

     h_electronp_pt->Fill(vp.Perp(),w0);
     h_electronp_eta->Fill(vp.Eta(),w0);
     h_electronp_phi->Fill(vp.Phi(),w0);

     h_electronm_pt->Fill(vm.Perp(),w0);
     h_electronm_eta->Fill(vm.Eta(),w0);
     h_electronm_phi->Fill(vm.Phi(),w0);
   }

   


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
DYanalyzerLHE::beginJob()
{
  cout<<"begin job"<<endl;

  edm::Service<TFileService> fs;
  h_dimuon_pt=fs->make<TH1F>("dimuon_pt","pT(#mu#mu)",50,0,100);
  h_dimuon_eta=fs->make<TH1F>("dimuon_eta","#eta(#mu#mu)",50,-5,5);
  h_dimuon_phi=fs->make<TH1F>("dimuon_phi","#phi(#mu#mu)",50,-4,4);
  h_dimuon_mass=fs->make<TH1F>("dimuon_mass","M(#mu#mu)",50,40,140);

  h_muonp_pt=fs->make<TH1F>("muonp_pt","pT(#mu+)",50,0,100);
  h_muonp_eta=fs->make<TH1F>("muonp_eta","#eta(#mu+)",50,-5,5);
  h_muonp_phi=fs->make<TH1F>("muonp_phi","#phi(#mu+)",50,-4,4);


  h_muonm_pt=fs->make<TH1F>("muonm_pt","pT(#mu-)",50,0,100);
  h_muonm_eta=fs->make<TH1F>("muonm_eta","#eta(#mu-)",50,-5,5);
  h_muonm_phi=fs->make<TH1F>("muonm_phi","#phi(#mu-)",50,-4,4);


  h_dielectron_pt=fs->make<TH1F>("dielectron_pt","pT(ee)",50,0,100);
  h_dielectron_eta=fs->make<TH1F>("dielectron_eta","#eta(ee)",50,-5,5);
  h_dielectron_phi=fs->make<TH1F>("dielectron_phi","#phi(ee)",50,-4,4);
  h_dielectron_mass=fs->make<TH1F>("dielectron_mass","M(ee)",50,40,140);


  h_electronp_pt=fs->make<TH1F>("electronp_pt","pT(e+)",50,0,100);
  h_electronp_eta=fs->make<TH1F>("electronp_eta","#eta(e+)",50,-5,5);
  h_electronp_phi=fs->make<TH1F>("electronp_phi","#phi(e+)",50,-4,4);


  h_electronm_pt=fs->make<TH1F>("electronm_pt","pT(e-)",50,0,100);
  h_electronm_eta=fs->make<TH1F>("electronm_eta","#eta(e-)",50,-5,5);
  h_electronm_phi=fs->make<TH1F>("electronm_phi","#phi(e-)",50,-4,4);

  cout<<"end of beginjob"<<endl;


}

// ------------ method called once each job just after ending the event loop  ------------
void 
DYanalyzerLHE::endJob() 
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
DYanalyzerLHE::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(DYanalyzerLHE);
