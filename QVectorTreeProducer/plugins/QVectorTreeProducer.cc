// -*- C++ -*-
//
// Package:    FlowCorr/QVectorTreeProducer
// Class:      QVectorTreeProducer
// 
/**\class QVectorTreeProducer QVectorTreeProducer.cc FlowCorr/QVectorTreeProducer/plugins/QVectorTreeProducer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Maxime Guilbaud
//         Created:  Fri, 23 Sep 2016 12:48:03 GMT
//
//


// system include files

// user include files
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "FlowCorr/QVectorTreeProducer/interface/QVectorTreeProducer.h"

QVectorTreeProducer::QVectorTreeProducer(const edm::ParameterSet& iConfig) :
   trackToken_(consumes<reco::TrackCollection>(iConfig.getParameter<edm::InputTag>("tracks"))),
   trackQualityTag_(iConfig.getUntrackedParameter<std::string>("trackQuality","highPurity")),
   vtxToken_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertex"))),
   centralityToken_(consumes<reco::Centrality>(iConfig.getParameter<edm::InputTag>("centralitySrc"))),
   centralityBinToken_(consumes<int>(iConfig.getUntrackedParameter<edm::InputTag>("centralityBinSrc"))),
   caloTowerToken_(consumes<CaloTowerCollection>(iConfig.getParameter<edm::InputTag>("calotower")))
{
   //now do what ever initialization is needed
   //usesResource("TFileService");
   edm::Service<TFileService> fs;

   flowTree_ = fs->make<TTree>("flowTree", "flowTree");
   flowTree_->Branch("noff", noff_, "noff/I");
   flowTree_->Branch("nref", nref_, "nref/I");
}


QVectorTreeProducer::~QVectorTreeProducer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
QVectorTreeProducer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   //----- Vertex selection -----
   edm::Handle< reco::VertexCollection > vertices;
   iEvent.getByToken(vtxToken_, vertices);
   if(!vertices->size())
   {
      edm::LogWarning ("Missing Collection") <<"Invalid or empty vertex collection!";
      return;
   }

   //----- Tracks selection -----
   edm::Handle< reco::TrackCollection > tracks;
   
   iEvent.getByToken(trackToken_, tracks);
   if( !tracks->size() )
   {
       edm::LogWarning ("Missing Collection") <<"Invalid or empty track collection!";
       return;
   }
   
   //----- Calotower selection -----
   edm::Handle< CaloTowerCollection > calotowers;
   iEvent.getByToken(caloTowerToken_, calotowers);
   if(!calotowers->size()) 
   { 
       edm::LogWarning ("Missing Collection") <<"Invalid or empty caloTower collection!";
       return; 
   }

   //----- Centrality selection ----
   edm::Handle< int > cbin;
   iEvent.getByToken(centralityBinToken_,cbin);
   int hiBin = *cbin;
   if(hiBin < 0) 
   { 
       edm::LogWarning ("Invalid value") <<"Invalid centrality value";
       return; 
   }
   
   edm::Handle< reco::Centrality > centrality;
   iEvent.getByToken(centralityToken_, centrality);

   
   //----- Event selection -----
   if(!isEventSelected(iEvent)) return;
}


// ------------ method called once each job just before starting event loop  ------------
void 
QVectorTreeProducer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
QVectorTreeProducer::endJob() 
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
QVectorTreeProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

// ------------ Other methods ------------
bool
QVectorTreeProducer::isEventSelected(const edm::Event& iEvent) {

  return true;
}


//define this as a plug-in
DEFINE_FWK_MODULE(QVectorTreeProducer);
