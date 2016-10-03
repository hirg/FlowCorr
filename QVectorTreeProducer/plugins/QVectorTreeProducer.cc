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
   flowTree_->Branch("noff", nOff_, "noff/I");
   flowTree_->Branch("nref", nRef_, "nref/I");
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
   if(!isEventSelected(iEvent, vertices, tracks)) return;
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
QVectorTreeProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) 
{
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

// ------------ Other methods ------------
bool
QVectorTreeProducer::isEventSelected(const edm::Event& iEvent, 
                                     edm::Handle< reco::VertexCollection > vertices, 
                                     edm::Handle< reco::TrackCollection >  tracks) 
{
  // ## Vertex ##
  //--- Vertex info 
  nVtx_ = 0;

  const reco::Vertex & bestVtx = (*vertices)[0];
   
  if(!bestVtx.isFake() && bestVtx.tracksSize()>=2)
  {
     xVtx_ = bestVtx.x();
     yVtx_ = bestVtx.y();
     zVtx_ = bestVtx.z();
     xVtxError_ = bestVtx.xError();
     yVtxError_ = bestVtx.yError();
     zVtxError_ = bestVtx.zError();
  }
  double rho = sqrt(xVtx_*xVtx_ + yVtx_*yVtx_);

  for(reco::VertexCollection::const_iterator itVtx = vertices->begin(); itVtx != vertices->end(); ++itVtx)
  {
      if(!itVtx->isFake() && itVtx->tracksSize() >= 2) nVtx_++;
  }

  //--- Vertex selection 
  if(nVtx_ > nVtxMax_)                 return false;
  if(zVtx_ < MinVz_ || zVtx_ > MaxVz_) return false;
  if(rho > MaxRho_ )                   return false;


  // ## Multiplicity ##
  //--- Multiplicity info (nOff)
  nOff_ = 0; 

  for(reco::TrackCollection::const_iterator itTrk = tracks->begin(); itTrk != tracks->end(); ++itTrk)
  {
     if(itTrk->pt() < 0.0001) continue;

     // Select tracks based on proximity to best vertex
     math::XYZPoint bestVtxPoint(xVtx_,yVtx_,zVtx_);

     double dzvtx   = itTrk->dz(bestVtxPoint);
     double d0vtx   = itTrk->dxy(bestVtxPoint);
     double dzerror = sqrt(itTrk->dzError()*itTrk->dzError()+zVtxError_*zVtxError_);
     double d0error = sqrt(itTrk->d0Error()*itTrk->d0Error()+xVtxError_*yVtxError_);
     double eta     = itTrk->eta();
     double pt      = itTrk->pt();
     double pterror = itTrk->ptError();
     double chi2n   = itTrk->normalizedChi2();
     double nlayers = itTrk->hitPattern().trackerLayersWithMeasurement();
     chi2n          = chi2n/nlayers;
     int nhits      = itTrk->numberOfValidHits();
     //int algo       = itTrk->originalAlgo();
     int charge     = itTrk->charge(); 

  //--- Multiplicity selection
     //if(!itTrk->quality(reco::TrackBase::qualityByName(trackQualityTag_))) continue;
     if(!itTrk->quality(reco::TrackBase::highPurity)) continue;
     if(fabs(dzvtx/dzerror) > 3.0)                    continue;
     if(fabs(d0vtx/d0error) > 3.0)                    continue;
     if(fabs(eta) > 2.4)                              continue;
     if(pt < 0.4)                                     continue;
     if(fabs(pterror)/pt > 0.1)                       continue;
     if(chi2n > 9999.)                                continue;
     if(nhits < 0)                                    continue;
     //if(find(trkAlgo_.begin(), trkAlgo_.end(), algo) == trkAlgo_.end() &&
     //   !trkAlgo_.empty())                            continue;
     if(charge == 0)                                  continue; 

     nOff_++;
  }

  return true;
}


//define this as a plug-in
DEFINE_FWK_MODULE(QVectorTreeProducer);
