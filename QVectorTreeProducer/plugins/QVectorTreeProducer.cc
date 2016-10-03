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
// #General info and tags
   trackToken_(consumes<reco::TrackCollection>(iConfig.getParameter<edm::InputTag>("tracks"))),
   trackQualityTag_(iConfig.getUntrackedParameter<std::string>("trackQuality","highPurity")),
   vtxToken_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertex"))),
   centralityToken_(consumes<reco::Centrality>(iConfig.getParameter<edm::InputTag>("centralitySrc"))),
   centralityBinToken_(consumes<int>(iConfig.getUntrackedParameter<edm::InputTag>("centralityBinSrc"))),
   caloTowerToken_(consumes<CaloTowerCollection>(iConfig.getParameter<edm::InputTag>("calotower"))),
// #Vertex
   MinVz_(iConfig.getUntrackedParameter<double>("minVz",-15.)),
   MaxVz_(iConfig.getUntrackedParameter<double>("maxVz", 15.)),
   MaxRho_(iConfig.getUntrackedParameter<double>("maxRho", 0.2)),
   nVtxMax_(iConfig.getUntrackedParameter<int>("nVtxMax", 9999)),
// #Offline tracks (noff)
   MinEtaOff_(iConfig.getUntrackedParameter<double>("minEtaOff",-2.4)),
   MaxEtaOff_(iConfig.getUntrackedParameter<double>("maxEtaOff", 2.4)),
   MinPtOff_(iConfig.getUntrackedParameter<double>("minPtOff",0.4)),
   MaxPtOff_(iConfig.getUntrackedParameter<double>("maxPtOff", 9999.)),
   ChargeOff_(iConfig.getUntrackedParameter< std::vector<int> >("chargeOff")),
   isPixelTrackingOff_(iConfig.getUntrackedParameter<bool>("isPixTrkOff", 9999.)),
   dzdzErrorOff_(iConfig.getUntrackedParameter<double>("dzdzErrorOff", 3.)),
   d0d0ErrorOff_(iConfig.getUntrackedParameter<double>("d0d0ErrorOff", 3.)),
   PtErrorPtOff_(iConfig.getUntrackedParameter<double>("PtErrorPtOff", 1.)),
   Chi2nOff_(iConfig.getUntrackedParameter<double>("Chi2nOff", 0.)),
   nHitsOff_(iConfig.getUntrackedParameter<int>("nHitsOff", 9999)),
   trkAlgoOff_(iConfig.getUntrackedParameter< std::vector<int> >("trkAlgoOff")),
// #Reference tracks (nref)
   MinEtaRef_(iConfig.getUntrackedParameter<double>("minEtaRef",-2.4)),
   MaxEtaRef_(iConfig.getUntrackedParameter<double>("maxEtaRef", 2.4)),
   MinPtRef_(iConfig.getUntrackedParameter<double>("minPtRef",0.4)),
   MaxPtRef_(iConfig.getUntrackedParameter<double>("maxPtRef", 9999.)),
   ChargeRef_(iConfig.getUntrackedParameter< std::vector<int> >("chargeRef")),
   isPixelTrackingRef_(iConfig.getUntrackedParameter<bool>("isPixTrkRef", 9999.)),
   dzdzErrorRef_(iConfig.getUntrackedParameter<double>("dzdzErrorRef", 3.)),
   d0d0ErrorRef_(iConfig.getUntrackedParameter<double>("d0d0ErrorRef", 3.)),
   PtErrorPtRef_(iConfig.getUntrackedParameter<double>("PtErrorPtRef", 1.)),
   Chi2nRef_(iConfig.getUntrackedParameter<double>("Chi2nRef", 0.)),
   nHitsRef_(iConfig.getUntrackedParameter<int>("nHitsRef", 9999)),
   trkAlgoRef_(iConfig.getUntrackedParameter< std::vector<int> >("trkAlgoRef"))
{
   //now do what ever initialization is needed
   edm::Service<TFileService> fs;

   flowTree_ = fs->make<TTree>("flowTree", "flowTree");
   flowTree_->Branch("cent", &Cent_, "cent/D");
   flowTree_->Branch("noff", &nOff_, "noff/I");
   flowTree_->Branch("nref", &nRef_, "nref/I");
   flowTree_->Branch("nVtx", &nVtx_, "nVtx/I");
   flowTree_->Branch("xVtx", &xVtx_, "xVtx/D");
   flowTree_->Branch("yVtx", &yVtx_, "yVtx/D");
   flowTree_->Branch("zVtx", &zVtx_, "zVtx/D");
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
   if(!isEventSelected(iEvent, vertices, hiBin, tracks)) return;

   nRef_ = 0;
   for(reco::TrackCollection::const_iterator itTrk = tracks->begin(); itTrk != tracks->end(); ++itTrk)
   {
      if(!isGoodTrack(*itTrk)) continue;
   }

   flowTree_->Fill();
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
                                     int hiBin, 
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
 
 
   // ## Centrality selection ##
   if((hiBin < 2*CentMin_ || hiBin >= 2*CentMax_) &&
      (CentMin_ != -1 && CentMax_ != -1)) 
                                        return false;
 
 
   // ## Multiplicity ##
   getNoff(tracks); 
   if(nOff_ < NoffMin_ || nOff_ > NoffMax_) return false;
 
 
   return true;
}

void 
QVectorTreeProducer::getNoff(edm::Handle< reco::TrackCollection > tracks)
{
   nOff_ = 0;
 
   for(reco::TrackCollection::const_iterator itTrk = tracks->begin(); itTrk != tracks->end(); ++itTrk)
   {
      //--- Multiplicity info (nOff)
      // Select tracks based on proximity to best vertex
      math::XYZPoint bestVtxPoint(xVtx_,yVtx_,zVtx_);
 
      double  dzvtx   = itTrk->dz(bestVtxPoint);
      double  d0vtx   = itTrk->dxy(bestVtxPoint);
      double  dzerror = sqrt(itTrk->dzError()*itTrk->dzError()+zVtxError_*zVtxError_);
      double  d0error = sqrt(itTrk->d0Error()*itTrk->d0Error()+xVtxError_*yVtxError_);
      double  eta     = itTrk->eta();
      double  pt      = itTrk->pt();
      double  pterror = itTrk->ptError();
      double chi2n   = itTrk->normalizedChi2();
      double nlayers = itTrk->hitPattern().trackerLayersWithMeasurement();
      if(nlayers > 0.) chi2n = chi2n/nlayers;
      else             chi2n = 999.;
      int nhits      = itTrk->numberOfValidHits();
      int algo       = itTrk->originalAlgo();
      int charge     = itTrk->charge(); 
 
      //--- Multiplicity selection
      //if(!itTrk->quality(reco::TrackBase::qualityByName(trackQualityTag_))) continue;
      if(!itTrk->quality(reco::TrackBase::highPurity)) continue;
      if(eta < MinEtaOff_ || eta > MaxEtaOff_)         continue;
      if(pt < MinPtOff_ || pt > MaxPtOff_)             continue;
      if(find(ChargeOff_.begin(), ChargeOff_.end(), charge) == ChargeOff_.end() &&
         !ChargeOff_.empty())                          continue; 
 
      if(isPixelTrackingOff_)
      {
         if(pt < 2.4 && (nhits < 3  || nhits > 6))     continue;
         else
         {
            if(fabs(dzvtx/dzerror) > dzdzErrorOff_)    continue;
            if(fabs(d0vtx/d0error) > d0d0ErrorOff_)    continue;
            if(fabs(pterror)/pt > PtErrorPtOff_)       continue;
            if(chi2n > Chi2nOff_)                      continue;
            if(nhits < nHitsOff_)                      continue;
            if(find(trkAlgoOff_.begin(), trkAlgoOff_.end(), algo) == trkAlgoOff_.end() &&
               !trkAlgoOff_.empty())                   continue;
         }
      }
      else
      {
         if(fabs(dzvtx/dzerror) > dzdzErrorOff_)       continue;
         if(fabs(d0vtx/d0error) > d0d0ErrorOff_)       continue;
         if(fabs(pterror)/pt > PtErrorPtOff_)          continue;
         if(chi2n > Chi2nOff_)                         continue;
         if(nhits < nHitsOff_)                         continue;
         if(find(trkAlgoOff_.begin(), trkAlgoOff_.end(), algo) == trkAlgoOff_.end() &&
            !trkAlgoOff_.empty())                      continue;
      }
 
      nOff_++;
   }
}

bool 
QVectorTreeProducer::isGoodTrack(const reco::Track & trk)
{
   //--- Multiplicity info (nRef)
   // Select tracks based on proximity to best vertex
   math::XYZPoint bestVtxPoint(xVtx_,yVtx_,zVtx_);
 
   double  dzvtx   = trk.dz(bestVtxPoint);
   double  d0vtx   = trk.dxy(bestVtxPoint);
   double  dzerror = sqrt(trk.dzError()*trk.dzError()+zVtxError_*zVtxError_);
   double  d0error = sqrt(trk.d0Error()*trk.d0Error()+xVtxError_*yVtxError_);
   double  eta     = trk.eta();
   double  pt      = trk.pt();
   double  pterror = trk.ptError();
   double chi2n   = trk.normalizedChi2();
   double nlayers = trk.hitPattern().trackerLayersWithMeasurement();
   if(nlayers > 0.) chi2n = chi2n/nlayers;
   else             chi2n = 999.;
   int nhits      = trk.numberOfValidHits();
   int algo       = trk.originalAlgo();
   int charge     = trk.charge(); 
 
   //--- Multiplicity selection
   if(trk.quality(reco::TrackBase::qualityByName(trackQualityTag_))) return false;
   if(eta < MinEtaRef_ || eta > MaxEtaRef_)                          return false;
   if(pt < MinPtRef_ || pt > MaxPtRef_)                              return false;
   if(find(ChargeRef_.begin(), ChargeRef_.end(), charge) == ChargeRef_.end() &&
      !ChargeRef_.empty())                                           return false; 
 
   if(isPixelTrackingRef_)
   {
      if(pt < 2.4 && (nhits < 3  || nhits > 6))                      return false;
      else
      {
         if(fabs(dzvtx/dzerror) > dzdzErrorRef_)                     return false;
         if(fabs(d0vtx/d0error) > d0d0ErrorRef_)                     return false;
         if(fabs(pterror)/pt > PtErrorPtRef_)                        return false;
         if(chi2n > Chi2nRef_)                                       return false;
         if(nhits < nHitsRef_)                                       return false;
         if(find(trkAlgoRef_.begin(), trkAlgoRef_.end(), algo) == trkAlgoRef_.end() &&
            !trkAlgoRef_.empty())                                    return false;
      }
   }
   else
   {
      if(fabs(dzvtx/dzerror) > dzdzErrorRef_)                        return false;
      if(fabs(d0vtx/d0error) > d0d0ErrorRef_)                        return false;
      if(fabs(pterror)/pt > PtErrorPtRef_)                           return false;
      if(chi2n > Chi2nRef_)                                          return false;
      if(nhits < nHitsRef_)                                          return false;
      if(find(trkAlgoRef_.begin(), trkAlgoRef_.end(), algo) == trkAlgoRef_.end() &&
         !trkAlgoRef_.empty())                                       return false;
   }
 
   return true;
}

//define this as a plug-in
DEFINE_FWK_MODULE(QVectorTreeProducer);
