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

#include <TString.h>

QVectorTreeProducer::QVectorTreeProducer(const edm::ParameterSet& iConfig) :
// #Tracks
   trackToken_(consumes<reco::TrackCollection>(iConfig.getParameter<edm::InputTag>("tracks"))),
   trackQualityTag_(iConfig.getUntrackedParameter<std::string>("trackQuality","highPurity")),
   MinNoff_(iConfig.getUntrackedParameter<int>("minNoff",0)),
   MaxNoff_(iConfig.getUntrackedParameter<int>("maxNoff", 9999)),
// -- #Offline tracks (noff)
   MinEtaOff_(iConfig.getUntrackedParameter<double>("minEtaOff",-2.4)),
   MaxEtaOff_(iConfig.getUntrackedParameter<double>("maxEtaOff", 2.4)),
   MinPtOff_(iConfig.getUntrackedParameter<double>("minPtOff",0.4)),
   MaxPtOff_(iConfig.getUntrackedParameter<double>("maxPtOff", 9999.)),
   ChargeOff_(iConfig.getUntrackedParameter< std::vector<int> >("chargeOff")),
   isPixelTrackingOff_(iConfig.getUntrackedParameter<bool>("isPixTrkOff", 9999.)),
   dzdzErrorOff_(iConfig.getUntrackedParameter<double>("dzdzErrorOff", 3.)),
   d0d0ErrorOff_(iConfig.getUntrackedParameter<double>("d0d0ErrorOff", 3.)),
   PtErrorPtOff_(iConfig.getUntrackedParameter<double>("ptErrorPtOff", 1.)),
   Chi2nOff_(iConfig.getUntrackedParameter<double>("chi2nOff", 0.)),
   nHitsOff_(iConfig.getUntrackedParameter<int>("nHitsOff", 9999)),
   trkAlgoOff_(iConfig.getUntrackedParameter< std::vector<int> >("trkAlgoOff")),
// -- #Reference tracks (nref)
   MinEtaRef_(iConfig.getUntrackedParameter<double>("minEtaRef",-2.4)),
   MaxEtaRef_(iConfig.getUntrackedParameter<double>("maxEtaRef", 2.4)),
   MinPtRef_(iConfig.getUntrackedParameter<double>("minPtRef",0.4)),
   MaxPtRef_(iConfig.getUntrackedParameter<double>("maxPtRef", 9999.)),
   ChargeRef_(iConfig.getUntrackedParameter< std::vector<int> >("chargeRef")),
   isPixelTrackingRef_(iConfig.getUntrackedParameter<bool>("isPixTrkRef", 9999.)),
   dzdzErrorRef_(iConfig.getUntrackedParameter<double>("dzdzErrorRef", 3.)),
   d0d0ErrorRef_(iConfig.getUntrackedParameter<double>("d0d0ErrorRef", 3.)),
   PtErrorPtRef_(iConfig.getUntrackedParameter<double>("ptErrorPtRef", 1.)),
   Chi2nRef_(iConfig.getUntrackedParameter<double>("chi2nRef", 0.)),
   nHitsRef_(iConfig.getUntrackedParameter<int>("nHitsRef", 9999)),
   trkAlgoRef_(iConfig.getUntrackedParameter< std::vector<int> >("trkAlgoRef")),
// #Vertex
   vtxToken_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertex"))),
   MinVz_(iConfig.getUntrackedParameter<double>("minVz",-15.)),
   MaxVz_(iConfig.getUntrackedParameter<double>("maxVz", 15.)),
   MaxRho_(iConfig.getUntrackedParameter<double>("maxRho", 0.2)),
   nVtxMax_(iConfig.getUntrackedParameter<int>("nVtxMax", 9999)),
// #CaloTower
   caloTowerToken_(consumes<CaloTowerCollection>(iConfig.getParameter<edm::InputTag>("caloTower"))),
// #Centrality
   centralityToken_(consumes<reco::Centrality>(iConfig.getParameter<edm::InputTag>("centralitySrc"))),
   centralityBinToken_(consumes<int>(iConfig.getUntrackedParameter<edm::InputTag>("centralityBinSrc"))),
   MinCent_(iConfig.getUntrackedParameter<int>("minCent", -1)),
   MaxCent_(iConfig.getUntrackedParameter<int>("maxCent", -1)),
// #Event plane
   evtPlaneTag_(consumes<reco::EvtPlaneCollection>(iConfig.getParameter<edm::InputTag>("evtPlane"))),
   epLvl_(iConfig.getUntrackedParameter<int>("epLvl", 2)),
// #cumulants
   cMode_(iConfig.getUntrackedParameter<int>("cMode")),
   cWeight_(iConfig.getUntrackedParameter<bool>("cWeight")),
// file acc & eff & fake
   fName_(iConfig.getUntrackedParameter<edm::InputTag>("fName"))
{
   //now do what ever initialization is needed
   //file acc & eff
   TString filename(fName_.label().c_str());
   fEff_ = 0x0;

   if(cWeight_ && !filename.IsNull())
   {
      edm::FileInPath fip(Form("FlowCorr/QVectorTreeProducer/data/%s",filename.Data()));
      fEff_ = new TFile(fip.fullPath().c_str(),"READ");

      hEff_.resize(fEff_->GetNkeys());
      for(int ihist = 0; ihist < fEff_->GetNkeys(); ++ihist)
      {
         hEff_.push_back((TH2D*) fEff_->Get(fEff_->GetListOfKeys()->At(ihist)->GetName()));
      }
      edm::LogInfo("Input file for accXeff corrections") <<"Using file " << fEff_->GetName();
   }


   //TTree
   edm::Service<TFileService> fs;

   flowTree_ = fs->make<TTree>("flowTree", "flowTree");
   flowTree_->Branch("cent",       &Cent_,   "cent/I");
   flowTree_->Branch("noff",       &nOff_,   "noff/I");
   flowTree_->Branch("nref",       &nRef_,   "nref/I");
   flowTree_->Branch("Vertex",     &Vtx_,    "xVtx/D:yVtx:zVtx:nVtx/I");
   flowTree_->Branch("HFmEPangle", &HFmPsi_, "psi1/D:psi2:psi3:psi4");
   flowTree_->Branch("HFpEPangle", &HFpPsi_, "psi1/D:psi2:psi3:psi4");
   flowTree_->Branch("HFEPangle",  &HFPsi_,  "psi1/D:psi2:psi3:psi4");
}

QVectorTreeProducer::~QVectorTreeProducer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
   if(fEff_) fEff_->Close();
   hEff_.clear();
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
   Cent_ = -1;
   edm::Handle< reco::Centrality > centrality;
   iEvent.getByToken(centralityToken_, centrality);

   edm::Handle< int > cbin;
   iEvent.getByToken(centralityBinToken_,cbin);
   Cent_ = *cbin;
   if(Cent_ < 0) 
   { 
       edm::LogWarning ("Invalid value") <<"Invalid centrality value";
   }
   
   //----- Event Plane selection ----
   //Index     Name   Detector Order hmin1 hmax1 hmin2 hmax2 minpt maxpt
   //    0      HFm1        HF     1 -5.00 -3.00  0.00  0.00  0.01 30.00
   //    1      HFp1        HF     1  3.00  5.00  0.00  0.00  0.01 30.00
   //    2       HF1        HF     1 -5.00 -3.00  3.00  5.00  0.01 30.00
   //    6      HFm2        HF     2 -5.00 -3.00  0.00  0.00  0.01 30.00
   //    7      HFp2        HF     2  3.00  5.00  0.00  0.00  0.01 30.00
   //    8       HF2        HF     2 -5.00 -3.00  3.00  5.00  0.01 30.00
   //   13      HFm3        HF     3 -5.00 -3.00  0.00  0.00  0.01 30.00
   //   14      HFp3        HF     3  3.00  5.00  0.00  0.00  0.01 30.00
   //   15       HF3        HF     3 -5.00 -3.00  3.00  5.00  0.01 30.00
   //   19      HFm4        HF     4 -5.00 -3.00  0.00  0.00  0.01 30.00
   //   20      HFp4        HF     4  3.00  5.00  0.00  0.00  0.01 30.00
   //   21       HF4        HF     4 -5.00 -3.00  3.00  5.00  0.01 30.00
   //   25    HFm1mc        HF     1 -5.00 -3.00  0.00  0.00  0.01 30.00
   //   26    HFp1mc        HF     1  3.00  5.00  0.00  0.00  0.01 30.00
   unsigned int HFm[4] = {0, 6, 13, 19};
   unsigned int HFp[4] = {1, 7, 14, 20};
   unsigned int HF[4]  = {2, 8, 15, 21};


   edm::Handle<reco::EvtPlaneCollection> evtPlanes;
   iEvent.getByToken(evtPlaneTag_, evtPlanes);
   if(evtPlanes.isValid())
   {
      //psi1
      if(HFm[0] < evtPlanes->size()) HFmPsi_.psi1 = (*evtPlanes)[HFm[0]].angle(epLvl_);
      if(HFp[0] < evtPlanes->size()) HFpPsi_.psi1 = (*evtPlanes)[HFp[0]].angle(epLvl_);
      if(HF[0]  < evtPlanes->size()) HFPsi_.psi1  = (*evtPlanes)[HF[0]].angle(epLvl_);

      //psi2
      if(HFm[1] < evtPlanes->size()) HFmPsi_.psi2 = (*evtPlanes)[HFm[1]].angle(epLvl_);
      if(HFp[1] < evtPlanes->size()) HFpPsi_.psi2 = (*evtPlanes)[HFp[1]].angle(epLvl_);
      if(HF[1]  < evtPlanes->size()) HFPsi_.psi2  = (*evtPlanes)[HF[1]].angle(epLvl_);

      //psi3
      if(HFm[2] < evtPlanes->size()) HFmPsi_.psi3 = (*evtPlanes)[HFm[2]].angle(epLvl_);
      if(HFp[2] < evtPlanes->size()) HFpPsi_.psi3 = (*evtPlanes)[HFp[2]].angle(epLvl_);
      if(HF[2]  < evtPlanes->size()) HFPsi_.psi3  = (*evtPlanes)[HF[2]].angle(epLvl_);

      //psi4
      if(HFm[3] < evtPlanes->size()) HFmPsi_.psi4 = (*evtPlanes)[HFm[3]].angle(epLvl_);
      if(HFp[3] < evtPlanes->size()) HFpPsi_.psi4 = (*evtPlanes)[HFp[3]].angle(epLvl_);
      if(HF[3]  < evtPlanes->size()) HFPsi_.psi4  = (*evtPlanes)[HF[3]].angle(epLvl_);
   }
   else
   {
       edm::LogWarning ("Invalid value") << "Invalid EP value";
   }

   //----- Event selection -----
   if(!isEventSelected(iEvent, vertices, tracks)) return;
  
   nRef_ = 0;
   for(reco::TrackCollection::const_iterator itTrk = tracks->begin(); itTrk != tracks->end(); ++itTrk)
   {
      if(!isGoodTrack(*itTrk)) continue;
      nRef_++;
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
                                     edm::Handle< reco::TrackCollection >  tracks) 
{
   // ## Vertex ##
   //--- Vertex info 
   const reco::Vertex & bestVtx = *(vertices->begin());

   Vtx_.xVtx = 0.;
   Vtx_.yVtx = 0.;
   Vtx_.zVtx = 0.;
   Vtx_.nVtx = 0;
    
   if(!bestVtx.isFake() && bestVtx.tracksSize()>=2)
   {
      Vtx_.xVtx = bestVtx.x();
      Vtx_.yVtx = bestVtx.y();
      Vtx_.zVtx = bestVtx.z();
      xVtxError_ = bestVtx.xError();
      yVtxError_ = bestVtx.yError();
      zVtxError_ = bestVtx.zError();
   }
   double rho = sqrt(Vtx_.xVtx*Vtx_.xVtx + Vtx_.yVtx*Vtx_.yVtx);
 
   for(reco::VertexCollection::const_iterator itVtx = vertices->begin(); itVtx != vertices->end(); ++itVtx)
   {
       if(!itVtx->isFake() && itVtx->tracksSize() >= 2) Vtx_.nVtx++;
   }
 
   //--- Vertex selection 
   if(Vtx_.nVtx > nVtxMax_)                     return false;
   if(Vtx_.zVtx < MinVz_ || Vtx_.zVtx > MaxVz_) return false;
   if(rho > MaxRho_ )                           return false;
 
 
   // ## Centrality selection ##
   if((Cent_ < 2*MinCent_ || Cent_ >= 2*MaxCent_) &&
      (MaxCent_ != -1 && MaxCent_ != -1)) 
                                        return false;
 
 
   // ## Multiplicity ##
   getNoff(tracks); 
   if(nOff_ < MinNoff_ || nOff_ > MaxNoff_) return false;
 
 
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
      math::XYZPoint bestVtxPoint(Vtx_.xVtx,Vtx_.yVtx,Vtx_.zVtx);
 
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
   math::XYZPoint bestVtxPoint(Vtx_.xVtx,Vtx_.yVtx,Vtx_.zVtx);
 
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
