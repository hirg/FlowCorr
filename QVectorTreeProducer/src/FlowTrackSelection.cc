// -*- C++ -*-
//
// Package:    FlowCorr/QVectorTreeProducer
// Class:      FlowTrackSelection
// 
/**\class FlowTrackSelection FlowTrackSelection.cc FlowCorr/FlowTrackSelection/src/FlowTrackSelection.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Maxime Guilbaud
//         Created:  Fri, 4 Oct 2016 14:12:00 GMT
//
//


// system include files
#include "FlowCorr/QVectorTreeProducer/interface/FlowTrackSelection.h"

FlowTrackSelection::FlowTrackSelection()
{
   trkCuts_.dzvtxdzerror =  3.0; 
   trkCuts_.d0vtxd0error =  3.0; 
   trkCuts_.etamin       = -2.4; 
   trkCuts_.etamax       =  2.4; 
   trkCuts_.ptmin        =  0.3; 
   trkCuts_.ptmax        =  3.0; 
   trkCuts_.pterrorpt    =  0.1; 
   trkCuts_.chi2nnlayers =  0.15; 
   trkCuts_.nhits        =  11; 
   trkCuts_.algo         =  std::vector<int>(); 
   trkCuts_.charge       =  std::vector<int>(); 
}

FlowTrackSelection::FlowTrackSelection(sTrkCut trkCuts) :
trkCuts_(trkCuts)
{

}

void
FlowTrackSelection::fillVar(const reco::Track & trk,  sVertex vtx, std::string trkQualityTag)
{
   // Select tracks based on proximity to best vertex
   math::XYZPoint bestVtxPoint(vtx.xVtx,vtx.yVtx,vtx.zVtx);

   double dzvtx   = trk.dz(bestVtxPoint);
   double d0vtx   = trk.dxy(bestVtxPoint);
   double dzerror = sqrt(trk.dzError()*trk.dzError()+vtx.zVtxError*vtx.zVtxError);
   double d0error = sqrt(trk.d0Error()*trk.d0Error()+vtx.xVtxError*vtx.yVtxError);
   double eta     = trk.eta();
   double pt      = trk.pt();
   double phi     = trk.phi();
   double pterror = trk.ptError();
   double chi2n   = trk.normalizedChi2();
   double nlayers = trk.hitPattern().trackerLayersWithMeasurement();
   int nhits      = trk.numberOfValidHits();
   int algo       = trk.originalAlgo();
   int charge     = trk.charge();

   trkVar_.quality = trk.quality(reco::TrackBase::qualityByName(trkQualityTag)); 
   trkVar_.dzvtx   = dzvtx;
   trkVar_.d0vtx   = d0vtx;
   trkVar_.dzerror = dzerror;
   trkVar_.d0error = d0error;
   trkVar_.eta     = eta;
   trkVar_.pt      = pt;
   trkVar_.phi     = phi;
   trkVar_.pterror = pterror;
   trkVar_.chi2n   = chi2n;
   trkVar_.nlayers = nlayers;
   trkVar_.nhits   = nhits;
   trkVar_.algo    = algo;
   trkVar_.charge  = charge;
}

void
FlowTrackSelection::resetVar()
{
   trkVar_.quality =  false; 
   trkVar_.dzvtx   = -999.;
   trkVar_.d0vtx   = -999.;
   trkVar_.dzerror =  1.;
   trkVar_.d0error =  1.;
   trkVar_.eta     = -999.;
   trkVar_.pt      = -999.;
   trkVar_.phi     = -999.;
   trkVar_.pterror =  1.;
   trkVar_.chi2n   =  999.;
   trkVar_.nlayers =  1.;
   trkVar_.nhits   =  0;
   trkVar_.algo    =  999;
   trkVar_.charge  = -999;
}

bool
FlowTrackSelection::isTrkPassCuts(bool isPixTrk)
{
   //Sel 
   double dzvtxdzerror = 999.;
   if(trkVar_.dzerror != 0.) dzvtxdzerror = trkVar_.dzvtx/trkVar_.dzerror; 
   double d0vtxd0error = 999.;
   if(trkVar_.d0error != 0.) d0vtxd0error = trkVar_.d0vtx/trkVar_.d0error; 
   double pterrorpt    = 999.; 
   if(trkVar_.pt != 0.)      pterrorpt    = trkVar_.pterror/trkVar_.pt; 
   double chi2n        = 999.;
   if(trkVar_.nlayers != 0.) chi2n        = trkVar_.chi2n/trkVar_.nlayers;

   if(!trkVar_.quality)                                               return false;
   if(trkVar_.eta < trkCuts_.etamin || trkVar_.eta > trkCuts_.etamax) return false;
   if(trkVar_.pt  < trkCuts_.ptmin  || trkVar_.pt  > trkCuts_.ptmax)  return false;
   if(find(trkCuts_.charge.begin(), trkCuts_.charge.end(), trkVar_.charge) == trkCuts_.charge.end() &&
      !trkCuts_.charge.empty())                                       return false;

   if(isPixTrk && trkVar_.pt < 2.4 && (trkVar_.nhits < 3  || trkVar_.nhits > 6))
   {
      return false;
   }
   else
   {
      if(fabs(dzvtxdzerror) > trkCuts_.dzvtxdzerror) return false;
      if(fabs(d0vtxd0error) > trkCuts_.d0vtxd0error) return false;
      if(fabs(pterrorpt) > trkCuts_.pterrorpt)       return false;
      if(chi2n > trkCuts_.chi2nnlayers)              return false;
      if(trkVar_.nhits < trkCuts_.nhits)             return false;
      if(find(trkCuts_.algo.begin(), trkCuts_.algo.end(), trkVar_.algo) == trkCuts_.algo.end() &&
         !trkCuts_.algo.empty())                     return false;
   }

   return true;
}
