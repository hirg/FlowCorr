#ifndef FLOWTRACKSELECTION_H
#define FLOWTRACKSELECTION_H

// -*- C++ -*-
//
// Package:    FlowCorr/QVectorTreeProducer
// Class:      FlowTrackSelection
// 
/**\class FlowTrackSelection FlowTrackSelection.h FlowCorr/FlowTrackSelection/interface/FlowTrackSelection.h

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
#include <string>
#include <vector>

// user include files
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

#include "FlowCorr/QVectorTreeProducer/interface/FlowUtils.h"

//
// Type def
//
typedef struct {
   double dzvtxdzerror;
   double d0vtxd0error;
   double etamin;
   double etamax;
   double ptmin;
   double ptmax;
   double pterrorpt;
   double chi2nnlayers;
   int    nhits;
   std::vector<int> algo;
   std::vector<int> charge;
} sTrkCut;

typedef struct {
   bool   quality; 
   double dzvtx;
   double d0vtx;
   double dzerror;
   double d0error;
   double eta;
   double pt;
   double phi;
   double pterror;
   double chi2n;
   double nlayers;
   int    nhits;
   int    algo;
   int    charge;
} sTrkSel;

class FlowTrackSelection {
   public:

      explicit FlowTrackSelection();
      explicit FlowTrackSelection(sTrkCut trkCuts);
      ~FlowTrackSelection(){};

      void fillVar(const reco::Track & trk, sVertex vtx, std::string trkQualityTag);
      void resetVar();
      bool isTrkPassCuts(bool isPixTrk);
      sTrkSel getTrk() {return trkVar_;};
   private:
      sTrkCut trkCuts_;
      sTrkSel trkVar_;
};

#endif

