#ifndef FLOWTOWERSELECTION_H
#define FLOWTOWERSELECTION_H

// -*- C++ -*-
//
// Package:    FlowCorr/QVectorTreeProducer
// Class:      FlowTowerSelection
// 
/**\class FlowTowerSelection FlowTowerSelection.h FlowCorr/FlowTowerSelection/interface/FlowTowerSelection.h

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
#include "DataFormats/CaloTowers/interface/CaloTower.h"
#include "DataFormats/CaloTowers/interface/CaloTowerFwd.h"

#include "FlowCorr/QVectorTreeProducer/interface/FlowUtils.h"

//
// Type def
//
typedef struct {
   double etamin;
   double etamax;
   double etmin;
} sTowCut;

typedef struct {
   double eta;
   double et;
   double phi;
   double weight;
} sTowSel;

class FlowTowerSelection {
   public:

      explicit FlowTowerSelection();
      explicit FlowTowerSelection(sTowCut towCuts);
      ~FlowTowerSelection(){};

      void fillVar(const CaloTower & tow);
      void resetVar();
      bool isTowerPassCuts();
      sTowSel getTower() {return towVar_;};
   private:
      sTowCut towCuts_;
      sTowSel towVar_;
};

#endif

