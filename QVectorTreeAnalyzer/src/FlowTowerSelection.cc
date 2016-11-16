// -*- C++ -*-
//
// Package:    FlowCorr/QVectorTreeProducer
// Class:      FlowTowerSelection
// 
/**\class FlowTowerSelection FlowTowerSelection.cc FlowCorr/FlowTowerSelection/src/FlowTowerSelection.cc

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
#include "FlowCorr/QVectorTreeProducer/interface/FlowTowerSelection.h"

FlowTowerSelection::FlowTowerSelection()
{
   towCuts_.etamin =  3; 
   towCuts_.etamax =  5; 
   towCuts_.etmin  =  0.; 
}

FlowTowerSelection::FlowTowerSelection(sTowCut towCuts) :
towCuts_(towCuts)
{

}

void
FlowTowerSelection::fillVar(const CaloTower & tow)
{
   // Select tower
   double eta  = tow.eta();
   double et   = tow.et();
   double phi  = tow.phi();
   double w    = 1.;
   if(et != 0.) w /= et;

   towVar_.eta    = eta;
   towVar_.et     = et;
   towVar_.phi    = phi;
   towVar_.weight = w;
}

void
FlowTowerSelection::resetVar()
{
   towVar_.eta    = -999.;
   towVar_.et     = -999.;
   towVar_.phi    = -999.;
   towVar_.weight = 1.;
}

bool
FlowTowerSelection::isTowerPassCuts()
{
   //Sel 
   if(fabs(towVar_.eta) < towCuts_.etamin || fabs(towVar_.eta) > towCuts_.etamax) return false;
   if(towVar_.et  < towCuts_.etmin)                                               return false;

   return true;
}
