// -*- C++ -*-
//
// Package:    FlowCorr/QVectorTreeProducer
// Class:      FlowEP
// 
/**\class FlowEPangle FlowEPangle.cc FlowCorr/FlowEP/src/FlowEPangle.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Maxime Guilbaud
//         Created:  Fri, 4 Oct 2016 17:12:00 GMT
//
//


// system include files
#include "FlowCorr/QVectorTreeProducer/interface/FlowEPangle.h"

FlowEPangle::FlowEPangle() :
psi1(0.),
psi2(0.),
psi3(0.),
psi4(0.)
{
}

FlowEPangle::~FlowEPangle()
{
}

void
FlowEPangle::fillEPangle(double angle, int harm)
{
   switch(harm)
   {
      case 1:
         psi1 = angle;
         break;
      case 2:
         psi2 = angle; 
         break;
      case 3:
         psi3 = angle;
         break;
      case 4:
         psi4 = angle;
         break;
      default:
         break;
   }
}
