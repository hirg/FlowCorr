// -*- C++ -*-
//
// Package:    FlowCorr/QVectorTreeProducer
// Class:      FlowQVector
// 
/**\class FlowQVector FlowQVector.cc FlowCorr/QVectorTreeProducer/src/FlowQVector.cc

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
#include <iostream> 
#include <algorithm> 

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FlowCorr/QVectorTreeProducer/interface/FlowQVector.h"

FlowQVector::FlowQVector() :
qv(std::vector< std::complex<double> >(1, std::complex<double>(0.,0.)))
{

}

FlowQVector::FlowQVector(std::vector<int> vHarm)
{
   if(find(vHarm.begin(), vHarm.end(), 0) == vHarm.end() && !vHarm.empty())
   {
     edm::LogInfo("Missing crucial harmonic") << "Harmonic 0 is not present in the harmonic vector \n" 
                                              << "We will insert it in the vector"; 
     vHarm.insert(vHarm.begin(), 0);
   }
   qv.resize(vHarm.size());

   for(unsigned int iharm = 0; iharm < qv.size(); ++iharm)
   {  
      qv[iharm] = std::complex<double>(0, 0);
   } 
}

FlowQVector::~FlowQVector()
{
   qv.clear();
}

void
FlowQVector::fill(double phi, double weight, std::vector<int> vHarm)
{

   for(unsigned int iharm = 0; iharm < qv.size(); iharm++)
   {
      if(iharm == 0) qv[iharm] += std::complex<double>(weight * cos(0. * phi), 
                                                       weight * sin(0. * phi));
      else           qv[iharm] += std::complex<double>(weight * cos(vHarm[iharm-1] * phi), 
                                                       weight * sin(vHarm[iharm-1] * phi));
   }
}

void
FlowQVector::reset()
{
   qv.clear();
   for(unsigned int iharm = 0; iharm < qv.size(); iharm++)
   {
      qv[iharm] = std::complex<double>(0., 0.);
   }
}
