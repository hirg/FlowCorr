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
corr(std::vector< std::complex<double> >(1, std::complex<double>(0.,0.))),
vHarm_(std::vector<int>(1,0))
{

}

FlowQVector::FlowQVector(std::vector<int> vHarm) :
vHarm_(vHarm)
{
   if(find(vHarm_.begin(), vHarm_.end(), 0) == vHarm_.end() && !vHarm_.empty())
   {
     edm::LogInfo("Missing crucial harmonic") << "Harmonic 0 is not present in the harmonic vector \n" 
                                              << "We will insert it in the vector"; 
     vHarm_.insert(vHarm_.begin(), 0);
   }
   corr.resize(vHarm_.size());

   for(unsigned int iharm = 0; iharm < corr.size(); ++iharm)
   {  
      corr[iharm] = std::complex<double>(0, 0);
   } 
}

FlowQVector::~FlowQVector()
{
   corr.clear();
}

void
FlowQVector::fill(double phi, double weight)
{

   for(unsigned int iharm = 0; iharm < vHarm_.size(); iharm++)
   {
      corr[iharm] += std::complex<double>(weight * cos(vHarm_[iharm] * phi), weight * sin(vHarm_[iharm] * phi));
   }
}

void
FlowQVector::reset()
{
   corr.clear();
   for(unsigned int iharm = 0; iharm < vHarm_.size(); iharm++)
   {
      corr[iharm] = std::complex<double>(0., 0.);
   }
}
