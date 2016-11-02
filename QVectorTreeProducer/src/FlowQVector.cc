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
#include "FlowCorr/QVectorTreeProducer/interface/FlowQVector.h"

FlowQVector::FlowQVector() :
corr(std::vector< std::complex<double> >(1, std::complex<double>(0.,0.))),
vHarm_(std::vector<int>(1,0))
{

}

FlowQVector::FlowQVector(std::vector<int> vHarm) :
vHarm_(vHarm)
{
   corr.resize(vHarm_.size());
}

FlowQVector::~FlowQVector()
{
   corr.clear();
}
