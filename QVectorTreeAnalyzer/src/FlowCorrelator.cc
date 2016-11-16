// -*- C++ -*-
//
// Package:    FlowCorr/QVectorTreeProducer
// Class:      FlowCorrelator
// 
/**\class FlowCorrelator FlowCorrelator.cc FlowCorr/QVectorTreeProducer/src/FlowCorrelator.cc

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
#include "FlowCorr/QVectorTreeAnalyzer/interface/FlowCorrelator.h"

FlowCorrelator::FlowCorrelator() :
corr(std::vector< std::complex<double> >(1, std::complex<double>(0.,0.))),
weight(1.)
{

}

FlowCorrelator::FlowCorrelator(unsigned int nHarm) :
weight(1.)
{
   corr.resize(nHarm);
}

FlowCorrelator::~FlowCorrelator()
{
   corr.clear();
}
