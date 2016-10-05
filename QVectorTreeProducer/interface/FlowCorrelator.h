#ifndef FLOWEP_H
#define FLOWEP_H
// -*- C++ -*-
//
// Package:    FlowCorr/QVectorTreeProducer
// Class:      FlowCorrelator
// 
/**\class FlowCorrelator FlowCorrelator.h FlowCorr/QVectorTreeProducer/interface/FlowCorrelator.h

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
#include "TObject.h"
#include <complex>
#include <vector>

// user include files

class FlowCorrelator {
   public:

      FlowCorrelator();
      FlowCorrelator(unsigned int nHarm);
      ~FlowCorrelator();

      std::vector< std::complex<double> > corr;
      double weight;
      
   private:
};

#endif
