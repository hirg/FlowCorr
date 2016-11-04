#ifndef FLOWQVECTOR_H
#define FLOWQVECTOR_H
// -*- C++ -*-
//
// Package:    FlowCorr/QVectorTreeProducer
// Class:      FlowQVector
// 
/**\class FlowQVector FlowQVector.h FlowCorr/QVectorTreeProducer/interface/FlowQVector.h

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

class FlowQVector {
   public:

      FlowQVector();
      FlowQVector(std::vector<int> vHarm);
      ~FlowQVector();

      void fill(double phi, double weight);
      void reset();

      std::vector< std::complex<double> > corr;
      
   private:
      std::vector<int> vHarm_;
};

typedef std::vector<FlowQVector*> FlowQVectorCollection;

#endif
