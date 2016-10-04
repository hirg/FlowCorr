// -*- C++ -*-
//
// Package:    FlowCorr/QVectorTreeProducer
// Class:      FlowEP
// 
/**\class FlowEPangle FlowEPangle.h FlowCorr/FlowEP/interface/FlowEPangle.h

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
#include <vector>

// user include files

class FlowEPangle {
   public:

      FlowEPangle();
      ~FlowEPangle();

      void fillEPangle(double angle, int harm);

      double psi1;
      double psi2;
      double psi3;
      double psi4;      

   private:
};
