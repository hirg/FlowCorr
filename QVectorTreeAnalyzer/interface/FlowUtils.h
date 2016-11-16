#ifndef FLOWUTILS_H
#define FLOWUTILS_H
// -*- C++ -*-
//
// Package:    FlowCorr/QVectorTreeAnalyzer
// Class:      FlowUtils
// 
/**\class FlowUtils FlowUtils.h FlowCorr/FlowEP/interface/FlowUtils.h

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
typedef struct {
   double noffcorr;
   double nrefcorr;
   int run;
   int event;
   int lumi;
   int cent;
   int noff;
   int nref;
} sEvent;

typedef struct {
  double psi1;
  double psi2;
  double psi3;
  double psi4;
} sFlowEPangle;

typedef struct {
  double xVtx;
  double yVtx;
  double zVtx;
  double xVtxError;
  double yVtxError;
  double zVtxError;
  int nVtx;
} sVertex;


#endif
