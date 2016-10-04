// -*- C++ -*-
//
// Package:    FlowCorr/QVectorTreeProducer
// Class:      QVectorTreeProducer
// 
/**\class QVectorTreeProducer QVectorTreeProducer.h FlowCorr/QVectorTreeProducer/interface/QVectorTreeProducer.h

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Maxime Guilbaud
//         Created:  Fri, 23 Sep 2016 12:48:03 GMT
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

//#include "DataFormats/HeavyIonEvent/interface/CentralityBins.h"
#include "DataFormats/HeavyIonEvent/interface/Centrality.h"
#include "DataFormats/HeavyIonEvent/interface/EvtPlane.h"

#include "DataFormats/CaloTowers/interface/CaloTower.h"
#include "DataFormats/CaloTowers/interface/CaloTowerFwd.h"

#include <TComplex.h>
#include <TTree.h>

//
// typedef
//
typedef struct {
  double xVtx;
  double yVtx;
  double zVtx;
  int nVtx;
} sVertex;

typedef struct {
  double psi1;  
  double psi2;
  double psi3;
  double psi4;
} sFlowEPangle;

//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<> and also remove the line from
// constructor "usesResource("TFileService");"
// This will improve performance in multithreaded jobs.

class QVectorTreeProducer : public edm::one::EDAnalyzer<>  {
   public:
      explicit QVectorTreeProducer(const edm::ParameterSet&);
      ~QVectorTreeProducer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      // ---------- member data -----------
      // ~~~> Config parameters <~~~ 
      // ## tracks ##
      edm::EDGetTokenT<reco::TrackCollection> trackToken_;
      std::string                             trackQualityTag_;
      int MinNoff_;
      int MaxNoff_;
      // -- ## offline tracks ##
      double MinEtaOff_;
      double MaxEtaOff_;
      double MinPtOff_;
      double MaxPtOff_;
      std::vector<int> ChargeOff_;
      bool isPixelTrackingOff_;
      double dzdzErrorOff_;
      double d0d0ErrorOff_;
      double PtErrorPtOff_;
      double Chi2nOff_;
      int nHitsOff_;
      std::vector<int> trkAlgoOff_;
      // -- ## ref tracks ##
      double MinEtaRef_;
      double MaxEtaRef_;
      double MinPtRef_;
      double MaxPtRef_;
      std::vector<int> ChargeRef_;
      bool isPixelTrackingRef_;
      double dzdzErrorRef_;
      double d0d0ErrorRef_;
      double PtErrorPtRef_;
      double Chi2nRef_;
      int nHitsRef_;
      std::vector<int> trkAlgoRef_;

      // ## vertex ##
      edm::EDGetTokenT<reco::VertexCollection>   vtxToken_;
      double MinVz_;
      double MaxVz_;
      double MaxRho_;
      int nVtxMax_;

      // ## calotower ##
      edm::EDGetTokenT<CaloTowerCollection> caloTowerToken_;

      // ## centrality ##	
      edm::EDGetTokenT<reco::Centrality> centralityToken_;
      edm::EDGetTokenT<int>              centralityBinToken_;
      int MinCent_;
      int MaxCent_;

      // ## event plane ##
      edm::EDGetTokenT<reco::EvtPlaneCollection> evtPlaneTag_;
      edm::EDGetTokenT<reco::EvtPlaneCollection> evtPlaneFlatTag_;
      int epLvl_;


      // ~~~> Class parameters <~~~ 
      // ## tracks ##
      int nOff_;
      int nRef_;

      // ## vertex ##
      sVertex Vtx_;
      double xVtxError_;
      double yVtxError_;
      double zVtxError_;

      // ## calotower ##
      double Etmin_;

      // ## centrality ##
      int Cent_;
      
      // ## event plane ##
      sFlowEPangle HFmPsi_;
      sFlowEPangle HFpPsi_;
      sFlowEPangle HFPsi_;

      // ## cumulant ##
      int cMode_;
      bool cWeight_;

      // ## gap ##
      double Gap_;

      // ## harmonics ##
      unsigned int nHarm_;
      std::vector<int> vHarm_;

      // ## output tree and histograms ##
      TTree* flowTree_;

      // ---------- public methods ------------
      bool isEventSelected(const edm::Event& iEvent, 
                           edm::Handle< reco::VertexCollection > vertices,
                           edm::Handle< reco::TrackCollection > tracks);
      void getNoff(edm::Handle< reco::TrackCollection > tracks);
      bool isGoodTrack(const reco::Track & trk);
};

