// -*- C++ -*-
//
// Package:    FlowCorr/QVectorTreeProducer
// Class:      QVectorTreeProducer
// 
/**\class QVectorTreeProducer QVectorTreeProducer.cc FlowCorr/QVectorTreeProducer/plugins/QVectorTreeProducer.cc

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

#include "DataFormats/HeavyIonEvent/interface/CentralityBins.h"
#include "DataFormats/HeavyIonEvent/interface/Centrality.h"

#include "DataFormats/CaloTowers/interface/CaloTower.h"
#include "DataFormats/CaloTowers/interface/CaloTowerFwd.h"

//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<> and also remove the line from
// constructor "usesResource("TFileService");"
// This will improve performance in multithreaded jobs.

class QVectorTreeProducer : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit QVectorTreeProducer(const edm::ParameterSet&);
      ~QVectorTreeProducer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      // ----------member data ---------------------------
      edm::EDGetTokenT<reco::TrackCollection>  trackToken_;
      std::string                              trackQualityTag_;
      edm::EDGetTokenT<reco::VertexCollection> vtxToken_;
      edm::EDGetTokenT<reco::Centrality>       centralityToken_;
      edm::EDGetTokenT<int>                    centralityBinToken_;
      edm::EDGetTokenT<CaloTowerCollection>    caloTowerToken_;

      // ## global for evt sel ##
      int NoffMin_;
      int NoffMax_;
      int CentMin_;
      int CentMax_;
      // ## vertex ##
      float MinVz_;
      float MaxVz_;
      float MaxRho_;
      float MinEta_;
      float MaxEta_;
      int nVtxMax_;	
      float xVtx_;
      float yVtx_;
      float zVtx_;
      float xVtxError_;
      float yVtxError_;
      float zVtxError_;
      int nVtx_;	
      // ## tracks ##
      float dzdzError_;
      float d0d0Error_;
      //double Chi2_;
      float PtErrorPt_;
      //int Charge_;
      float PtMin_;
      float PtMax_;
      // ## calotower ##
      float Etmin_;
      // ## cumulant ##
      int cMode_;
      bool cWeight_;
      // ## gap ##
      float Gap_;

      unsigned int nHarm_;
      std::vector<int> vHarm_;
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
