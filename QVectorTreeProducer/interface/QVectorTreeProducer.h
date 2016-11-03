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
#include <TFile.h>
#include <TH2D.h>

#include <FlowCorr/QVectorTreeProducer/interface/FlowCorrelator.h>
#include <FlowCorr/QVectorTreeProducer/interface/FlowUtils.h>
#include <FlowCorr/QVectorTreeProducer/interface/FlowQVector.h>
#include <FlowCorr/QVectorTreeProducer/interface/FlowTrackSelection.h>
#include <FlowCorr/QVectorTreeProducer/interface/correlations/Types.hh>
#include <FlowCorr/QVectorTreeProducer/interface/correlations/Result.hh>
#include <FlowCorr/QVectorTreeProducer/interface/correlations/QVector.hh>
#include <FlowCorr/QVectorTreeProducer/interface/correlations/FromQVector.hh>
#include <FlowCorr/QVectorTreeProducer/interface/correlations/recursive/FromQVector.hh>
#include <FlowCorr/QVectorTreeProducer/interface/correlations/recurrence/FromQVector.hh>

//
// typedef
//

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
      double MinEtHF_;
      double MinEtaHF_;
      double MaxEtaHF_;

      // ## centrality ##	
      edm::EDGetTokenT<reco::Centrality> centralityToken_;
      edm::EDGetTokenT<int>              centralityBinToken_;
      int MinCent_;
      int MaxCent_;

      // ## event plane ##
      edm::EDGetTokenT<reco::EvtPlaneCollection> evtPlaneTag_;
      edm::EDGetTokenT<reco::EvtPlaneCollection> evtPlaneFlatTag_;
      int epLvl_;

      // ## harmonics ##
      std::vector<int> vHarmTrk_;
      std::vector<int> vHarmHF_;

      // ## cumulant ##
      int  cMode_;
      bool cWeight_;

      // ## Acc x Eff correction ##
      edm::InputTag    fName_;
      std::vector<int> effCentBin_;


      // ~~~> Class parameters <~~~ 
      // ## tracks ##
      sEvent Evt_;
      int nOff_;
      int nRef_;
      double nOffcorr_;
      double nRefcorr_;
      FlowTrackSelection trkSelectorOff_;
      FlowTrackSelection trkSelectorRef_;

      // ## vertex ##
      sVertex Vtx_;

      // ## calotower ##
      double Etmin_;

      // ## centrality ##
      int Cent_;
      
      // ## event plane ##
      sFlowEPangle HFmPsi_;
      sFlowEPangle HFpPsi_;
      sFlowEPangle HFPsi_;

      // ## harmonics ##
      unsigned int nHarmTrk_;
      unsigned int nHarmHF_;

      // ## cumulant ##
      //   --- diagonal terms
      std::vector<correlations::QVector>        qNvec_;
      std::vector<correlations::HarmonicVector> hcNvec_;
      std::vector<correlations::FromQVector*>   cqNvec_;
   
      FlowCorrelator *qN2_;
      FlowCorrelator *qN4_;
      FlowCorrelator *qN6_;
      FlowCorrelator *qN8_;
      //   --- non-diagonal terms
      std::vector< std::vector<correlations::QVector> >        qNMvec_;
      std::vector< std::vector<correlations::HarmonicVector> > hcNMvec_;
      std::vector< std::vector<correlations::FromQVector*> >   cqNMvec_;

      std::vector<FlowCorrelator*> qNM4_;
      //   --- diagonal terms w/ gap

      // ## Qvector eta binning && gap ##
      double EtaBinWidthTrk_;
      unsigned int nEtaBinTrk_;
      double EtaBinWidthHF_;
      unsigned int nEtaBinHF_;

      std::vector<FlowQVector*> vQn_trkP_; 
      std::vector<FlowQVector*> vQn_trkN_; 
      std::vector<FlowQVector*> vQn_hfP_; 
      std::vector<FlowQVector*> vQn_hfM_; 

      // ## Acc x Eff correction ##
      TFile * fEff_;
      std::vector<TH2D*> hEff_;

      // ## output tree and histograms ##
      TH1F* htrk_eta_ref_;
      TH1F* htrk_phi_ref_;
      TH1F* htrk_pt_ref_ ;

      TH1F* htrk_eta_corr_ref_; 
      TH1F* htrk_phi_corr_ref_;
      TH1F* htrk_pt_corr_ref_;

      TH1F* htrk_dzdzerr_ref_; 
      TH1F* htrk_d0d0err_ref_;
      TH1F* htrk_ptpterr_ref_;
      TH1F* htrk_chi2nlayers_ref_;
      TH1I* htrk_nhits_ref_;
      TH1I* htrk_algo_ref_; 
      //--
      TH1F* htrk_eta_off_;
      TH1F* htrk_phi_off_;
      TH1F* htrk_pt_off_ ;

      TH1F* htrk_eta_corr_off_; 
      TH1F* htrk_phi_corr_off_;
      TH1F* htrk_pt_corr_off_;

      TH1F* htrk_dzdzerr_off_; 
      TH1F* htrk_d0d0err_off_;
      TH1F* htrk_ptpterr_off_;
      TH1F* htrk_chi2nlayers_off_;
      TH1I* htrk_nhits_off_;
      TH1I* htrk_algo_off_; 

      TFileDirectory flowHistListRef_;
      TFileDirectory flowHistListOff_;

      TTree* globalTree_;
      TTree* correlatorTree_;
      TTree* qvectorTree_;

      // ---------- public methods ------------
      bool isEventSelected(const edm::Event& iEvent, 
                           edm::Handle< reco::VertexCollection > vertices,
                           edm::Handle< reco::TrackCollection > tracks);
      void getNoff(edm::Handle< reco::TrackCollection > tracks);
      void initDiagQ();
      void initNonDiagQ();
      void doneDiagQ();
      void doneNonDiagQ();
      void resetQvectors();
      double getAccEffWeight(double eta, double pt);
      void fillHisto(double weight);
      void fillQVectorCorr(double weight);
      void fillQVectorTrk(double weight);
      void fillQVectorHF(double weight);
      void fillCorrelators();
      void fillEPangle(const edm::Handle<reco::EvtPlaneCollection> & ep);
};

