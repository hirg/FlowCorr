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
#include <cmath>

// user include files
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "FlowCorr/QVectorTreeProducer/interface/QVectorTreeProducer.h"

#include <TString.h>

QVectorTreeProducer::QVectorTreeProducer(const edm::ParameterSet& iConfig) :
// #Tracks
   trackToken_(consumes<reco::TrackCollection>(iConfig.getParameter<edm::InputTag>("tracks"))),
   trackQualityTag_(iConfig.getUntrackedParameter<std::string>("trackQuality","highPurity")),
   MinNoff_(iConfig.getUntrackedParameter<int>("minNoff",0)),
   MaxNoff_(iConfig.getUntrackedParameter<int>("maxNoff", 9999)),
// -- #Offline tracks (noff)
   MinEtaOff_(iConfig.getUntrackedParameter<double>("minEtaOff",-2.4)),
   MaxEtaOff_(iConfig.getUntrackedParameter<double>("maxEtaOff", 2.4)),
   MinPtOff_(iConfig.getUntrackedParameter<double>("minPtOff",0.4)),
   MaxPtOff_(iConfig.getUntrackedParameter<double>("maxPtOff", 9999.)),
   ChargeOff_(iConfig.getUntrackedParameter< std::vector<int> >("chargeOff")),
   isPixelTrackingOff_(iConfig.getUntrackedParameter<bool>("isPixTrkOff", 9999.)),
   dzdzErrorOff_(iConfig.getUntrackedParameter<double>("dzdzErrorOff", 3.)),
   d0d0ErrorOff_(iConfig.getUntrackedParameter<double>("d0d0ErrorOff", 3.)),
   PtErrorPtOff_(iConfig.getUntrackedParameter<double>("ptErrorPtOff", 1.)),
   Chi2nOff_(iConfig.getUntrackedParameter<double>("chi2nOff", 0.)),
   nHitsOff_(iConfig.getUntrackedParameter<int>("nHitsOff", 9999)),
   trkAlgoOff_(iConfig.getUntrackedParameter< std::vector<int> >("trkAlgoOff")),
// -- #Reference tracks (nref)
   MinEtaRef_(iConfig.getUntrackedParameter<double>("minEtaRef",-2.4)),
   MaxEtaRef_(iConfig.getUntrackedParameter<double>("maxEtaRef", 2.4)),
   MinPtRef_(iConfig.getUntrackedParameter<double>("minPtRef",0.4)),
   MaxPtRef_(iConfig.getUntrackedParameter<double>("maxPtRef", 9999.)),
   ChargeRef_(iConfig.getUntrackedParameter< std::vector<int> >("chargeRef")),
   isPixelTrackingRef_(iConfig.getUntrackedParameter<bool>("isPixTrkRef", 9999.)),
   dzdzErrorRef_(iConfig.getUntrackedParameter<double>("dzdzErrorRef", 3.)),
   d0d0ErrorRef_(iConfig.getUntrackedParameter<double>("d0d0ErrorRef", 3.)),
   PtErrorPtRef_(iConfig.getUntrackedParameter<double>("ptErrorPtRef", 1.)),
   Chi2nRef_(iConfig.getUntrackedParameter<double>("chi2nRef", 0.)),
   nHitsRef_(iConfig.getUntrackedParameter<int>("nHitsRef", 9999)),
   trkAlgoRef_(iConfig.getUntrackedParameter< std::vector<int> >("trkAlgoRef")),
// #Vertex
   vtxToken_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertex"))),
   MinVz_(iConfig.getUntrackedParameter<double>("minVz",-15.)),
   MaxVz_(iConfig.getUntrackedParameter<double>("maxVz", 15.)),
   MaxRho_(iConfig.getUntrackedParameter<double>("maxRho", 0.2)),
   nVtxMax_(iConfig.getUntrackedParameter<int>("nVtxMax", 9999)),
// #CaloTower
   caloTowerToken_(consumes<CaloTowerCollection>(iConfig.getParameter<edm::InputTag>("caloTower"))),
   MinEtHF_(iConfig.getUntrackedParameter<double>("minEt",0.)),
   MinEtaHF_(iConfig.getUntrackedParameter<double>("minEtaHF",2.9)),
   MaxEtaHF_(iConfig.getUntrackedParameter<double>("maxEtaHF",5.1)),
// #Centrality
   centralityToken_(consumes<reco::Centrality>(iConfig.getParameter<edm::InputTag>("centralitySrc"))),
   centralityBinToken_(consumes<int>(iConfig.getUntrackedParameter<edm::InputTag>("centralityBinSrc"))),
   MinCent_(iConfig.getUntrackedParameter<int>("minCent", -1)),
   MaxCent_(iConfig.getUntrackedParameter<int>("maxCent", -1)),
// #Event plane
   evtPlaneTag_(consumes<reco::EvtPlaneCollection>(iConfig.getParameter<edm::InputTag>("evtPlane"))),
   epLvl_(iConfig.getUntrackedParameter<int>("epLvl", 2)),
// #harmonics
   vHarmTrk_(iConfig.getUntrackedParameter< std::vector<int> >("vHarmTrk")),
   vHarmHF_(iConfig.getUntrackedParameter< std::vector<int> >("vHarmHF")),
// #cumulants
   cMode_(iConfig.getUntrackedParameter<int>("cMode")),
   cWeight_(iConfig.getUntrackedParameter<bool>("cWeight")),
// file acc & eff & fake
   fName_(iConfig.getUntrackedParameter<edm::InputTag>("fName")),
   effCentBin_(iConfig.getUntrackedParameter< std::vector<int> >("effCentBin")),
// Qvector eta binning
   EtaBinWidthTrk_(iConfig.getUntrackedParameter<double>("etaBinWidthTrk", 0.5)),
   EtaBinWidthHF_(iConfig.getUntrackedParameter<double>("etaBinWidthHF"  , 0.5))
{
   //now do what ever initialization is needed
   //Qn vectors trk
   nHarmTrk_ = vHarmTrk_.size();
   qNMvec_.resize(nHarmTrk_-1);
   for(unsigned int iharm = 0; iharm < nHarmTrk_; ++iharm)
   {
      qNvec_.push_back(correlations::QVector(0, 0, cWeight_));   

      for(unsigned int jharm = iharm+1; jharm < nHarmTrk_; ++jharm)
      {
         qNMvec_[iharm].push_back(correlations::QVector(0, 0, cWeight_)); 
      }
   }

   initDiagQ();
   qN2_ = new FlowCorrelator(nHarmTrk_);
   qN4_ = new FlowCorrelator(nHarmTrk_);
   qN6_ = new FlowCorrelator(nHarmTrk_);
   qN8_ = new FlowCorrelator(nHarmTrk_);

   initNonDiagQ();
   qNM4_.resize(nHarmTrk_-1);
   for(unsigned int iharm = 0; iharm < nHarmTrk_-1; ++iharm)
   {
      qNM4_[iharm] = new FlowCorrelator(nHarmTrk_-iharm-1);
   } 

   //file acc & eff
   TString filename(fName_.label().c_str());
   fEff_ = 0x0;

   if(cWeight_ && !filename.IsNull())
   {
      edm::FileInPath fip(Form("FlowCorr/QVectorTreeProducer/data/%s",filename.Data()));
      fEff_ = new TFile(fip.fullPath().c_str(),"READ");

      for(int ihist = 0; ihist < fEff_->GetNkeys(); ++ihist)
      {
         hEff_.push_back((TH2D*) fEff_->Get(fEff_->GetListOfKeys()->At(ihist)->GetName()));
      }
      edm::LogInfo("Input file for accXeff corrections") <<"Using file " << fEff_->GetName();
   }

   //qvector eta binning
   if(EtaBinWidthTrk_ <= 0) EtaBinWidthTrk_ = 0.1;
   nEtaBinTrk_ = ceil((MaxEtaRef_ - MinEtaRef_) / EtaBinWidthTrk_);
   for(unsigned int ieta = 0; ieta < nEtaBinTrk_; ++ieta)
   {
      vQn_trkP_.push_back(new FlowQVector(vHarmTrk_));
      vQn_trkN_.push_back(new FlowQVector(vHarmTrk_));
   }

   if(EtaBinWidthHF_ <= 0) EtaBinWidthTrk_ = 0.1;
   nEtaBinHF_ = ceil((MaxEtaHF_ - MinEtaHF_) / EtaBinWidthHF_);

   for(unsigned int ieta = 0; ieta < nEtaBinHF_; ++ieta)
   {
      vQn_hfP_.push_back(new FlowQVector(vHarmHF_));
      vQn_hfM_.push_back(new FlowQVector(vHarmHF_));
   }

   //Track selector
   sTrkCut cutsOff;
   cutsOff.dzvtxdzerror = dzdzErrorOff_; 
   cutsOff.d0vtxd0error = d0d0ErrorOff_; 
   cutsOff.etamin       = MinEtaOff_; 
   cutsOff.etamax       = MaxEtaOff_; 
   cutsOff.ptmin        = MinPtOff_; 
   cutsOff.ptmax        = MaxPtOff_; 
   cutsOff.pterrorpt    = PtErrorPtOff_; 
   cutsOff.chi2nnlayers = Chi2nOff_; 
   cutsOff.nhits        = nHitsOff_; 
   cutsOff.algo         = trkAlgoOff_; 
   cutsOff.charge       = ChargeOff_;
   trkSelectorOff_ = FlowTrackSelection(cutsOff); 
  
   sTrkCut cutsRef;
   cutsRef.dzvtxdzerror = dzdzErrorRef_; 
   cutsRef.d0vtxd0error = d0d0ErrorRef_; 
   cutsRef.etamin       = MinEtaRef_; 
   cutsRef.etamax       = MaxEtaRef_; 
   cutsRef.ptmin        = MinPtRef_; 
   cutsRef.ptmax        = MaxPtRef_; 
   cutsRef.pterrorpt    = PtErrorPtRef_; 
   cutsRef.chi2nnlayers = Chi2nRef_; 
   cutsRef.nhits        = nHitsRef_; 
   cutsRef.algo         = trkAlgoRef_; 
   cutsRef.charge       = ChargeRef_;
   trkSelectorRef_ = FlowTrackSelection(cutsRef);

   //Tower selector
   sTowCut cutsTow; 
   cutsTow.etamin = MinEtaHF_; 
   cutsTow.etamax = MaxEtaHF_; 
   cutsTow.etmin  = MinEtHF_;
   towSelector_ = FlowTowerSelection(cutsTow); 
  
   edm::Service<TFileService> fs;
   //Histo && TFileDirectory
   flowHistListRef_ = fs->mkdir("nRef_TrkHists");
   htrk_eta_ref_      = flowHistListRef_.make<TH1F>("htrk_eta_ref",      "", 60,   -3., 3.);
   htrk_phi_ref_      = flowHistListRef_.make<TH1F>("htrk_phi_ref",      "", 80,   -4., 4.);
   htrk_pt_ref_       = flowHistListRef_.make<TH1F>("htrk_pt_ref" ,      "", 1000,  0,  100.);
   htrk_eta_corr_ref_ = flowHistListRef_.make<TH1F>("htrk_eta_corr_ref", "", 60,   -3., 3.);
   htrk_phi_corr_ref_ = flowHistListRef_.make<TH1F>("htrk_phi_corr_ref", "", 80,   -4., 4.);
   htrk_pt_corr_ref_  = flowHistListRef_.make<TH1F>("htrk_pt_corr_ref" , "", 1000,  0,  100.);
   htrk_dzdzerr_ref_  = flowHistListRef_.make<TH1F>("htrk_dzdzerr_ref" , "", 200, -10., 10.);
   htrk_d0d0err_ref_  = flowHistListRef_.make<TH1F>("htrk_d0d0err_ref" , "", 200, -10., 10.);
   htrk_ptpterr_ref_  = flowHistListRef_.make<TH1F>("htrk_ptpterr_ref" , "",1200, -0.2, 1.);
   htrk_chi2nlayers_ref_ = flowHistListRef_.make<TH1F>("htrk_chi2nlayers_ref" , "", 200, 0., 1.);
   htrk_nhits_ref_    = flowHistListRef_.make<TH1I>("htrk_nhits_ref"   , "", 100,   0,  100);
   htrk_algo_ref_     = flowHistListRef_.make<TH1I>("htrk_algo_ref"    , "", 20,    0,  20);
   //--
   flowHistListOff_ = fs->mkdir("nOff_TrkHists");
   htrk_eta_off_      = flowHistListOff_.make<TH1F>("htrk_eta_off",      "", 60,   -3., 3.);
   htrk_phi_off_      = flowHistListOff_.make<TH1F>("htrk_phi_off",      "", 80,   -4., 4.);
   htrk_pt_off_       = flowHistListOff_.make<TH1F>("htrk_pt_off" ,      "", 1000,  0,  100.);
   htrk_eta_corr_off_ = flowHistListOff_.make<TH1F>("htrk_eta_corr_off", "", 60,   -3., 3.);
   htrk_phi_corr_off_ = flowHistListOff_.make<TH1F>("htrk_phi_corr_off", "", 80,   -4., 4.);
   htrk_pt_corr_off_  = flowHistListOff_.make<TH1F>("htrk_pt_corr_off" , "", 1000,  0,  100.);
   htrk_dzdzerr_off_  = flowHistListOff_.make<TH1F>("htrk_dzdzerr_off" , "", 200, -10., 10.);
   htrk_d0d0err_off_  = flowHistListOff_.make<TH1F>("htrk_d0d0err_off" , "", 200, -10., 10.);
   htrk_ptpterr_off_  = flowHistListOff_.make<TH1F>("htrk_ptpterr_off" , "",1200, -0.2, 1.);
   htrk_chi2nlayers_off_ = flowHistListOff_.make<TH1F>("htrk_chi2nlayers_off" , "", 200, 0., 1.);
   htrk_nhits_off_    = flowHistListOff_.make<TH1I>("htrk_nhits_off"   , "", 100,   0,  100);
   htrk_algo_off_     = flowHistListOff_.make<TH1I>("htrk_algo_off"    , "", 20,    0,  20);
   //--
   flowHistListHF_ = fs->mkdir("tow_HFHists");
   hHF_eta_tow_   = flowHistListHF_.make<TH1F>("hHF_eta_tow",      "", 120,  -6., 6.);
   hHF_et_tow_    = flowHistListHF_.make<TH1F>("hHF_et_tow" ,      "", 1000,  0,  100.);
   hHF_et_tow_p_  = flowHistListHF_.make<TH1F>("hHF_et_tow_plus",  "", 1000,  0,  100.);
   hHF_et_tow_m_  = flowHistListHF_.make<TH1F>("hHF_et_tow_minus", "", 1000,  0,  100.);
   hHF_phi_tow_   = flowHistListHF_.make<TH1F>("hHF_phi_tow",      "", 40,   -4., 4.);
   //TTree
   globalTree_ = fs->make<TTree>("globalTree", "globalTree");
   globalTree_->Branch("Event",    &Evt_,    "noff_corr/D:nref_corr:run/I:event:lumi:cent:noff:nref");
   globalTree_->Branch("Vertex",   &Vtx_,    "xVtx/D:yVtx:zVtx:xVtxError:yVtxError:zVtxError:nVtx/I");
   globalTree_->Branch("PsiMinus", &HFmPsi_, "psi1/D:psi2:psi3:psi4");
   globalTree_->Branch("PsiPlus",  &HFpPsi_, "psi1/D:psi2:psi3:psi4");
   globalTree_->Branch("Psi",      &HFPsi_,  "psi1/D:psi2:psi3:psi4");
   //--
   correlatorTree_ = fs->make<TTree>("correlatorTree", "correlatorTree");
   correlatorTree_->Branch("<2>_nn", "FlowCorrelator", &qN2_, 32000, 3);
   correlatorTree_->Branch("<4>_nn", "FlowCorrelator", &qN4_, 32000, 3);
   correlatorTree_->Branch("<6>_nn", "FlowCorrelator", &qN6_, 32000, 3);
   correlatorTree_->Branch("<8>_nn", "FlowCorrelator", &qN8_, 32000, 3);
   for(unsigned int iharm = 0; iharm < nHarmTrk_-1; ++iharm)
   {
      correlatorTree_->Branch(Form("<4>_%dm",vHarmTrk_[iharm]), "FlowCorrelator", &qNM4_[iharm], 32000, 3);
   }
   //--
   qvectorTreeTrk_ = fs->make<TTree>("qvectorTree_Trk","qvectorTree_Trk");
   for(unsigned int ieta = 0; ieta < nEtaBinTrk_; ++ieta)
   {
      if(find(ChargeRef_.begin(), ChargeRef_.end(),  1) != ChargeRef_.end())
         qvectorTreeTrk_->Branch(Form("Qtrk_p_%d",ieta), "FlowQVector", &vQn_trkP_[ieta], 32000, 3);
      if(find(ChargeRef_.begin(), ChargeRef_.end(), -1) != ChargeRef_.end()) 
         qvectorTreeTrk_->Branch(Form("Qtrk_n_%d",ieta), "FlowQVector", &vQn_trkN_[ieta], 32000, 3);
   }
   //--
   qvectorTreeHF_p_ = fs->make<TTree>("qvectorTree_HFp","qvectorTree_HFp");
   qvectorTreeHF_m_ = fs->make<TTree>("qvectorTree_HFm","qvectorTree_HFm");
   for(unsigned int ieta = 0; ieta < nEtaBinHF_; ++ieta)
   {
      qvectorTreeHF_p_->Branch(Form("QHF_p_%d",ieta), "FlowQVector", &vQn_hfP_[ieta], 32000, 3);
      qvectorTreeHF_m_->Branch(Form("QHF_m_%d",ieta), "FlowQVector", &vQn_hfM_[ieta], 32000, 3);
   }
}

QVectorTreeProducer::~QVectorTreeProducer()
{
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
   if(fEff_) fEff_->Close();
   hEff_.clear();

   qNvec_.clear();
   hcNvec_.clear();
   for(unsigned icq = 0; icq < cqNvec_.size(); ++icq)
   {
      delete cqNvec_[icq];
   }
   cqNvec_.clear();

   delete qN2_;
   delete qN4_;
   delete qN6_;
   delete qN8_;


   for(unsigned int iharm = 0; iharm < nHarmTrk_-1; ++iharm)
   {
      qNMvec_[iharm].clear();
      hcNMvec_[iharm].clear();
      for(unsigned int icq = 0; icq < cqNMvec_[iharm].size(); ++icq)
      {
         delete cqNMvec_[iharm][icq];
      }
      cqNMvec_[iharm].clear();
   } 
   qNMvec_.clear();
   hcNMvec_.clear();
   cqNMvec_.clear();

   for(unsigned int iqnm = 0; iqnm < qNM4_.size(); ++iqnm)
   {
      delete qNM4_[iqnm];
   }
   qNM4_.clear();


   for(unsigned int ieta = 0; ieta < nEtaBinTrk_; ++ieta)
   {
      delete vQn_trkP_[ieta];
      delete vQn_trkN_[ieta];
   }
   for(unsigned int ieta = 0; ieta < nEtaBinHF_; ++ieta)
   {
      delete vQn_hfP_[ieta];
      delete vQn_hfM_[ieta];
   }
   vQn_trkP_.clear(); 
   vQn_trkN_.clear();
   vQn_hfP_.clear();
   vQn_hfM_.clear();
}


//
// member functions
//

// ------------ method called for each event  ------------
void
QVectorTreeProducer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   //----- General informations -----
   Evt_.run   = iEvent.id().run();
   Evt_.event = iEvent.id().event();
   Evt_.lumi  = iEvent.luminosityBlock();
   Evt_.cent  = -1;
   Evt_.noff  = -1;
   Evt_.nref  = -1;
   Evt_.noffcorr = -1;
   Evt_.noffcorr = -1;

 
   //----- Vertex selection -----
   edm::Handle< reco::VertexCollection > vertices;
   iEvent.getByToken(vtxToken_, vertices);
   if(!vertices->size())
   {
      edm::LogWarning ("Missing Collection") <<"Invalid or empty vertex collection!";
      return;
   }

   //----- Tracks selection -----
   edm::Handle< reco::TrackCollection > tracks;
   
   iEvent.getByToken(trackToken_, tracks);
   if( !tracks->size() )
   {
       edm::LogWarning ("Missing Collection") <<"Invalid or empty track collection!";
       return;
   }
   
   //----- Calotower selection -----
   edm::Handle< CaloTowerCollection > calotowers;
   iEvent.getByToken(caloTowerToken_, calotowers);
   if(!calotowers->size()) 
   { 
       edm::LogWarning ("Missing Collection") <<"Invalid or empty caloTower collection!";
       return; 
   }

   //----- Centrality selection ----
   Cent_ = -1;
   edm::Handle< reco::Centrality > centrality;
   iEvent.getByToken(centralityToken_, centrality);

   edm::Handle< int > cbin;
   iEvent.getByToken(centralityBinToken_,cbin);
   Cent_ = *cbin;
   if(Cent_ < 0) 
   { 
       edm::LogWarning ("Invalid value") <<"Invalid centrality value";
   }
   Evt_.cent = Cent_;
   
   //----- Event Plane selection ----
   edm::Handle<reco::EvtPlaneCollection> evtPlanes;
   iEvent.getByToken(evtPlaneTag_, evtPlanes);
   fillEPangle(evtPlanes);

   //----- Event selection -----
   if(!isEventSelected(iEvent, vertices, tracks)) return;
  
   nRef_     = 0;
   nRefcorr_ = 0.;

   //Loop on tracks
   for(reco::TrackCollection::const_iterator itTrk = tracks->begin(); 
       itTrk != tracks->end(); 
       ++itTrk)
   {
      trkSelectorRef_.resetVar(); 
      trkSelectorRef_.fillVar(*itTrk, Vtx_, trackQualityTag_);
      if(!trkSelectorRef_.isTrkPassCuts(isPixelTrackingRef_)) continue; 

      //Weight calculation from histograms
      double weight = 1.;
      if(cWeight_)
      {
         weight = getAccEffWeight(trkSelectorRef_.getTrk().eta, trkSelectorRef_.getTrk().pt);
      }

      //Mult calculation
      nRef_++;
      nRefcorr_ += weight;

      //Fill histos
      fillHistoTrk(weight);

      //Filling Qvectors for correlators
      fillQVectorCorr(weight);

      //Filling trk QVectors for the QVector tree
      fillQVectorTrk(weight);
   }

   //Loop on calo towers
   for (CaloTowerCollection::const_iterator itHF = calotowers->begin();
        itHF != calotowers->end();
        ++itHF)
   {
      towSelector_.resetVar(); 
      towSelector_.fillVar(*itHF);
      if(!towSelector_.isTowerPassCuts()) continue; 

      //Fill histos
      fillHistoHF();

      //Filling HF QVectors for the QVector tree
      fillQVectorHF(towSelector_.getTower().weight, 
                    towSelector_.getTower().eta, 
                    towSelector_.getTower().phi);
   }

   //qNvec_[0].print();
   //std::cout << "qvector from me" << std::endl;
   //std::cout << vQn_trkP_[0]->corr[0]+vQn_trkN_[0]->corr[0] << std::endl; 
   //std::cout << vQn_trkP_[0]->corr[1]+vQn_trkN_[0]->corr[1] << std::endl; 

   if(nOff_ != 0) Evt_.noff = nOff_; 
   if(nRef_ != 0) Evt_.nref = nRef_; 

   if(nOffcorr_ != 0) Evt_.noffcorr = nOffcorr_; 
   if(nRefcorr_ != 0) Evt_.nrefcorr = nRefcorr_; 

   //Calculate correlators
   fillCorrelators();

   //Fill TTrees
   globalTree_->Fill();
   correlatorTree_->Fill();
   qvectorTreeTrk_->Fill();
   qvectorTreeHF_p_->Fill();
   qvectorTreeHF_m_->Fill();

   //Reset correlators and their Q vectors
   doneDiagQ();
   doneNonDiagQ();
   resetQvectors();
}


// ------------ method called once each job just before starting event loop  ------------
void 
QVectorTreeProducer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
QVectorTreeProducer::endJob() 
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
QVectorTreeProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) 
{
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

// ------------ Other methods ------------
bool
QVectorTreeProducer::isEventSelected(const edm::Event& iEvent, 
                                     edm::Handle< reco::VertexCollection > vertices,
                                     edm::Handle< reco::TrackCollection >  tracks) 
{
   // ## Vertex ##
   //--- Vertex info 
   const reco::Vertex & bestVtx = *(vertices->begin());

   Vtx_.xVtx = 0.;
   Vtx_.yVtx = 0.;
   Vtx_.zVtx = 0.;
   Vtx_.nVtx = 0;
    
   if(!bestVtx.isFake() && bestVtx.tracksSize()>=2)
   {
      Vtx_.xVtx = bestVtx.x();
      Vtx_.yVtx = bestVtx.y();
      Vtx_.zVtx = bestVtx.z();
      Vtx_.xVtxError = bestVtx.xError();
      Vtx_.yVtxError = bestVtx.yError();
      Vtx_.zVtxError = bestVtx.zError();
   }
   double rho = sqrt(Vtx_.xVtx*Vtx_.xVtx + Vtx_.yVtx*Vtx_.yVtx);
 
   for(reco::VertexCollection::const_iterator itVtx = vertices->begin(); itVtx != vertices->end(); ++itVtx)
   {
       if(!itVtx->isFake() && itVtx->tracksSize() >= 2) Vtx_.nVtx++;
   }
 
   //--- Vertex selection 
   if(Vtx_.nVtx > nVtxMax_)                     return false;
   if(Vtx_.zVtx < MinVz_ || Vtx_.zVtx > MaxVz_) return false;
   if(rho > MaxRho_ )                           return false;
 
 
   // ## Centrality selection ##
   if((Cent_ < 2*MinCent_ || Cent_ >= 2*MaxCent_) &&
      (MaxCent_ != -1 && MaxCent_ != -1)) 
                                        return false;
 
 
   // ## Multiplicity ##
   getNoff(tracks); 
   if(nOff_ < MinNoff_ || nOff_ > MaxNoff_) return false;
 
 
   return true;
}

void 
QVectorTreeProducer::getNoff(edm::Handle< reco::TrackCollection > tracks)
{
   nOff_ = 0;
   nOffcorr_ = 0.;
 
   for(reco::TrackCollection::const_iterator itTrk = tracks->begin(); itTrk != tracks->end(); ++itTrk)
   {
      trkSelectorOff_.resetVar(); 
      trkSelectorOff_.fillVar(*itTrk, Vtx_, "highPurity");

      if(!trkSelectorOff_.isTrkPassCuts(isPixelTrackingOff_)) continue; 
 
      //Weight calculation from histograms
      double weight = 1.;
      if(cWeight_)
      {
         weight = getAccEffWeight(trkSelectorOff_.getTrk().eta, trkSelectorOff_.getTrk().pt);
      }

      //Ntrk offline calculation
      nOff_++;
      nOffcorr_ += weight; 
   }
}

void
QVectorTreeProducer::initDiagQ()
{
   cqNvec_.resize(nHarmTrk_);
   for(unsigned int iharm = 0; iharm < nHarmTrk_; ++iharm)
   {
       hcNvec_.push_back(correlations::HarmonicVector(8));
       for(int ipart = 0; ipart < 8; ++ipart)
       {
           if(ipart%2 == 0) hcNvec_[iharm][ipart] =    vHarmTrk_[iharm];
           else             hcNvec_[iharm][ipart] = -1*vHarmTrk_[iharm];
       }

       qNvec_[iharm].resize(hcNvec_[iharm]);
       switch ( cMode_ ) {
          case 1:
             cqNvec_[iharm] = new correlations::recurrence::FromQVector(qNvec_[iharm]);
             break;
          default:
      	     cqNvec_[iharm] = new correlations::recursive::FromQVector(qNvec_[iharm]);
      	     break;
       }
   }
}

void
QVectorTreeProducer::initNonDiagQ()
{
   cqNMvec_.resize(nHarmTrk_-1);
   hcNMvec_.resize(nHarmTrk_-1);
   for(unsigned int iharm = 0; iharm < nHarmTrk_-1; ++iharm)
   {
       for(unsigned int jharm = iharm+1; jharm < nHarmTrk_; ++jharm)
       {
          hcNMvec_[iharm].push_back(correlations::HarmonicVector(4));
          hcNMvec_[iharm][jharm-iharm-1][0] =    vHarmTrk_[iharm];
          hcNMvec_[iharm][jharm-iharm-1][1] = -1*vHarmTrk_[iharm];
          hcNMvec_[iharm][jharm-iharm-1][2] =    vHarmTrk_[jharm];
          hcNMvec_[iharm][jharm-iharm-1][3] = -1*vHarmTrk_[jharm];

          qNMvec_[iharm][jharm-iharm-1].resize(hcNMvec_[iharm][jharm-iharm-1]);
          switch ( cMode_ ) {
             case 1:
                cqNMvec_[iharm].push_back(new correlations::recurrence::FromQVector(qNMvec_[iharm][jharm-iharm-1]));
                break;
             default:
                cqNMvec_[iharm].push_back(new correlations::recursive::FromQVector(qNMvec_[iharm][jharm-iharm-1]));
                break;
          }
       }
   }

}

void
QVectorTreeProducer::doneDiagQ()
{
   for(unsigned int iharm = 0; iharm < nHarmTrk_; ++iharm)
   {
      qNvec_[iharm].reset();
   }
 
//   qNvec_.clear();
}

void
QVectorTreeProducer::doneNonDiagQ()
{
   for(unsigned int iharm = 0; iharm < nHarmTrk_-1; ++iharm)
   {
      for(unsigned int jharm = iharm+1; jharm < nHarmTrk_; ++jharm)
      {
         qNMvec_[iharm][jharm-iharm-1].reset();
      }
//      qNMvec_[iharm].clear();
   }
 
//   qNMvec_.clear();
}

void
QVectorTreeProducer::resetQvectors()
{
   for(unsigned int ieta = 0; ieta < nEtaBinTrk_; ++ieta)
   {
      vQn_trkP_[ieta]->reset();
      vQn_trkN_[ieta]->reset();
   }
   for(unsigned int ieta = 0; ieta < nEtaBinHF_; ++ieta)
   {
      vQn_hfP_[ieta]->reset();
      vQn_hfM_[ieta]->reset();
   } 
}

double
QVectorTreeProducer::getAccEffWeight(double eta, double pt)
{
   double weight = 1.;
   if(hEff_.size() != effCentBin_.size() - 1 || effCentBin_.empty())
   {
      if(hEff_[0]->GetBinContent(hEff_[0]->FindBin(eta,pt)) != 0.) 
         weight = 1./hEff_[0]->GetBinContent(hEff_[0]->FindBin(eta,pt));
      else weight = 0.;
   }
   else
   {
      //weight = 0.;
      for(unsigned int icent = 0; icent < effCentBin_.size() - 1; ++icent)
      {
         if(Cent_ > 2*effCentBin_[icent] && Cent_ <= 2*effCentBin_[icent+1])
         { 
            if(hEff_[icent]->GetBinContent(hEff_[icent]->FindBin(eta,pt)) != 0.) 
               weight = 1./hEff_[icent]->GetBinContent(hEff_[icent]->FindBin(eta,pt));
            else weight = 0.;
            break;
         }
      }
   }
   return weight;
}

void
QVectorTreeProducer::fillHistoTrk(double weight)
{
   htrk_eta_ref_->Fill(trkSelectorRef_.getTrk().eta);
   htrk_phi_ref_->Fill(trkSelectorRef_.getTrk().phi);
   htrk_pt_ref_ ->Fill(trkSelectorRef_.getTrk().pt);

   htrk_eta_corr_ref_->Fill(trkSelectorRef_.getTrk().eta, weight);
   htrk_phi_corr_ref_->Fill(trkSelectorRef_.getTrk().phi, weight);
   htrk_pt_corr_ref_ ->Fill(trkSelectorRef_.getTrk().pt,  weight);

   if(trkSelectorRef_.getTrk().dzerror != 0.) 
      htrk_dzdzerr_ref_    ->Fill(trkSelectorRef_.getTrk().dzvtx/trkSelectorRef_.getTrk().dzerror);
   if(trkSelectorRef_.getTrk().d0error != 0.) 
      htrk_d0d0err_ref_    ->Fill(trkSelectorRef_.getTrk().d0vtx/trkSelectorRef_.getTrk().d0error);
   if(trkSelectorRef_.getTrk().pt      != 0.) 
      htrk_ptpterr_ref_    ->Fill(trkSelectorRef_.getTrk().pterror/trkSelectorRef_.getTrk().pt);
   if(trkSelectorRef_.getTrk().nlayers != 0.) 
      htrk_chi2nlayers_ref_->Fill(trkSelectorRef_.getTrk().chi2n/trkSelectorRef_.getTrk().nlayers);
   htrk_nhits_ref_->Fill(trkSelectorRef_.getTrk().nhits);
   htrk_algo_ref_ ->Fill(trkSelectorRef_.getTrk().algo);
   //--
   htrk_eta_off_->Fill(trkSelectorOff_.getTrk().eta);
   htrk_phi_off_->Fill(trkSelectorOff_.getTrk().phi);
   htrk_pt_off_ ->Fill(trkSelectorOff_.getTrk().pt);

   htrk_eta_corr_off_->Fill(trkSelectorOff_.getTrk().eta, weight);
   htrk_phi_corr_off_->Fill(trkSelectorOff_.getTrk().phi, weight);
   htrk_pt_corr_off_ ->Fill(trkSelectorOff_.getTrk().pt,  weight);

   if(trkSelectorOff_.getTrk().dzerror != 0.) 
      htrk_dzdzerr_off_    ->Fill(trkSelectorOff_.getTrk().dzvtx/trkSelectorOff_.getTrk().dzerror);
   if(trkSelectorOff_.getTrk().d0error != 0.) 
      htrk_d0d0err_off_    ->Fill(trkSelectorOff_.getTrk().d0vtx/trkSelectorOff_.getTrk().d0error);
   if(trkSelectorOff_.getTrk().pt      != 0.) 
      htrk_ptpterr_off_    ->Fill(trkSelectorOff_.getTrk().pterror/trkSelectorOff_.getTrk().pt);
   if(trkSelectorOff_.getTrk().nlayers != 0.) 
      htrk_chi2nlayers_off_->Fill(trkSelectorOff_.getTrk().chi2n/trkSelectorOff_.getTrk().nlayers);
   htrk_nhits_off_->Fill(trkSelectorOff_.getTrk().nhits);
   htrk_algo_off_ ->Fill(trkSelectorOff_.getTrk().algo);
}

void
QVectorTreeProducer::fillHistoHF()
{
   hHF_eta_tow_ ->Fill(towSelector_.getTower().eta);
   hHF_et_tow_  ->Fill(towSelector_.getTower().et);
   if(towSelector_.getTower().eta > 0.) 
      hHF_et_tow_p_->Fill(towSelector_.getTower().et);
   if(towSelector_.getTower().eta < 0.) 
      hHF_et_tow_m_->Fill(towSelector_.getTower().et);
   hHF_phi_tow_ ->Fill(towSelector_.getTower().phi);
}

void
QVectorTreeProducer::fillQVectorCorr(double weight)
{
      for(unsigned int iharm = 0; iharm < nHarmTrk_; ++iharm)
      {
      //-- Diagonal
         qNvec_[iharm].fill(trkSelectorRef_.getTrk().phi, weight);
      //-- Non diagonal
         for(unsigned int jharm = iharm+1; jharm < nHarmTrk_; ++jharm)
         {
             qNMvec_[iharm][jharm-iharm-1].fill(trkSelectorRef_.getTrk().phi, weight);
         }
      }
}

void
QVectorTreeProducer::fillQVectorTrk(double weight)
{
      int etaIdx = -1;
      for(unsigned int ieta = 0; ieta < nEtaBinTrk_; ++ieta)
      {
         if(trkSelectorRef_.getTrk().eta >= (MinEtaRef_ + ieta    *EtaBinWidthTrk_) && 
            trkSelectorRef_.getTrk().eta <  (MinEtaRef_ + (ieta+1)*EtaBinWidthTrk_) &&
            ieta != nEtaBinTrk_-1)
         {
              etaIdx = ieta;
         }
         else if(trkSelectorRef_.getTrk().eta >= (MinEtaRef_ + ieta    *EtaBinWidthTrk_) && 
                 trkSelectorRef_.getTrk().eta <= (MinEtaRef_ + (ieta+1)*EtaBinWidthTrk_) &&
                 ieta == nEtaBinTrk_-1)
         {
              etaIdx = ieta;
         }
      }
      if(etaIdx < 0)
      {
         edm::LogInfo("Track out of Qvector eta range") << "The track with eta = " 
                                                        << trkSelectorRef_.getTrk().eta 
                                                        << " is not filled in the Qvectors";
         return;
      }
      else
      {
         if(trkSelectorRef_.getTrk().charge > 0) vQn_trkP_[etaIdx]->fill(trkSelectorRef_.getTrk().phi, weight);
         if(trkSelectorRef_.getTrk().charge < 0) vQn_trkN_[etaIdx]->fill(trkSelectorRef_.getTrk().phi, weight);
      }
}

void
QVectorTreeProducer::fillQVectorHF(double weight, double eta, double phi)
{
      int etaIdx = -1;
      for(unsigned int ieta = 0; ieta < nEtaBinHF_; ++ieta)
      {
         if(fabs(eta) >= (MinEtaHF_ + ieta    *EtaBinWidthHF_) && 
            fabs(eta) <  (MinEtaHF_ + (ieta+1)*EtaBinWidthHF_) &&
            ieta != nEtaBinHF_-1)
         {
              etaIdx = ieta;
         }
         else if(fabs(eta) >= (MinEtaHF_ + ieta    *EtaBinWidthHF_) && 
                 fabs(eta) <= (MinEtaHF_ + (ieta+1)*EtaBinWidthHF_) &&
                 ieta == nEtaBinHF_-1)
         {
              etaIdx = ieta;
         }
      }
      if(etaIdx < 0)
      {
         edm::LogInfo("Tower out of Qvector eta range") << "The tower with eta = " 
                                                        << eta 
                                                        << " is not filled in the Qvectors";
         return;
      }
      else
      {
         if(eta > 0) vQn_hfP_[etaIdx]->fill(phi, weight);
         if(eta < 0) vQn_hfM_[etaIdx]->fill(phi, weight);
      }
}

void
QVectorTreeProducer::fillCorrelators()
{
   //-- Diagonal
   correlations::Result rN2; 
   correlations::Result rN4;
   correlations::Result rN6;
   correlations::Result rN8;
   //-- Non diagonal
   correlations::Result rNM4;
   if(nHarmTrk_ != 0)
   {
      for(unsigned int iharm = 0; iharm < nHarmTrk_; ++iharm)
      {
         //-- Diagonal
         rN2 = cqNvec_[iharm]->calculate(2, hcNvec_[iharm]);
         rN4 = cqNvec_[iharm]->calculate(4, hcNvec_[iharm]);
         rN6 = cqNvec_[iharm]->calculate(6, hcNvec_[iharm]);
         rN8 = cqNvec_[iharm]->calculate(8, hcNvec_[iharm]);
         
         qN2_->corr[iharm] = TComplex(rN2._sum.real(), rN2._sum.imag());
         qN4_->corr[iharm] = TComplex(rN4._sum.real(), rN4._sum.imag());
         qN6_->corr[iharm] = TComplex(rN6._sum.real(), rN6._sum.imag());
         qN8_->corr[iharm] = TComplex(rN8._sum.real(), rN8._sum.imag());


         //-- Non diagonal
         for(unsigned int jharm = iharm+1; jharm < nHarmTrk_; ++jharm)
         {
            rNM4 = cqNMvec_[iharm][jharm-iharm-1]->calculate(4, hcNMvec_[iharm][jharm-iharm-1]);

            qNM4_[iharm]->corr[jharm-iharm-1] = TComplex(rNM4._sum.real(), rNM4._sum.imag());
         }

         if(iharm < nHarmTrk_ -1)
         {
            qNM4_[iharm]->weight = rNM4._weights;
         }
      }
      qN2_->weight = rN2._weights; 
      qN4_->weight = rN4._weights; 
      qN6_->weight = rN6._weights; 
      qN8_->weight = rN8._weights;
      //std::cout << "<2> = " << TComplex(rN2._sum.real(), rN2._sum.imag()) << std::endl;
      //std::cout << "<2> weight = " << qN2_->weight << std::endl;  
   }
}

void 
QVectorTreeProducer::fillEPangle(const edm::Handle<reco::EvtPlaneCollection> & ep)
{
   //----- Event Plane selection ----
   //Index     Name   Detector Order hmin1 hmax1 hmin2 hmax2 minpt maxpt
   //    0      HFm1        HF     1 -5.00 -3.00  0.00  0.00  0.01 30.00
   //    1      HFp1        HF     1  3.00  5.00  0.00  0.00  0.01 30.00
   //    2       HF1        HF     1 -5.00 -3.00  3.00  5.00  0.01 30.00
   //    6      HFm2        HF     2 -5.00 -3.00  0.00  0.00  0.01 30.00
   //    7      HFp2        HF     2  3.00  5.00  0.00  0.00  0.01 30.00
   //    8       HF2        HF     2 -5.00 -3.00  3.00  5.00  0.01 30.00
   //   13      HFm3        HF     3 -5.00 -3.00  0.00  0.00  0.01 30.00
   //   14      HFp3        HF     3  3.00  5.00  0.00  0.00  0.01 30.00
   //   15       HF3        HF     3 -5.00 -3.00  3.00  5.00  0.01 30.00
   //   19      HFm4        HF     4 -5.00 -3.00  0.00  0.00  0.01 30.00
   //   20      HFp4        HF     4  3.00  5.00  0.00  0.00  0.01 30.00
   //   21       HF4        HF     4 -5.00 -3.00  3.00  5.00  0.01 30.00
   //   25    HFm1mc        HF     1 -5.00 -3.00  0.00  0.00  0.01 30.00
   //   26    HFp1mc        HF     1  3.00  5.00  0.00  0.00  0.01 30.00
   unsigned int HFm[4] = {0, 6, 13, 19};
   unsigned int HFp[4] = {1, 7, 14, 20};
   unsigned int HF[4]  = {2, 8, 15, 21};

   if(ep.isValid())
   {
      //psi1
      if(HFm[0] < ep->size()) HFmPsi_.psi1 = (*ep)[HFm[0]].angle(epLvl_);
      if(HFp[0] < ep->size()) HFpPsi_.psi1 = (*ep)[HFp[0]].angle(epLvl_);
      if(HF[0]  < ep->size()) HFPsi_.psi1  = (*ep)[HF[0]].angle(epLvl_);

      //psi2
      if(HFm[1] < ep->size()) HFmPsi_.psi2 = (*ep)[HFm[1]].angle(epLvl_);
      if(HFp[1] < ep->size()) HFpPsi_.psi2 = (*ep)[HFp[1]].angle(epLvl_);
      if(HF[1]  < ep->size()) HFPsi_.psi2  = (*ep)[HF[1]].angle(epLvl_);

      //psi3
      if(HFm[2] < ep->size()) HFmPsi_.psi3 = (*ep)[HFm[2]].angle(epLvl_);
      if(HFp[2] < ep->size()) HFpPsi_.psi3 = (*ep)[HFp[2]].angle(epLvl_);
      if(HF[2]  < ep->size()) HFPsi_.psi3  = (*ep)[HF[2]].angle(epLvl_);

      //psi4
      if(HFm[3] < ep->size()) HFmPsi_.psi4 = (*ep)[HFm[3]].angle(epLvl_);
      if(HFp[3] < ep->size()) HFpPsi_.psi4 = (*ep)[HFp[3]].angle(epLvl_);
      if(HF[3]  < ep->size()) HFPsi_.psi4  = (*ep)[HF[3]].angle(epLvl_);
   }
   else
   {
       edm::LogWarning ("Invalid value") << "Invalid EP value";
   }

}


//define this as a plug-in
DEFINE_FWK_MODULE(QVectorTreeProducer);
