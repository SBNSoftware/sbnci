////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Class:       PDSValidation                                                                             //
// Modlue Type: Analyser                                                                                  //
// File         PDSValidation_module.cc                                                                //
// Author:      Chris Hilgenberg chilgenb@fnal.gov                                                        //
//                                                                                                        //
// Usage:       Simple workflow to validate PDS reconstruction                                            //
//09/10/2021 Adopted for  (Sergey Martynenko smartynen@bnl.gov)                                         //
////////////////////////////////////////////////////////////////////////////////////////////////////////////


// Framework includes
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/View.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/PtrVector.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileService.h"
#include "art_root_io/TFileDirectory.h"
#include "canvas/Persistency/Common/FindMany.h"
#include "canvas/Persistency/Common/FindOneP.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

//Larsoft includes
#include "larcore/Geometry/Geometry.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "lardataobj/RecoBase/OpFlash.h"
#include "lardataobj/RecoBase/OpHit.h"

//Root Includes
#include "TMath.h"
#include "TTree.h"
#include "TSystem.h"
#include "TH1F.h"

//C++ Includes
#include <vector>
#include <iostream>
#include <fstream>

using namespace std;

namespace ana {
  class PDSValidation;
}

class ana::PDSValidation : public art::EDAnalyzer {
  public:

    explicit PDSValidation(const fhicl::ParameterSet& pset);
    // The compiler-generated destructor is fine for non-base
    // classes without bare pointers or other resource use.

    // Plugins should not be copied or assigned.
    PDSValidation(PDSValidation const&) = delete;
    PDSValidation(PDSValidation&&) = delete;
    PDSValidation& operator=(PDSValidation const&) = delete;
    PDSValidation& operator=(PDSValidation&&) = delete;

    // required abstract methods we need to impliment 
    void analyze(const art::Event& evt) override;
    void endJob() override;
    void beginJob() override;

  private:

    //fcl parameters
    string fGenieGenModuleLabel;
    string fLArGeantModuleLabel;
    string fOpFlashTPC1Label;
    string fOpHitPMTLabel;
    string fOpHitArapucaLabel;
    bool   fVerbose;

    //TTree
    TTree* fTree;
//    TH1F* fTimeHist;

    //Service handlesacktracker
    art::ServiceHandle<cheat::BackTrackerService> backtracker;
    art::ServiceHandle<cheat::ParticleInventoryService> particleInventory;
    art::ServiceHandle<geo::Geometry> geom;
    art::ServiceHandle<art::TFileService> tfs;

    int numevents;

    // tree variables
    int fRun;    ///< art run number
    int fSubRun; ///< art subrun number
    int fEvent;  ///< art event number

    int fNFlash; ///< number of OpFlashes in the event
    vector<double> fFlashTotalPE; ///< total PE associated with single OpFlash
    vector<double> fFlashTime; ///< time of flash
    vector<double> fFlashTimeWidth;
    vector<double> fYCenter;
    vector<double> fZCenter;
    vector<double> fYWidth;
    vector<double> fZWidth;

    int fPMTNHit;
    vector<int> fPMTOpChannel;
    vector<double> fPMTPeakTimeAbs;
    vector<double> fPMTPeakTime;
    vector<double> fPMTWidth;
    vector<double> fPMTArea;
    vector<double> fPMTAmplitude;
    vector<double> fPMTPE;
    vector<double> fPMTFastToTotal;

    int fArapucaNHit;
    vector<int> fArapucaOpChannel;
    vector<double> fArapucaPeakTimeAbs;
    vector<double> fArapucaPeakTime;
    vector<double> fArapucaWidth;
    vector<double> fArapucaArea;
    vector<double> fArapucaAmplitude;
    vector<double> fArapucaPE;
    vector<double> fArapucaFastToTotal;

}; //end def class PDSValidation

//////////////////////////////////////////////////////////////////////////////////////
ana::PDSValidation::PDSValidation(const fhicl::ParameterSet& pset) : 
  EDAnalyzer(pset),
  numevents(0)
{

  fGenieGenModuleLabel = pset.get<string>("GenieGenModuleLabel","generator");
  fLArGeantModuleLabel = pset.get<string>("LArGeantModuleLabel","largeant");
  fOpFlashTPC1Label    = pset.get<string>("OpFlashTPC1Label","opflashTPC1");
  fOpHitPMTLabel       = pset.get<string>("OpHitPMTLabel", "ophit");
  fVerbose             = pset.get<bool>("Verbose",false); 

} //end constructor

///////////////////////////////////////////////////////////////////////////////
void ana::PDSValidation::beginJob() {
  
  fTree = tfs->make<TTree>("pdsTree", "Tree with PDS validation information");
  //gInterpreter->GenerateDictionary("vector<vector<float> > ","vector");

  fTree->Branch("Run",          &fRun,          "Run/I");
  fTree->Branch("SubRun",       &fSubRun,       "SubRun/I");
  fTree->Branch("Event",        &fEvent,        "Event/I");
  fTree->Branch("NFlash",       &fNFlash,       "NFlash/I"); 
  fTree->Branch("PMTNHit",         &fPMTNHit,         "PMTNHit/I");

  fTree->Branch("TotalFlashPE", &fFlashTotalPE);
  fTree->Branch("FlashTime", &fFlashTime);
  fTree->Branch("FlashTimeWidth", &fFlashTimeWidth);
  fTree->Branch("YCenter", &fYCenter);
  fTree->Branch("ZCenter", &fZCenter);
  fTree->Branch("YWidth", &fYWidth);
  fTree->Branch("ZWidth", &fZWidth);
 
  fTree->Branch("PMTOpChannel",&fPMTOpChannel);
  fTree->Branch("PMTPeakTime",&fPMTPeakTime);
  fTree->Branch("PMTPeakTimeAbs",&fPMTPeakTimeAbs);
  fTree->Branch("PMTWidth",&fPMTWidth);
  fTree->Branch("PMTArea",&fPMTArea);
  fTree->Branch("PMTAmplitude",&fPMTAmplitude);
  fTree->Branch("PMTPE",&fPMTPE);
  fTree->Branch("PMTFastToTotal",&fPMTFastToTotal);



//  fTimeHist=tfs->make<TH1F>("fTimeHist", ";FlashTime (ns)", 20, 0, 100);
}// end beginJob

////////////////////////////////////////////////////////////////////
void ana::PDSValidation::analyze(const art::Event& evt) {

  fRun    = evt.run();
  fSubRun = evt.subRun();
  fEvent  = evt.event();

  if(fVerbose) {
      cout << "Analysing run " << fRun << ", subrun " << fSubRun << ", event: " << fEvent << endl;
  }
  numevents++;

  //Getting  MC truth information
  art::Handle< vector<simb::MCTruth> > mctruthListHandle;
  vector<art::Ptr<simb::MCTruth> > mclist;
  if(evt.getByLabel(fGenieGenModuleLabel,mctruthListHandle)){
      art::fill_ptr_vector(mclist, mctruthListHandle);
  }

  //###############################################
  //### Get the Truth information for the event ###
  //###############################################

  //List the particles in the event
  const sim::ParticleList& particles = particleInventory->ParticleList();

  //Loop over the particles
  map<int,const simb::MCParticle*> trueParticles;
  map<int,bool> mcparticlescontained;
  map<int,float> trueParticleEnergy;

  //auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(evt);

  //Make a map of Track id and pdgcode
  for (const auto& particleIt: particles) {

      const simb::MCParticle* particle = particleIt.second;
      trueParticles[particle->TrackId()] = particle;

      if(fVerbose){
          cout << "True Particle with track ID: " << particle->TrackId() << " Has code of: " 
               << particle->PdgCode() << " and Energy of: " << particle->E() << " With Mother: " 
               << particle->Mother() << " Proccess: " << particle->Process() << " End Process: "  
               << particle->EndProcess() 
          << endl;
      }// if verbose

  } // for MCParticles


  //###############################################
  //### Get PDS reco info for the event         ###
  //###############################################
  //Getting  OpFlash information
  art::Handle< vector<recob::OpFlash> > opFlashListHandle;
  vector<art::Ptr<recob::OpFlash> > opFlashList;
  if(evt.getByLabel(fOpFlashTPC1Label, opFlashListHandle)){
      art::fill_ptr_vector(opFlashList, opFlashListHandle);
  }

  fNFlash = opFlashList.size();
std::cout<<"Number of opFlash = "<<fNFlash<<std::endl;
  for(auto const& flash : opFlashList) {
    fFlashTotalPE.push_back(flash->TotalPE());
    fFlashTime.push_back(flash->Time());
    fFlashTimeWidth.push_back(flash->TimeWidth());
    fYCenter.push_back(flash->YCenter());
    fZCenter.push_back(flash->ZCenter());
    fYWidth.push_back(flash->YWidth());
    fZWidth.push_back(flash->ZWidth());
//    fTimeHist->Fill(flash->Time());
  }

  art::Handle< vector<recob::OpHit> > opHitPMTListHandle;
  vector<art::Ptr<recob::OpHit> > opHitPMTList;
  if(evt.getByLabel(fOpHitPMTLabel, opHitPMTListHandle)){
      art::fill_ptr_vector(opHitPMTList, opHitPMTListHandle);
  }

  fPMTNHit =opHitPMTList.size();
std::cout<<"Numerb of opHits="<<fPMTNHit<<std::endl;
  for(auto const& hit : opHitPMTList) {
    fPMTOpChannel.push_back(hit->OpChannel());
    fPMTPeakTimeAbs.push_back(hit->PeakTimeAbs());
    fPMTPeakTime.push_back(hit->PeakTime());
    fPMTWidth.push_back(hit->Width());
    fPMTArea.push_back(hit->Area());
    fPMTAmplitude.push_back(hit->Amplitude());
    fPMTPE.push_back(hit->PE());
    fPMTFastToTotal.push_back(hit->FastToTotal());
}




  //Fill the tree
  fTree->Fill();
  fFlashTotalPE.clear();
  fFlashTime.clear();
  fFlashTimeWidth.clear();
  fYCenter.clear();
  fZCenter.clear();
  fYWidth.clear();
  fZWidth.clear();

  fPMTOpChannel.clear();
  fPMTPeakTimeAbs.clear();
  fPMTPeakTime.clear();
  fPMTWidth.clear();
  fPMTArea.clear();
  fPMTAmplitude.clear();
  fPMTPE.clear();
  fPMTFastToTotal.clear();


  return;
}// end analyze

//////////////////////////////////////////////////////////////////////////
void ana::PDSValidation::endJob() {

  cout << "Number of events processed: " <<  numevents  << endl;

}// end endJob


DEFINE_ART_MODULE(ana::PDSValidation)
