////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Class:       CRTValidation                                                                             //
// Modlue Type: Analyser                                                                                  //
// File         CRTValidation_module.cc                                                                //
// Author:      Chris Hilgenberg chilgenb@fnal.gov                                                        //
//                                                                                                        //
// Usage:       The Validation module takes an array of shower modules and perfoms truth matching on      //
//              reconstructed EM showers. It calculates the truth information from geant information      //
//              such as the diretion and starting position and energy deposited. Metrics are then defined //
//              to compare the efficiency of the shower reconstruction modules. These are presented       //
//              as histograms and in terms of energy. Energy plots look at the metrics for specific       //
//              energies between +- fEnergyWidth of the values described in the fcl table                 //
//                                                                                                        //
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
#include "sbnobj/Common/CRT/CRTHit.hh"
#include "sbndcode/RecoUtils/RecoUtils.h"

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
  class CRTValidation;
}

class ana::CRTValidation : public art::EDAnalyzer {
  public:

    explicit CRTValidation(const fhicl::ParameterSet& pset);
    // The compiler-generated destructor is fine for non-base
    // classes without bare pointers or other resource use.

    // Plugins should not be copied or assigned.
    CRTValidation(CRTValidation const&) = delete;
    CRTValidation(CRTValidation&&) = delete;
    CRTValidation& operator=(CRTValidation const&) = delete;
    CRTValidation& operator=(CRTValidation&&) = delete;

    // required abstract methods we need to impliment 
    void analyze(const art::Event& evt) override;
    void endJob() override;
    void beginJob() override;

  private:

    //fcl parameters
    string fGenieGenModuleLabel;
    string fLArGeantModuleLabel;
    string fCRTHitLabel;
    bool   fVerbose;

    //TTree
    TTree* fTree;

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

    int fNHit; ///< Number of CRTHits in event
    vector<float> fPESHit; ///<Number of photoelectrons associated with hit
    vector<float> fx_pos;
    vector<float> fx_err;
    vector<float> fy_pos;
    vector<float> fy_err;
    vector<float> fz_pos;
    vector<float> fz_err;

}; //end def class CRTValidation

//////////////////////////////////////////////////////////////////////////////////////
ana::CRTValidation::CRTValidation(const fhicl::ParameterSet& pset) : 
  EDAnalyzer(pset),
  numevents(0)
{

  fGenieGenModuleLabel = pset.get<string>("GenieGenModuleLabel","generator");
  fLArGeantModuleLabel = pset.get<string>("LArGeantModuleLabel","largeant");
  fCRTHitLabel = pset.get<string>("CRTHitLabel","crthit");
  fVerbose             = pset.get<bool>("Verbose",false); 

} //end constructor

///////////////////////////////////////////////////////////////////////////////
void ana::CRTValidation::beginJob() {
  
  fTree = tfs->make<TTree>("crtTree", "Tree with CRT validation information");
  //gInterpreter->GenerateDictionary("vector<vector<float> > ","vector");

  fTree->Branch("Run",          &fRun,          "Run/I");
  fTree->Branch("SubRun",       &fSubRun,       "SubRun/I");
  fTree->Branch("Event",        &fEvent,        "Event/I");
  fTree->Branch("NHit",       &fNHit,       "NHit/I");
  fTree->Branch("PESHit",&fPESHit);
  fTree->Branch("x_pos", &fx_pos);
  fTree->Branch("x_err", &fx_err);
  fTree->Branch("y_pos", &fy_pos);
  fTree->Branch("y_err", &fy_err);
  fTree->Branch("z_pos", &fz_pos);
  fTree->Branch("z_err", &fz_err);

}// end beginJob

////////////////////////////////////////////////////////////////////
void ana::CRTValidation::analyze(const art::Event& evt) {

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
  //### Get CRT info for the event         ###
  //###############################################
  //Getting  CRTHit information
  art::Handle< vector<sbn::crt::CRTHit> > CRTHitListHandle;
  vector<art::Ptr<sbn::crt::CRTHit> > CRTHitList;
  if(evt.getByLabel(fCRTHitLabel, CRTHitListHandle)){
      art::fill_ptr_vector(CRTHitList, CRTHitListHandle);
  }

  fNHit = CRTHitList.size();
  for(auto const& hit : CRTHitList) {
    fPESHit.push_back(hit->peshit); 
    fx_pos.push_back(hit->x_pos);
    fx_err.push_back(hit->x_err);
    fy_pos.push_back(hit->y_pos);
    fy_err.push_back(hit->y_err);
    fz_pos.push_back(hit->z_pos);
    fz_err.push_back(hit->z_err);
  }

  //Fill the tree
  fTree->Fill();
  fPESHit.clear();
  fx_pos.clear();
  fx_err.clear();
  fy_pos.clear();
  fy_err.clear();
  fz_pos.clear();
  fz_err.clear();

  return;
}// end analyze

//////////////////////////////////////////////////////////////////////////
void ana::CRTValidation::endJob() {

  cout << "Number of events processed: " <<  numevents  << endl;

}// end endJob


DEFINE_ART_MODULE(ana::CRTValidation)
