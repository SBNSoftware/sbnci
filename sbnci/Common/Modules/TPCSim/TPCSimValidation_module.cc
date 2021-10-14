////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Class:       TPCSimValidation                                                                             //
// Modlue Type: Analyser                                                                                  //
// File         TPCSimValidation_module.cc                                                                //
// Author:      Chris Hilgenberg chilgenb@fnal.gov                                                        //
//                                                                                                        //
// Usage:       The Validation module for TPC simulation                                                   //
//09/10/2021 Adopted for  (Sergey Martynenko smartynen@bnl.gov)                                      //
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
#include "lardataobj/RawData/RawDigit.h"

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
  class TPCSimValidation;
}

class ana::TPCSimValidation : public art::EDAnalyzer {
  public:

    explicit TPCSimValidation(const fhicl::ParameterSet& pset);
    // The compiler-generated destructor is fine for non-base
    // classes without bare pointers or other resource use.

    // Plugins should not be copied or assigned.
    TPCSimValidation(TPCSimValidation const&) = delete;
    TPCSimValidation(TPCSimValidation&&) = delete;
    TPCSimValidation& operator=(TPCSimValidation const&) = delete;
    TPCSimValidation& operator=(TPCSimValidation&&) = delete;

    // required abstract methods we need to impliment 
    void analyze(const art::Event& evt) override;
    void endJob() override;
    void beginJob() override;

  private:

    //fcl parameters
    string fGenieGenModuleLabel;
    string fLArGeantModuleLabel;
    std::vector<art::InputTag> fRawDigitLabelVect;
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
	
    int fNRawDigits; ///<Number of RawDigits in event
//    vector<short> fADCs; ///<Reference to the compressed ADC count vector. not working right now.
    vector<size_t> fNADC; ///<Number of elements in the compressed ADC sample vector
    vector<short> fADC; ///<ADC vector element number i; no decompression is applied.
    vector<raw::ChannelID_t> fChannel; ///<DAQ channel this raw data was read from.
    vector<ULong64_t> fSamples; ///<Number of samples in the uncompressed ADC data.
    vector<float> fPedestal; 
    vector<float> fSigma;
//    vector<raw::Compress_t> fCompression; ///<Compression algorithm used to store the ADC counts.

}; //end def class TPCSimValidation

//////////////////////////////////////////////////////////////////////////////////////
ana::TPCSimValidation::TPCSimValidation(const fhicl::ParameterSet& pset) : 
  EDAnalyzer(pset),
  numevents(0)
{

  fGenieGenModuleLabel = pset.get<string>("GenieGenModuleLabel","generator");
  fLArGeantModuleLabel = pset.get<string>("LArGeantModuleLabel","largeant");
  fRawDigitLabelVect = pset.get<std::vector<art::InputTag>>("RawDigitLabel",std::vector<art::InputTag>()={"daq"});
  fVerbose             = pset.get<bool>("Verbose",false); 

} //end constructor

///////////////////////////////////////////////////////////////////////////////
void ana::TPCSimValidation::beginJob() {
  
  fTree = tfs->make<TTree>("tpcsimTree", "Tree with TPC validation information");
  //gInterpreter->GenerateDictionary("vector<vector<float> > ","vector");
  

  fTree->Branch("Run",          &fRun,          "Run/I");
  fTree->Branch("SubRun",       &fSubRun,       "SubRun/I");
  fTree->Branch("Event",        &fEvent,        "Event/I");
  fTree->Branch("NRawDigits",       &fNRawDigits,       "NRawDigits/I");
//  fTree->Branch("ADCs",&fADCs);
  fTree->Branch("NADC",&fNADC);
  fTree->Branch("ADC", &fADC);
  fTree->Branch("Channel", &fChannel);
  fTree->Branch("Samples", &fSamples);
  fTree->Branch("Pedestal", &fPedestal);
  fTree->Branch("Sigma", &fSigma);
//  fTree->Branch("Compression", &fCompression);
}// end beginJob

////////////////////////////////////////////////////////////////////
void ana::TPCSimValidation::analyze(const art::Event& evt) {

  fRun    = evt.run();
  fSubRun = evt.subRun();
  fEvent  = evt.event();
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
          MF_LOG_INFO("TPCSimValidation")
               << "True Particle with track ID: " << particle->TrackId() << " Has code of: " 
               << particle->PdgCode() << " and Energy of: " << particle->E() << " With Mother: " 
               << particle->Mother() << " Proccess: " << particle->Process() << " End Process: "  
               << particle->EndProcess();
      }// if verbose

  } // for MCParticles


  //###############################################
  //### Get TPC info for the event         ###
  //###############################################
  //Getting  RawDigit information
  for(const auto& rawDigitLabel : fRawDigitLabelVect)
  { 
     art::Handle< vector<raw::RawDigit> > RawDigitListHandle;
     vector<art::Ptr<raw::RawDigit> > RawDigitList;
     if(evt.getByLabel(rawDigitLabel, RawDigitListHandle)){
         art::fill_ptr_vector(RawDigitList, RawDigitListHandle);
     }

     fNRawDigits = RawDigitList.size();
     MF_LOG_DEBUG("TPCSimValidation") << "RawDigit label: " << rawDigitLabel;

     for(auto const& digit : RawDigitList) {
         MF_LOG_DEBUG("TPCSimValidation") << "NADC ="<< digit->NADC();
         for(size_t iadc=0; iadc<digit->NADC(); iadc++) {
           fADC.push_back(digit->ADC(iadc));
         }
         
//         fADCs.push_back(digit->ADCs());
         fNADC.push_back(digit->NADC()); 
         fChannel.push_back(digit->Channel());
         fSamples.push_back(digit->Samples());
         fPedestal.push_back(digit->GetPedestal());
         fSigma.push_back(digit->GetSigma());
//         fCompression.push_back(digit->Compression());
     }  
  }

  //Fill the tree
  fTree->Fill();
//  fADCs.clear();
  fNADC.clear();
  fADC.clear();
  fSamples.clear();
  fPedestal.clear();
  fSigma.clear();
//  fCompression.clear();

  return;
}// end analyze

//////////////////////////////////////////////////////////////////////////
void ana::TPCSimValidation::endJob() {

  MF_LOG_INFO("TPCSimValidation") << "Number of events processed: " <<  numevents;

}// end endJob


DEFINE_ART_MODULE(ana::TPCSimValidation)
