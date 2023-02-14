////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Class:       PDSSimValidation                                                                             //
// Modlue Type: Analyser                                                                                  //
// File         PDSSimValidation_module.cc                                                                   //
// Author:      Chris Hilgenberg chilgenb@fnal.gov                                                        //
// Usage:       Simple workflow to validate PDS reconstruction                                            //
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
#include "art/Utilities/make_tool.h"

//Larsoft includes
#include "larcore/CoreUtils/ServiceUtil.h"
#include "larcore/Geometry/Geometry.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardataobj/Simulation/SimPhotons.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCParticle.h"

//Root Includes
#include "TMath.h"
#include "TTree.h"
#include "TSystem.h"
#include "TH1F.h"

//C++ Includes
#include <iostream>
#include <fstream>
#include <limits>
#include <map>
#include <numeric>
#include <vector>

//sbncode includes
#include "sbncode/OpDet/PDMapAlg.h"

using namespace std;


namespace ana {
  class PDSSimValidation;
}

class ana::PDSSimValidation : public art::EDAnalyzer {
  public:

    explicit PDSSimValidation(const fhicl::ParameterSet& pset);
    // The compiler-generated destructor is fine for non-base
    // classes without bare pointers or other resource use.

    // Plugins should not be copied or assigned.
    PDSSimValidation(PDSSimValidation const&) = delete;
    PDSSimValidation(PDSSimValidation&&) = delete;
    PDSSimValidation& operator=(PDSSimValidation const&) = delete;
    PDSSimValidation& operator=(PDSSimValidation&&) = delete;

    // required abstract methods we need to impliment
    void analyze(const art::Event& evt) override;
    void endJob() override;
    void beginJob() override;

  private:
    void ResetTreeVariables();
    void FillPERecoveringResolution();

    // fcl parameters
    const string fNuGeneratorLabel;
    //std::vector<std::string> fLArGeantModuleLabel;
    const bool fVerbose;
    const bool fUseReflectedLight;
    const bool fSaveByPDType;
    const double fPECut;

    // PDS map
    const size_t fNOpChannels;
    const std::unique_ptr<opdet::PDMapAlg> fPDSMapPtr;
    const std::vector<std::string> fPMTMapLabel, fArapucaMapLabel;

    // Service handles backtracker
    art::ServiceHandle<art::TFileService> fTFS;
    detinfo::DetectorClocksData fClockData = art::ServiceHandle<detinfo::DetectorClocksService>()->DataForJob();

    int fNumEvents;

    // tree variables
    TTree* fTree;
    int fRun;    ///< art run number
    int fSubRun; ///< art subrun number
    int fEvent;  ///< art event number
    std::vector<int> fNPhotonsVec;
    std::map<int,int> fNTimePhotonsMap;
    int fNPhotons;
    int fTimePhotons;
    int fNPhotonsPMTCoated;
    int fNPhotonsPMTUncoated;
    int fNPhotonsXARAPUCAVuv;
    int fNPhotonsXARAPUCAVis;

}; //end def class PDSSimValidation

//////////////////////////////////////////////////////////////////////////////////////
ana::PDSSimValidation::PDSSimValidation(const fhicl::ParameterSet& pset) :
  EDAnalyzer(pset),
  fNuGeneratorLabel (pset.get<string>("NuGeneratorLabel","generator")),
  //fLArGeantModuleLabel = pset.get<std::vector<std::string>>("LArGeantModuleLabel",{"pdfastsim", "pdfastsimout"});
  fVerbose (pset.get<bool>("Verbose",false)),
  fUseReflectedLight (pset.get<bool>("UseReflectedLight",true)),
  fSaveByPDType (pset.get<bool>("SaveByPDType",true)),
  fPECut (pset.get<double>("PECut", 10)),
  fNOpChannels ((lar::providerFrom<geo::Geometry>())->NOpChannels()),
  fPDSMapPtr (art::make_tool<opdet::PDMapAlg>(pset.get<fhicl::ParameterSet>("PDSMapTool"))),
  fPMTMapLabel (pset.get<std::vector<std::string>>("PMTMapLabel",{"pmt_coated", "pmt_uncoated"})),
  fArapucaMapLabel (pset.get<std::vector<std::string>>("ArapucaMapLabel",{"xarapuca_vuv", "xarapuca_vis"})),
  fNumEvents(0)
{
} //end constructor

///////////////////////////////////////////////////////////////////////////////
void ana::PDSSimValidation::beginJob() {

  fTree = fTFS->make<TTree>("pdsSimTree", "Tree with PDS simulation validation information");

  fTree->Branch("Run",          &fRun,          "Run/I");
  fTree->Branch("SubRun",       &fSubRun,       "SubRun/I");
  fTree->Branch("Event",        &fEvent,        "Event/I");

  fTree->Branch("NPhotonsVec", &fNPhotonsVec);
  fTree->Branch("NTimePhotonsMap", &fNTimePhotonsMap);
  fTree->Branch("NPhotons", &fNPhotons, "NPhotons/I");
  if(fSaveByPDType){
    fTree->Branch("NPhotonsPMTCoated", &fNPhotonsPMTCoated, "NPhotonsPMTCoated/I");
    fTree->Branch("NPhotonsPMTUncoated", &fNPhotonsPMTUncoated, "NPhotonsPMTUncoated/I");
    fTree->Branch("NPhotonsXARAPUCAVuv", &fNPhotonsXARAPUCAVuv, "NPhotonsXARAPUCAVuv/I");
    fTree->Branch("NPhotonsXARAPUCAVis", &fNPhotonsXARAPUCAVis, "NPhotonsXARAPUCAVis/I");
  }

}// end beginJob

////////////////////////////////////////////////////////////////////
void ana::PDSSimValidation::analyze(const art::Event& evt) {
  //Reset variables
  ResetTreeVariables();

  fRun    = evt.run();
  fSubRun = evt.subRun();
  fEvent  = evt.event();

  if(fVerbose) {
    cout << "Analysing run " << fRun << ", subrun " << fSubRun << ", event: " << fEvent << endl;
  }
  fNumEvents++;

  //Get truth info
  art::Handle< std::vector<simb::MCTruth> > MCTruthHandle;
  evt.getByLabel(fNuGeneratorLabel, MCTruthHandle);
  std::vector<art::Ptr<simb::MCTruth> > mctruth_v;
  art::fill_ptr_vector(mctruth_v, MCTruthHandle);
  int nuT0 = std::numeric_limits<int>::max();
  for (size_t n = 0; n < mctruth_v.size(); n++) {
    art::Ptr<simb::MCTruth> evtTruth = mctruth_v[n];
    for (int  p= 0; p < evtTruth->NParticles(); p++){
      simb::MCParticle const& mcpar = evtTruth->GetParticle(p);
      if( (std::abs(mcpar.PdgCode()) == 14 || std::abs(mcpar.PdgCode()) == 12)
          && mcpar.Mother() == -1) {
        if(mcpar.T() < nuT0){
          nuT0 = mcpar.T(); // ns
          // nuT0 = clockData.G4ToElecTime(mcpar.T());
          // std::cout<<clockData.G4ToElecTime(mcpar.T())<<" "<<mcpar.T()<<std::endl;
          if(fVerbose) std::cout << "T0: "<< nuT0 << " NMCTruth= " << n << std::endl;
        }
      }
    }
  }
  if(nuT0 >= std::numeric_limits<int>::max()) nuT0 = 0;

  //###############################################
  //### Get PDS sim info for the event         ###
  //###############################################
  std::vector<art::Handle<std::vector<sim::SimPhotonsLite> >> litePhotonHandle_list;
  litePhotonHandle_list = evt.getMany<std::vector<sim::SimPhotonsLite>>();
  std::vector<int> simPhotonsPE_v(fNOpChannels, 0);
  std::vector<std::map<int,int>> simTimePhotonsPE_v(fNOpChannels);

  for ( const auto& litePhotonHandle: (litePhotonHandle_list) ){
    vector<art::Ptr<sim::SimPhotonsLite> > litePhotonList;
    art::fill_ptr_vector(litePhotonList, litePhotonHandle);
    bool reflected = (litePhotonHandle.provenance()->productInstanceName()== "Reflected");
    if(fVerbose){
      std::cout << "*********SimPhotonLite List: "
                << litePhotonHandle.provenance()->moduleLabel()
                << " Provenance:" << litePhotonHandle.provenance()->productInstanceName()
                << std::endl;
    }
    for ( auto const& litePhotons : (*litePhotonHandle) ){
      std::string pd_type = fPDSMapPtr->pdType(litePhotons.OpChannel);
      if(reflected && pd_type=="xarapuca_vuv") continue;
      if(!reflected && (pd_type=="xarapuca_vis" || pd_type=="pmt_uncoated")) continue;

      std::map<int, int> litePhotons_map = litePhotons.DetectedPhotons;
      int nphotons = std::accumulate(
        litePhotons_map.begin(), litePhotons_map.end(), 0,
        [] (int sum, const auto& phmap){ return sum+phmap.second; } );
      simPhotonsPE_v[litePhotons.OpChannel] += nphotons;

      for(const auto& lp : litePhotons_map) {
        simTimePhotonsPE_v[litePhotons.OpChannel][lp.first-nuT0] = lp.second;
      }
    }
  }
  for(size_t opch=0; opch<simPhotonsPE_v.size(); opch++){
    if(simPhotonsPE_v[opch]>0){
      fNPhotonsVec.push_back(simPhotonsPE_v[opch]);
      fNPhotons += simPhotonsPE_v[opch];
      if(fSaveByPDType){
        std::string pd_type=fPDSMapPtr->pdType(opch);
        if(pd_type=="pmt_coated") {
          fNPhotonsPMTCoated += simPhotonsPE_v[opch];
        }
        else if(pd_type=="pmt_uncoated"){
          fNPhotonsPMTUncoated += simPhotonsPE_v[opch];
        }
        else if(pd_type=="xarapuca_vuv"){
          fNPhotonsXARAPUCAVuv += simPhotonsPE_v[opch];
        }
        else if(pd_type=="xarapuca_vis"){
          fNPhotonsXARAPUCAVis += simPhotonsPE_v[opch];
        }
      }
    }
  }
  for(size_t opch=0; opch<simTimePhotonsPE_v.size(); opch++){
    const auto& lpm = simTimePhotonsPE_v.at(opch);
    for(const auto& lp : lpm) {
      if(lp.second > 0){
        if ( auto it{fNTimePhotonsMap.find(lp.first)}; it!=fNTimePhotonsMap.end() )
          it->second += lp.second;
        else fNTimePhotonsMap[lp.first] = lp.second;
      }
    }
  }
  if(fVerbose) {
    std::cout << "NPhPMTCoat: " << fNPhotonsPMTCoated   << ", NPhPMTUnco: " << fNPhotonsPMTUncoated << ", "
              << "NPhXAraVUV: " << fNPhotonsXARAPUCAVuv << ", NPhXAraVIS: " << fNPhotonsXARAPUCAVis << std::endl;
  }


  //Fill the tree
  fTree->Fill();

  return;
}// end analyze

//////////////////////////////////////////////////////////////////////////
void ana::PDSSimValidation::endJob(){
  cout << "Number of events processed: " <<  fNumEvents  << endl;
}


void ana::PDSSimValidation::ResetTreeVariables() {
  fNPhotons = 0;
  fNPhotonsVec.clear();
  fNTimePhotonsMap.clear();
  if(fSaveByPDType){
    fNPhotonsPMTCoated = 0;
    fNPhotonsPMTUncoated = 0;
    fNPhotonsXARAPUCAVuv = 0;
    fNPhotonsXARAPUCAVis = 0;
  }

}

DEFINE_ART_MODULE(ana::PDSSimValidation)
