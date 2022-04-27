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
#include "larcore/Geometry/Geometry.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "lardataobj/Simulation/SimPhotons.h"
#include "lardataobj/RecoBase/OpFlash.h"
#include "lardataobj/RecoBase/OpHit.h"

//Root Includes
#include "TMath.h"
#include "TTree.h"
#include "TSystem.h"
#include "TH1F.h"

//C++ Includes
#include <vector>
#include <numeric>
#include <iostream>
#include <fstream>

//sbncode includes
#include "sbncode/OpDet/PDMapAlg.h"

using namespace std;

const double DEFAULT_T0VALUE = 1e9;

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

    //fcl parameters
    string fNuGeneratorLabel;
    //std::vector<std::string> fLArGeantModuleLabel;

    //PDS map
    std::unique_ptr<opdet::PDMapAlg> fPDSMapPtr;
    std::vector<std::string> fPMTMapLabel, fArapucaMapLabel;

    bool fVerbose;
    bool fUseReflectedLight;
    bool fSaveByPDType;
    double fPECut;
    double fPDSResponseDelay;

    //TTree
    TTree* fTree;

    //Service handlesacktracker
    art::ServiceHandle<geo::Geometry> geom;
    art::ServiceHandle<art::TFileService> tfs;
    detinfo::DetectorClocksData clockData = art::ServiceHandle<detinfo::DetectorClocksService>()->DataForJob();

    int numevents;

    // tree variables
    int fRun;    ///< art run number
    int fSubRun; ///< art subrun number
    int fEvent;  ///< art event number

    std::vector<int> fSimPhotonsPE_v;
    double fNuT0;

    //True Variable
    std::vector<int> fNPhotonsVec;
    int fNPhotons;
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
  fPDSMapPtr (art::make_tool<opdet::PDMapAlg>(pset.get<fhicl::ParameterSet>("PDSMapTool"))),
  fPMTMapLabel (pset.get<std::vector<std::string>>("PMTMapLabel",{"pmt_coated", "pmt_uncoated"})),
  fArapucaMapLabel (pset.get<std::vector<std::string>>("ArapucaMapLabel",{"xarapuca_vuv", "xarapuca_vis"})),
  fVerbose (pset.get<bool>("Verbose",false)),
  fUseReflectedLight (pset.get<bool>("UseReflectedLight",true)),
  fSaveByPDType (pset.get<bool>("SaveByPDType",true)),
  fPECut (pset.get<double>("PECut", 10)),
  fPDSResponseDelay (pset.get<double>("PDSResponseDelay", 0.2)),
  numevents(0),
  fSimPhotonsPE_v (geom->NOpChannels(), 0)
{
} //end constructor

///////////////////////////////////////////////////////////////////////////////
void ana::PDSSimValidation::beginJob() {

  fTree = tfs->make<TTree>("pdsSimTree", "Tree with PDS simulation validation information");

  fTree->Branch("Run",          &fRun,          "Run/I");
  fTree->Branch("SubRun",       &fSubRun,       "SubRun/I");
  fTree->Branch("Event",        &fEvent,        "Event/I");

  fTree->Branch("NPhotonsVec", &fNPhotonsVec);
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
  numevents++;


  //###############################################
  //### Get PDS sim info for the event         ###
  //###############################################
  std::vector<art::Handle<std::vector<sim::SimPhotonsLite> >> fLitePhotonHandle_list;
  fLitePhotonHandle_list = evt.getMany<std::vector<sim::SimPhotonsLite>>();
  for ( const art::Handle<std::vector<sim::SimPhotonsLite>>& fLitePhotonHandle: (fLitePhotonHandle_list) ){
    vector<art::Ptr<sim::SimPhotonsLite> > fLitePhotonList;
    art::fill_ptr_vector(fLitePhotonList, fLitePhotonHandle);
    bool reflected = (fLitePhotonHandle.provenance()->productInstanceName()== "Reflected");
    if(fVerbose){
      std::cout<<"*********SimPhotonLite List: "<<fLitePhotonHandle.provenance()->moduleLabel()<<" Provenance:"<<fLitePhotonHandle.provenance()->productInstanceName()<<std::endl;
    }
    for ( auto const& fLitePhotons : (*fLitePhotonHandle) ){
      std::string pd_type=fPDSMapPtr->pdType(fLitePhotons.OpChannel);
      if(reflected && pd_type=="xarapuca_vuv") continue;
      if(!reflected && (pd_type=="xarapuca_vis" || pd_type=="pmt_uncoated")) continue;

      std::map<int, int> fLitePhotons_map = fLitePhotons.DetectedPhotons;

      int nphotons=std::accumulate( std::begin(fLitePhotons_map), std::end(fLitePhotons_map), 0
            , [] (int sum, const std::map<int, int>::value_type& phmap){ return sum+phmap.second; } );

      fSimPhotonsPE_v[fLitePhotons.OpChannel]+=nphotons;
    }
  }
  for(size_t opch=0; opch<fSimPhotonsPE_v.size(); opch++){
    if(fSimPhotonsPE_v[opch]>0){
      fNPhotonsVec.push_back(fSimPhotonsPE_v[opch]);
      fNPhotons+=fSimPhotonsPE_v[opch];
      if(fSaveByPDType){
        std::string pd_type=fPDSMapPtr->pdType(opch);
        if(pd_type=="pmt_coated") {
          fNPhotonsPMTCoated+=fSimPhotonsPE_v[opch];
        }
        else if(pd_type=="pmt_uncoated"){
          fNPhotonsPMTUncoated+=fSimPhotonsPE_v[opch];
        }
        else if(pd_type=="xarapuca_vuv"){
          fNPhotonsXARAPUCAVuv+=fSimPhotonsPE_v[opch];
        }
        else if(pd_type=="xarapuca_vis"){
          fNPhotonsXARAPUCAVis+=fSimPhotonsPE_v[opch];
        }
      }
    }
  }
  if(fVerbose) std::cout<<"NPH: "<<fNPhotonsPMTCoated<<" "<<fNPhotonsPMTUncoated<<" "<<fNPhotonsXARAPUCAVuv<<" "<<fNPhotonsXARAPUCAVis<<"\n";


  //Get truth info
  art::Handle< std::vector<simb::MCTruth> > MCTruthHandle;
  evt.getByLabel(fNuGeneratorLabel, MCTruthHandle);
  std::vector<art::Ptr<simb::MCTruth> > mctruth_v;
  art::fill_ptr_vector(mctruth_v, MCTruthHandle);

  for (size_t n = 0; n < mctruth_v.size(); n++) {
    art::Ptr<simb::MCTruth> evtTruth = mctruth_v[n];
    for (int  p= 0; p < evtTruth->NParticles(); p++){
      simb::MCParticle const& mcpar = evtTruth->GetParticle(p);
      if( (std::abs(mcpar.PdgCode()) == 14 || std::abs(mcpar.PdgCode()) == 12) && mcpar.Mother()==-1) {
        if(mcpar.T()<fNuT0){
          fNuT0 = mcpar.T()/1000;
          //fNuT0 = clockData.G4ToElecTime(mcpar.T());
          //std::cout<<clockData.G4ToElecTime(mcpar.T())<<" "<<mcpar.T()<<std::endl;
          if(fVerbose) std::cout<<"T0: "<<fNuT0<<" NMCTruth="<<n<<std::endl;
        }
      }
    }
  }
  if(fNuT0==DEFAULT_T0VALUE) fNuT0=0;

  //Fill the tree
  fTree->Fill();

  return;
}// end analyze

//////////////////////////////////////////////////////////////////////////
void ana::PDSSimValidation::endJob(){
  cout << "Number of events processed: " <<  numevents  << endl;
}


void ana::PDSSimValidation::ResetTreeVariables() {
  std::fill(fSimPhotonsPE_v.begin(), fSimPhotonsPE_v.end(), 0);
  fNuT0=DEFAULT_T0VALUE;

  fNPhotons=0;
  fNPhotonsVec.clear();
  if(fSaveByPDType){
    fNPhotonsPMTCoated=0;
    fNPhotonsPMTUncoated=0;
    fNPhotonsXARAPUCAVuv=0;
    fNPhotonsXARAPUCAVis=0;
  }

}

DEFINE_ART_MODULE(ana::PDSSimValidation)
