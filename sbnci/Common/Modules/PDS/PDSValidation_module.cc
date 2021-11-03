////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Class:       PDSValidation                                                                             //
// Modlue Type: Analyser                                                                                  //
// File         PDSValidation_module.cc                                                                   //
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
    void ResetTreeVariables();
    void FillPERecoveringResolution();

    //fcl parameters
    string fNuGeneratorLabel;
    //std::vector<std::string> fLArGeantModuleLabel;
    string fOpHitPMTLabel;
    string fOpHitArapucaLabel;
    std::vector<string> fOpFlashLabel;

    //PDS map
    std::unique_ptr<opdet::PDMapAlg> fPDSMapPtr;
    std::vector<std::string> fPMTMapLabel, fArapucaMapLabel;

    bool fVerbose;
    bool fUseReflectedLight;
    bool fSaveByPDType;
    bool fUseArapucas;
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
    std::vector<double> fOpHitPE_v;
    double fNuT0;

    //True Variable
    std::vector<int> fNPhotonsVec;
    int fNPhotons;
    int fNPhotonsPMTCoated;
    int fNPhotonsPMTUncoated;
    int fNPhotonsXARAPUCAVuv;
    int fNPhotonsXARAPUCAVis;

    //Reco variables
    int fPMTNHit;
    vector<double> fPMTPeakTime;
    vector<double> fPMTWidth;
    vector<double> fPMTAmplitude;
    vector<double> fPMTPE;
    std::vector<double> fPMTPEResolutionVec;

    int fArapucaNHit;
    vector<double> fArapucaPeakTime;
    vector<double> fArapucaWidth;
    vector<double> fArapucaAmplitude;
    vector<double> fArapucaPE;
    std::vector<double> fArapucaPEResolutionVec;

    int fNFlash; ///< number of OpFlashes in the event
    vector<double> fFlashTotalPE; ///< total PE associated with single OpFlash
    vector<double> fFlashTime; ///< time of flash
    vector<double> fYCenter;
    vector<double> fZCenter;
    std::vector<double> fFlashTimeResolution;
}; //end def class PDSValidation

//////////////////////////////////////////////////////////////////////////////////////
ana::PDSValidation::PDSValidation(const fhicl::ParameterSet& pset) :
  EDAnalyzer(pset),
  fNuGeneratorLabel (pset.get<string>("NuGeneratorLabel","generator")),
  //fLArGeantModuleLabel = pset.get<std::vector<std::string>>("LArGeantModuleLabel",{"pdfastsim", "pdfastsimout"});
  fOpHitPMTLabel (pset.get<string>("OpHitPMTLabel", "ophitpmt")),
  fOpHitArapucaLabel (pset.get<string>("OpHitArapucaLabel", "ophitarapuca")),
  fOpFlashLabel (pset.get<std::vector<std::string>>("OpFlashLabel",{"opflashtpc0", "opflashtpc1"})),
  fPDSMapPtr (art::make_tool<opdet::PDMapAlg>(pset.get<fhicl::ParameterSet>("PDSMapTool"))),
  fPMTMapLabel (pset.get<std::vector<std::string>>("PMTMapLabel",{"pmt_coated", "pmt_uncoated"})),
  fArapucaMapLabel (pset.get<std::vector<std::string>>("ArapucaMapLabel",{"xarapuca_vuv", "xarapuca_vis"})),
  fVerbose (pset.get<bool>("Verbose",false)),
  fUseReflectedLight (pset.get<bool>("UseReflectedLight",true)),
  fSaveByPDType (pset.get<bool>("SaveByPDType",true)),
  fUseArapucas (pset.get<bool>("UseArapucas", true)),
  fPECut (pset.get<double>("PECut", 10)),
  fPDSResponseDelay (pset.get<double>("PDSResponseDelay", 0.2)),
  numevents(0),
  fSimPhotonsPE_v (geom->NOpChannels(), 0),
  fOpHitPE_v (geom->NOpChannels(), 0)
{
} //end constructor

///////////////////////////////////////////////////////////////////////////////
void ana::PDSValidation::beginJob() {

  fTree = tfs->make<TTree>("pdsTree", "Tree with PDS validation information");

  fTree->Branch("Run",          &fRun,          "Run/I");
  fTree->Branch("SubRun",       &fSubRun,       "SubRun/I");
  fTree->Branch("Event",        &fEvent,        "Event/I");
  fTree->Branch("NFlash",       &fNFlash,       "NFlash/I");
  fTree->Branch("PMTNHit",         &fPMTNHit,         "PMTNHit/I");
  fTree->Branch("ArapucaNHit",         &fArapucaNHit,         "ArapucaNHit/I");

  fTree->Branch("NPhotonsVec", &fNPhotonsVec);
  fTree->Branch("NPhotons", &fNPhotons, "NPhotons/I");
  if(fSaveByPDType){
    fTree->Branch("NPhotonsPMTCoated", &fNPhotonsPMTCoated, "NPhotonsPMTCoated/I");
    fTree->Branch("NPhotonsPMTUncoated", &fNPhotonsPMTUncoated, "NPhotonsPMTUncoated/I");
    fTree->Branch("NPhotonsXARAPUCAVuv", &fNPhotonsXARAPUCAVuv, "NPhotonsXARAPUCAVuv/I");
    fTree->Branch("NPhotonsXARAPUCAVis", &fNPhotonsXARAPUCAVis, "NPhotonsXARAPUCAVis/I");
  }

  fTree->Branch("PMTPeakTime",&fPMTPeakTime);
  fTree->Branch("PMTWidth",&fPMTWidth);
  fTree->Branch("PMTAmplitude",&fPMTAmplitude);
  fTree->Branch("PMTPE",&fPMTPE);
  fTree->Branch("PMTPEResolutionVec",&fPMTPEResolutionVec);

  if(fUseArapucas){
    fTree->Branch("ArapucaPeakTime",&fArapucaPeakTime);
    fTree->Branch("ArapucaWidth",&fArapucaWidth);
    fTree->Branch("ArapucaAmplitude",&fArapucaAmplitude);
    fTree->Branch("ArapucaPE",&fArapucaPE);
    fTree->Branch("ArapucaPEResolutionVec",&fArapucaPEResolutionVec);
  }

  fTree->Branch("TotalFlashPE", &fFlashTotalPE);
  fTree->Branch("FlashTime", &fFlashTime);
  fTree->Branch("YCenter", &fYCenter);
  fTree->Branch("ZCenter", &fZCenter);
  fTree->Branch("FlashTimeResolution", &fFlashTimeResolution);

}// end beginJob

////////////////////////////////////////////////////////////////////
void ana::PDSValidation::analyze(const art::Event& evt) {
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

  //###############################################
  //### Get PDS reco info for the event         ###
  //###############################################
  art::Handle< vector<recob::OpHit> > opHitListHandle;

  //### Get OpHits PMT ###
  vector<art::Ptr<recob::OpHit> > opHitPMTList;
  if(evt.getByLabel(fOpHitPMTLabel, opHitListHandle)){
      art::fill_ptr_vector(opHitPMTList, opHitListHandle);
  }

  fPMTNHit =opHitPMTList.size();
  for(auto const& hit : opHitPMTList) {
    fPMTPeakTime.push_back(hit->PeakTime()-fNuT0-fPDSResponseDelay);
    fPMTWidth.push_back(hit->Width());
    fPMTAmplitude.push_back(hit->Amplitude());
    fPMTPE.push_back(hit->PE());
    fOpHitPE_v[hit->OpChannel()]+=hit->PE();
  }

  //### Get OpHits XARAPUCAs ###
  if(fUseArapucas){
    vector<art::Ptr<recob::OpHit> > opHitArapucaList;
    if(evt.getByLabel(fOpHitArapucaLabel, opHitListHandle)){
        art::fill_ptr_vector(opHitArapucaList, opHitListHandle);
    }

    fArapucaNHit =opHitArapucaList.size();
    for(auto const& hit : opHitArapucaList) {
      fArapucaPeakTime.push_back(hit->PeakTime());
      fArapucaWidth.push_back(hit->Width());
      fArapucaAmplitude.push_back(hit->Amplitude());
      fArapucaPE.push_back(hit->PE());
      fOpHitPE_v[hit->OpChannel()]+=hit->PE();
    }
  }

  //Fill PE recovering resolution for OpHits
  FillPERecoveringResolution();


  //### Get OpFlash   ###
  art::Handle< vector<recob::OpFlash> > opFlashListHandle;
  for (size_t s = 0; s < fOpFlashLabel.size(); s++) {
    vector<art::Ptr<recob::OpFlash> > opFlashList;
    if(evt.getByLabel(fOpFlashLabel[s], opFlashListHandle)){
      art::fill_ptr_vector(opFlashList, opFlashListHandle);
    }

    fNFlash +=opFlashList.size();
    for(auto const& flash : opFlashList) {
      fFlashTotalPE.push_back(flash->TotalPE());
      fFlashTime.push_back(flash->Time());
      fYCenter.push_back(flash->YCenter());
      fZCenter.push_back(flash->ZCenter());
      fFlashTimeResolution.push_back(flash->Time()-fNuT0-fPDSResponseDelay);
    }
  }

  //Fill the tree
  fTree->Fill();

  return;
}// end analyze

//////////////////////////////////////////////////////////////////////////
void ana::PDSValidation::endJob(){
  cout << "Number of events processed: " <<  numevents  << endl;
}


void ana::PDSValidation::ResetTreeVariables() {
  std::fill(fSimPhotonsPE_v.begin(), fSimPhotonsPE_v.end(), 0);
  std::fill(fOpHitPE_v.begin(), fOpHitPE_v.end(), 0);
  fNuT0=DEFAULT_T0VALUE;

  fNPhotons=0;
  fNPhotonsVec.clear();
  if(fSaveByPDType){
    fNPhotonsPMTCoated=0;
    fNPhotonsPMTUncoated=0;
    fNPhotonsXARAPUCAVuv=0;
    fNPhotonsXARAPUCAVis=0;
  }

  fPMTPeakTime.clear();
  fPMTWidth.clear();
  fPMTAmplitude.clear();
  fPMTPE.clear();
  fPMTPEResolutionVec.clear();

  if(fUseArapucas){
    fArapucaPeakTime.clear();
    fArapucaWidth.clear();
    fArapucaAmplitude.clear();
    fArapucaPE.clear();
    fArapucaPEResolutionVec.clear();
  }

  fNFlash=0;
  fFlashTotalPE.clear();
  fFlashTime.clear();
  fYCenter.clear();
  fZCenter.clear();
  fFlashTimeResolution.clear();
}

void ana::PDSValidation::FillPERecoveringResolution(){
  //double NPhotonsPMT=fNPhotonsPMTCoated+fNPhotonsPMTUncoated;
  //fPMTPEResolution=(std::accumulate( std::begin(fPMTPE), std::end(fPMTPE), 0 )-NPhotonsPMT)/NPhotonsPMT;
  //double NPhotonsArapuca=fNPhotonsXARAPUCAVuv+fNPhotonsXARAPUCAVuv;
  //fArapucaPEResolution=(std::accumulate( std::begin(fArapucaPE), std::end(fArapucaPE), 0 )-NPhotonsArapuca)/NPhotonsArapuca;
  for(size_t opch=0; opch<fOpHitPE_v.size(); opch++){
    if( std::find(fPMTMapLabel.begin(), fPMTMapLabel.end(), fPDSMapPtr->pdType(opch))!=fPMTMapLabel.end() ) {
      if(fVerbose) std::cout<<" Reco PE:"<<fOpHitPE_v[opch]<<" True PE:"<<fSimPhotonsPE_v[opch]<<"\n";
      if(fSimPhotonsPE_v[opch]>=fPECut) fPMTPEResolutionVec.push_back( (fOpHitPE_v[opch]-fSimPhotonsPE_v[opch])/fSimPhotonsPE_v[opch] );
    }
    else if(std::find(fArapucaMapLabel.begin(), fArapucaMapLabel.end(), fPDSMapPtr->pdType(opch))!=fArapucaMapLabel.end() && fUseArapucas ){
      if(fVerbose) std::cout<<"Reco PE:"<<fOpHitPE_v[opch]<<" True PE:"<<fSimPhotonsPE_v[opch]<<"\n";
      if(fSimPhotonsPE_v[opch]>=fPECut) fArapucaPEResolutionVec.push_back( (fOpHitPE_v[opch]-fSimPhotonsPE_v[opch])/fSimPhotonsPE_v[opch] );
    }
  }
}




DEFINE_ART_MODULE(ana::PDSValidation)
