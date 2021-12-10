///////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Class:       TrackValidation                                                                              //
// Plugin Type: analyzer (art v3_06_03)                                                                      //
// File:        TrackValidation_module.cc                                                                    //
// Author:      Ala Zglam aazglam1@sheffield.ac.uk                                                           //
//                                                                                                           //
// Usage:       The Validation module takes an array of track modules and perfoms truth matching on          //
//              reconstructed pandora tracks. It calculates the truth information from geant information     //
//              such as the diretion and starting position and energy deposited. Metrics are then defined    //
//              to compare the efficiency of the track reconstruction modules. These are presented           //
//              as histograms and in terms of energy.                                                        //
//                                                                                                           //
// Generated at Fri Aug 13 03:16:57 2021 by Ala Zglam using cetskelgen                                       //
// from cetlib version v3_11_01.                                                                             //
//                                                                                                           //
// Updates:                                                                                                  //
//                                                                                                           //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////

//Larsoft includes
#include "larcore/Geometry/Geometry.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Wire.h"
#include "lardataobj/RecoBase/Track.h"
//#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "larcoreobj/SimpleTypesAndConstants/PhysicalConstants.h"
#include "lardataalg/DetectorInfo/DetectorClocksData.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "larsim/MCCheater/BackTrackerService.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
//LAr data objects includes tests
#include "lardataobj/RecoBase/MCSFitResult.h"
#include "sbnobj/Common/Reco/RangeP.h"

// Additional framework includes
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileService.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "sbnci/Common/Modules/MCRecoUtils/RecoUtils.h"
#include "sbnci/Common/Modules/MCRecoUtils/TrackUtils.h"

// Root includes
#include "TTree.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "TVector3.h"
#include "TSystem.h"

// C++ includes
#include <vector>
#include <string>
#include <utility>
#include <map>
#include <iostream>
#include <fstream>
//#include <array>

namespace ana {
  class TrackValidation;
}


class ana::TrackValidation : public art::EDAnalyzer {
public:
  explicit TrackValidation(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  TrackValidation(TrackValidation const&) = delete;
  TrackValidation(TrackValidation&&) = delete;
  TrackValidation& operator=(TrackValidation const&) = delete;
  TrackValidation& operator=(TrackValidation&&) = delete;

  // Required functions.
  void analyze(art::Event const& evt) override;

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;

  //
  void initTree(TTree* Tree, std::string branchName, float& Metric);
  void initTree(TTree* Tree, std::string branchName, std::vector<float>& Metric);
  void initTree(TTree* Tree, std::string branchName, std::map<std::string,std::vector<float> >& Metric, std::vector<std::string> fTrackModuleLabel);
  void initTree(TTree* Tree, std::string branchName, std::map<std::string,std::vector<std::vector<float> > >& Metric, std::vector<std::string> fTrackModuleLabel); 
  void initClusterTree(TTree* Tree,  std::string branchName,  std::map<std::string,std::vector<std::vector<std::vector<float> > > >& Metric, std::vector<std::string> fTrackModuleLabels);
  
 
private:

  // Declare member data here.
  
  //fcl parameters
  std::string fGenieGenModuleLabel;
  std::string fLArGeantModuleLabel;
  std::string fHitsModuleLabel;
  //std::string fTrkModuleLabel;
  std::string fCalorimetryLabel;
  //std::string fShowerModuleLabel;
  std::string fPFParticleLabel;
  std::string fpandoraTrackMCS;
  std::string fpandoraTrackRange;
  std::vector<std::string> fTrackModuleLabels;
  std::vector<std::string> fHitModuleLabels;

  bool  fRemoveNonContainedParticles;
  bool  fUseNeutrinoTracks;
  bool  fHitValidation;
  bool  fMatchTracksInTruth;
  int   fVerbose;
  int   fMinHitSize;
  float fSimEnergyCut;
  float fLengthCut;
  float fMaxSimEnergy;
  float fMinRecoEnergy;
  float fMatchedCut;
  bool  fUseMedian;
  std::vector<unsigned int> fTrackPdgCodeVect;
  int fUseBestPlane;
  std::vector<std::string> fPIDnames;
  std::vector<std::string> frangePIDnames;

  // Output tree declaration
  //TTree *fRecoTrkTree;
  //TTree *fMCTrkTree;
  TTree *fTree;
  
  //Service handlesacktracker
  art::ServiceHandle<cheat::BackTrackerService> bt_serv;
  art::ServiceHandle<cheat::ParticleInventoryService> particleInventory;
  art::ServiceHandle<geo::Geometry> geom;
  art::ServiceHandle<art::TFileService> tfs;
  //calo::CalorimetryAlg          fCalorAlg;

  //Declaring Global Variables
  float fEventRun_TreeVal;
  float fEventSubrun_TreeVal;
  float fEventNumber_TreeVal;
  float fNumTrueTrks_TreeVal;//
  float fNumTrueTrksviaECut_TreeVal;//
  float fNumTrueTrksviaDisCut_TreeVal;
  float fNoShowerTrack_TreeVal;         //Number of showers clustered as track.
  //
  std::vector<std::string> fMCTrackList;
  std::map<std::string, std::vector<float> > TrueHitNum_TreeVal;
  std::map<std::string, std::vector<float> > TrueTrkE_TreeValMap;
  std::map<std::string, std::vector<float> > TrueTrkEviaECut_TreeValMap;
  std::map<std::string, std::vector<float> > TrueTrkEviaDisCut_TreeValMap;
  //Other declaration for tree
  std::map<std::string,std::vector<std::vector<float> > > hEnergyCompl_TreeVal;
  
  //PFParticle and Neutrino information
  std::map<std::string,std::vector<float> > PFPNeutrinos_TreeVal;
  std::map<std::string,std::vector<float> > PFPAllTracks_TreeVal;
  std::map<std::string,std::vector<float> > PFPAllShowers_TreeVal;
  std::map<std::string,std::vector<float> > PFPrimaryTracks_TreeVal;
  std::map<std::string,std::vector<float> > PFPrimaryShowers_TreeVal;
 
  //Vertex information
  std::map<std::string,std::vector<float> > PFPVertexDistX_TreeVal;
  std::map<std::string,std::vector<float> > PFPVertexDistY_TreeVal;
  std::map<std::string,std::vector<float> > PFPVertexDistZ_TreeVal;
  std::map<std::string,std::vector<float> > PFPVertexDistMag_TreeVal;

  std::map<std::string,std::vector<float> > PFPLongestTrkLength_TreeVal;
  std::map<std::string,std::vector<float> > PFPLongestTrkHits_TreeVal;
  //Track information
  std::map<std::string,std::vector<float> > trkAllHitsPurity_TreeVal;
  std::map<std::string,std::vector<float> > trkAllHitsComp_TreeVal;
  std::map<std::string,std::vector<float> > trkEnergyPurityBestPlane_TreeVal;
  std::map<std::string,std::vector<float> > trkEnergyCompBestPlane_TreeVal;
  std::map<std::string,std::vector<float> > trkTrackPDG_MostHit_TreeVal;
  std::map<std::string,std::vector<float> > trkPDG_TreeVal;
  std::map<std::string,std::vector<float> > trkMother_TreeVal;
  std::map<std::string,std::vector<float> > trkTrueEnergy_TreeVal;
  std::map<std::string,std::vector<float> > trkTrueKEnergy_TreeVal;

  std::map<std::string,std::vector<float> > trkSignalToNoiseRatio_TreeVal;
  //std::map<std::string,std::vector<float> > trkSignalToNoiseRatioU_TreeVal;
  //std::map<std::string,std::vector<float> > trkSignalToNoiseRatioV_TreeVal;
  //std::map<std::string,std::vector<float> > trkSignalToNoiseRatioY_TreeVal;
  
  std::map<std::string,std::vector<std::string> > trkStartEndProcess_TreeVal;
  std::map<std::string,std::vector<std::string> > trkStartProcess_TreeVal;
  std::map<std::string,std::vector<float> > trkTrueTrackLength_TreeVal;
  std::map<std::string,std::vector<float> > trkTrueTrackSpread_TreeVal;
  std::map<std::string,std::vector<float> > trkTrueTrackResidual_TreeVal;

  //Track  direction information
  std::map<std::string,std::vector<float> > trkDirX_TreeVal;
  std::map<std::string,std::vector<float> > trkDirY_TreeVal;
  std::map<std::string,std::vector<float> > trkDirZ_TreeVal;
  std::map<std::string,std::vector<float> > trkTrueDirX_TreeVal;
  std::map<std::string,std::vector<float> > trkTrueDirY_TreeVal;
  std::map<std::string,std::vector<float> > trkTrueDirZ_TreeVal;
  std::map<std::string,std::vector<float> > trkDirDiff_TreeVal;
  //Track start direction information
  std::map<std::string,std::vector<float> > trkStartDirX_TreeVal;
  std::map<std::string,std::vector<float> > trkStartDirY_TreeVal;
  std::map<std::string,std::vector<float> > trkStartDirZ_TreeVal;
  std::map<std::string,std::vector<float> > trkStartTrueDirX_TreeVal;
  std::map<std::string,std::vector<float> > trkStartTrueDirY_TreeVal;
  std::map<std::string,std::vector<float> > trkStartTrueDirZ_TreeVal;
  std::map<std::string,std::vector<float> > trkStratDirDiff_TreeVal;
  //Track end direction information
  std::map<std::string,std::vector<float> > trkEndDirX_TreeVal;
  std::map<std::string,std::vector<float> > trkEndDirY_TreeVal;
  std::map<std::string,std::vector<float> > trkEndDirZ_TreeVal;
  std::map<std::string,std::vector<float> > trkEndTrueDirX_TreeVal;
  std::map<std::string,std::vector<float> > trkEndTrueDirY_TreeVal;
  std::map<std::string,std::vector<float> > trkEndTrueDirZ_TreeVal;
  std::map<std::string,std::vector<float> > trkEndDirDiff_TreeVal;
  //Start position information
  std::map<std::string,std::vector<float> > trkStartX_TreeVal;
  std::map<std::string,std::vector<float> > trkStartY_TreeVal;
  std::map<std::string,std::vector<float> > trkStartZ_TreeVal;
  std::map<std::string,std::vector<float> > trkStartDist_TreeVal; 
  //
  std::map<std::string,std::vector<float> > trkBestPlane_TreeVal;
  std::map<std::string,std::vector<float> > trkNumHits_TreeVal;
  std::map<std::string,std::vector<float> > trkID_TreeVal;
  std::map<std::string,std::vector<float> > trkParticleId_TreeVal;
  std::map<std::string,std::vector<float> > trkNdof_TreeVal;
  std::map<std::string,std::vector<float> > trkLength_TreeVal;
  std::map<std::string,std::vector<float> > trkLengthDiff_TreeVal;
  std::map<std::string,std::vector<float> > trkGeoProjectionMatched_TreeVal;
  std::map<std::string,std::vector<float> > trkChi2_TreeVal;
  //Track momentum information
  std::map<std::string,std::vector<float> > trkStartMomentum_TreeVal;
  std::map<std::string,std::vector<float> > trkEndMomentum_TreeVal;
  std::map<std::string,std::vector<float> > trkVertexMomentum_TreeVal;
  //Track Angles information
  std::map<std::string,std::vector<float> > trkAzimuthAngle_TreeVal;
  std::map<std::string,std::vector<float> > trkPhi_TreeVal;
  std::map<std::string,std::vector<float> > trkTheta_TreeVal;
  std::map<std::string,std::vector<float> > trkZenithAngle_TreeVal;
  //Track calorimetry information
  std::map<std::string,std::vector<float> > calotrkPlaneId_TreeVal;
  
  std::map<std::string,std::vector<float> > trkKEBP_TreeVal;
  std::map<std::string,std::vector<float> > trkKEU_TreeVal;
  std::map<std::string,std::vector<float> > trkKEV_TreeVal;
  std::map<std::string,std::vector<float> > trkKEY_TreeVal;
  
  std::map<std::string,std::vector<float> > trkRange_TreeVal;
  std::map<std::string,std::vector<float> > trkPitchY_TreeVal;
  
  std::map<std::string,std::vector<std::vector<float> > > trkResRangeBP_TreeVal;
  std::map<std::string,std::vector<std::vector<float> > > trkResRangeU_TreeVal;
  std::map<std::string,std::vector<std::vector<float> > > trkResRangeV_TreeVal;
  std::map<std::string,std::vector<std::vector<float> > > trkResRangeY_TreeVal;
  
  std::map<std::string,std::vector<std::vector<float> > > trkdEdxVecBP_TreeVal;
  std::map<std::string,std::vector<std::vector<float> > > trkdEdxVecU_TreeVal;
  std::map<std::string,std::vector<std::vector<float> > > trkdEdxVecV_TreeVal;
  std::map<std::string,std::vector<std::vector<float> > > trkdEdxVecY_TreeVal;
  
  std::map<std::string,std::vector<std::vector<float> > > trkdQdxVecBP_TreeVal;
  std::map<std::string,std::vector<std::vector<float> > > trkdQdxVecU_TreeVal;
  std::map<std::string,std::vector<std::vector<float> > > trkdQdxVecV_TreeVal;
  std::map<std::string,std::vector<std::vector<float> > > trkdQdxVecY_TreeVal;

  //Track energy information
  std::map<std::string,std::vector<float> > trkKinEnergy_TreeVal;
  std::map<std::string,std::vector<float> > trkKinEnergyRat_TreeVal;
  std::map<std::string,std::vector<float> > trkKinEnergyDiff_TreeVal;
  //
  std::map<std::string,std::vector<float> > recotrk_range_p_TreeVal;
  std::map<std::string,std::vector<float> > recotrk_bestLogLikelihood_p_TreeVal;
  std::map<std::string,std::vector<float> > recotrk_bestMomUncertainty_p_TreeVal;
  std::map<std::string,std::vector<float> > recotrk_bestMomentum_p_TreeVal;
  std::map<std::string,std::vector<float> > recotrk_bwdLogLikelihood_p_TreeVal;
  std::map<std::string,std::vector<float> > recotrk_bwdMomentum_p_TreeVal;
  std::map<std::string,std::vector<float> > recotrk_bwdMomUncertainty_p_TreeVal;
  std::map<std::string,std::vector<float> > recotrk_deltaLogLikelihood_p_TreeVal;
  std::map<std::string,std::vector<float> > recotrk_fwdLogLikelihood_p_TreeVal;
  std::map<std::string,std::vector<float> > recotrk_fwdMomentum_p_TreeVal;
  std::map<std::string,std::vector<float> > recotrk_fwdMomUncertainty_p_TreeVal;
  //
  std::map<std::string,std::vector<float> > trkdEdxBP_TreeVal; //For Best Plane
  std::map<std::string,std::vector<float> > trkdEdxU_TreeVal;
  std::map<std::string,std::vector<float> > trkdEdxV_TreeVal;
  std::map<std::string,std::vector<float> > trkdEdxY_TreeVal;
  //std::map<std::string,std::vector<double> > ;
  //
  std::map<std::string,std::vector<float> > eSegmentation_TreeVal;
  std::map<std::string,std::vector<float> > eNumRecoTracks_TreeVal;

  //
  int fnumevents;
  int fcontainmentCutNum;
  int fnumtracks;
  int fnumtrackspassTPC;
  //int numtrackspassdensity;
  int fnumtrackspassenergy;
  int fnumrecotracks;
  int fnumrecotracksana;
  //int fNoShowerTrack;         //Number of showers clustered as track.
  int fTrueParticleShowers;
  int fTrueParticleTracks;


};


ana::TrackValidation::TrackValidation(fhicl::ParameterSet const& p) : EDAnalyzer(p){
  // More initializers here.
  // Call appropriate consumes<>() for any products to be retrieved by this module.
  fGenieGenModuleLabel = p.get<std::string>("GenieGenModuleLabel" );//
  fLArGeantModuleLabel = p.get<std::string>("LArGeantModuleLabel" );//
  fHitsModuleLabel     = p.get<std::string>("HitsModuleLabel"     );//
  //fTrkModuleLabel      = p.get<std::string>("TrkModuleLabel"      );//
  fCalorimetryLabel    = p.get<std::string>("CalorimetryLabel"    );//
  //fShowerModuleLabel   = p.get<std::string>("ShowerModuleLabel"   );//
  fPFParticleLabel     = p.get<std::string>("PFParticleLabel"     );//
  fpandoraTrackMCS     = p.get<std::string>("pandoraTrackMCS"    );//
  fpandoraTrackRange   = p.get<std::string>("pandoraTrackRange"  );//
  fTrackModuleLabels   = p.get<std::vector<std::string> >("TrackModuleLabels");//
  fHitModuleLabels     = p.get<std::vector<std::string> >("HitModuleLabels"  );//

  fRemoveNonContainedParticles = p.get<bool>("RemoveNonContainedParticles" );//
  fUseNeutrinoTracks           = p.get<bool>("UseNeutrinoTracks"           );//
  fHitValidation               = p.get<bool>("HitValidation"               );//
  fMatchTracksInTruth          = p.get<bool>("MatchTracksInTruth"          );//
  fVerbose                     = p.get<int>("Verbose"                      );//
  fMinHitSize                  = p.get<int>("MinHitSize"                   );//
  fSimEnergyCut                = p.get<float>("SimEnergyCut"               );//
  fLengthCut                   = p.get<float>("LengthCut"                  );//
  fMaxSimEnergy                = p.get<float>("MaxSimEnergy"               );//
  fMinRecoEnergy               = p.get<float>("MinRecoEnergy"              );//
  fMatchedCut                  = p.get<float>("MatchedCut"                 );//
  fUseMedian                   = p.get<bool>("UseMedian"                   );//
  fUseBestPlane                = p.get<int>("UseBestPlane"                 );//
  fTrackPdgCodeVect = p.get<std::vector<unsigned int> >("TrackPdgCodeVect" );//
  fPIDnames         =p.get<std::vector<std::string> >("PIDnames"           );//
  frangePIDnames    =p.get<std::vector<std::string> >("rangePIDnames"      );//


}//Close ParameterSet block

void ana::TrackValidation::initTree(TTree* Tree, std::string branchName, std::vector<float>& Metric){
  const char* brname = branchName.c_str();
  Tree->Branch(brname,"std::vector<float>", &Metric, 32000, 0);
}

void ana::TrackValidation::initTree(TTree* Tree, std::string branchName, float& Metric){
  const char* brname = branchName.c_str();
  Tree->Branch(brname, &Metric, 32000, 0);
}

void ana::TrackValidation::initTree(TTree* Tree, std::string branchName, 
                                    std::map<std::string,std::vector<float> >& Metric, 
                                    std::vector<std::string> fTrackModuleLabels){
  for(auto const& fTrackModuleLabel : fTrackModuleLabels){
    std::string branchString = branchName + "_" + fTrackModuleLabel;
    const char* brname = branchString.c_str();
    Tree->Branch(brname,"std::vector<float>", &Metric[fTrackModuleLabel], 32000, 0);
  }
}

void ana::TrackValidation::initTree(TTree* Tree, std::string branchName, 
                                    std::map<std::string,std::vector<std::vector<float> > >& Metric, 
                                    std::vector<std::string> fTrackModuleLabels){
  for(auto const& fTrackModuleLabel : fTrackModuleLabels){
    std::string branchString = branchName + "_" + fTrackModuleLabel;
    const char* brname = branchString.c_str();
    Tree->Branch(brname,"std::vector<std::vector<float> >", &Metric[fTrackModuleLabel], 32000, 0);
  }
}

void ana::TrackValidation::initClusterTree(TTree* Tree, std::string branchName,
                           std::map<std::string,std::vector<std::vector<std::vector<float> > > >& Metric,
                           std::vector<std::string> fTrackModuleLabels){
  for(auto const& fTrackModuleLabel : fTrackModuleLabels){
    unsigned int numberPlanes = geom->Nplanes();
    Metric[fTrackModuleLabel].resize(numberPlanes);
    for(unsigned int plane=0; plane<numberPlanes; plane++){
      std::string branchString = "Plane_" + std::to_string(plane) + "_" + branchName + "_" + fTrackModuleLabel;
      const char* brname = branchString.c_str();
      Tree->Branch(brname,"std::vector<std::vector<float> >", &Metric[fTrackModuleLabel].at(plane), 32000, 0);
    }
  }
}


void ana::TrackValidation::beginJob()
{
  // Implementation of optional member function here.
  fnumevents             = 0;
  fcontainmentCutNum     = 0;
  fnumtracks             = 0;
  fnumtrackspassTPC      = 0;
  //fnumtrackspassdensity;
  fnumtrackspassenergy   = 0;
  fnumrecotracks         = 0;
  fnumrecotracksana      = 0;
  fTrueParticleShowers   = 0;
  fTrueParticleTracks    = 0;
  fNoShowerTrack_TreeVal = 0;
  fMCTrackList.clear();

  //
  //
  fTree = tfs->make<TTree>("MetricTree", "Tree Holding all metric information");
  gInterpreter->GenerateDictionary("vector<vector<float> > ","vector");

  fTree->Branch("EventRun",    &fEventRun_TreeVal, 32000, 0);
  fTree->Branch("EventSubrun", &fEventSubrun_TreeVal, 32000, 0);
  fTree->Branch("EventNumber", &fEventNumber_TreeVal, 32000, 0);
  fTree->Branch("NoShowerRecoAsTrack", &fNoShowerTrack_TreeVal, 32000, 0);

  initTree(fTree,"TrueTrkE",TrueTrkE_TreeValMap,fMCTrackList);//TODO add it to mc tree and add pdgcode
  initTree(fTree,"TrueTrkEviaECut",TrueTrkEviaECut_TreeValMap,fMCTrackList);
  initTree(fTree,"TrueTrkEviaDisCut",TrueTrkEviaDisCut_TreeValMap,fMCTrackList);

  initTree(fTree,"trkAllHitsPurity",trkAllHitsPurity_TreeVal,fTrackModuleLabels);
  initTree(fTree,"trkAllHitsComp",trkAllHitsComp_TreeVal,fTrackModuleLabels);
  initTree(fTree,"trkEnergyPurityBestPlane",trkEnergyPurityBestPlane_TreeVal,fTrackModuleLabels);
  initTree(fTree,"trkEnergyCompBestPlane",trkEnergyCompBestPlane_TreeVal,fTrackModuleLabels);

  initTree(fTree,"trkTrackPDG_MostHit",trkTrackPDG_MostHit_TreeVal,fTrackModuleLabels);
  initTree(fTree,"trkPDG",trkPDG_TreeVal,fTrackModuleLabels);
  initTree(fTree,"trkMother",trkMother_TreeVal,fTrackModuleLabels);
  initTree(fTree,"trkTrueEnergy",trkTrueEnergy_TreeVal,fTrackModuleLabels);
  initTree(fTree,"trkTrueKinEnergy",trkTrueKEnergy_TreeVal,fTrackModuleLabels);

  initTree(fTree,"trkSignalToNoiseRatioBP",trkSignalToNoiseRatio_TreeVal,fTrackModuleLabels);
  //initTree(fTree,"trkSignalToNoiseRatioU",trkSignalToNoiseRatioU_TreeVal,fTrackModuleLabels);
  //initTree(fTree,"trkSignalToNoiseRatioV",trkSignalToNoiseRatioV_TreeVal,fTrackModuleLabels);
  //initTree(fTree,"trkSignalToNoiseRatioY",trkSignalToNoiseRatioY_TreeVal,fTrackModuleLabels);

  initTree(fTree,"trkTrueTrackLength",trkTrueTrackLength_TreeVal,fTrackModuleLabels);
  initTree(fTree,"trkTrueTrackSpread",trkTrueTrackSpread_TreeVal,fTrackModuleLabels);
  initTree(fTree,"trkTrueTrackResidual",trkTrueTrackResidual_TreeVal,fTrackModuleLabels);
  initTree(fTree,"trkDirX",trkDirX_TreeVal,fTrackModuleLabels);
  initTree(fTree,"trkDirY",trkDirY_TreeVal,fTrackModuleLabels);
  initTree(fTree,"trkDirZ",trkDirZ_TreeVal,fTrackModuleLabels);
  initTree(fTree,"trkTrueDirX",trkTrueDirX_TreeVal,fTrackModuleLabels);
  initTree(fTree,"trkTrueDirY",trkTrueDirY_TreeVal,fTrackModuleLabels);
  initTree(fTree,"trkTrueDirZ",trkTrueDirZ_TreeVal,fTrackModuleLabels);
  initTree(fTree,"trkDirDiff",trkDirDiff_TreeVal,fTrackModuleLabels);
  
  initTree(fTree,"trkStartDirX",trkStartDirX_TreeVal,fTrackModuleLabels);
  initTree(fTree,"trkStartDirY",trkStartDirY_TreeVal,fTrackModuleLabels);
  initTree(fTree,"trkStartDirZ",trkStartDirZ_TreeVal,fTrackModuleLabels);
  initTree(fTree,"trkStartTrueDirX",trkStartTrueDirX_TreeVal,fTrackModuleLabels);
  initTree(fTree,"trkStartTrueDirY",trkStartTrueDirY_TreeVal,fTrackModuleLabels);
  initTree(fTree,"trkStartTrueDirZ",trkStartTrueDirZ_TreeVal,fTrackModuleLabels);
  initTree(fTree,"trkStratDirDiff",trkStratDirDiff_TreeVal,fTrackModuleLabels);
  
  initTree(fTree,"trkEndDirX",trkEndDirX_TreeVal,fTrackModuleLabels);
  initTree(fTree,"trkEndDirY",trkEndDirY_TreeVal,fTrackModuleLabels);
  initTree(fTree,"trkEndDirZ",trkEndDirZ_TreeVal,fTrackModuleLabels);
  initTree(fTree,"trkEndTrueDirX",trkEndTrueDirX_TreeVal,fTrackModuleLabels);
  initTree(fTree,"trkEndTrueDirY",trkEndTrueDirY_TreeVal,fTrackModuleLabels);
  initTree(fTree,"trkEndTrueDirZ",trkEndTrueDirZ_TreeVal,fTrackModuleLabels);
  initTree(fTree,"trkEndDirDiff",trkEndDirDiff_TreeVal,fTrackModuleLabels);
  
  initTree(fTree,"trkStartX",trkStartX_TreeVal,fTrackModuleLabels);
  initTree(fTree,"trkStartY",trkStartY_TreeVal,fTrackModuleLabels);
  initTree(fTree,"trkStartZ",trkStartZ_TreeVal,fTrackModuleLabels);
  initTree(fTree,"trkStartDist",trkStartDist_TreeVal,fTrackModuleLabels);
  
  initTree(fTree,"trkBestPlane",trkBestPlane_TreeVal,fTrackModuleLabels);
  initTree(fTree,"trkNumHits",trkNumHits_TreeVal,fTrackModuleLabels);
  initTree(fTree,"trkID",trkID_TreeVal,fTrackModuleLabels);
  initTree(fTree,"trkParticleId",trkParticleId_TreeVal,fTrackModuleLabels);
  initTree(fTree,"trkNdof",trkNdof_TreeVal,fTrackModuleLabels);
  initTree(fTree,"trkLength",trkLength_TreeVal,fTrackModuleLabels);
  initTree(fTree,"trkLengthDiff",trkLengthDiff_TreeVal,fTrackModuleLabels);
  initTree(fTree,"trkGeoProjectionMatched",trkGeoProjectionMatched_TreeVal,fTrackModuleLabels);
  initTree(fTree,"trkChi2",trkChi2_TreeVal,fTrackModuleLabels);
  
  initTree(fTree,"trkStartMomentum",trkStartMomentum_TreeVal,fTrackModuleLabels);
  initTree(fTree,"trkEndMomentum",trkEndMomentum_TreeVal,fTrackModuleLabels);
  initTree(fTree,"trkVertexMomentum",trkVertexMomentum_TreeVal,fTrackModuleLabels);
  
  initTree(fTree,"trkAzimuthAngle",trkAzimuthAngle_TreeVal,fTrackModuleLabels);
  initTree(fTree,"trkPhi",trkPhi_TreeVal,fTrackModuleLabels);
  initTree(fTree,"trkTheta",trkTheta_TreeVal,fTrackModuleLabels);
  initTree(fTree,"trkZenithAngle",trkZenithAngle_TreeVal,fTrackModuleLabels);
  
  initTree(fTree,"calotrkPlaneId",calotrkPlaneId_TreeVal,fTrackModuleLabels);
  
  initTree(fTree,"trkKEBP",trkKEBP_TreeVal,fTrackModuleLabels);
  initTree(fTree,"trkKEU",trkKEU_TreeVal,fTrackModuleLabels);
  initTree(fTree,"trkKEV",trkKEV_TreeVal,fTrackModuleLabels);
  initTree(fTree,"trkKEY",trkKEY_TreeVal,fTrackModuleLabels);

  initTree(fTree,"trkRange",trkRange_TreeVal,fTrackModuleLabels);
  initTree(fTree,"trkPitchY",trkPitchY_TreeVal,fTrackModuleLabels);
  
  initTree(fTree,"trkResRangeBP",trkResRangeBP_TreeVal,fTrackModuleLabels);
  initTree(fTree,"trkResRangeU",trkResRangeU_TreeVal,fTrackModuleLabels);
  initTree(fTree,"trkResRangeV",trkResRangeV_TreeVal,fTrackModuleLabels);
  initTree(fTree,"trkResRangeY",trkResRangeY_TreeVal,fTrackModuleLabels);

  initTree(fTree,"trkdEdxVecBP",trkdEdxVecBP_TreeVal,fTrackModuleLabels);
  initTree(fTree,"trkdEdxVecU",trkdEdxVecU_TreeVal,fTrackModuleLabels);
  initTree(fTree,"trkdEdxVecV",trkdEdxVecV_TreeVal,fTrackModuleLabels);
  initTree(fTree,"trkdEdxVecY",trkdEdxVecY_TreeVal,fTrackModuleLabels);

  initTree(fTree,"trkdQdxVecBP",trkdQdxVecBP_TreeVal,fTrackModuleLabels);
  initTree(fTree,"trkdQdxVecU",trkdQdxVecU_TreeVal,fTrackModuleLabels);
  initTree(fTree,"trkdQdxVecV",trkdQdxVecV_TreeVal,fTrackModuleLabels);
  initTree(fTree,"trkdQdxVecY",trkdQdxVecY_TreeVal,fTrackModuleLabels);
  
  initTree(fTree,"recoTrkKinEnergy",trkKinEnergy_TreeVal,fTrackModuleLabels);
  initTree(fTree,"recoTrkKinEnergyRat",trkKinEnergyRat_TreeVal,fTrackModuleLabels);
  initTree(fTree,"recoTrkKinEnergyDiff",trkKinEnergyDiff_TreeVal,fTrackModuleLabels);

  initTree(fTree,"recotrk_range_p",recotrk_range_p_TreeVal,fTrackModuleLabels);
  initTree(fTree,"recotrk_bestLogLikelihood_p",recotrk_bestLogLikelihood_p_TreeVal,fTrackModuleLabels);
  initTree(fTree,"recotrk_bestMomUncertainty_p",recotrk_bestMomUncertainty_p_TreeVal,fTrackModuleLabels);
  initTree(fTree,"recotrk_bestMomentum_p",recotrk_bestMomentum_p_TreeVal,fTrackModuleLabels);
  initTree(fTree,"recotrk_bwdLogLikelihood_p",recotrk_bwdLogLikelihood_p_TreeVal,fTrackModuleLabels);
  initTree(fTree,"recotrk_bwdMomentum_p",recotrk_bwdMomentum_p_TreeVal,fTrackModuleLabels);
  initTree(fTree,"recotrk_bwdMomUncertainty_p",recotrk_bwdMomUncertainty_p_TreeVal,fTrackModuleLabels);
  initTree(fTree,"recotrk_deltaLogLikelihood_p",recotrk_deltaLogLikelihood_p_TreeVal,fTrackModuleLabels);
  initTree(fTree,"recotrk_fwdLogLikelihood_p",recotrk_fwdLogLikelihood_p_TreeVal,fTrackModuleLabels);
  initTree(fTree,"recotrk_fwdMomentum_p",recotrk_fwdMomentum_p_TreeVal,fTrackModuleLabels);
  initTree(fTree,"recotrk_fwdMomUncertainty_p",recotrk_fwdMomUncertainty_p_TreeVal,fTrackModuleLabels);
  
  initTree(fTree,"trkdEdx_bp",trkdEdxBP_TreeVal,fTrackModuleLabels);
  initTree(fTree,"trkdEdx_u",trkdEdxU_TreeVal,fTrackModuleLabels);
  initTree(fTree,"trkdEdx_v",trkdEdxV_TreeVal,fTrackModuleLabels);
  initTree(fTree,"trkdEdx_y",trkdEdxY_TreeVal,fTrackModuleLabels);
  
  initTree(fTree,"PFPNeutrinos",PFPNeutrinos_TreeVal,fTrackModuleLabels);
  initTree(fTree,"PFPAllTracks",PFPAllTracks_TreeVal,fTrackModuleLabels);
  initTree(fTree,"PFPAllShowers",PFPAllShowers_TreeVal,fTrackModuleLabels);
  initTree(fTree,"PFPrimaryTracks",PFPrimaryTracks_TreeVal,fTrackModuleLabels);
  initTree(fTree,"PFPrimaryShowers",PFPrimaryShowers_TreeVal,fTrackModuleLabels);
  initTree(fTree,"PFPVertexDistX",PFPVertexDistX_TreeVal,fTrackModuleLabels);
  initTree(fTree,"PFPVertexDistY",PFPVertexDistY_TreeVal,fTrackModuleLabels);
  initTree(fTree,"PFPVertexDistZ",PFPVertexDistZ_TreeVal,fTrackModuleLabels);
  initTree(fTree,"PFPVertexDistMag",PFPVertexDistMag_TreeVal,fTrackModuleLabels);
  initTree(fTree,"PFPLongestTrkLength",PFPLongestTrkLength_TreeVal,fTrackModuleLabels);
  initTree(fTree,"PFPLongestTrkHits",PFPLongestTrkHits_TreeVal,fTrackModuleLabels);
  //initTree(fTree,"",,fTrackModuleLabels);  eLongestTrkLength_TreeVal eLongestTrkHits_TreeVal
  
  initTree(fTree,"eSegmentation",eSegmentation_TreeVal,fTrackModuleLabels);
  initTree(fTree,"eNumRecoTracks",eNumRecoTracks_TreeVal,fTrackModuleLabels);
 
  initTree(fTree,"fNumTrueTrks",fNumTrueTrks_TreeVal);
  initTree(fTree,"fNumTrueTrksviaECut",fNumTrueTrksviaECut_TreeVal);
  initTree(fTree,"fNumTrueTrksviaDisCut",fNumTrueTrksviaDisCut_TreeVal);

  for(auto const& fHitModuleLabel : fHitModuleLabels){
    hEnergyCompl_TreeVal[fHitModuleLabel].resize(geom->Nplanes());

    std::string hitString = "hEnergyComp_" + fHitModuleLabel;
    const char* hitChar   = hitString.c_str();
    fTree->Branch(hitChar,"std::vector<std::vector<float> > >", &hEnergyCompl_TreeVal[fHitModuleLabel], 32000, 0);
  }

  for(auto const& fTrackModuleLabel : fTrackModuleLabels){
    
    std::string processString = "trkStartEndProcess_" + fTrackModuleLabel;
    const char* processChar   = processString.c_str();
    
    fTree->Branch(processChar,"<std::string,std::vector<std::string>", &trkStartEndProcess_TreeVal[fTrackModuleLabel], 32000, 0);

    std::string sprocessString = "trkStartProcess_" + fTrackModuleLabel;
    const char* sprocessChar   = sprocessString.c_str();

    fTree->Branch(sprocessChar,"<std::string,std::vector<std::string>", &trkStartProcess_TreeVal[fTrackModuleLabel], 32000, 0);
  }

  fTree->Print();

}//Close beginJob block

void ana::TrackValidation::analyze(art::Event const& evt)
{
  // Implementation of required member function here.
  //###### Declaration and Initiate variables #######
  //
  //services 
  auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(evt);
  //Access Event Information
  fEventRun_TreeVal = evt.run();
  fEventSubrun_TreeVal = evt.subRun();
  fEventNumber_TreeVal = evt.event();
  //int fEventID = evt.id().event(); //it gives same output to the evt.event()
  
  ++fnumevents;

  //-------------------------- Getting the SimWire information -------------------------- 
  //Get the SimChannels so that we can find the energy deposited on them.
  art::Handle<std::vector<sim::SimChannel> > simChannelHandle;
  std::vector<art::Ptr<sim::SimChannel> > simchannels;
  if(evt.getByLabel(fLArGeantModuleLabel,simChannelHandle))// Make sure the handle is valid
  {art::fill_ptr_vector(simchannels, simChannelHandle);}   // Fill the vector with art::Ptr hits

  //-------------------------- Getting all reconstructed hits --------------------------
  art::Handle<std::vector<recob::Hit> > hitListHandle;
  std::vector<art::Ptr<recob::Hit> > allHits;
  if(evt.getByLabel(fHitsModuleLabel, hitListHandle)) // Make sure the handle is valid
  {art::fill_ptr_vector(allHits, hitListHandle);}   // Fill the vector with art::Ptr hits
  if(fVerbose > 1){
    const size_t NallHits = allHits.size();
    std::cout << "The all Hits vector size = " << NallHits << std::endl;
    if(NallHits ==0){std::cout << "There is no Hits in this event" << std::endl;}
  }

  //Get other information by getManyByType which initalises the handles giving every particle product id.
  //Get all the hits
  std::vector<art::Handle<std::vector<recob::Hit> > > hitHandles;
  evt.getManyByType(hitHandles);

  //Get all the clusters
  std::vector<art::Handle<std::vector<recob::Cluster> > > clusterHandles;
  evt.getManyByType(clusterHandles);

  //Get all the pfparticles
  std::vector<art::Handle<std::vector<recob::PFParticle> > > pfpHandles;
  evt.getManyByType(pfpHandles);

  //
  ////////////////////////////////////////////////////////////////////////////////////////////////
  ////                      Getting the truth information from MCParticles                       //
  ////////////////////////////////////////////////////////////////////////////////////////////////
  //###### Declaration and Initiate variables #######
  //Create a map of trackIDs to MCParticles
  std::map<int, const simb::MCParticle*> trueparticlesMap;// Map of all true mc particles
  std::map<int, const simb::MCParticle*> truePrimaryShowerMap;//Map of primary showr-like particles
  std::map<int, const simb::MCParticle*> truePrimaryTrackMap;//Map of primary track-like particles
  std::vector<int> EventTrksTrackIdVec;   //Vector of trackids for track particels
  std::map<int, float> MCTrackEnergyMap;  //Map of the MC Energy deposited for each MC track
  std::map<art::ProductID,std::map<int,std::map<geo::PlaneID, int> > > MCTrack_hit_map; //
  //Variables
  int numTrueTraks = 0;
  float trueTrkEnergy = -999;
  const sim::ParticleList& particles = particleInventory->ParticleList();
  //
  for(auto const& particleIter : particles){
    const simb::MCParticle *particle = particleIter.second;
    //Fill the map with track id and particles
    trueparticlesMap[particle->TrackId()] = particle;
    //
    if(fVerbose > 2){std::cout << "True Particle with track ID: " << particle->TrackId() << " Has code of: " << particle->PdgCode() << " and Energy of: " << particle->E() << " With Mother: " << particle->Mother() << " Proccess: " << particle->Process() << " End Process: "  << particle->EndProcess() << std::endl;}
    //Find the shower particle and fill the map with track id and particles
    if(particles.IsPrimary(particle->TrackId())){
      int ParticlePdgCode = std::abs(particle->PdgCode());
      if(ParticlePdgCode == 11 || ParticlePdgCode == 22){
        truePrimaryShowerMap[particle->TrackId()] = particle;
      }
      int number_of_tracks = 0;
      if(std::find(fTrackPdgCodeVect.begin(),fTrackPdgCodeVect.end(), ParticlePdgCode) != fTrackPdgCodeVect.end()){
        truePrimaryTrackMap[particle->TrackId()] = particle;
        ++number_of_tracks;
      }//
    }
  }// close the loop over the ParticleList

  if(fMatchTracksInTruth){//open truth information block
    EventTrksTrackIdVec = TrackUtils::EventMCTrkIdsFun(evt, trueparticlesMap, fTrackPdgCodeVect);
    //Getting the number of paricle shower-like and track-like per event
    fTrueParticleShowers = truePrimaryShowerMap.size();
    fTrueParticleTracks += truePrimaryTrackMap.size();
    //Get the MC Energy deposited for each MC track.
    MCTrackEnergyMap = RecoUtils::TrueEnergyDepositedFromMCTracks(simchannels);
    //std::cout<< " MCTrackEnergyMap: "<< MCTrackEnergyMap[1] <<std::endl;
    //
    //Find number of true tracks per events
    fNumTrueTrks_TreeVal = EventTrksTrackIdVec.size();
    //Loop over Tracks mother and get their energy
    for(auto const& mctrks : EventTrksTrackIdVec){
      const simb::MCParticle* trackMCParticle = trueparticlesMap[(mctrks)];
      //Get a vector of track mothers 
      fMCTrackList.push_back((std::to_string(mctrks)));
      //Save the energy of tracks with their track 
      TrueTrkE_TreeValMap[(std::to_string(mctrks))].push_back(trackMCParticle->E());
    }//close the loop over trkmothers map
    if(fVerbose > 1){
      std::cout<<"Inital track candidates: "<< EventTrksTrackIdVec.size() <<std::endl;
    }//Close if statement
    if(fRemoveNonContainedParticles){
      TrackUtils::RemoveNoneContainedParticles(EventTrksTrackIdVec, trueparticlesMap, MCTrackEnergyMap);
    }//Close if statement 
    if(fVerbose > 1){
      std::cout<<" Track candidates after containment cut are: "<<EventTrksTrackIdVec.size()<<std::endl;
    }//Close if statement

    //Cut on the Track mothers with energy
    TrackUtils::CutTrkMothersByE(EventTrksTrackIdVec, trueparticlesMap, fSimEnergyCut);
    fNumTrueTrksviaECut_TreeVal = EventTrksTrackIdVec.size();
    for(auto const& mctrks : EventTrksTrackIdVec){
      const simb::MCParticle* TracksParticle = trueparticlesMap.at(mctrks);
      TrueTrkEviaECut_TreeValMap[(std::to_string(mctrks))].push_back(TracksParticle->E());
    }//close the loop over trkmothers after energy cut
    if(fVerbose > 1){
        std::cout<<" Track candidates after E cut are: "<<EventTrksTrackIdVec.size()<<std::endl;
    }//Close if statement
    //Cut on track mothers with length and uncontained tracks
    TrackUtils::CalculateTrackLengthCut(EventTrksTrackIdVec, trueparticlesMap, fLengthCut);
    //
    fNumTrueTrksviaDisCut_TreeVal = EventTrksTrackIdVec.size();
    for(auto const& mctrks : EventTrksTrackIdVec){
      const simb::MCParticle* TracksParticle = trueparticlesMap.at(mctrks);
      //Select first traj point where the track loses energy, last be default
      ////Get the number of Traj points to loop over
      unsigned int MCNumTrajPoints     = TracksParticle->NumberTrajectoryPoints();
      TLorentzVector PositionTrajStart = TracksParticle->Position(0);
      TLorentzVector PositionTrajEnd   = TracksParticle->Position(MCNumTrajPoints-1);
      //Get the track length of the track.
      double MCTrkLength = (PositionTrajStart - PositionTrajEnd).Vect().Mag(); //in cm
      //
      TrueTrkEviaDisCut_TreeValMap[(std::to_string(mctrks))].push_back(MCTrkLength);
    }
    if(fVerbose > 1){
      std::cout<< "Track candidates after length cut are: " << EventTrksTrackIdVec.size() <<std::endl;
    }

    //If there are no true tracks we can't validate
    if(EventTrksTrackIdVec.size() == 0){
      if(fVerbose > 0){std::cout << "No track; finishing" << std::endl;}
      return;
    }
    //
    if(fVerbose > 2){
      for(auto const& mctrks : EventTrksTrackIdVec){
        const simb::MCParticle* TracksMCParticle = trueparticlesMap.at(mctrks);//test
        std::cout << " A Track has track id: " << mctrks << " with Energy: "<< TracksMCParticle->E() << std::endl;
      }//close the loop over tracks mothers 
    }//close fVerbose > 2 

    numTrueTraks = EventTrksTrackIdVec.size();
    fnumtracks += numTrueTraks;

    //Get the number of hits associated wit each Track. This is done for every hits handle in the event.
    for(auto& handle : hitHandles){
      if(!handle.isValid()){
        mf::LogError("TrackValidation") << "Bad hit handle from the all the hit handles" << std::endl;
        continue;
      }//close if statement
      //Getting the Hit Information
      art::Handle<std::vector<recob::Hit> > hitHandle;
      std::vector<art::Ptr<recob::Hit> > hits_fromhandle;
      if(evt.getByLabel(handle.provenance()->moduleLabel(),hitHandle))
      {art::fill_ptr_vector(hits_fromhandle, hitHandle);}


      //Get a map of the number of hits per plane each track has deposited.
      MCTrack_hit_map[handle.id()] = RecoUtils::NumberofPlaneHitsPerTrack(clockData, hits_fromhandle);
    }//close the loop over hit handle
  }//close true information block

  //======================
  //=== Hit Validation ===
  //======================
  
  if(fHitValidation){
    //Loop over the hit handles
    for(auto const& fHitModuleLabel : fHitModuleLabels){
      
      //Getting the Hit Information
      art::Handle<std::vector<recob::Hit> > hitValHandle;
      std::vector<art::Ptr<recob::Hit> > hits_fromValhandle;
      if(evt.getByLabel(fHitModuleLabel, hitValHandle))
      {art::fill_ptr_vector(hits_fromValhandle, hitValHandle);}

      //Calculate the total energy deposited
      float TotalEnergyDeposited = 0.f;
      for(auto const& track_iter: MCTrackEnergyMap){
        TotalEnergyDeposited += track_iter.second;
      }//close the loop over MCTrackEnergyMap

      //Calculate how much energy was deposited in the hits
      for(geo::PlaneID plane_id: geom->IteratePlaneIDs()){

        float TotalEnergyDepinHits = RecoUtils::TotalEnergyDepinHits(clockData, allHits, plane_id.Plane);
        float hitCompleteness = TotalEnergyDepinHits!=0 ? TotalEnergyDepinHits/TotalEnergyDeposited : -99999;
        hEnergyCompl_TreeVal[fHitModuleLabel][plane_id.Plane].push_back(hitCompleteness);

      }//close loop over planes

    }//close the loop over Hit Module
    if(fVerbose > 1){std::cout << "Hit Validation Complete" << std::endl;}
  }//Close the hit validation statement

  //====================================================
  //=== Get the reconstructed info and Match with MC ===
  //====================================================

  for(auto const& fTrackModuleLabel : fTrackModuleLabels){

    if(EventTrksTrackIdVec.size() == 0 && fMatchTracksInTruth){
      if(fVerbose > 0){std::cout << "No track to Match to the true tracks" << std::endl;}
      continue;
    }//Close if statement

    //Getting the Track Information
    art::Handle<std::vector<recob::Track> > trackListHandle;
    std::vector<art::Ptr<recob::Track> > recoTracks;
    if(evt.getByLabel(fTrackModuleLabel,trackListHandle))
    {art::fill_ptr_vector(recoTracks, trackListHandle);}

    if(recoTracks.size() == 0){
      if(fVerbose > 0){std::cout << "No Track in the Event" << std::endl;}
      continue;
    }//close fVerbose statement

    art::Handle<std::vector<recob::PFParticle> > pfpListHandle;
    std::vector<art::Ptr<recob::PFParticle> > pfps;
    if(evt.getByLabel(fPFParticleLabel,pfpListHandle))
    {art::fill_ptr_vector(pfps,pfpListHandle);}

    //Getting the Track Information
    // ----------------- Association bettween: --- Tracks  <---> 2d Hits --------
    art::FindManyP<recob::Hit> fmtrkh(trackListHandle, evt, fTrackModuleLabel);

    // ----------------- Association bettween: --- SpacdPoints  <---> Tracks ----
    art::FindManyP<recob::SpacePoint> fmsptrk(trackListHandle, evt, fTrackModuleLabel);

    // ----------------- Association bettween: --- Tracks  <---> Clusters -------
    art::FindManyP<recob::Cluster> fmtrkch(trackListHandle, evt, fTrackModuleLabel);

    // ----------------- Association bettween: --- PFParticle  <---> Clusters ---
    art::FindManyP<recob::Cluster> fmpfc(pfpListHandle, evt, fPFParticleLabel);

    // ----------------- Association bettween: --- PFParticle  <---> Vertex -----
    art::FindManyP<recob::Vertex> fmpfv(pfpListHandle, evt, fPFParticleLabel);

    // ----------------- Association bettween: --- Tracks  <---> PFParticle -----
    art::FindManyP<recob::PFParticle> fmtrkpf(trackListHandle, evt, fTrackModuleLabel);

    // ----------------- Association bettween: --- Tracks  <---> Calorimetery ---
    art::FindManyP<anab::Calorimetry> fmcalotrk(trackListHandle, evt, fCalorimetryLabel);

    //======================
    //=== PFP Validation ===
    //======================

    std::vector< art::Ptr<recob::Vertex> > neutrinoVertices;
    std::map<int, art::Ptr<recob::PFParticle> > pfpsMap;
    for (auto const& pfp: pfps){
      pfpsMap[pfp->Self()] = pfp;
    }//close the loop over PFParticle


    if(fVerbose > 1){std::cout<<"Doing PFP Validation"<<std::endl;}

    std::vector<int> PFPrimaries;
    //Create a map between PFParticles and their IDs
    int pfpAllTrackCounter  = 0;
    int pfpAllShowerCounter = 0;
    // Only the primary pfps (Daughters of the pfp netrinos)
    int pfpTrackCounter     = 0;
    int pfpShowerCounter    = 0;
    int pfpNeutrinoCounter  = 0;

    if(fVerbose > 1){std::cout<<"Number of PFP: "<<pfps.size()<<std::endl;}

    for (auto const& pfp: pfps){
      if(fVerbose > 2){
        std::cout<<"PFParticle: "<<pfp->Self()<<" with: PDG: "<<pfp->PdgCode()
          <<" Parent: "<<pfp->Parent()<<std::endl;

      }//close fVerbose statement
      if(pfp->PdgCode()==11) {++pfpAllShowerCounter;}
      if(pfp->PdgCode()==13) {++pfpAllTrackCounter;}
      if((pfp->PdgCode()==12) || (pfp->PdgCode()==14)){ //Find Neutrino and primary daughters
        // Gives a vector of Duaghter particle ID's, get Daughters using PFParticlesMap
        const std::vector<size_t> Daughters = pfp->Daughters();
        for (auto const& daughterIter : Daughters){
          const art::Ptr<recob::PFParticle>& daughter = pfpsMap[daughterIter];

          // Split into shower like, 11, and track like, 13.
          if(daughter->PdgCode() == 11){
            ++pfpShowerCounter;
          }else if(daughter->PdgCode() == 13){
            ++pfpTrackCounter;
          }else{
            std::cout<<"Something has gone horribly wrong, PFP PDG != 11||13"<<std::endl;
            return;
          }//Close else statement 
          PFPrimaries.push_back(daughter->Self());
        }//Close the loop over the neutrino daughters
      
        //Get the PFP neutrino vertex
        if(fmpfv.isValid() && fmpfv.size()!=0 && fmpfv.at(0).size()!=0){
          art::Handle<std::vector<recob::Vertex > > vertexHandle;
          evt.get(fmpfv.at(0).front().id(),vertexHandle);

          if(vertexHandle.isValid()){
            std::vector< art::Ptr<recob::Vertex> > pfpVertexVector = fmpfv.at(pfp.key());
            art::Ptr<recob::Vertex> pfpVertex = pfpVertexVector.at(0);
            neutrinoVertices.push_back(pfpVertex);
          }//Close if statement that check the vertex handle
        }//close getting the pfp neutrino vertex
        ++pfpNeutrinoCounter;
      }//close the looping over pfps
    }//Close the if statement for finding neutrino and primary daughters
  
    if(fVerbose > 1){std::cout<<"Primary Tracks: "<<pfpTrackCounter
      <<" and Primary Showers: "<<pfpShowerCounter<<std::endl;}
  //
    PFPNeutrinos_TreeVal[fTrackModuleLabel].push_back(pfpNeutrinoCounter);
    PFPAllTracks_TreeVal[fTrackModuleLabel].push_back(pfpAllTrackCounter);
    PFPAllShowers_TreeVal[fTrackModuleLabel].push_back(pfpAllShowerCounter);
    PFPrimaryTracks_TreeVal[fTrackModuleLabel].push_back(pfpTrackCounter); 
    PFPrimaryShowers_TreeVal[fTrackModuleLabel].push_back(pfpShowerCounter);
    if(fVerbose > 1){std::cout<<"Number of PFP Neutrinos:"<<pfpNeutrinoCounter
      <<" and vertex: "<<neutrinoVertices.size()<<std::endl;}
    
    //Looping over PFParticle vertices
    for(auto const& pfpVertex : neutrinoVertices){

      //get the pfp neutrino vertex
      double pfpvtx[3];
      pfpVertex->XYZ(pfpvtx);

      //Take the first primary track start position as the true start vertex
      std::map<int,const simb::MCParticle*>::iterator primaryTrkParticleIt = truePrimaryTrackMap.begin();
      const simb::MCParticle* primaryTrkParticle = (*primaryTrkParticleIt).second;

      const TLorentzVector PositionTrajP = primaryTrkParticle->Position();
      double truevtx[3] = {PositionTrajP.X() ,PositionTrajP.Y(), PositionTrajP.Z()};

      double PFP_Start_Diff_X = pfpvtx[0] - truevtx[0];
      double PFP_Start_Diff_Y = pfpvtx[1] - truevtx[1];
      double PFP_Start_Diff_Z = pfpvtx[2] - truevtx[2];

      double PFP_Start_Diff = std::sqrt((std::pow((PFP_Start_Diff_X),2))+(std::pow((PFP_Start_Diff_Y),2))+(std::pow((PFP_Start_Diff_Z),2)));
      //
      PFPVertexDistX_TreeVal[fTrackModuleLabel].push_back(PFP_Start_Diff_X);
      PFPVertexDistY_TreeVal[fTrackModuleLabel].push_back(PFP_Start_Diff_Y);
      PFPVertexDistZ_TreeVal[fTrackModuleLabel].push_back(PFP_Start_Diff_Z);
      PFPVertexDistMag_TreeVal[fTrackModuleLabel].push_back(PFP_Start_Diff);

    }//Close the loop over pfp vertex

    //Check if the hits associated to track valid
    if(fmtrkh.size() == 0 || fmtrkh.at(0).size() == 0){
      std::cout << " No hits in a recob Track. Association is made incorrectly. Bailing" << std::endl;
      continue;
    }//close checking the hit associated to tracks

    //Get the ID of the track hit module
    art::ProductID trkhit_productid = fmtrkh.at(0).front().id();

    //--------------------------------------------------------------------------------------------
    unsigned int max_hitnum       = 0;
    unsigned int longestTrkNoHits = 0;
    float        longestrk        = 0.0;
    art::Ptr<recob::Track> longestTrack;
    art::Ptr<recob::Track> MostHitsTrack;

    //Find the longest track (by number of hits)
    for(auto const& trk : recoTracks){
      ++fnumrecotracks;
      double TrkLength = std::abs(trk->Length());
      //Get the hits vector from the track
      std::vector<art::Ptr<recob::Hit> > trkhits = fmtrkh.at(trk.key());
      if(trkhits.size() > max_hitnum){
        max_hitnum = trkhits.size();
        MostHitsTrack = trk;
      }//close if statement
      if(TrkLength > longestrk){
        longestrk = TrkLength;
        longestTrkNoHits = trkhits.size();
        longestTrack = trk;
      }//Close if statement
    }//close the loop over tracks
    PFPLongestTrkLength_TreeVal[fTrackModuleLabel].push_back(longestrk);
    PFPLongestTrkHits_TreeVal[fTrackModuleLabel].push_back(longestTrkNoHits);
    //--------------------------------------------------------------------------------------------

    //Loop over tracks to get reco information
    for(auto const& trk : recoTracks){
      //Get the information for the track
      //Getting the Azimuth Angle (The angle is measured on the plane orthogonal to the z axis,
      //The angle is measured counterclockwise from the x axis) basiclly the angle between x and z.
      double trkAzimuthAngle = trk->AzimuthAngle();
      //Getting the Chi^2 
      float trkChi2 = trk->Chi2();
      //Track start position 
      TVector3 trkStartPosition = trk->Start<TVector3>();
      //Getting the track end position
      TVector3 trkEndPosition = trk->End<TVector3>();
      //getting track start direction
      TVector3 trkstartdir = trk->StartDirection<TVector3>();
      //Getting the track end directin using EndDirection function
      TVector3 trkenddir = trk->EndDirection<TVector3>();
      //Getting Track Start Momentum
      double trkstartmom = trk->StartMomentum();
      //Access to track end momentum.
      double trkEndMomentum = trk->EndMomentum();
      //Getting Track Start Momentum Vector
      TVector3 trkstartmomVec = trk->StartMomentumVector<TVector3>();
      //Momentum vector at end point.
      TVector3 trkEndMomentumVec = trk->EndMomentumVector<TVector3>();
      //Get the track Id 
      int trkID = trk->ID();
      //Getting length of the track
      double trkLength = trk->Length();
      //Access the track postion at different points
      //TVector3 trkposition = trk->LocationAtPoint<TVector3>(p);
      //getting the number of degrees of freedom of the fit
      int trkNdof = trk->Ndof();
      //Getting the number of valid trajectory points
      //size_t trkNPoints = trk->NPoints();
      //Getting number of trajectory points
      //size_t trkNoTrajPoints = trk->NPoints();
      //Getting ParticleId 
      int trkParticleId = trk->ParticleId();
      //Access the phi angle for track (the angle to y axis)
      double trkPhi = trk->Phi();
      //Access the Theta angle for track (the angle between z and x axises respect to z axis)
      double trkTheta = trk->Theta();
      //Access the zenith angle for track (the zenith angle is defined by the angle between the track starting
      // direction and the y axis) The y component of the starting direction is in fact the cosine of that angle.
      double trkZenithAngle = trk->ZenithAngle();
      //Access the track vertex should be similar to start position
      TVector3 trkVertex = trk->Vertex<TVector3>();
      //Get the direction of the track vertex
      TVector3 trkVertexdir = trk->VertexDirection<TVector3>();
      //Get the track vertex momentum
      double trkVertexMom = trk->VertexMomentum();
      //Get the track vertex momentum vector
      TVector3 trkVertexMomVec = trk->VertexMomentumVector<TVector3>();
      //
      //Getting track calorimetry information
      std::vector<art::Ptr<anab::Calorimetry>> trkCalos = fmcalotrk.at(trk.key());
      if(trkCalos.empty()){
        mf::LogError("TrackValidation") << "No Calorimetry information in a recob Track. Association is made incorrectly." << std::endl;
      }//close if statement
      if(trkCalos.size() > 3){
        mf::LogError("TrackValidation") << 
          "Calorimetry information is more than one for this Track. Association is made incorrectly." <<
          " EventID: " <<evt.event()<<" TrackID: "<< trkID << std::endl;
      }//close the if statement that check trkcalo size
      
      //Loop over calorimetry 
      //-------
      //double trk_max_energy = -99999;
      //
      int Nocaloloop = 0;
      std::map<int, int> caloPlaneIdMap;
      std::map<int, float> trkPitchYMap;
      std::map<int, float> trkKEMap;
      std::map<int, float> trkRangeMap;
      std::map<int, std::vector<float> > trkdEdxMap;
      std::map<int, std::vector<float> > trkdQdxMap;
      std::map<int, std::vector<float> > trkResRangeMap;
      std::map<int, std::vector<TVector3> > trkXYZYMap;

      for(auto const& trkcalo : trkCalos){
        //Get the plane Id of the trkcalo
        int caloPlaneId = trkcalo->PlaneID().Plane;
        //if(caloPlaneId != 2){continue;} // I change it to the bestplane
        caloPlaneIdMap[caloPlaneId] = caloPlaneId;
        //Getting the dEdx for the track
        trkdEdxMap[caloPlaneId] = trkcalo->dEdx();
        //Getting the dQdx for the track
        trkdQdxMap[caloPlaneId] = trkcalo->dQdx();
        //Getting the kinetic Energy for the track.
        trkKEMap[caloPlaneId] = trkcalo->KineticEnergy();
        //Getting the range of the track
        trkRangeMap[caloPlaneId] = trkcalo->Range();
        //Getting the Residual Range
        //The residual range is defined as the distance from a given hit to the last hit in the track trajectory.
        trkResRangeMap[caloPlaneId] = trkcalo->ResidualRange();
        //Getting track pitch collection
        trkPitchYMap[caloPlaneId] = trkcalo->TrkPitchC();
        //Getting the track XYZ 
        //trkXYZYMap[caloPlaneId] = trkcalo->XYZ<std::vector<TVector3>>();
        ++Nocaloloop;
        //std::cout<<" No of calorimetry loop: "<<Nocaloloop<<std::endl;
        //
      }//close the loop over calorimetry
      //std::cout<<" "<<trkdEdxMap.size()<<" "<<trkdQdxMap.size()<<" "<<trkKEMap.size()<<" "<<trkRangeMap.size()<<std::endl;

      //Check if user wants neutrino tracks.
      if(fmtrkpf.isValid() && fUseNeutrinoTracks){
        if(fmtrkpf.at(trk.key()).size() > 0){
          if(pfpsMap.find(fmtrkpf.at(trk.key()).at(0)->Parent()) != pfpsMap.end()){
            if(pfpsMap[fmtrkpf.at(trk.key()).at(0)->Parent()]->PdgCode() != 12 && 
                pfpsMap[fmtrkpf.at(trk.key()).at(0)->Parent()]->PdgCode() != 14){
              if(fVerbose > 0){
                std::cout << "Removing the track as it is not daughter of the neutrino" << std::endl;
              }
              continue;
            }
          }
        }
      }//Close the check if user wants neutrino tracks 

      //========================
      //=== Track Validation ===
      //========================

      //Get the hits vector from the track
      std::vector<art::Ptr<recob::Hit>> trkhits = fmtrkh.at(trk.key());
      int NumberofHitsinRecoTrk = trkhits.size();
      if(fVerbose > 1){std::cout << "track hits size: " <<NumberofHitsinRecoTrk<<std::endl;}
      if(NumberofHitsinRecoTrk < fMinHitSize){
        if(fVerbose > 0){
          std::cout << "Removing track as it doesn't have enough hits" << std::endl;
        }
        continue;
      }

      //If fUseBestPlane is 2 Y plane, 3 best-plane from number of hits or 4 the best plane from energy deposited.
      int BestPlaneID = TrackUtils::FindingTrackBeastPlane(clockData, trkhits, fUseBestPlane);
      if(BestPlaneID > 2){
        mf::LogError("TrackValidation")
          << " There is an error to determine the best plane for a track; please check this event might be a bug in the 'FindingTrackBeastPlane' function" <<evt.event()<<" TrackID: "<< trkID << std::endl;
      } 
      
      //================================================
      //=== Hit Completeness and purity Calculations ===
      //================================================
      
      std::pair<int,float> TrktrackIdInfo = TrackUtils::pairTrackTrackIdToEnergyDep(clockData, EventTrksTrackIdVec, trkhits, BestPlaneID);

      int TrktrackID = TrktrackIdInfo.first;
      double TrueEnergyDepWithinTrack_FromTrueTrack = TrktrackIdInfo.second;
      double TrueEnergyDep_FromTrack = 0.0;
      double TrueEnergyDep_WithinRecoTrack = 0.0;
      int TrueHitDep_FromTrueTrack = 0;
      int TrueHitsDep_WithinRecoTrack = 0;

      //Check if the track was correctly matched.
      if(TrueEnergyDepWithinTrack_FromTrueTrack == -99999 && fMatchTracksInTruth){
        if(fVerbose > 0){
          std::cout << "Reco track not matched to true track. Think of reducing the energy threshold" << std::endl;
        }
        continue;
      }//Close the if statement to check if the track was correctly matched.

      //Tracks that have been matched to the truth information
      if(fMatchTracksInTruth){
        //Calculate the true Energy deposited By Track
        TrueEnergyDep_FromTrack = MCTrackEnergyMap[TrktrackID];
        for(geo::PlaneID planeId : geom->IteratePlaneIDs()){
          TrueHitDep_FromTrueTrack += MCTrack_hit_map[trkhit_productid][TrktrackID][planeId];
        }//close the loop over plane id from geom 

        //Get the number of hits in the reco track from the true track.
        std::map<int,std::map<geo::PlaneID,int> > MCTrack_trkhit_map = RecoUtils::NumberofPlaneHitsPerTrack(clockData, trkhits);
        for(auto const& planehitMap : MCTrack_trkhit_map){
          if(TrktrackID != planehitMap.first){continue;}
          for(auto const& NoHit : planehitMap.second){
            TrueHitsDep_WithinRecoTrack += NoHit.second;
          }//close the loop over planehitMap.second
        }//close the loop over MCTrack_trkhit_map 
      }//close if statement of checking if we need to match the tracks to truth (fMatchTracksInTruth).

      //Loop over the hits and find the IDEs
      int NoTrueHits  = 0;
      int NoAllHits = 0;
      for(auto const& trkhitIt : trkhits){
        //Get the plane ID
        geo::WireID wireid = trkhitIt->WireID();
        int PlaneID = wireid.Plane;
        ++NoAllHits;

        //if(PlaneID != BestPlaneID){continue;}

        //Split the Hit into its IDE for each track it associates with.
        std::vector<sim::TrackIDE> trackIDEs = bt_serv->HitToTrackIDEs(clockData, trkhitIt);
        if(trackIDEs.size() != 0){++NoTrueHits;}

        if(PlaneID != BestPlaneID){continue;}
        for(auto const& trackIDE: trackIDEs){
          
          //Find the true total energy deposited in a set of hits.
          TrueEnergyDep_WithinRecoTrack += trackIDE.energy;
        }//Close the loop over trackIDEs
      }//Close the loop over the reco track hits
      //std::cout<<" NoTrueHits: "<<NoTrueHits<<" NoAllHits: "<<NoAllHits<<std::endl;

      double hitcompl  = -99999;     //Hit completeness
      double hitpurity = -99999;     //Hit purity
      double energycompl  = -99999;  // energy completeness
      double energypurity = -99999;  // energy purity

      if(TrueHitDep_FromTrueTrack != 0){
        hitcompl = ((double)TrueHitsDep_WithinRecoTrack / (double)TrueHitDep_FromTrueTrack);
      }//

      if(NumberofHitsinRecoTrk != 0){
        hitpurity = ((double)TrueHitsDep_WithinRecoTrack / (double)NumberofHitsinRecoTrk);
      }

      if(TrueEnergyDep_FromTrack != 0){
        energycompl = ((double)TrueEnergyDepWithinTrack_FromTrueTrack / (double)TrueEnergyDep_FromTrack);
      }

      if(TrueEnergyDep_WithinRecoTrack != 0){
        energypurity = ((double)TrueEnergyDepWithinTrack_FromTrueTrack / (double)TrueEnergyDep_WithinRecoTrack);
      }

      trkAllHitsPurity_TreeVal[fTrackModuleLabel].push_back(hitcompl);
      trkAllHitsComp_TreeVal[fTrackModuleLabel].push_back(hitpurity);
      trkEnergyPurityBestPlane_TreeVal[fTrackModuleLabel].push_back(energycompl);
      trkEnergyCompBestPlane_TreeVal[fTrackModuleLabel].push_back(energypurity);

      if(NoTrueHits != 0){
        trkSignalToNoiseRatio_TreeVal[fTrackModuleLabel].push_back(((NoTrueHits - (NoAllHits - NoTrueHits)) / NoTrueHits));
      }

      //Calculate metrics
      TVector3 TrueTrackStartDirection = {-99999,-99999,-99999};
      TVector3 TrueTrackEndDirection   = {-99999,-99999,-99999};
      TVector3 TrueTrackDirection      = {-99999,-99999,-99999};
      TVector3 TrueTrackStart          = {-99999,-99999, -99999};
      float MCKEnergy                  = -99999.9; //in GeV
      double TrueTrackLength           = 0; // in cm

      if(fMatchTracksInTruth){
        TLorentzVector PositionTrajStart = {-99999,-99999,-99999, -99999};
        TLorentzVector PositionTrajEnd   = {-99999,-99999,-99999, -99999};

        const simb::MCParticle* MCTrackParticle = trueparticlesMap.at(TrktrackID);

        int track_pdgcode = MCTrackParticle->PdgCode();

        //Check if another particle actually gave more energy.
        int MostHitsTrackID = std::abs(RecoUtils::TrueParticleIDFromTotalRecoHits(clockData, trkhits));
        int trk_most_hit_pdg = -99999;
        if(trueparticlesMap.find(MostHitsTrackID) != trueparticlesMap.end()){
          trk_most_hit_pdg = trueparticlesMap.at(MostHitsTrackID)->PdgCode();
        }//close if statement that look if MostHitsTrackId is a true particle
        trkTrackPDG_MostHit_TreeVal[fTrackModuleLabel].push_back(trk_most_hit_pdg);

        trkPDG_TreeVal[fTrackModuleLabel].push_back((std::abs(track_pdgcode)));
        trkMother_TreeVal[fTrackModuleLabel].push_back(MCTrackParticle->Mother());

        //Find the Energy of the particle:
        float mctrkEnergy = MCTrackParticle->E();
        trkTrueEnergy_TreeVal[fTrackModuleLabel].push_back((mctrkEnergy*1000)); //Energy in MeV
        //kinetic energy
        MCKEnergy = (mctrkEnergy - (MCTrackParticle->Mass())); //in GeV
        trkTrueKEnergy_TreeVal[fTrackModuleLabel].push_back((MCKEnergy*1000)); //Energy in MeV

        trkStartEndProcess_TreeVal[fTrackModuleLabel].push_back(MCTrackParticle->EndProcess());
        trkStartProcess_TreeVal[fTrackModuleLabel].push_back(MCTrackParticle->Process());

        //Get the number of Traj points to loop over
        unsigned int trkTrajPoints = MCTrackParticle->NumberTrajectoryPoints();

        //Select first traj point where the track loses energy, last be default
        PositionTrajStart = MCTrackParticle->Position(0);
        PositionTrajEnd   = MCTrackParticle->Position(trkTrajPoints-1);

        //Get the track length of the track.
        TrueTrackLength = (PositionTrajStart - PositionTrajEnd).Vect().Mag(); //in cm

        //How straight are the IDE points. First Define the line.
        TVector3 track_line = (PositionTrajStart - PositionTrajEnd).Vect().Unit();
        //
        double TrackResidual = 0;
        double TrackResidual_iter = 0;
        
        //Get the energy depositions
        std::vector<const sim::IDE*> IDEs = bt_serv->TrackIdToSimIDEs_Ps(TrktrackID);
        for(auto const& ide : IDEs){
           //Calculate the perpendicualr distance from the straight line
           if(ide->trackID == TrktrackID){
             TVector3 pos = {ide->x, ide->y,ide->z};
             TVector3 rel_pos = pos - PositionTrajStart.Vect();
             double len = rel_pos.Dot(track_line);
             if(len > TrueTrackLength){continue;}
             double perp = (rel_pos - (len * track_line)).Mag();
             TrackResidual += perp;
             ++TrackResidual_iter;
           }//Close the if statement
        }//Close the loop over IDEs

        trkTrueTrackLength_TreeVal[fTrackModuleLabel].push_back(TrueTrackLength);
        trkTrueTrackSpread_TreeVal[fTrackModuleLabel].push_back(TrackResidual/TrackResidual_iter);
        trkTrueTrackResidual_TreeVal[fTrackModuleLabel].push_back(TrackResidual);
        //TODO calculate the true residual range by using the trajactory points or
        //true space points of track using truth hits of TrktrackID to do so check Calorimetry_module.cc

        //The three vector for track length is the track direction
        TrueTrackDirection      = MCTrackParticle->Momentum().Vect();
        //The track start direction and end direction are done in the same way that have been in reco. algs.
        TrueTrackStartDirection = ((MCTrackParticle->Momentum(0) - MCTrackParticle->Momentum(20)).Vect());
        TrueTrackEndDirection   = (((MCTrackParticle->Momentum(0) - MCTrackParticle->Momentum(trkTrajPoints-21)) - MCTrackParticle->Momentum(trkTrajPoints-1)).Vect());
        TrueTrackStart          = PositionTrajStart.Vect();
        //
      }//close fMatchTracksInTruth

      //===========================
      //=== Projection Matching ===
      //===========================
      
      float geoprojectionmatched_score = -99999;
      if(fmsptrk.isValid()){
        //Get the spacepoints
        std::vector<art::Ptr<recob::SpacePoint> > trksps = fmsptrk.at(trk.key());

        if(trksps.size() > 0){
          art::Handle<std::vector<recob::SpacePoint> > spHandle;
          evt.get(trksps.front().id(),spHandle);

          if(spHandle.isValid()){
            art::FindManyP<recob::Hit> fmsph(spHandle, evt, spHandle.provenance()->moduleLabel());

            float geomatched = 0;
            for(auto const& sp: trksps){

              //Get the the hit
              std::vector< art::Ptr<recob::Hit> > hit = fmsph.at(sp.key());
              if(hit.size() < 1){continue;}
              
              //Get the true position
              std::vector<double> hitXYZ = {-999,-999,-999};
              try{
                hitXYZ = bt_serv->HitToXYZ(clockData, hit[0]);
              }
              catch(...){
                if(fVerbose>2){std::cout << "Noise Hit" << std::endl;}
                continue;
              }

              TVector3 trueposition = {hitXYZ[0],hitXYZ[1],hitXYZ[2]};
              //Get the spacepoint xyz.
              const Double32_t* sp_xyz = sp->XYZ();

              TVector3 sp_postiion = {sp_xyz[0], sp_xyz[1], sp_xyz[2]};
              //If the spacepoint was within 0.5 cm of the ide then its a good match.
              if((trueposition-sp_postiion).Mag() < fMatchedCut){
                ++geomatched;
              }
            }//close the loop over space points
            geoprojectionmatched_score = geomatched/(float) trksps.size();
          }//close the statement checking if spHandle is Valid
        }//Close the statement checking the size of trksps size
      }//Close statement checking if association is valid
      // std::cout<< "finished to spacepoint" << std::endl;
      
      //=============================
      //=== Calculate the metrics ===
      //=============================

      //Get the angles between the direction
      float TrackDirection_Xdiff      = -99999;
      float TrackDirection_Ydiff      = -99999;
      float TrackDirection_Zdiff      = -99999;
      float TrackStartDirection_Xdiff = -99999;
      float TrackStartDirection_Ydiff = -99999;
      float TrackStartDirection_Zdiff = -99999;
      float TrackEndDirection_Xdiff   = -99999;
      float TrackEndDirection_Ydiff   = -99999;
      float TrackEndDirection_Zdiff   = -99999;
      float TrackDirection_XTrue      = -99999;
      float TrackDirection_YTrue      = -99999;
      float TrackDirection_ZTrue      = -99999;
      float TrackStartDirection_XTrue = -99999;
      float TrackStartDirection_YTrue = -99999;
      float TrackStartDirection_ZTrue = -99999;
      float TrackEndDirection_XTrue   = -99999;
      float TrackEndDirection_YTrue   = -99999;
      float TrackEndDirection_ZTrue   = -99999;
      float TrackDirection_diff       = -99999;
      float TrackStartDirection_diff  = -99999;
      float TrackEndDirection_diff    = -99999;
      double Start_diff               = -99999;
      double Start_diff_X             = -99999;
      double Start_diff_Y             = -99999;
      double Start_diff_Z             = -99999;
      double trkdEdx                  = -99999;
      std::map<int, double> PlanetrkdEdx;
      //double recotrkEnergy            = -99999;
      //double recotrkEnergyRat         = -99999;
      //double recotrkEnergyDiff        = -99999;
      //For kinetic energy
      float recotrkKE                = -99999;
      float recotrkKERat             = -99999;
      float recotrkKEDiff            = -99999;
      float recotrk_range_p          = -99999;
      
      float recotrk_bestLogLikelihood_p  = -99999;
      float recotrk_bestMomUncertainty_p = -99999;
      float recotrk_bestMomentum_p       = -99999;
      float recotrk_bwdLogLikelihood_p   = -99999;
      float recotrk_bwdMomentum_p        = -99999;
      float recotrk_bwdMomUncertainty_p  = -99999;
      float recotrk_deltaLogLikelihood_p = -99999;
      float recotrk_fwdLogLikelihood_p   = -99999;
      float recotrk_fwdMomentum_p        = -99999;
      float recotrk_fwdMomUncertainty_p  = -99999;

      
      if(fMatchTracksInTruth){
        
        //The difference for the track start direction
        if(std::sqrt(((TrueTrackDirection.Y()*TrueTrackDirection.Y()) + (TrueTrackDirection.Z()*TrueTrackDirection.Z())))* std::sqrt(((trkstartdir.Y()*trkstartdir.Y()) + (trkstartdir.Z()*trkstartdir.Z()))) !=0){
          //
          TrackDirection_Xdiff = (TrueTrackDirection.Y()*trkstartdir.Y() + TrueTrackDirection.Z()*trkstartdir.Z()) / (std::sqrt(((TrueTrackDirection.Y()*TrueTrackDirection.Y()) + (TrueTrackDirection.Z()*TrueTrackDirection.Z())))* std::sqrt(((trkstartdir.Y()*trkstartdir.Y()) + (trkstartdir.Z()*trkstartdir.Z()))));
          
          TrackDirection_XTrue      = TrueTrackDirection.X();
        }//close if statement for x direction difference 
        
        if(std::sqrt(((TrueTrackDirection.X()*TrueTrackDirection.X()) + (TrueTrackDirection.Z()*TrueTrackDirection.Z())))* std::sqrt(((trkstartdir.X()*trkstartdir.X()) + (trkstartdir.Z()*trkstartdir.Z()))) !=0){
          //
          TrackDirection_Ydiff = (TrueTrackDirection.X()*trkstartdir.X() + TrueTrackDirection.Z()*trkstartdir.Z()) / (std::sqrt(((TrueTrackDirection.X()*TrueTrackDirection.X()) + (TrueTrackDirection.Z()*TrueTrackDirection.Z())))* std::sqrt(((trkstartdir.X()*trkstartdir.X()) + (trkstartdir.Z()*trkstartdir.Z()))));

          TrackDirection_YTrue = TrueTrackDirection.Y();
        }//close if statement for y direction difference 
        
        if(std::sqrt(((TrueTrackDirection.X()*TrueTrackDirection.X()) + (TrueTrackDirection.Y()*TrueTrackDirection.Y())))* std::sqrt(((trkstartdir.X()*trkstartdir.X()) + (trkstartdir.Y()*trkstartdir.Y()))) !=0){
          //
          TrackDirection_Zdiff = (TrueTrackDirection.X()*trkstartdir.X() + TrueTrackDirection.Y()*trkstartdir.Y()) / (std::sqrt(((TrueTrackDirection.X()*TrueTrackDirection.X()) + (TrueTrackDirection.Y()*TrueTrackDirection.Y())))* std::sqrt(((trkstartdir.X()*trkstartdir.X()) + (trkstartdir.Y()*trkstartdir.Y()))));

          TrackDirection_ZTrue      = TrueTrackDirection.Z();
        }//close if statement for z direction difference 


        //The difference for the track start direction
        if(std::sqrt(((TrueTrackStartDirection.Y()*TrueTrackStartDirection.Y()) + (TrueTrackStartDirection.Z()*TrueTrackStartDirection.Z())))* std::sqrt(((trkstartdir.Y()*trkstartdir.Y()) + (trkstartdir.Z()*trkstartdir.Z()))) !=0){
          //
          TrackStartDirection_Xdiff = (TrueTrackStartDirection.Y()*trkstartdir.Y() + TrueTrackStartDirection.Z()*trkstartdir.Z()) / (std::sqrt(((TrueTrackStartDirection.Y()*TrueTrackStartDirection.Y()) + (TrueTrackStartDirection.Z()*TrueTrackStartDirection.Z())))* std::sqrt(((trkstartdir.Y()*trkstartdir.Y()) + (trkstartdir.Z()*trkstartdir.Z()))));
          TrackStartDirection_XTrue = TrueTrackStartDirection.X();

        }//close if statement for x start direction difference 

        if(std::sqrt(((TrueTrackStartDirection.X()*TrueTrackStartDirection.X()) + (TrueTrackStartDirection.Z()*TrueTrackStartDirection.Z())))* std::sqrt(((trkstartdir.X()*trkstartdir.X()) + (trkstartdir.Z()*trkstartdir.Z()))) !=0){
          //
          TrackStartDirection_Ydiff = (TrueTrackStartDirection.X()*trkstartdir.X() + TrueTrackStartDirection.Z()*trkstartdir.Z()) / (std::sqrt(((TrueTrackStartDirection.X()*TrueTrackStartDirection.X()) + (TrueTrackStartDirection.Z()*TrueTrackStartDirection.Z())))* std::sqrt(((trkstartdir.X()*trkstartdir.X()) + (trkstartdir.Z()*trkstartdir.Z()))));
          TrackStartDirection_YTrue = TrueTrackStartDirection.Y();

        }//close if statement for y start direction difference 

        if(std::sqrt(((TrueTrackStartDirection.X()*TrueTrackStartDirection.X()) + (TrueTrackStartDirection.Y()*TrueTrackStartDirection.Y())))* std::sqrt(((trkstartdir.X()*trkstartdir.X()) + (trkstartdir.Y()*trkstartdir.Y()))) !=0){
          // 
          TrackStartDirection_Zdiff = (TrueTrackStartDirection.X()*trkstartdir.X() + TrueTrackStartDirection.Y()*trkstartdir.Y()) / (std::sqrt(((TrueTrackStartDirection.X()*TrueTrackStartDirection.X()) + (TrueTrackStartDirection.Y()*TrueTrackStartDirection.Y())))* std::sqrt(((trkstartdir.X()*trkstartdir.X()) + (trkstartdir.Y()*trkstartdir.Y()))));
          TrackStartDirection_ZTrue = TrueTrackStartDirection.Z();

        }//close if statement for z start direction difference 


        //The difference for the track end direction
        if(std::sqrt(((TrueTrackEndDirection.X()*TrueTrackEndDirection.X()) + (TrueTrackEndDirection.Y()*TrueTrackEndDirection.Y())))* std::sqrt(((trkenddir.X()*trkenddir.X()) + (trkenddir.Y()*trkenddir.Y()))) !=0){
          //
         
          TrackEndDirection_Xdiff = (TrueTrackEndDirection.Y()*trkenddir.Y() + TrueTrackEndDirection.Z()*trkenddir.Z()) / (std::sqrt(((TrueTrackEndDirection.Y()*TrueTrackEndDirection.Y()) + (TrueTrackEndDirection.Z()*TrueTrackEndDirection.Z())))* std::sqrt(((trkenddir.Y()*trkenddir.Y()) + (trkenddir.Z()*trkenddir.Z()))));
 
          TrackEndDirection_XTrue   = TrueTrackEndDirection.X();
        }//close if statement for x end direction difference 
        
        if(std::sqrt(((TrueTrackEndDirection.X()*TrueTrackEndDirection.X()) + (TrueTrackEndDirection.Y()*TrueTrackEndDirection.Y())))* std::sqrt(((trkenddir.X()*trkenddir.X()) + (trkenddir.Y()*trkenddir.Y()))) !=0){
          //
          TrackEndDirection_Ydiff = (TrueTrackEndDirection.X()*trkenddir.X() + TrueTrackEndDirection.Z()*trkenddir.Z()) / (std::sqrt(((TrueTrackEndDirection.X()*TrueTrackEndDirection.X()) + (TrueTrackEndDirection.Z()*TrueTrackEndDirection.Z())))* std::sqrt(((trkenddir.X()*trkenddir.X()) + (trkenddir.Z()*trkenddir.Z()))));
 
          TrackEndDirection_YTrue   = TrueTrackEndDirection.Y();
        }//close if statement for y end direction difference
        
        if(std::sqrt(((TrueTrackEndDirection.X()*TrueTrackEndDirection.X()) + (TrueTrackEndDirection.Y()*TrueTrackEndDirection.Y())))* std::sqrt(((trkenddir.X()*trkenddir.X()) + (trkenddir.Y()*trkenddir.Y()))) !=0){
          //
         
          TrackEndDirection_Zdiff = (TrueTrackEndDirection.X()*trkenddir.X() + TrueTrackEndDirection.Y()*trkenddir.Y()) / (std::sqrt(((TrueTrackEndDirection.X()*TrueTrackEndDirection.X()) + (TrueTrackEndDirection.Y()*TrueTrackEndDirection.Y())))* std::sqrt(((trkenddir.X()*trkenddir.X()) + (trkenddir.Y()*trkenddir.Y()))));
 
          TrackEndDirection_ZTrue   = TrueTrackEndDirection.Z();
        }//close if statement for z end direction difference


        if(TrueTrackDirection.Mag() != 0 || trkstartdir.Mag() !=0){
          TrackDirection_diff = (TrueTrackDirection.Dot(trkstartdir)/(TrueTrackDirection.Mag()*trkstartdir.Mag()));
        }

        if(TrueTrackStartDirection.Mag() != 0 || trkstartdir.Mag() !=0){
          TrackStartDirection_diff = (TrueTrackStartDirection.Dot(trkstartdir)/(TrueTrackStartDirection.Mag()*trkstartdir.Mag()));
        }

        if(TrueTrackEndDirection.Mag() != 0 || trkenddir.Mag() !=0){
          TrackEndDirection_diff = (TrueTrackEndDirection.Dot(trkenddir)/(TrueTrackEndDirection.Mag()*trkenddir.Mag()));
        }

        //Get the Error in the position. intialised as 0,0,0 this is a problem here.
        Start_diff   = (TrueTrackStart - trkStartPosition).Mag();
        Start_diff_X = (TrueTrackStart.X() - trkStartPosition.X());
        Start_diff_Y = (TrueTrackStart.Y() - trkStartPosition.Y());
        Start_diff_Z = (TrueTrackStart.Z() - trkStartPosition.Z());

      }//Close fMatchTracksInTruth

      //Getting median or mean of track's dEdx
      for(int PlaneID = 0; PlaneID < 3; ++PlaneID){
        if(trkdEdxMap[PlaneID].size() != 0){
          if(fUseMedian){
            trkdEdx = TMath::Median((trkdEdxMap[PlaneID].size()), &(trkdEdxMap[PlaneID][0]));
            PlanetrkdEdx[PlaneID] = trkdEdx;
          }//close if statement to check if we want to use the median or mean
          else{
            //Else calculate the mean value.
            double trkdEdx_mean = 0;
            int    numdEdx_mean = 0;
            for(auto const& dEdx : trkdEdxMap[BestPlaneID]){
              //if(dEdx > 10 || dEdx < 0){continue;}
              ++numdEdx_mean;
              trkdEdx_mean += dEdx;
            }//close the loop over track dEdx
            double trkdEdx = (trkdEdx_mean/(float)(numdEdx_mean));
            PlanetrkdEdx[PlaneID] = trkdEdx;
          }//Close else statement
        }//close if statement to get the track dEdx mean or median
      }

      //Getting the kinetic energy ratio and difference
      if(trkKEMap.size() != 0){
        recotrkKE = trkKEMap[BestPlaneID];
      }

      if(MCKEnergy != 0 && recotrkKE > 0){
        recotrkKERat = (recotrkKE / (MCKEnergy * 1000));
        recotrkKEDiff = ((recotrkKE - (MCKEnergy * 1000)) / (MCKEnergy * 1000));
      }

      {//####################### range_p and mcs_p #######################
        //calculating momentum range and momentum coulomb scattering.
        const simb::MCParticle* MCTrackParticle = trueparticlesMap.at(TrktrackID);
        int track_pdgcode = MCTrackParticle->PdgCode();
        if(std::abs(track_pdgcode) == 11 || std::abs(track_pdgcode) == 22){
          ++fNoShowerTrack_TreeVal;
        }
        
        //Other initialisations
        std::map<int, int> PdgMCSIter = {{13, 0},{211, 1},{321, 2},{2212, 3}};
        std::map<int, int> PdgRangeIter = {{13, 0},{211, 1},{2212, 2}};

        // Make sure that we are looking at track particles and daughter
        if(std::find(fTrackPdgCodeVect.begin(),fTrackPdgCodeVect.end(),track_pdgcode) != fTrackPdgCodeVect.end()){
          std::vector<art::FindManyP<recob::MCSFitResult>> fmMCSs;
          //static const std::vector<std::string> PIDnames {"muon", "pion", "kaon", "proton"};
          for(std::string pid: fPIDnames){
            art::InputTag tag(fpandoraTrackMCS, pid);
            fmMCSs.push_back(art::FindManyP<recob::MCSFitResult>(trackListHandle, evt, tag));
          }

          std::vector<art::FindManyP<sbn::RangeP>> fmRanges;
          //static const std::vector<std::string> rangePIDnames {"muon", "pion", "proton"};
          for(std::string pid: frangePIDnames){
            art::InputTag tag(fpandoraTrackRange, pid);
            fmRanges.push_back(art::FindManyP<sbn::RangeP>(trackListHandle, evt, tag));
          }
          //trkf::TrajectoryMCSFitter fMCSCalculator;
          if(fmMCSs.size()!=0 && fmMCSs.at(PdgMCSIter[(std::abs(track_pdgcode))]).size()!=0){
            std::vector<art::Ptr<recob::MCSFitResult>> mcscol = fmMCSs[(PdgMCSIter[(std::abs(track_pdgcode))])].at(trk.key());//

            if(!mcscol.empty()){
              //std::cout<<" mcscol: "<<mcscol.size()<<" : "<< (PdgMCSIter[(std::abs(track_pdgcode))])<<std::endl;
              recotrk_bestLogLikelihood_p  = mcscol[0]->bestLogLikelihood();
              recotrk_bestMomUncertainty_p = mcscol[0]->bestMomUncertainty();
              recotrk_bestMomentum_p       = mcscol[0]->bestMomentum();
              recotrk_bwdLogLikelihood_p   = mcscol[0]->bwdLogLikelihood();
              recotrk_bwdMomentum_p        = mcscol[0]->bwdMomentum();
              recotrk_bwdMomUncertainty_p  = mcscol[0]->bwdMomUncertainty();
              recotrk_deltaLogLikelihood_p = mcscol[0]->deltaLogLikelihood();
              recotrk_fwdLogLikelihood_p   = mcscol[0]->fwdLogLikelihood();
              recotrk_fwdMomentum_p        = mcscol[0]->fwdMomentum();
              recotrk_fwdMomUncertainty_p  = mcscol[0]->fwdMomUncertainty();

            }
          }
          
          if(fmRanges.size()!=0 && fmRanges.at(PdgRangeIter[(std::abs(track_pdgcode))]).size()!=0){
            std::vector<art::Ptr<sbn::RangeP>> mrange = fmRanges[(PdgRangeIter[(std::abs(track_pdgcode))])].at(trk.key());
           
            if(!mrange.empty()){
              recotrk_range_p = mrange[0]->range_p;
            }

          }
          //std::cout<<" recotrk_mcs_p: "<<recotrk_mcs_p<<" recotrk_range_p: "<<recotrk_range_p<<std::endl; 
          //std::cout<<" fmMCSs: "<<fmMCSs.size()<<" fmRanges: "<<fmRanges.size()<<std::endl;
        }//Close the if statement for checking the track pdgcode is for a track particle
        //
      }//close the block
 

      // ****************** Fill the histograms ******************
      trkDirX_TreeVal[fTrackModuleLabel].push_back(TrackDirection_Xdiff);
      trkDirY_TreeVal[fTrackModuleLabel].push_back(TrackDirection_Ydiff);
      trkDirZ_TreeVal[fTrackModuleLabel].push_back(TrackDirection_Zdiff);
      trkTrueDirX_TreeVal[fTrackModuleLabel].push_back(TrackDirection_XTrue);
      trkTrueDirY_TreeVal[fTrackModuleLabel].push_back(TrackDirection_YTrue);
      trkTrueDirZ_TreeVal[fTrackModuleLabel].push_back(TrackDirection_ZTrue);
      trkDirDiff_TreeVal[fTrackModuleLabel].push_back(TrackDirection_diff);
      //For the start of the track
      trkStartDirX_TreeVal[fTrackModuleLabel].push_back(TrackStartDirection_Xdiff);
      trkStartDirY_TreeVal[fTrackModuleLabel].push_back(TrackStartDirection_Ydiff);
      trkStartDirZ_TreeVal[fTrackModuleLabel].push_back(TrackStartDirection_Zdiff);
      trkStartTrueDirX_TreeVal[fTrackModuleLabel].push_back(TrackStartDirection_XTrue);
      trkStartTrueDirY_TreeVal[fTrackModuleLabel].push_back(TrackStartDirection_YTrue);
      trkStartTrueDirZ_TreeVal[fTrackModuleLabel].push_back(TrackStartDirection_ZTrue);
      trkStratDirDiff_TreeVal[fTrackModuleLabel].push_back(TrackStartDirection_diff);
      //For the end of the track
      trkEndDirX_TreeVal[fTrackModuleLabel].push_back(TrackEndDirection_Xdiff);
      trkEndDirY_TreeVal[fTrackModuleLabel].push_back(TrackEndDirection_Ydiff);
      trkEndDirZ_TreeVal[fTrackModuleLabel].push_back(TrackEndDirection_Zdiff);
      trkEndTrueDirX_TreeVal[fTrackModuleLabel].push_back(TrackEndDirection_XTrue);
      trkEndTrueDirY_TreeVal[fTrackModuleLabel].push_back(TrackEndDirection_YTrue);
      trkEndTrueDirZ_TreeVal[fTrackModuleLabel].push_back(TrackEndDirection_ZTrue);
      trkEndDirDiff_TreeVal[fTrackModuleLabel].push_back(TrackEndDirection_diff);
      //
      trkStartX_TreeVal[fTrackModuleLabel].push_back(Start_diff_X);
      trkStartY_TreeVal[fTrackModuleLabel].push_back(Start_diff_Y);
      trkStartZ_TreeVal[fTrackModuleLabel].push_back(Start_diff_Z);
      trkStartDist_TreeVal[fTrackModuleLabel].push_back(Start_diff);
      //
      trkBestPlane_TreeVal[fTrackModuleLabel].push_back(BestPlaneID);//
      trkLength_TreeVal[fTrackModuleLabel].push_back(trkLength);//
      trkLengthDiff_TreeVal[fTrackModuleLabel].push_back((TrueTrackLength - trkLength));
      trkNumHits_TreeVal[fTrackModuleLabel].push_back(NumberofHitsinRecoTrk);
      trkGeoProjectionMatched_TreeVal[fTrackModuleLabel].push_back(geoprojectionmatched_score);
      trkID_TreeVal[fTrackModuleLabel].push_back(trkID);
      trkParticleId_TreeVal[fTrackModuleLabel].push_back(trkParticleId);
      trkChi2_TreeVal[fTrackModuleLabel].push_back(trkChi2);
      trkNdof_TreeVal[fTrackModuleLabel].push_back(trkNdof);
      trkStartMomentum_TreeVal[fTrackModuleLabel].push_back(trkstartmom);
      trkEndMomentum_TreeVal[fTrackModuleLabel].push_back(trkEndMomentum);
      trkVertexMomentum_TreeVal[fTrackModuleLabel].push_back(trkVertexMom);
      trkAzimuthAngle_TreeVal[fTrackModuleLabel].push_back(trkAzimuthAngle);
      trkPhi_TreeVal[fTrackModuleLabel].push_back(trkPhi);
      trkTheta_TreeVal[fTrackModuleLabel].push_back(trkTheta);
      trkZenithAngle_TreeVal[fTrackModuleLabel].push_back(trkZenithAngle);
      //Track calorimetry infrormation
      calotrkPlaneId_TreeVal[fTrackModuleLabel].push_back(caloPlaneIdMap[BestPlaneID]);
      trkdEdxVecBP_TreeVal[fTrackModuleLabel].push_back(trkdEdxMap[BestPlaneID]);
      trkdEdxVecU_TreeVal[fTrackModuleLabel].push_back(trkdEdxMap[0]);
      trkdEdxVecV_TreeVal[fTrackModuleLabel].push_back(trkdEdxMap[1]);
      trkdEdxVecY_TreeVal[fTrackModuleLabel].push_back(trkdEdxMap[2]);

      trkdQdxVecBP_TreeVal[fTrackModuleLabel].push_back(trkdQdxMap[BestPlaneID]);
      trkdQdxVecU_TreeVal[fTrackModuleLabel].push_back(trkdQdxMap[0]);
      trkdQdxVecV_TreeVal[fTrackModuleLabel].push_back(trkdQdxMap[1]);
      trkdQdxVecY_TreeVal[fTrackModuleLabel].push_back(trkdQdxMap[2]);
      
      trkKEBP_TreeVal[fTrackModuleLabel].push_back(trkKEMap[BestPlaneID]); //
      trkKEU_TreeVal[fTrackModuleLabel].push_back(trkKEMap[0]); //
      trkKEV_TreeVal[fTrackModuleLabel].push_back(trkKEMap[1]); //
      trkKEY_TreeVal[fTrackModuleLabel].push_back(trkKEMap[2]); //

      trkRange_TreeVal[fTrackModuleLabel].push_back(trkRangeMap[BestPlaneID]);
      
      trkResRangeBP_TreeVal[fTrackModuleLabel].push_back(trkResRangeMap[BestPlaneID]);
      trkResRangeU_TreeVal[fTrackModuleLabel].push_back(trkResRangeMap[0]);
      trkResRangeV_TreeVal[fTrackModuleLabel].push_back(trkResRangeMap[1]);
      trkResRangeY_TreeVal[fTrackModuleLabel].push_back(trkResRangeMap[2]);

      trkPitchY_TreeVal[fTrackModuleLabel].push_back(trkPitchYMap[BestPlaneID]);
      //
      //trkEnergy_TreeVal[fTrackModuleLabel].push_back(recotrkEnergy);
      //trkEnergyRat_TreeVal[fTrackModuleLabel].push_back(recotrkEnergyRat);
      //trkEnergyDiff_TreeVal[fTrackModuleLabel].push_back(recotrkEnergyDiff);
      //
      trkKinEnergy_TreeVal[fTrackModuleLabel].push_back(recotrkKE);
      trkKinEnergyRat_TreeVal[fTrackModuleLabel].push_back(recotrkKERat);
      trkKinEnergyDiff_TreeVal[fTrackModuleLabel].push_back(recotrkKEDiff);
      //
      //Momentum values range_p -> momentum range and mcs_p -> momentum coulomb scattering.
      recotrk_range_p_TreeVal[fTrackModuleLabel].push_back(recotrk_range_p);
      recotrk_bestLogLikelihood_p_TreeVal[fTrackModuleLabel].push_back(recotrk_bestLogLikelihood_p);
      recotrk_bestMomUncertainty_p_TreeVal[fTrackModuleLabel].push_back(recotrk_bestMomUncertainty_p);
      recotrk_bestMomentum_p_TreeVal[fTrackModuleLabel].push_back(recotrk_bestMomentum_p);
      recotrk_bwdLogLikelihood_p_TreeVal[fTrackModuleLabel].push_back(recotrk_bwdLogLikelihood_p);
      recotrk_bwdMomentum_p_TreeVal[fTrackModuleLabel].push_back(recotrk_bwdMomentum_p);
      recotrk_bwdMomUncertainty_p_TreeVal[fTrackModuleLabel].push_back(recotrk_bwdMomUncertainty_p);
      recotrk_deltaLogLikelihood_p_TreeVal[fTrackModuleLabel].push_back(recotrk_deltaLogLikelihood_p);
      recotrk_fwdLogLikelihood_p_TreeVal[fTrackModuleLabel].push_back(recotrk_fwdLogLikelihood_p);
      recotrk_fwdMomentum_p_TreeVal[fTrackModuleLabel].push_back(recotrk_fwdMomentum_p);
      recotrk_fwdMomUncertainty_p_TreeVal[fTrackModuleLabel].push_back(recotrk_fwdMomUncertainty_p);
      
      trkdEdxBP_TreeVal[fTrackModuleLabel].push_back(PlanetrkdEdx[BestPlaneID]);
      trkdEdxU_TreeVal[fTrackModuleLabel].push_back(PlanetrkdEdx[0]);
      trkdEdxV_TreeVal[fTrackModuleLabel].push_back(PlanetrkdEdx[1]);
      trkdEdxY_TreeVal[fTrackModuleLabel].push_back(PlanetrkdEdx[2]);

      if(fVerbose > 2){
        std::cout<<"#################################################" <<std::endl;
        std::cout<<"            Global Event Information             " <<std::endl;
        std::cout<<"#################################################" <<std::endl;
        
        std::cout<<"Track Label:--------------------------"<<fTrackModuleLabel<<std::endl;
        std::cout<<"Number of Reco Tracks: ---------------"<<recoTracks.size()<<std::endl;
        std::cout<<"Energy Simulated: --------------------"<<trueTrkEnergy<<std::endl;
        std::cout<<"Number of True Tracks Before Cuts: ---"<<numTrueTraks<<std::endl;
        std::cout<<"Number of True Tracks in the Event:---"<<EventTrksTrackIdVec.size()<<std::endl;

        std::cout<<"#################################################" <<std::endl;
        std::cout<<"   Reconstructed/Associated Truth Information    " <<std::endl;
        std::cout<<"#################################################" <<std::endl;
        
        std::cout<<"Hit Size:--------------- "<<NumberofHitsinRecoTrk<<std::endl;
        std::cout<<"TrackBest_Plane:-------- "<<BestPlaneID <<std::endl;
        std::cout<<"True Start:------------- "<<TrueTrackStart.X()<<" Track Start: "<<trkStartPosition.X()<<"\n";
        std::cout<<"Reco. Poisition (x,y,z): ("<<trkStartPosition.X()<<","<<trkStartPosition.Y()<<","<<trkStartPosition.Z()<<") \n";
        //std::cout<<"TrueTrackLength: "<<TrueTrackLength<<" Reco Track Length: "<<trkLength<<std::endl;
        std::cout<<"Track Phi Angle:-------- "<<trkPhi<<std::endl;
        std::cout<<"Track Theta Angle:------ "<<trkTheta<<std::endl;
        std::cout<<"Best Plane:------------- "<<BestPlaneID<<std::endl;
        std::cout<<"Best Plane Reco Trk KE:- "<<recotrkKE<<std::endl;
        std::cout<<"Best Plane Reco Trk dEdx:"<<PlanetrkdEdx[BestPlaneID]<<std::endl;
        
        std::cout<<"True Energy Deposited from true track within track: "<<TrueEnergyDepWithinTrack_FromTrueTrack<<std::endl;
        std::cout<<"True Energy Deposited From true associated track:-- "<<TrueEnergyDep_FromTrack<<std::endl;
        std::cout<<"True Energy Deposited by hits in track:------------ "<<TrueEnergyDep_WithinRecoTrack<<"\n";
        std::cout<<"Energy Purity:--- "<<energypurity<<" Energy completeness: "<<energycompl<<std::endl;
        std::cout<<"Hit Purity:------ "<<hitpurity<<   " Hit Completeness:--- "<<hitcompl<<std::endl;
        std::cout<<"#################################################"<<std::endl;
        //
      }//Close printing the track information

      if(fVerbose > 0){std::cout << "Track Validation Complete" << std::endl;}

      ++fnumrecotracksana;
      //
    }//close the loop over tracks

    //========================
    //===   Event Metric   ===
    //========================

    //Calculate the True Hit number
    for(auto const& track_id : EventTrksTrackIdVec){
      int TrueHitDep_FromTrueTracks = 0;
      for(geo::PlaneID plane_id: geom->IteratePlaneIDs()){
        TrueHitDep_FromTrueTracks += MCTrack_hit_map[trkhit_productid][track_id][plane_id];
      }
      TrueHitNum_TreeVal[(std::to_string(track_id))].push_back(TrueHitDep_FromTrueTracks);
    }//Close the loop over track mothers

    //This is rather pandora track specific for tuning pandora so lets add it if the track module 
    //used pandora and the track module Pandora track is used.

    //Whats the segementyness of the event.
    eSegmentation_TreeVal[fTrackModuleLabel].push_back((float)recoTracks.size() / (float)numTrueTraks);
    eNumRecoTracks_TreeVal[fTrackModuleLabel].push_back(recoTracks.size());

    if(fVerbose>2){
      std::cout<<"finished track module: "<<fTrackModuleLabel<<std::endl;
    }
    //********************************
  }//Close the loop over Track Module

  //Fill the tree
  fTree->Fill();

  for(auto const& fTrackModuleLabel : fTrackModuleLabels){
    PFPNeutrinos_TreeVal[fTrackModuleLabel].clear();
    PFPAllTracks_TreeVal[fTrackModuleLabel].clear();
    PFPAllShowers_TreeVal[fTrackModuleLabel].clear();
    PFPrimaryTracks_TreeVal[fTrackModuleLabel].clear();
    PFPrimaryShowers_TreeVal[fTrackModuleLabel].clear();
    PFPVertexDistX_TreeVal[fTrackModuleLabel].clear();
    PFPVertexDistY_TreeVal[fTrackModuleLabel].clear();
    PFPVertexDistZ_TreeVal[fTrackModuleLabel].clear();
    PFPVertexDistMag_TreeVal[fTrackModuleLabel].clear();

    PFPLongestTrkLength_TreeVal[fTrackModuleLabel].clear();
    PFPLongestTrkHits_TreeVal[fTrackModuleLabel].clear();
    
    trkAllHitsPurity_TreeVal[fTrackModuleLabel].clear();
    trkAllHitsComp_TreeVal[fTrackModuleLabel].clear();
    trkEnergyPurityBestPlane_TreeVal[fTrackModuleLabel].clear();
    trkEnergyCompBestPlane_TreeVal[fTrackModuleLabel].clear();
    trkTrackPDG_MostHit_TreeVal[fTrackModuleLabel].clear();
    trkPDG_TreeVal[fTrackModuleLabel].clear();
    trkMother_TreeVal[fTrackModuleLabel].clear();
    trkTrueEnergy_TreeVal[fTrackModuleLabel].clear();
 
    trkSignalToNoiseRatio_TreeVal[fTrackModuleLabel].clear();
    //trkSignalToNoiseRatioU_TreeVal[fTrackModuleLabel].clear();
    //trkSignalToNoiseRatioV_TreeVal[fTrackModuleLabel].clear();
    //trkSignalToNoiseRatioY_TreeVal[fTrackModuleLabel].clear();

    trkTrueKEnergy_TreeVal[fTrackModuleLabel].clear();
    trkStartEndProcess_TreeVal[fTrackModuleLabel].clear();
    trkStartProcess_TreeVal[fTrackModuleLabel].clear();
    trkTrueTrackLength_TreeVal[fTrackModuleLabel].clear();
    trkTrueTrackSpread_TreeVal[fTrackModuleLabel].clear();
    trkTrueTrackResidual_TreeVal[fTrackModuleLabel].clear();
    
    trkDirX_TreeVal[fTrackModuleLabel].clear();
    trkDirY_TreeVal[fTrackModuleLabel].clear();
    trkDirZ_TreeVal[fTrackModuleLabel].clear();
    trkTrueDirX_TreeVal[fTrackModuleLabel].clear();
    trkTrueDirY_TreeVal[fTrackModuleLabel].clear();
    trkTrueDirZ_TreeVal[fTrackModuleLabel].clear();
    trkDirDiff_TreeVal[fTrackModuleLabel].clear();
    trkStartDirX_TreeVal[fTrackModuleLabel].clear();
    trkStartDirY_TreeVal[fTrackModuleLabel].clear();
    trkStartDirZ_TreeVal[fTrackModuleLabel].clear();
    trkStartTrueDirX_TreeVal[fTrackModuleLabel].clear();
    trkStartTrueDirY_TreeVal[fTrackModuleLabel].clear();
    trkStartTrueDirZ_TreeVal[fTrackModuleLabel].clear();
    trkStratDirDiff_TreeVal[fTrackModuleLabel].clear();
    trkEndDirX_TreeVal[fTrackModuleLabel].clear();
    trkEndDirY_TreeVal[fTrackModuleLabel].clear();
    trkEndDirZ_TreeVal[fTrackModuleLabel].clear();
    trkEndTrueDirX_TreeVal[fTrackModuleLabel].clear();
    trkEndTrueDirY_TreeVal[fTrackModuleLabel].clear();
    trkEndTrueDirZ_TreeVal[fTrackModuleLabel].clear();
    trkEndDirDiff_TreeVal[fTrackModuleLabel].clear();
    trkStartX_TreeVal[fTrackModuleLabel].clear();
    trkStartY_TreeVal[fTrackModuleLabel].clear();
    trkStartZ_TreeVal[fTrackModuleLabel].clear();
    trkStartDist_TreeVal[fTrackModuleLabel].clear();
    trkBestPlane_TreeVal[fTrackModuleLabel].clear();
    trkNumHits_TreeVal[fTrackModuleLabel].clear();
    trkID_TreeVal[fTrackModuleLabel].clear();
    trkParticleId_TreeVal[fTrackModuleLabel].clear();
    trkNdof_TreeVal[fTrackModuleLabel].clear();
    trkLength_TreeVal[fTrackModuleLabel].clear();
    trkLengthDiff_TreeVal[fTrackModuleLabel].clear();
    trkGeoProjectionMatched_TreeVal[fTrackModuleLabel].clear();
    trkChi2_TreeVal[fTrackModuleLabel].clear();
    trkStartMomentum_TreeVal[fTrackModuleLabel].clear();
    trkEndMomentum_TreeVal[fTrackModuleLabel].clear();
    trkVertexMomentum_TreeVal[fTrackModuleLabel].clear();
    trkAzimuthAngle_TreeVal[fTrackModuleLabel].clear();
    trkPhi_TreeVal[fTrackModuleLabel].clear();
    trkTheta_TreeVal[fTrackModuleLabel].clear();
    trkZenithAngle_TreeVal[fTrackModuleLabel].clear();
    
    calotrkPlaneId_TreeVal[fTrackModuleLabel].clear();
    
    trkKEBP_TreeVal[fTrackModuleLabel].clear();
    trkKEU_TreeVal[fTrackModuleLabel].clear();
    trkKEV_TreeVal[fTrackModuleLabel].clear();
    trkKEY_TreeVal[fTrackModuleLabel].clear();

    trkRange_TreeVal[fTrackModuleLabel].clear();
    trkPitchY_TreeVal[fTrackModuleLabel].clear();
    
    trkResRangeBP_TreeVal[fTrackModuleLabel].clear();
    trkResRangeU_TreeVal[fTrackModuleLabel].clear();
    trkResRangeV_TreeVal[fTrackModuleLabel].clear();
    trkResRangeY_TreeVal[fTrackModuleLabel].clear();

    trkdEdxVecBP_TreeVal[fTrackModuleLabel].clear();
    trkdEdxVecU_TreeVal[fTrackModuleLabel].clear();
    trkdEdxVecV_TreeVal[fTrackModuleLabel].clear();
    trkdEdxVecY_TreeVal[fTrackModuleLabel].clear();

    trkdQdxVecBP_TreeVal[fTrackModuleLabel].clear();
    trkdQdxVecU_TreeVal[fTrackModuleLabel].clear();
    trkdQdxVecV_TreeVal[fTrackModuleLabel].clear();
    trkdQdxVecY_TreeVal[fTrackModuleLabel].clear();
    
    //trkEnergy_TreeVal[fTrackModuleLabel].clear();
    //trkEnergyRat_TreeVal[fTrackModuleLabel].clear();
    //trkEnergyDiff_TreeVal[fTrackModuleLabel].clear();
    
    trkKinEnergy_TreeVal[fTrackModuleLabel].clear();
    trkKinEnergyRat_TreeVal[fTrackModuleLabel].clear();
    trkKinEnergyDiff_TreeVal[fTrackModuleLabel].clear();

    recotrk_range_p_TreeVal[fTrackModuleLabel].clear();
    recotrk_bestLogLikelihood_p_TreeVal[fTrackModuleLabel].clear();
    recotrk_bestMomUncertainty_p_TreeVal[fTrackModuleLabel].clear();
    recotrk_bestMomentum_p_TreeVal[fTrackModuleLabel].clear();
    recotrk_bwdLogLikelihood_p_TreeVal[fTrackModuleLabel].clear();
    recotrk_bwdMomentum_p_TreeVal[fTrackModuleLabel].clear();
    recotrk_bwdMomUncertainty_p_TreeVal[fTrackModuleLabel].clear();
    recotrk_deltaLogLikelihood_p_TreeVal[fTrackModuleLabel].clear();
    recotrk_fwdLogLikelihood_p_TreeVal[fTrackModuleLabel].clear();
    recotrk_fwdMomentum_p_TreeVal[fTrackModuleLabel].clear();
    recotrk_fwdMomUncertainty_p_TreeVal[fTrackModuleLabel].clear();

    trkdEdxBP_TreeVal[fTrackModuleLabel].clear();
    trkdEdxU_TreeVal[fTrackModuleLabel].clear();
    trkdEdxV_TreeVal[fTrackModuleLabel].clear();
    trkdEdxY_TreeVal[fTrackModuleLabel].clear();
    
    eSegmentation_TreeVal[fTrackModuleLabel].clear();
    eNumRecoTracks_TreeVal[fTrackModuleLabel].clear();

  }// close fTrackModuleLabels loop

  fNumTrueTrks_TreeVal        = -99999;
  fNumTrueTrksviaECut_TreeVal = -99999;
  fNumTrueTrksviaDisCut_TreeVal = -99999;

  for(auto const& trkMother : fMCTrackList){
    TrueHitNum_TreeVal[trkMother].clear();
    TrueTrkE_TreeValMap[trkMother].clear();
    TrueTrkEviaECut_TreeValMap[trkMother].clear();
    TrueTrkEviaDisCut_TreeValMap[trkMother].clear();
  }//close the loop over 

  for(unsigned int hitlab_it=0; hitlab_it<fHitModuleLabels.size(); ++hitlab_it){
    for(unsigned int plane_it=0; plane_it<geom->Nplanes(); ++plane_it){
      hEnergyCompl_TreeVal[fHitModuleLabels[hitlab_it]][plane_it].clear();
    }//close the loop over planes
  }//Close the loop over fHitModuleLabels

  //return;
  //
}//Close the events block

//#################################################################################################
//
//=========================================================
//======= These function can be moved to TrackUtils =======
//=========================================================

void ana::TrackValidation::endJob()
{
  // Implementation of optional member function here.
  std::cout << "Number of events ran over:      " << fnumevents           << std::endl;
  std::cout << "Number of initial MC Tracks:    " << fnumtracks           << std::endl;
  std::cout << "Number of reco tracks:          " << fnumrecotracks       << std::endl;
  std::cout << "Number of reco tracks analysed: " << fnumrecotracksana    << std::endl;
  std::cout << "Number of True Tracks:          " << fTrueParticleTracks  << std::endl;
  std::cout << "Number of True Showers:         " << fTrueParticleShowers << std::endl;
  //fTrueParticleShowers
  //fTrueParticleTracks 

  if(fnumrecotracksana == 0){
    std::ofstream myfile;
    myfile.open ("trackhistcomp.txt");
    myfile << "0";
    myfile.close();

  }
}



DEFINE_ART_MODULE(ana::TrackValidation)
