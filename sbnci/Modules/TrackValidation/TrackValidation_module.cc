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
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "larcoreobj/SimpleTypesAndConstants/PhysicalConstants.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardataalg/DetectorInfo/DetectorClocksData.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "larsim/MCCheater/BackTrackerService.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
//#include "larreco/Calorimetry/CalorimetryAlg.h"
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
#include "sbndcode/RecoUtils/RecoUtils.h"

// Root includes
#include "TTree.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "TVector3.h"
#include "TSystem.h"
#include "TSpline.h"
//#include "TSpline3.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TDatabasePDG.h"

// C++ includes
#include <vector>
#include <string>
#include <utility>
#include <map>
#include <iostream>
#include <fstream>
#include <array>
#include <cassert>
// #include <cmath>  //acos

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

  //Cluster Validation Function
  void ClusterValidation(std::vector< art::Ptr<recob::Cluster> >& clusters,
                       const art::Event& evt,
                       const detinfo::DetectorClocksData& clockData,
                       art::Handle<std::vector<recob::Cluster> >& clusterHandle,
                       std::map<int,std::vector<int> >& TrackMotherTrackIDs,
                       std::map<int,float>& MCTrack_Energy_map,
                       std::map<art::ProductID,std::map<int,std::map<geo::PlaneID, int> > >& MCTrack_hit_map,
                       int& TrueTrackID,
                       float& trueTrackEnergy,
                       const std::string& fTrackModuleLabel);
 
  //PFPValidation Function
  void PFPValidation(std::vector<art::Ptr<recob::Cluster> >& clusters, art::Ptr<recob::PFParticle> pfp,
                     const art::Event& evt,
                     std::map<int, const simb::MCParticle*> TrueParticlesMap,
                     const detinfo::DetectorClocksData& clockData,
                     art::Handle<std::vector<recob::Cluster> >& clusterHandle,
                     std::map<int,std::vector<int> >& TrackMotherTrackIDs,
                     std::map<int,float>& MCTrack_Energy_map,
                     std::map<art::ProductID,std::map<int,std::map<geo::PlaneID, int> > >& MCTrack_hit_map,
                     const std::string& fTrackModuleLabel);

  //Function retunr a map of ShowerMother trackid as a key and a vector of daughters ids
  void ShowerMotherTrackIdsFun(art::Event const& evt,
                               std::map<int, const simb::MCParticle*> TrueParticlesMap,
                               std::map<int, std::vector<int> >& ShwrMothTrkIdWithDauIdsMap);
  
  //Function retunr a map of Track's Mother trackid as a key and a vector of daughters ids
  void MotherTrkIdsTodauIdsFun(art::Event const& evt,
                               std::map<int, const simb::MCParticle*> TrueParticlesMap,
                               std::map<int, std::vector<int> >& MothTrkIdWithDauIdsMap);
  
  // Function retrun the true track id
  int FunTrueParticleID(detinfo::DetectorClocksData const& clockData,
                        const art::Ptr<recob::Hit>& hit) const;
  
  //Function return all true hits of a particle
  void TruthTrkIdsToHits(art::Event const& evt,
                         detinfo::DetectorClocksData const& clockData,
                         int MotherTrkIds,
                         std::vector<int>& DaugthersIdsVec,
                         std::vector<art::Ptr<recob::Hit>> allrecohits,
                         std::vector<art::Ptr<recob::Hit>>& hitsOfTrkIds);
  
  // Function retrun the true track id of hits that deposit the most energy in a shower or track
  int TureParticleIdDepMostEFun(detinfo::DetectorClocksData const& clockData,
                                std::vector<art::Ptr<recob::Hit> > ParticleHits) const;
  //
  //Function return a vector of true hits of a reco particle
  std::vector<art::Ptr<recob::Hit> > TrueHitsOfParticle(detinfo::DetectorClocksData const& clockData,
                                                        int MotherTrkIds,
                                                        std::vector<int>& DaugthersIdsVec,
                                                        std::vector<art::Ptr<recob::Hit> > RecoParticleHits);

  //Function return a map of 
  //Gets the total energy deposited by a track  by looping through the MC info on each wire.
  std::map<int,float> TrueEnergyDepositedFromMCTracks(const std::vector<art::Ptr<sim::SimChannel>> &simchannels);

  //Function remove non contained particles
  void RemoveNoneContainedParticles(std::map<int,std::vector<int> >&  TracksMothers,
                                    std::map<int,const simb::MCParticle*>& trueParticlesMap,
                                    std::map<int,float>& MCTrackEnergyMap);
  //
  //Function cut on the Track mothers with energy
  void CutTrkMothersByE(std::map<int,std::vector<int> >&  TracksMothers,
                                    std::map<int,const simb::MCParticle*>& trueParticlesMap,
                                    float& EnergyCut);
  //
  //Function for cut on track via length 
  void CalculateTrackLengthCut(std::map<int,std::vector<int> >&  TracksMothers,
                               std::map<int,const simb::MCParticle*>& trueParticlesMap,
                               float& LengthCut);
  //Function returns a bool if the track or the shower inside the tpc
  bool IsInsideTPC(TVector3 position, double distance_buffer);
  //Function return a map of number of hits per plane per track
  std::map<int,std::map<geo::PlaneID,int>> NumberofPlaneHitsPerTrack(detinfo::DetectorClocksData const& clockData,
                                                               const std::vector<art::Ptr<recob::Hit> >& hits);

  //Function for finding the track Best plane
  int FindingTrackBeastPlane(detinfo::DetectorClocksData const& clockData,
                             const std::vector<art::Ptr<recob::Hit> >& trkhits);

  //Function to finds the most probable track id associated with the set of hits from 
  //there true energy depositions. The pair returns the energy as well.
  std::pair<int,double> pairTrackTrackIdToEnergyDep(detinfo::DetectorClocksData const& clockData,
                                                    std::map<int,std::vector<int>> &TracksMothers,
                                                    const std::vector<art::Ptr<recob::Hit> >& hits,
                                                    int planeid);

  //Function return the track id of particle that has the most hits in a track.
  int TrueParticleIDFromTotalRecoHits(detinfo::DetectorClocksData const& clockData,
                                      const std::vector<art::Ptr<recob::Hit> >& hits);
  
  //Function return the track id of particle that deposite the most energy in the track
  int TrueParticleIDFromTotalTrueEnergy(detinfo::DetectorClocksData const& clockData,
                                        const std::vector<art::Ptr<recob::Hit> >& hits);

  //Function returns a float value of total energy deposited from a vector of hits for a specific plane
  float TotalEnergyDepinHitsFun(detinfo::DetectorClocksData const& clockData,
                             const std::vector<art::Ptr<recob::Hit> >& hits,
                             int Plane);

private:

  // Declare member data here.
  
  //fcl parameters
  std::string fGenieGenModuleLabel;
  std::string fLArGeantModuleLabel;
  std::string fHitsModuleLabel;
  //std::string fTrkModuleLabel;
  std::string fCalorimetryLabel;
  std::string fShowerModuleLabel;
  std::string fPFParticleLabel;
  std::string fpandoraTrackMCS;
  std::string fpandoraTrackRange;
  std::vector<std::string> fTrackModuleLabels;
  std::vector<std::string> fHitModuleLabels;

  bool  fRemoveNonContainedParticles;
  bool  fPFPValidation;
  bool  fUseNeutrinoTracks;
  bool  fClusterValidation;
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
  std::vector<std::string> fTrackMotherList;
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
  //PFParticle Validation, tracks and showers
  std::map<std::string,std::vector<float> > PFProjectionMatched_TreeVal;
  std::map<std::string,std::vector<float> > PFPShowerProjectionMatched_TreeVal;
  std::map<std::string,std::vector<float> > PFPTrackProjectionMatched_TreeVal;

  std::map<std::string,std::vector<float> > PFPHitsCompl_TreeVal;
  std::map<std::string,std::vector<float> > PFPHitsPurity_TreeVal;
  std::map<std::string,std::vector<float> > PFPEnergyCompl_TreeVal;
  std::map<std::string,std::vector<float> > PFPEnergyPurity_TreeVal;

  std::map<std::string,std::vector<float> > PFPTrackHitsCompl_TreeVal;
  std::map<std::string,std::vector<float> > PFPTrackHitsPurity_TreeVal;
  std::map<std::string,std::vector<float> > PFPTrackEnergyCompl_TreeVal;
  std::map<std::string,std::vector<float> > PFPTrackEnergyPurity_TreeVal;
  
  std::map<std::string,std::vector<float> > PFPShowerHitsCompl_TreeVal;
  std::map<std::string,std::vector<float> > PFPShowerHitsPurity_TreeVal;
  std::map<std::string,std::vector<float> > PFPShowerEnergyCompl_TreeVal;
  std::map<std::string,std::vector<float> > PFPShowerEnergyPurity_TreeVal;

  //Vertex information
  std::map<std::string,std::vector<float> > PFPVertexDistX_TreeVal;
  std::map<std::string,std::vector<float> > PFPVertexDistY_TreeVal;
  std::map<std::string,std::vector<float> > PFPVertexDistZ_TreeVal;
  std::map<std::string,std::vector<float> > PFPVertexDistMag_TreeVal;
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
  //std::map<std::string,std::vector<float> > trkEnergy_TreeVal;
  //std::map<std::string,std::vector<float> > trkEnergyRat_TreeVal;
  //std::map<std::string,std::vector<float> > trkEnergyDiff_TreeVal;
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
  std::map<std::string,std::vector<float> > eNumRecoShowers_TreeVal;
  std::map<std::string,std::vector<float> > eNumRecoTracks_TreeVal;

  //Cluster informations
  //std::map<std::string,std::vector<std::vector<std::vector<float> > > > cProjectionMatchedEnergy_TreeVal;
  std::map<std::string,std::vector<std::vector<std::vector<float> > > > cEnergyCompl_TreeVal;
  std::map<std::string,std::vector<std::vector<std::vector<float> > > > cEnergyPurity_TreeVal;
  std::map<std::string,std::vector<std::vector<std::vector<float> > > > cHitsCompl_TreeVal;
  std::map<std::string,std::vector<std::vector<std::vector<float> > > > cHitsPurity_TreeVal;
  
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
  fShowerModuleLabel   = p.get<std::string>("ShowerModuleLabel"   );//
  fPFParticleLabel     = p.get<std::string>("PFParticleLabel"     );//
  fpandoraTrackMCS     = p.get<std::string>("pandoraTrackMCS"    );//
  fpandoraTrackRange   = p.get<std::string>("pandoraTrackRange"  );//
  fTrackModuleLabels   = p.get<std::vector<std::string> >("TrackModuleLabels");//
  fHitModuleLabels     = p.get<std::vector<std::string> >("HitModuleLabels"  );//

  fRemoveNonContainedParticles = p.get<bool>("RemoveNonContainedParticles" );//
  fPFPValidation               = p.get<bool>("PFPValidation"               );//
  fUseNeutrinoTracks           = p.get<bool>("UseNeutrinoTracks"           );//
  fClusterValidation           = p.get<bool>("ClusterValidation"           );//
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
  fTrackPdgCodeVect = p.get<std::vector<unsigned int> >("TrackPdgCodeVect" );//
  fUseBestPlane                = p.get<int>("UseBestPlane"                 );//

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
  fTrackMotherList.clear();

  //
  //
  fTree = tfs->make<TTree>("MetricTree", "Tree Holding all metric information");
  gInterpreter->GenerateDictionary("vector<vector<float> > ","vector");

  fTree->Branch("EventRun",    &fEventRun_TreeVal, 32000, 0);
  fTree->Branch("EventSubrun", &fEventSubrun_TreeVal, 32000, 0);
  fTree->Branch("EventNumber", &fEventNumber_TreeVal, 32000, 0);
  fTree->Branch("NoShowerRecoAsTrack", &fNoShowerTrack_TreeVal, 32000, 0);

  initTree(fTree,"TrueTrkE",TrueTrkE_TreeValMap,fTrackMotherList);//TODO add it to mc tree and add pdgcode
  initTree(fTree,"TrueTrkEviaECut",TrueTrkEviaECut_TreeValMap,fTrackMotherList);
  initTree(fTree,"TrueTrkEviaDisCut",TrueTrkEviaDisCut_TreeValMap,fTrackMotherList);

  initTree(fTree,"trkAllHitsPurity",trkAllHitsPurity_TreeVal,fTrackModuleLabels);
  initTree(fTree,"trkAllHitsComp",trkAllHitsComp_TreeVal,fTrackModuleLabels);
  initTree(fTree,"trkEnergyPurityBestPlane",trkEnergyPurityBestPlane_TreeVal,fTrackModuleLabels);
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
  
  //initTree(fTree,"recoTrkEnergy",trkEnergy_TreeVal,fTrackModuleLabels);
  //initTree(fTree,"recoTrkEnergyRat",trkEnergyRat_TreeVal,fTrackModuleLabels);
  //initTree(fTree,"recoTrkEnergyDiff",trkEnergyDiff_TreeVal,fTrackModuleLabels);

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
  //initTree(fTree,"",,fTrackModuleLabels);

  if(fPFPValidation){
    initTree(fTree,"PFProjectionMatched",PFProjectionMatched_TreeVal,fTrackModuleLabels);
    initTree(fTree,"PFPHitsCompl",PFPHitsCompl_TreeVal,fTrackModuleLabels);
    initTree(fTree,"PFPHitsPurity",PFPHitsPurity_TreeVal,fTrackModuleLabels);
    initTree(fTree,"PFPEnergyCompl",PFPEnergyCompl_TreeVal,fTrackModuleLabels);
    initTree(fTree,"PFPEnergyPurity",PFPEnergyPurity_TreeVal,fTrackModuleLabels);
    initTree(fTree,"PFPTrackProjectionMatched",PFPTrackProjectionMatched_TreeVal,fTrackModuleLabels);
    initTree(fTree,"PFPTrackHitsCompl",PFPTrackHitsCompl_TreeVal,fTrackModuleLabels);
    initTree(fTree,"PFPTrackHitsPurity",PFPTrackHitsPurity_TreeVal,fTrackModuleLabels);
    initTree(fTree,"PFPTrackEnergyCompl",PFPTrackEnergyCompl_TreeVal,fTrackModuleLabels);
    initTree(fTree,"PFPTrackEnergyPurity",PFPTrackEnergyPurity_TreeVal,fTrackModuleLabels);
    initTree(fTree,"PFPShowerProjectionMatched",PFPShowerProjectionMatched_TreeVal,fTrackModuleLabels);
    initTree(fTree,"PFPShowerHitsCompl",PFPShowerHitsCompl_TreeVal,fTrackModuleLabels);
    initTree(fTree,"PFPShowerHitsPurity",PFPShowerHitsPurity_TreeVal,fTrackModuleLabels);
    initTree(fTree,"PFPShowerEnergyCompl",PFPShowerEnergyCompl_TreeVal,fTrackModuleLabels);
    initTree(fTree,"PFPShowerEnergyPurity",PFPShowerEnergyPurity_TreeVal,fTrackModuleLabels);
  }

  initTree(fTree,"eSegmentation",eSegmentation_TreeVal,fTrackModuleLabels);
  initTree(fTree,"eNumRecoShowers",eNumRecoShowers_TreeVal,fTrackModuleLabels);
  initTree(fTree,"eNumRecoTracks",eNumRecoTracks_TreeVal,fTrackModuleLabels);
 
  initTree(fTree,"fNumTrueTrks",fNumTrueTrks_TreeVal);
  initTree(fTree,"fNumTrueTrksviaECut",fNumTrueTrksviaECut_TreeVal);
  initTree(fTree,"fNumTrueTrksviaDisCut",fNumTrueTrksviaDisCut_TreeVal);

  if (fClusterValidation){
    initClusterTree(fTree,"cEnergyCompl",cEnergyCompl_TreeVal,fTrackModuleLabels);
    initClusterTree(fTree,"cEnergyPurity",cEnergyPurity_TreeVal,fTrackModuleLabels);
    initClusterTree(fTree,"cHitsCompl",cHitsCompl_TreeVal,fTrackModuleLabels);
    initClusterTree(fTree,"cHitsPurity",cHitsPurity_TreeVal,fTrackModuleLabels);
  }

  for(auto const& fHitModuleLabel : fHitModuleLabels){
    std::cout << geom->Nplanes() << std::endl;
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
  auto const detProp   = art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataFor(evt, clockData);
  //auto const detProp   = art::ServiceHandle<detinfo::DetectorPropertiesData const>()->DataFor(evt, clockData);
  //Access Event Information
  fEventRun_TreeVal = evt.run();
  fEventSubrun_TreeVal = evt.subRun();
  fEventNumber_TreeVal = evt.event();
  //int fEventID = evt.id().event(); //it gives same output to the evt.event()
  //std::cout << "Evnet ID  = " << fEventID << std::endl;
  std::cout<<"Analysing Event: "<<fEventNumber_TreeVal<<":"<<fEventSubrun_TreeVal<<":"<<fEventRun_TreeVal<<std::endl;
  
  ++fnumevents;

  //-------------------------- Getting  MC truth information -------------------------- 
  art::Handle< std::vector<simb::MCTruth> > mctruthHandle;
  std::vector<art::Ptr<simb::MCTruth> > mclist;
  if(evt.getByLabel(fGenieGenModuleLabel, mctruthHandle))// Make sure the handle is valid
  {art::fill_ptr_vector(mclist, mctruthHandle);}            // Fill the vector with art::Ptr hits

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

  //------------------------ getting the reconstructed showers -----------------------
  art::Handle<std::vector<recob::Shower> > showerHandle;
  std::vector<art::Ptr<recob::Shower> > showerlist;
  if(evt.getByLabel(fShowerModuleLabel, showerHandle)){// Make sure the handle is valid
    art::fill_ptr_vector(showerlist, showerHandle);}// Fill the vector with art::Ptr hits
  const size_t NShowers = showerlist.size();
  std::cout << "The shower vector size = " << NShowers << std::endl;
  if(NShowers ==0){std::cout << "There is no shower in this event" << std::endl;}
  //----------------------------------------------------------------------------------------------

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
  std::map<int, std::vector<int> > ShwrMotherTrackIdMap;//Map of mother to daughters trackids for shower particels
  std::map<int, std::vector<int> > TrkMotherTrackIdMap;//Map of mother to daughters trackids for track particels
  std::map<int, float> MCTrackEnergyMap;               //Map of the MC Energy deposited for each MC track
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
      //if(ParticlePdgCode ==13 || ParticlePdgCode == 211 || ParticlePdgCode == 2212 || ParticlePdgCode == 321){
      int number_of_tracks = 0;
      if(std::find(fTrackPdgCodeVect.begin(),fTrackPdgCodeVect.end(), ParticlePdgCode) != fTrackPdgCodeVect.end()){
        truePrimaryTrackMap[particle->TrackId()] = particle;
        ++number_of_tracks;
        std::cout<<"This working and push one more value" << truePrimaryTrackMap.size() << " : " << number_of_tracks <<std::endl;
      }//
    //
    }
  }// close the loop over the ParticleList

  if(fMatchTracksInTruth){//open truth information block
    ShowerMotherTrackIdsFun(evt, trueparticlesMap, ShwrMotherTrackIdMap);
    MotherTrkIdsTodauIdsFun(evt, trueparticlesMap, TrkMotherTrackIdMap);
    //Getting the number of paricle shower-like and track-like per event
    fTrueParticleShowers = truePrimaryShowerMap.size();
    fTrueParticleTracks = truePrimaryTrackMap.size();
    //Get the MC Energy deposited for each MC track.
    MCTrackEnergyMap = TrueEnergyDepositedFromMCTracks(simchannels);
    std::cout<< " MCTrackEnergyMap: "<< MCTrackEnergyMap[1] <<std::endl;
    //
    //Find number of true tracks per events
    fNumTrueTrks_TreeVal = TrkMotherTrackIdMap.size();
    //Loop over Tracks mother and get their energy
    for(auto const& trkmother : TrkMotherTrackIdMap){
      const simb::MCParticle* trackMotherParticle = trueparticlesMap[(trkmother.first)];
      //Get a vector of track mothers 
      fTrackMotherList.push_back((std::to_string(trkmother.first)));
      //Save the energy of tracks with their track 
      TrueTrkE_TreeValMap[(std::to_string(trkmother.first))].push_back(trackMotherParticle->E());
    }//close the loop over trkmothers map
    if(fVerbose > 1){
      std::cout<<"Inital track candidates: "<< TrkMotherTrackIdMap.size() <<std::endl;
    }//Close if statement
    if(fRemoveNonContainedParticles){
      RemoveNoneContainedParticles(TrkMotherTrackIdMap, trueparticlesMap, MCTrackEnergyMap);
    }//Close if statement 
    if(fVerbose > 1){
      std::cout<<" Track candidates after containment cut are: "<<TrkMotherTrackIdMap.size()<<std::endl;
    }//Close if statement

    //Cut on the Track mothers with energy
    CutTrkMothersByE(TrkMotherTrackIdMap, trueparticlesMap, fSimEnergyCut);
    fNumTrueTrksviaECut_TreeVal = TrkMotherTrackIdMap.size();
    for(auto const& trkmother : TrkMotherTrackIdMap){
      const simb::MCParticle* TrackMotherParticle = trueparticlesMap.at(trkmother.first);
      TrueTrkEviaECut_TreeValMap[(std::to_string(trkmother.first))].push_back(TrackMotherParticle->E());
    }//close the loop over trkmothers after energy cut
    if(fVerbose > 1){
        std::cout<<" Track candidates after E cut are: "<<TrkMotherTrackIdMap.size()<<std::endl;
    }//Close if statement
    //Cut on track mothers with length and uncontained tracks
    CalculateTrackLengthCut(TrkMotherTrackIdMap, trueparticlesMap, fLengthCut);
    //
    fNumTrueTrksviaDisCut_TreeVal = TrkMotherTrackIdMap.size();
    for(auto const& trkmother : TrkMotherTrackIdMap){
      const simb::MCParticle* TrackMotherParticle = trueparticlesMap.at(trkmother.first);
      //Select first traj point where the track loses energy, last be default
      ////Get the number of Traj points to loop over
      unsigned int MCNumTrajPoints     = TrackMotherParticle->NumberTrajectoryPoints();
      TLorentzVector PositionTrajStart = TrackMotherParticle->Position(0);
      TLorentzVector PositionTrajEnd   = TrackMotherParticle->Position(MCNumTrajPoints-1);
      //Get the track length of the track.
      double MotherTrkLength = (PositionTrajStart - PositionTrajEnd).Vect().Mag(); //in cm
      //
      TrueTrkEviaDisCut_TreeValMap[(std::to_string(trkmother.first))].push_back(MotherTrkLength);
    }
    if(fVerbose > 1){
      std::cout<< "Track candidates after length cut are: " << TrkMotherTrackIdMap.size() <<std::endl;
    }

    //If there are no true tracks we can't validate
    if(TrkMotherTrackIdMap.size() == 0){
      if(fVerbose > 0){std::cout << "No track mothers finishing" << std::endl;}
      return;
    }
    //
    if(fVerbose > 2){
      for(auto const& trkmother : TrkMotherTrackIdMap){
        const simb::MCParticle* TrackMotherParticle = trueparticlesMap.at(trkmother.first);//test
        std::cout << " A Mother has track id: " << trkmother.first << " with Energy: "<< TrackMotherParticle->E() << std::endl;
        for (auto const& daughterID : trkmother.second){
          std::cout << "  Daughter of id: " << daughterID << std::endl;
        }//close the loop over track daughter
      }//close the loop over tracks mothers 
    }//close fVerbose > 2 

    numTrueTraks = TrkMotherTrackIdMap.size();
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
      MCTrack_hit_map[handle.id()] = NumberofPlaneHitsPerTrack(clockData, hits_fromhandle);
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

        float TotalEnergyDepinHits = TotalEnergyDepinHitsFun(clockData, allHits, plane_id.Plane);
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

    if(TrkMotherTrackIdMap.size() == 0 && fMatchTracksInTruth){
      if(fVerbose > 0){std::cout << "No Mothers to Match to the Tracks" << std::endl;}
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
          if(fPFPValidation){
            if(fmpfc.isValid() && fmpfc.size()!=0){
              art::Handle<std::vector<recob::Cluster > > clusterHandle;
              evt.get(fmpfc.at(0).front().id(),clusterHandle);
              if(clusterHandle.isValid()){
                std::vector<art::Ptr<recob::Cluster> > pfpClusters = fmpfc.at(daughter.key());
                ana::TrackValidation::PFPValidation(pfpClusters, daughter, evt, trueparticlesMap, clockData, clusterHandle, TrkMotherTrackIdMap, MCTrackEnergyMap, MCTrack_hit_map, fTrackModuleLabel);
              }//close checking if cluster handle is valid
            }//close checking if cluseter is valid
          }//colse the fPFPValidation

          // Split into shower like, 11, and track like, 13
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
    unsigned int max_hitnum=0;
    art::Ptr<recob::Track> longestTrack;

    //Find the longest track (by number of hits)
    for(auto const& trk : recoTracks){
      ++fnumrecotracks;
      //Get the hits vector from the track
      std::vector<art::Ptr<recob::Hit> > trkhits = fmtrkh.at(trk.key());
      if(trkhits.size() > max_hitnum){
        max_hitnum = trkhits.size();
        longestTrack = trk;
      }//close if statement
    }//close the loop over tracks
    //--------------------------------------------------------------------------------------------

    //Loop over tracks to get reco information
    for(auto const& trk : recoTracks){
      //Get the information for the track
      //Getting the Azimuth Angle (The angle is measured on the plane orthogonal to the z axis,
      //The angle is measured counterclockwise from the x axis) basiclly the angle between x and z.
      double trkAzimuthAngle = trk->AzimuthAngle();
      //Getting the Chi^2 
      float trkChi2 = trk->Chi2();
      //Getting the track end position
      TVector3 trkEndPosition = trk->End<TVector3>();
      //Getting the track start direction and end direction
      TVector3 TrackStartDirection, TrackEndDirection;
      std::tie(TrackStartDirection, TrackEndDirection) = trk->Direction<TVector3>();
      //Getting the track end directin using EndDirection function
      TVector3 trkenddir = trk->EndDirection<TVector3>();
      //Access to track end momentum.
      double trkEndMomentum = trk->EndMomentum();
      //Momentum vector at end point.
      TVector3 trkEndMomentumVec = trk->EndMomentumVector<TVector3>();
      //Access position at start and end points
      TVector3 trkStartPosition, trkendpos;
      std::tie(trkStartPosition, trkendpos) = trk->Extent<TVector3>();
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
      //Track start position 
      TVector3 TrkStartPos = trk->Start<TVector3>();
      //getting track start direction
      TVector3 trkstartdir = trk->StartDirection<TVector3>();
      //Getting Track Start Momentum
      double trkstartmom = trk->StartMomentum();
      //Getting Track Start Momentum Vector
      TVector3 trkstartmomVec = trk->StartMomentumVector<TVector3>();
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
      int BestPlaneID = FindingTrackBeastPlane(clockData, trkhits);
      if(BestPlaneID > 2){
        BestPlaneID = 2;
      } 
      
      //================================================
      //=== Hit Completeness and purity Calculations ===
      //================================================
      
      std::pair<int,float> TrktrackIdInfo = pairTrackTrackIdToEnergyDep(clockData, TrkMotherTrackIdMap, trkhits, BestPlaneID);

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
        for(auto const& trkDaughter: TrkMotherTrackIdMap[TrktrackID]){
          TrueEnergyDep_FromTrack += MCTrackEnergyMap[trkDaughter];
          for(geo::PlaneID planeId : geom->IteratePlaneIDs()){
            TrueHitDep_FromTrueTrack += MCTrack_hit_map[trkhit_productid][trkDaughter][planeId];
          }//close the loop over plane id from geom
        }//Close the loop over the tracks daughters 

        //Get the number of hits in the reco track from the true track.
        std::map<int,std::map<geo::PlaneID,int> > MCTrack_trkhit_map = NumberofPlaneHitsPerTrack(clockData, trkhits);
        for(auto const& planehitMap : MCTrack_trkhit_map){
          if(std::find(TrkMotherTrackIdMap[TrktrackID].begin(), TrkMotherTrackIdMap[TrktrackID].end(), planehitMap.first) == TrkMotherTrackIdMap[TrktrackID].end()){continue;}
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
        //Getting the position at first 20 traj and last 20 traj
        TLorentzVector PositionTrajStartEnd = {-99999,-99999,-99999, -99999};
        TLorentzVector PositionTrajEndStart = {-99999,-99999,-99999, -99999};

        const simb::MCParticle* MCTrackParticle = trueparticlesMap.at(TrktrackID);

        int track_pdgcode = MCTrackParticle->PdgCode();

        //Check if another particle actually gave more energy.
        int MostHitsTrackID = std::abs(TrueParticleIDFromTotalRecoHits(clockData, trkhits));
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

        //Getting the position of trajectory point for first 20 and last 20
        PositionTrajStartEnd = MCTrackParticle->Position(20);
        PositionTrajEndStart = MCTrackParticle->Position(trkTrajPoints-21);

        //Get the track length of the track.
        TrueTrackLength = (PositionTrajStart - PositionTrajEnd).Vect().Mag(); //in cm
        //double TrackLengthStart = (PositionTrajStart - PositionTrajStartEnd).Vect().Mag(); TODO
        //double TrackLengthEnd   = (PositionTrajEndStart - PositionTrajEnd).Vect().Mag(); TODO

        //How straight are the IDE points. First Define the line.
        TVector3 track_line = (PositionTrajStart - PositionTrajEnd).Vect().Unit();
        TVector3 track_line_start = (PositionTrajStart - PositionTrajStartEnd).Vect().Unit();
        TVector3 track_line_end = (PositionTrajEndStart - PositionTrajEnd).Vect().Unit();
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
              if(dEdx > 10 || dEdx < 0){continue;}
              ++numdEdx_mean;
              trkdEdx_mean += dEdx;
            }//close the loop over track dEdx
            double trkdEdx = (trkdEdx_mean/(float)(numdEdx_mean));
            PlanetrkdEdx[PlaneID] = trkdEdx;
          }//Close else statement
        }//close if statement to get the track dEdx mean or median
      }

      //std::cout<<" Plane trkdEdx: "<<PlanetrkdEdx.size()<<std::endl;
      
      // the track energy
      /*if(TrackEnergyPlanes.size()!=0){
        recotrkEnergy = (TrackEnergyPlanes.at(BestPlaneID));
      }
      //recotrkEnergy = trk_max_energy;
      
      if(TrueEnergyDep_FromTrack != 0 && recotrkEnergy > 0){
        recotrkEnergyRat = recotrkEnergy / TrueEnergyDep_FromTrack;
        recotrkEnergyDiff = (recotrkEnergy - TrueEnergyDep_FromTrack)/TrueEnergyDep_FromTrack;
      } */

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
          static const std::vector<std::string> PIDnames {"muon", "pion", "kaon", "proton"};
          for(std::string pid: PIDnames){
            art::InputTag tag(fpandoraTrackMCS, pid);
            fmMCSs.push_back(art::FindManyP<recob::MCSFitResult>(trackListHandle, evt, tag));
          }

          std::vector<art::FindManyP<sbn::RangeP>> fmRanges;
          static const std::vector<std::string> rangePIDnames {"muon", "pion", "proton"};
          for(std::string pid: rangePIDnames){
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
        std::cout<<"Number of True Tracks in the Event:---"<<TrkMotherTrackIdMap.size()<<std::endl;

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


      //==========================
      //=== Cluster Validation ===
      //==========================
      
      if(fClusterValidation){
        //Get the clusters associated to the tracks.
        if(fmtrkch.isValid()){
          art::Handle<std::vector<recob::Cluster > > clusterHandle;
          evt.get(fmtrkch.at(trk.key()).front().id(),clusterHandle);
          if(clusterHandle.isValid()){
            std::vector<art::Ptr<recob::Cluster> > trackclusters = fmtrkch.at(trk.key());
            ana::TrackValidation::ClusterValidation(trackclusters,evt,clockData,clusterHandle,TrkMotherTrackIdMap,
                                                    MCTrackEnergyMap,MCTrack_hit_map,TrktrackID,trueTrkEnergy,
                                                    fTrackModuleLabel); //TODO
          }//close the statement that check if the cluster handle is valid or not

          else{
            mf::LogError("TrackValidation")<<"Cluster handle is stale. No clustering validation done"<<std::endl;
          }//close else statement
        }//close the statement that check if association between track and clusters working

        else if(fmtrkpf.isValid()){
          //Find the Clusters associated to PF particle. fmpfc
          art::Handle<std::vector<recob::PFParticle> > pfpHandle;
          evt.get(fmtrkpf.at(trk.key()).front().id(),pfpHandle);
          if(pfpHandle.isValid()){
            art::FindManyP<recob::Cluster> fmcpf(pfpHandle, evt, pfpHandle.provenance()->moduleLabel());
            if(fmcpf.isValid()){
              art::Handle<std::vector<recob::Cluster > > clusterHandle;
              evt.get(fmcpf.at(0).front().id(),clusterHandle);
              if(clusterHandle.isValid()){
                std::vector<art::Ptr<recob::Cluster> > trackclusters = fmcpf.at(trk.key());
                ana::TrackValidation::ClusterValidation(trackclusters,evt,clockData,clusterHandle,
                                                        TrkMotherTrackIdMap,MCTrackEnergyMap,MCTrack_hit_map,
                                                        TrktrackID,trueTrkEnergy,fTrackModuleLabel); //TODO
              }//close the statement of clusterHandle.isValid
              else{
                mf::LogError("TrackValidation")<<" Cluster handle is stale. No clustering validation done"<<std::endl;
              }//close else statement of fmcpf.isValid
            }//close the statement of fmcpf.isValid
            else{
              mf::LogError("TrackValidation")<<"No Assosoication between pf particles and clusters was found for track made from a pf particle. No clustering validation was done." <<std::endl;
            }//close else statement
          }// close the statment of pfpHandle.isValid
          else{
            mf::LogError("TrackValidation")<<"pf particle handle is stale" << std::endl;
          }//close the else statement of fmtrkpf.isValid
        }//close the statement of fmtrkpf.isValid
        else{
          mf::LogError("TrackValidation") << "No cluster or pandora association" << std::endl;
        }

        if(fVerbose > 1){std::cout << "Cluster Validation Complete" << std::endl;}
      
      }// close the statement of fClusterValidation

      ++fnumrecotracksana;
      //
    }//close the loop over tracks

    //========================
    //===   Event Metric   ===
    //========================

    //Calculate the True Hit number
    for(auto const& trackmother: TrkMotherTrackIdMap){
      int TrueHitDep_FromTrueTracks = 0;
       for(auto const& track_id: trackmother.second){
         for(geo::PlaneID plane_id: geom->IteratePlaneIDs()){
           TrueHitDep_FromTrueTracks += MCTrack_hit_map[trkhit_productid][track_id][plane_id];
         }
       }//Close the second loop
       //TODO make it as map with track id 
       TrueHitNum_TreeVal[(std::to_string(trackmother.first))].push_back(TrueHitDep_FromTrueTracks);
    }//Close the loop over track mothers

    //This is rather pandora track specific for tuning pandora so lets add it if the track module 
    //used pandora and the track module Pandora track is used.
    eNumRecoShowers_TreeVal[fTrackModuleLabel].push_back(showerlist.size());

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
    PFProjectionMatched_TreeVal[fTrackModuleLabel].clear();
    PFPShowerProjectionMatched_TreeVal[fTrackModuleLabel].clear();
    PFPTrackProjectionMatched_TreeVal[fTrackModuleLabel].clear();
    PFPHitsCompl_TreeVal[fTrackModuleLabel].clear();
    PFPHitsPurity_TreeVal[fTrackModuleLabel].clear();
    PFPEnergyCompl_TreeVal[fTrackModuleLabel].clear();
    PFPEnergyPurity_TreeVal[fTrackModuleLabel].clear();
    PFPTrackHitsCompl_TreeVal[fTrackModuleLabel].clear();
    PFPTrackHitsPurity_TreeVal[fTrackModuleLabel].clear();
    PFPTrackEnergyCompl_TreeVal[fTrackModuleLabel].clear();
    PFPTrackEnergyPurity_TreeVal[fTrackModuleLabel].clear();
    PFPShowerHitsCompl_TreeVal[fTrackModuleLabel].clear();
    PFPShowerHitsPurity_TreeVal[fTrackModuleLabel].clear();
    PFPShowerEnergyCompl_TreeVal[fTrackModuleLabel].clear();
    PFPShowerEnergyPurity_TreeVal[fTrackModuleLabel].clear();
    PFPVertexDistX_TreeVal[fTrackModuleLabel].clear();
    PFPVertexDistY_TreeVal[fTrackModuleLabel].clear();
    PFPVertexDistZ_TreeVal[fTrackModuleLabel].clear();
    PFPVertexDistMag_TreeVal[fTrackModuleLabel].clear();

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
    eNumRecoShowers_TreeVal[fTrackModuleLabel].clear();
    eNumRecoTracks_TreeVal[fTrackModuleLabel].clear();

    if(fClusterValidation){
      for(unsigned int plane=0; plane<geom->Nplanes(); plane++){
        cEnergyCompl_TreeVal[fTrackModuleLabel][plane].clear();
        cEnergyPurity_TreeVal[fTrackModuleLabel][plane].clear();
        cHitsCompl_TreeVal[fTrackModuleLabel][plane].clear();
        cHitsPurity_TreeVal[fTrackModuleLabel][plane].clear();
      }//close the loop over planes
    }//Close fClusterValidation
  }// close fTrackModuleLabels loop

  fNumTrueTrks_TreeVal        = -99999;
  fNumTrueTrksviaECut_TreeVal = -99999;
  fNumTrueTrksviaDisCut_TreeVal = -99999;

  for(auto const& trkMother : fTrackMotherList){
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
//======================================================
//======= These Function should be in the module =======
//======================================================

//Cluster Validation Function
void ana::TrackValidation::ClusterValidation(std::vector< art::Ptr<recob::Cluster> >& clusters,
                           const art::Event& evt,
                           const detinfo::DetectorClocksData& clockData,
                           art::Handle<std::vector<recob::Cluster> >& clusterHandle,
                           std::map<int,std::vector<int> >& TrackMotherTrackIDs,
                           std::map<int,float>& MCTrack_Energy_map,
                           std::map<art::ProductID,std::map<int,std::map<geo::PlaneID, int> > >& MCTrack_hit_map,
                           int& TrueTrackID,
                           float& trueTrackEnergy,
                           const std::string& fTrackModuleLabel){
  
  //
  //Initialise Trees
   unsigned int numPlanes = (int)geom->Nplanes();

   std::vector<std::vector<float> > cProjectionMatchedEnergy_TreeVec(numPlanes);
   std::vector<std::vector<float> > cEnergyCompl_TreeVec(numPlanes);
   std::vector<std::vector<float> > cEnergyPurity_TreeVec(numPlanes);
   std::vector<std::vector<float> > cHitsCompl_TreeVec(numPlanes);
   std::vector<std::vector<float> > cHitsPurity_TreeVec(numPlanes);

   //Get the associated hits
   art::FindManyP<recob::Hit> fmhc(clusterHandle, evt, clusterHandle.provenance()->moduleLabel());

   //Holder for cluster hits
   art::Handle<std::vector<recob::Hit > > hitHandle;
   // std::cout<<"Clusters size: "<<clusters.size()<<" and "<<fmhc.at(clusters.at(0).key()).size()
   // <<" and "<<clusterHandle.provenance()->moduleLabel()<<std::endl;
   
   if(fmhc.at(clusters.at(0).key()).size()==0){
     mf::LogError("TrackValidation")<<"Cluster has no hits. Trying next cluster."<<std::endl;
     return;
   }//close fmhc statement

   //Get the Hits Handle used for this cluster type WARNING
   evt.get(fmhc.at(clusters.at(0).key()).front().id(),hitHandle);

   if(!hitHandle.isValid()){
     mf::LogError("TrackValidation")<<"Hits handle is stale. No clustering validation done"<<std::endl;
     return;
   }

   //Get the hits vector from the track
   for(auto const& cluster : clusters){
     
     std::vector< art::Ptr<recob::Hit> > clusterhits = fmhc.at(cluster.key());
     //Function, finds the most probable track ID associated with the set of hits from there true energy 
     //depositons. The pair returns the energy in the hits as well.
     std::pair<int,float> TrktrackIdInfo = pairTrackTrackIdToEnergyDep(clockData, 
                                           TrackMotherTrackIDs, clusterhits, (cluster->Plane().Plane));

     //Make sure the cluster has been matched.
     if(TrktrackIdInfo.second == -99999){
       if(fVerbose > 0){
         std::cout << "Reco cluster not matched to a True track"<<std::endl;
       }//close fVerbose
       continue;
     }//close the statement of TrktrackIdInfo
     
     float TotalTrueEnergy = 0;
     float signalhits      = 0;
     float totalhits       = 0;

     std::map<int, std::map<geo::PlaneID, int> > hitPlaneMap = NumberofPlaneHitsPerTrack(clockData, clusterhits);

     //Close the loop over track daughter
     for(auto const& daughterID : TrackMotherTrackIDs[TrktrackIdInfo.first]){
       
       //Calculate the true Energy deposited By Track
       TotalTrueEnergy += MCTrack_Energy_map[daughterID];
       
       // Get the map of planeIDs to number of hits
       auto hitPlaneMapIter = hitPlaneMap.find(daughterID);

       if(hitPlaneMapIter == hitPlaneMap.end()){
         if(fVerbose>2){
           std::cout<<"daughter: "<<daughterID<<" not found in cluster hit map"<<std::endl;
         }
         continue;
       }
       
       for(const auto& hitPlaneIt: hitPlaneMapIter->second){
         //Count how many hits are from the true track.
         signalhits += hitPlaneIt.second;
         //Count how many hits are missed in the plane from the true track id.
         totalhits += MCTrack_hit_map[hitHandle.id()][daughterID][hitPlaneIt.first];
       }//close the loop over hit planes
     }//Close the loop over track daughters

     int projection_match = -99999;

     //Have we matched the 2D cluster to the correct track correctly. In terms of Energy depositions:
     if(TrktrackIdInfo.first == TrueTrackID){projection_match = 1;}
     else{projection_match = 0;}

     //Calculate the purity and completeness metrics
    float completeness_hits   = -99999;
    float purity_hits         = -99999;
    float completeness_energy = -99999;
    float purity_energy       = -99999;

    float TotalEnergyDepinHits = TotalEnergyDepinHitsFun(clockData, clusterhits,(cluster->Plane().Plane));

    cProjectionMatchedEnergy_TreeVec[(cluster->Plane().Plane)].push_back(projection_match);

    if(totalhits != 0){
      completeness_hits = ((signalhits)/totalhits);
    }//Close if statement for totalhits

    if(clusterhits.size() != 0){
      purity_hits = (signalhits/clusterhits.size());
    }//close clusterhits statement

    if(TotalTrueEnergy != 0){
      completeness_energy = ((TrktrackIdInfo.second)/TotalTrueEnergy);
    }

    if(TotalEnergyDepinHits != 0){
      purity_energy = (TrktrackIdInfo.second/TotalEnergyDepinHits);
    }
        
    cHitsCompl_TreeVec[cluster->Plane().Plane].push_back(completeness_hits);
    cHitsPurity_TreeVec[cluster->Plane().Plane].push_back(purity_hits);
    cEnergyCompl_TreeVec[cluster->Plane().Plane].push_back(completeness_energy);
    cEnergyPurity_TreeVec[cluster->Plane().Plane].push_back(purity_energy);

    if(fVerbose > 1){
      int clusterPlaneID = cluster->Plane().Plane;
      std::cout<<"#################################################"    <<std::endl;
      std::cout<<"               Cluster Metrics                   "    <<std::endl;
      std::cout<<"#################################################"    <<std::endl;
      std::cout<<"Cluster Plane ID :           " << clusterPlaneID      <<std::endl;
      std::cout<<"Projection matched:          " << projection_match    <<std::endl;
      std::cout<<"Cluster hit completeness:    " << completeness_hits   <<std::endl;
      std::cout<<"Cluster hit purity:          " << purity_hits         <<std::endl;
      std::cout<<"Cluster energy completeness: " << completeness_energy <<std::endl;
      std::cout<<"Cluster energy purity:       " << purity_energy       <<std::endl;
      std::cout<<"#################################################"    <<std::endl;

    }//Close fVerbose>1
   }//close the loop over clusters

   for(unsigned int plane=0; plane<numPlanes; plane++){
     //cProjectionMatchedEnergy_TreeVal[fModuleLabel][plane].push_back(cProjectionMatchedEnergy_TreeVec.at(plane));
     cEnergyCompl_TreeVal[fTrackModuleLabel][plane].push_back(cEnergyCompl_TreeVec.at(plane));
     cEnergyPurity_TreeVal[fTrackModuleLabel][plane].push_back(cEnergyPurity_TreeVec.at(plane));
     cHitsCompl_TreeVal[fTrackModuleLabel][plane].push_back(cHitsCompl_TreeVec.at(plane));
     cHitsPurity_TreeVal[fTrackModuleLabel][plane].push_back(cHitsPurity_TreeVec.at(plane));
   }//Close the loop over the planes
   return;
}//Close ClusterValidation function

//PFParticle Validation Function
//PFPValidation(pfpClusters,daughter,evt,clockData,clusterHandle,TrkMotherTrackIdMap,MCTrackEnergyMap,
//MCTrack_hit_map, fTrackModuleLabel);
void ana::TrackValidation::PFPValidation(std::vector<art::Ptr<recob::Cluster> >& clusters,
                           art::Ptr<recob::PFParticle> pfp,
                           const art::Event& evt,
                           std::map<int, const simb::MCParticle*> TrueParticlesMap,
                           const detinfo::DetectorClocksData& clockData,
                           art::Handle<std::vector<recob::Cluster> >& clusterHandle,
                           std::map<int,std::vector<int> >& TrackMotherTrackIDs,
                           std::map<int,float>& MCTrack_Energy_map,
                           std::map<art::ProductID,std::map<int,std::map<geo::PlaneID, int> > >& MCTrack_hit_map,
                           const std::string& fTrackModuleLabel){


  //Get the associated hits
  art::FindManyP<recob::Hit> fmhc(clusterHandle, evt, clusterHandle.provenance()->moduleLabel());
  
  //Holder for cluster hits
  //Get the Hits Handle used for this cluster type
  art::Handle<std::vector<recob::Hit > > hitHandle;
  evt.get(fmhc.at(clusters.at(0).key()).front().id(),hitHandle);

  if(!hitHandle.isValid()){
    mf::LogError("TrackValidation") << "Hits handle is stale. No pfp validation done" << std::endl;
    return;
  }//close the statement that check if hitHandle is valid
  
  float TotalTrueEnergy      = 0;
  float TotalEnergyDepinHits = 0;
  float signalhits           = 0;
  float totalhits            = 0;
  float pfphits              = 0;
  float clusterenergy        = 0;
  float pfpenergy            = 0;

  //Calculate the purity and completeness metrics  
  float completeness_hits   = -99999;
  float purity_hits         = -99999;
  float completeness_energy = -99999;
  float purity_energy       = -99999;

  int projectionMatched = 1;
  std::vector<int> projected_IDs;

  //Get the hits vector from the track
  for(auto const& cluster : clusters){
    std::vector< art::Ptr<recob::Hit> > clusterhits = fmhc.at(cluster.key());

    std::map<int, std::map<geo::PlaneID, int> > hitPlaneMap = NumberofPlaneHitsPerTrack(clockData, clusterhits);

    TotalEnergyDepinHits = TotalEnergyDepinHitsFun(clockData, clusterhits,(cluster->Plane().Plane));
    pfphits += clusterhits.size();

    std::pair<int,double> TrackShowerInfo;
    std::vector<int>  daughters;

    // if track-like use utils from here, else if shower-like use showerutils
    int pdg = abs(pfp->PdgCode()); // Track or shower
    if(pdg==13){
      //Function, finds the most probable track ID associated with the set of hits from there true energy 
      //depositons. The pair returns the energy in the hits as well.
      TrackShowerInfo = pairTrackTrackIdToEnergyDep(clockData,TrackMotherTrackIDs, clusterhits, 
                                                   (cluster->Plane().Plane));

      //Make sure the cluster has been matched.
      if(TrackShowerInfo.second == -99999){
        if(fVerbose>0){
          std::cout << "Reco cluster not matched to a true track" << std::endl;
        }
        continue;
      }
      
      daughters = TrackMotherTrackIDs[TrackShowerInfo.first];

    }else{
      std::map<int,std::vector<int> > ShowerMotherTrackIDs;
      ShowerMotherTrackIdsFun(evt, TrueParticlesMap, ShowerMotherTrackIDs);
      
      TrackShowerInfo = pairTrackTrackIdToEnergyDep(clockData, ShowerMotherTrackIDs, clusterhits,
                                                   (cluster->Plane().Plane));

      if(TrackShowerInfo.first == -99999){
        if(fVerbose>0){
          std::cout << "Reco cluster not matched to a True shower" << std::endl;
        }
        continue;
      }
      
      daughters = ShowerMotherTrackIDs[TrackShowerInfo.first];

    }//close else statement

    for(auto const& daughter: daughters){
      TotalTrueEnergy += MCTrack_Energy_map[daughter];

      auto hitPlaneMapIter = hitPlaneMap.find(daughter);
      if(hitPlaneMapIter==hitPlaneMap.end()){
        if(fVerbose>2){
          std::cout<<"daughter: "<<daughter<<" not found in cluster hit map"<<std::endl;
        }
        continue;
      }//Close the if statement
      for(const auto& hitPlaneIt: hitPlaneMapIter->second){
        //Count how many hits are from the true track or shower.
        signalhits += hitPlaneIt.second;
        //Count how many hits are missed in the plane from the true track id.
        totalhits += MCTrack_hit_map[hitHandle.id()][daughter][hitPlaneIt.first];
      }//close loop over the planes
      projected_IDs.push_back(TrackShowerInfo.second); //TODO check with Ed what he means by it.
    }//close the loop over daughters
    clusterenergy = TrackShowerInfo.second;
  }//Close cluster loop

  //Get the total reco energy of pfp particle
  pfpenergy += clusterenergy;

  //Have we matched the 2D cluster to the correct track or shower correctly. In terms of Energy depositions:
  //
  for(auto ID : projected_IDs){
    if(ID != projected_IDs.at(0)){projectionMatched=0;}
  }

  PFProjectionMatched_TreeVal[fTrackModuleLabel].push_back(projectionMatched);
  if(pfp->PdgCode() == 11){PFPShowerProjectionMatched_TreeVal[fTrackModuleLabel].push_back(projectionMatched);
  }else if(pfp->PdgCode() == 13){PFPTrackProjectionMatched_TreeVal[fTrackModuleLabel].push_back(projectionMatched);}


  if(totalhits != 0){
    completeness_hits = (signalhits)/totalhits;
  }//close if statement

  if(pfphits != 0){
    purity_hits = signalhits/pfphits;
  }

  if(TotalTrueEnergy != 0){
    completeness_energy = pfpenergy/TotalTrueEnergy;
  }

  if(TotalEnergyDepinHits != 0){
    purity_energy = pfpenergy/TotalEnergyDepinHits;
  }

  PFPHitsCompl_TreeVal[fTrackModuleLabel].push_back(completeness_hits);
  PFPHitsPurity_TreeVal[fTrackModuleLabel].push_back(purity_hits);
  PFPEnergyCompl_TreeVal[fTrackModuleLabel].push_back(completeness_energy);
  PFPEnergyPurity_TreeVal[fTrackModuleLabel].push_back(purity_energy);
  //For Tracks
  if(pfp->PdgCode() == 13){
    PFPTrackHitsCompl_TreeVal[fTrackModuleLabel].push_back(completeness_hits);
    PFPTrackHitsPurity_TreeVal[fTrackModuleLabel].push_back(purity_hits);
    PFPTrackEnergyCompl_TreeVal[fTrackModuleLabel].push_back(completeness_energy);
    PFPTrackEnergyPurity_TreeVal[fTrackModuleLabel].push_back(purity_energy);
  }else if(pfp->PdgCode() == 11){
    PFPShowerHitsCompl_TreeVal[fTrackModuleLabel].push_back(completeness_hits);
    PFPShowerHitsPurity_TreeVal[fTrackModuleLabel].push_back(purity_hits);
    PFPShowerEnergyCompl_TreeVal[fTrackModuleLabel].push_back(completeness_energy);
    PFPShowerEnergyPurity_TreeVal[fTrackModuleLabel].push_back(purity_energy);
  }

  if(fVerbose>1){
    std::cout << "#################################################" << std::endl;
    std::cout << "                 PFP Metrics                     " << std::endl;
    std::cout << "#################################################" << std::endl;
    std::cout << "Projection matched:      " << projectionMatched    << std::endl;
    std::cout << "PFP hit completeness:    " << completeness_hits    << std::endl;
    std::cout << "PFP hit purity:          " << purity_hits          << std::endl;
    std::cout << "PFP energy completeness: " << completeness_energy  << std::endl;
    std::cout << "PFP energy purity:       " << purity_energy        << std::endl;
    std::cout << "#################################################" << std::endl;

  }//close fverbose
  return;
}//close phparticel validation


//
//
//#################################################################################################
//
//=========================================================
//======= These function can be moved to TrackUtils =======
//=========================================================

//Function retunr a map of ShowerMother trackid as a key and a vector of daughters ids
void ana::TrackValidation::ShowerMotherTrackIdsFun(art::Event const& evt,
                                                   std::map<int, const simb::MCParticle*> TrueParticlesMap,
                                                   std::map<int, std::vector<int> >& ShwrMothTrkIdWithDauIdsMap){
  //###### Declaration and Initiate variables #######
  for(auto const& mapiter : TrueParticlesMap){
    //Declaration and Initiate variables
    const simb::MCParticle *particle = mapiter.second;
    // looking at the PdgCode of the particles
    unsigned int particlePdg = std::abs(particle->PdgCode());
    // Make sure that we are looking at shower particles and daughter
    if(particlePdg != 11 && particlePdg != 22){continue;}
    //Check the mother is not a shower particle
    if(TrueParticlesMap.find(particle->Mother()) != TrueParticlesMap.end()){
      if(std::abs(TrueParticlesMap[particle->Mother()]->PdgCode()) == 11
          || std::abs(TrueParticlesMap[particle->Mother()]->PdgCode()) == 22){
        //order the daughter based on it's mother.
        int MotherId = particle->Mother();
        int part_temp = MotherId;
        //Match particle with mother
        while(MotherId != 0){
          part_temp = TrueParticlesMap[MotherId]->Mother();
          if(TrueParticlesMap.find(part_temp) == TrueParticlesMap.end()){break;}
          if(std::abs(TrueParticlesMap[part_temp]->PdgCode()) != 11 && std::abs(TrueParticlesMap[part_temp]->PdgCode()) != 22){break;}
          MotherId = part_temp;
        }
        //Add to the mother chain.
        ShwrMothTrkIdWithDauIdsMap[MotherId].push_back(particle->TrackId());
        continue;
      }
    }//close the first condition block
    //
    ShwrMothTrkIdWithDauIdsMap[particle->TrackId()].push_back(particle->TrackId());
  }//close the loop over true mc particles
  return;
}//close ShowerMotherTrackIdsFun()

//Function retunr a map of Track's Mother trackid as a key and a vector of daughters ids
void ana::TrackValidation::MotherTrkIdsTodauIdsFun(art::Event const& evt,
                                                   std::map<int, const simb::MCParticle*> TrueParticlesMap,
                                                   std::map<int, std::vector<int> >& MothTrkIdWithDauIdsMap){
  //###### Declaration and Initiate variables #######
  for(auto const& mapiter : TrueParticlesMap){
    //Declaration and Initiate variables
    //
    std::vector<unsigned int> ParTrkPdgCode = fTrackPdgCodeVect;//{13,211,2212,321} PdgCode of track like particles
    //
    const simb::MCParticle *particle = mapiter.second;
    // looking at the PdgCode of the particles
    unsigned int particlePdg = std::abs(particle->PdgCode());
    // Make sure that we are looking at track particles and daughter
    if(std::find(ParTrkPdgCode.begin(), ParTrkPdgCode.end(), particlePdg) == ParTrkPdgCode.end()){continue;}
    //Check the mother is not a track particle
    if(TrueParticlesMap.find(particle->Mother()) != TrueParticlesMap.end()){
      if(std::find(ParTrkPdgCode.begin(), ParTrkPdgCode.end(), (std::abs(TrueParticlesMap[particle->Mother()]->PdgCode()))) != ParTrkPdgCode.end()){
        //order the daughter based on it's mother.
        int MotherId = particle->Mother();
        int part_temp = MotherId;
        //Match particle with mother
        while(MotherId != 0){
          part_temp = TrueParticlesMap[MotherId]->Mother();
          if(TrueParticlesMap.find(part_temp) == TrueParticlesMap.end()){break;}
          if(std::find(ParTrkPdgCode.begin(), ParTrkPdgCode.end(), (std::abs(TrueParticlesMap[part_temp]->PdgCode()))) == ParTrkPdgCode.end()){break;}
          MotherId = part_temp;
        }
        //Add to the mother chain.
        MothTrkIdWithDauIdsMap[MotherId].push_back(particle->TrackId());
        continue;
      }
    }//close the first condition block
    //
    MothTrkIdWithDauIdsMap[particle->TrackId()].push_back(particle->TrackId());
  }//close the loop over true mc particles
  return;
}//close MotherTrackIdsTodauIdsFun()

// Function retrun the true track id
int ana::TrackValidation::FunTrueParticleID(detinfo::DetectorClocksData const& clockData,
                                            const art::Ptr<recob::Hit>& hit) const{
  //###### Declaration and Initiate variables #######
  std::map<int,double> id_to_energy_map;
  double particleEnergy = -99999;
  int likelyTrackID = 0;
  std::vector<sim::TrackIDE> trackIDs = bt_serv->HitToTrackIDEs(clockData, hit);
  for(unsigned int idIt = 0; idIt < trackIDs.size(); ++idIt){
    //if(trackIDs.at(idIt).energy > particleEnergy){
    double Energy = trackIDs.at(idIt).energy;  
    int TrkID = std::abs(trackIDs.at(idIt).trackID);
    id_to_energy_map[TrkID] += Energy;
  }
  //
  for(auto const& mapIter : id_to_energy_map){
    double energy_map = mapIter.second;
    if(energy_map > particleEnergy){
      particleEnergy = energy_map;
      likelyTrackID = mapIter.first;
    }
  }//close the loop over the map
  return likelyTrackID;
}//Close the function retrun the true track id

//Function return all true hits of a particle
void ana::TrackValidation::TruthTrkIdsToHits(art::Event const& evt,
                                             detinfo::DetectorClocksData const& clockData,
                                             int MotherTrkIds,
                                             std::vector<int>& DaugthersIdsVec,
                                             std::vector<art::Ptr<recob::Hit>> allrecohits,
                                             std::vector<art::Ptr<recob::Hit>>& hitsOfTrkIds){
  //###### Declaration and Initiate variables #######
  for(auto const& hit : allrecohits){
    int HitsTrkIds = FunTrueParticleID(clockData, hit);
    if(HitsTrkIds == MotherTrkIds){
      hitsOfTrkIds.push_back(hit);
    }//colse condition statement
    else if(std::find(DaugthersIdsVec.begin(), DaugthersIdsVec.end(), HitsTrkIds) != DaugthersIdsVec.end()){
      hitsOfTrkIds.push_back(hit);
    }//Close else statement
  }//Close looping over hits
  return;
}//Close the function that return all ture hits of particle

// Function retrun the true track id of hits that deposit the most energy in a shower or track
int ana::TrackValidation::TureParticleIdDepMostEFun(detinfo::DetectorClocksData const& clockData,
                                                    std::vector<art::Ptr<recob::Hit> > ParticleHits) const{
  //###### Declaration and Initiate variables #######
  float max_energy = -999.9;
  int likelyParticleID = 0;
  std::map<int, float> TrkidToEnergyMap;
  for(auto const& hit : ParticleHits){
    double particleEnergy = 0;
    int likelyTrackID = 0;
    std::vector<sim::TrackIDE> trackIDs = bt_serv->HitToTrackIDEs(clockData, hit);
    for(unsigned int idIt = 0; idIt < trackIDs.size(); ++idIt){
      if(trackIDs.at(idIt).energy > particleEnergy){
        particleEnergy = trackIDs.at(idIt).energy;
        likelyTrackID = trackIDs.at(idIt).trackID;
      }
    }
    TrkidToEnergyMap[likelyTrackID] += particleEnergy;
    //return likelyTrackID;
  }// Close the loop over the hits
  for(auto const& [trkId , energyMap] : TrkidToEnergyMap){
    if(energyMap > max_energy){
      max_energy = energyMap;
      likelyParticleID = trkId;
    }
  }//Close the loop over TrkidToEnergy map
  if(likelyParticleID == 0){
    std::cout<< "The Particle trkid is zero check before continue "<<" Energt is:"<<max_energy<<std::endl;
  }
  return likelyParticleID;
}//Close funtion that return the true trackids that participate the most in showers or tracks  hits

//Function return a vector of true hits of a reco particle
std::vector<art::Ptr<recob::Hit> > ana::TrackValidation::TrueHitsOfParticle(detinfo::DetectorClocksData const& clockData,
                                                                            int MotherTrkIds,
                                                                            std::vector<int>& DaugthersIdsVec,
                                                                            std::vector<art::Ptr<recob::Hit> > RecoParticleHits){
  //###### Declaration and Initiate variables #######
  std::vector<art::Ptr<recob::Hit> >  hitsOfTrkIds;
  for(auto const& hit : RecoParticleHits){
    int TrkIdOfHit = FunTrueParticleID(clockData, hit);
    if(TrkIdOfHit == MotherTrkIds){
      hitsOfTrkIds.push_back(hit);
    }//colse condition statement
    else if(std::find(DaugthersIdsVec.begin(), DaugthersIdsVec.end(), TrkIdOfHit) != DaugthersIdsVec.end()){
      hitsOfTrkIds.push_back(hit);
    }//Close else statement
  }//close looping over all reco hits
  return hitsOfTrkIds;
}//Close the function that return a vector of ture hits of a reco particle

//Gets the total energy deposited by a track  by looping through the MC info on each wire.
std::map<int,float> ana::TrackValidation::TrueEnergyDepositedFromMCTracks(const std::vector<art::Ptr<sim::SimChannel> > &simchannels){
  //###### Declaration and Initiate variables #######
  art::ServiceHandle<geo::Geometry> geom;
  std::map<int,float> total_energies;
  //
  // Loop over SimChannel
  for(auto const& simchannel : simchannels){
    //
    //Get the plane.
    raw::ChannelID_t Channel = simchannel->Channel();
    std::vector<geo::WireID> Wire = geom->ChannelToWire(Channel);
    //sbnd is one channel per wire so the vector should be size one.
    int  PlaneID = (Wire[0]).Plane;

    //Get the TDCIDEMap (Charge vs Time)
    auto tdc_ide_map = simchannel->TDCIDEMap();

    //Loop through the map
    for(auto const& tdc_ide_pair : tdc_ide_map){
      //Get the IDEs associated to the TDC?
      auto const& ide_v = tdc_ide_pair.second;

      //Loop over the IDEs and add the energy. Only count from the collection plane.
      for(auto const& ide : ide_v){
        if(PlaneID == 2){
          total_energies[std::abs(ide.trackID)] +=  ide.energy;
        }
      }
    }
  }//Close the loop over the simchannels
  return total_energies;
}//Close the function that return the total energy deposited by a track.

//Function remove non contained particles
void ana::TrackValidation::RemoveNoneContainedParticles(std::map<int,std::vector<int> >&  TracksMothers,
                                                        std::map<int,const simb::MCParticle*>& trueParticlesMap, 
                                                        std::map<int,float>& MCTrackEnergyMap){
  //###### Declaration and Initiate variables #######
  for(auto const& trkmother : TracksMothers){
    
    //Loop over the daughter
    float deposited_energy = 0.f;
    float sim_energy       = (trueParticlesMap[trkmother.first]->E())*1000;
    for(auto const& daughter : trkmother.second){
      deposited_energy += MCTrackEnergyMap[daughter];
    }
    //If the energy of track is over 90% seen on the wires in truth, It is contained.
    if((deposited_energy/sim_energy) < 0.9){
      TracksMothers.erase(trkmother.first);
    }
  }
  return;
}//Close the function that remove non contained particles

//Function cut on the Track mothers with energy
void ana::TrackValidation::CutTrkMothersByE(std::map<int,std::vector<int> >&  TracksMothers,
                                            std::map<int,const simb::MCParticle*>& trueParticlesMap,
                                            float& EnergyCut){
  //Cut the true tracks and make sure they are a tracks.
  for(auto const& mapiter : TracksMothers){
    const simb::MCParticle* TrkMotherParticle = trueParticlesMap[mapiter.first];
    if(TrkMotherParticle->E() < EnergyCut){
      TracksMothers.erase(mapiter.first);
    }
  }//close the loop over mothers
  return;
}//close the function

//Function for cut on track via length 
void ana::TrackValidation::CalculateTrackLengthCut(std::map<int,std::vector<int> >&  TracksMothers,
                                                   std::map<int,const simb::MCParticle*>& trueParticlesMap,
                                                   float& LengthCut){
  //###### Declaration and Initiate variables #######
  for(auto const& mapiter : TracksMothers){
    float TrackLength = 0.f;
    bool startpc = false;
    bool endtpc  = false;
    const simb::MCParticle* TrkMotherParticle = trueParticlesMap[mapiter.first];
    //Get the number of Traj points to loop over
    unsigned int TrajPoints = TrkMotherParticle->NumberTrajectoryPoints();
    // Select first traj point where the track loses energy, last be default
    TLorentzVector PositionTrajStart =  TrkMotherParticle->Position(0);
    TLorentzVector PositionTrajEnd   =  TrkMotherParticle->Position(TrajPoints-1);
    //Check if the start and end position of track is inside the tpc
    startpc = IsInsideTPC((PositionTrajStart.Vect()), 0);
    endtpc  = IsInsideTPC((PositionTrajEnd.Vect()), 0);
    //Get the track length.
    TrackLength = (PositionTrajStart-PositionTrajEnd).Vect().Mag();
    if((!startpc) || (!endtpc) || (TrackLength < LengthCut)){
      TracksMothers.erase(mapiter.first);
    }
  }//close the loop over tracks mothers
  return;
}//close the function 

//Function return a bool; true if the track inside the tpc and false if not.
bool ana::TrackValidation::IsInsideTPC(TVector3 position, double distance_buffer){
  //###### Declaration and Initiate variables #######
  double vtx[3] = {position.X(), position.Y(), position.Z()};
  bool inside = false;
  geo::TPCID idtpc = geom->FindTPCAtPosition(vtx);
  if(geom->HasTPC(idtpc)){
    const geo::TPCGeo& tpcgeo = geom->GetElement(idtpc);
    double minx = tpcgeo.MinX(); double maxx = tpcgeo.MaxX();
    double miny = tpcgeo.MinY(); double maxy = tpcgeo.MaxY();
    double minz = tpcgeo.MinZ(); double maxz = tpcgeo.MaxZ();

    for(size_t a = 0; a < geom->Ncryostats(); a++){
      const geo::CryostatGeo& cryostat = geom->Cryostat(a);
      for(size_t b = 0; b < cryostat.NTPC(); b++){
        const geo::TPCGeo& tpcg = cryostat.TPC(b);
        if (tpcg.MinX() < minx) minx = tpcg.MinX();
        if (tpcg.MaxX() > maxx) maxx = tpcg.MaxX();
        if (tpcg.MinY() < miny) miny = tpcg.MinY();
        if (tpcg.MaxY() > maxy) maxy = tpcg.MaxY();
        if (tpcg.MinZ() < minz) minz = tpcg.MinZ();
        if (tpcg.MaxZ() > maxz) maxz = tpcg.MaxZ();
      }//close the second loop (b)
    }//close the first loop (a)

    //x
    double dista = fabs(minx - position.X());
    double distb = fabs(position.X() - maxx);
    if((position.X() > minx) && (position.X() < maxx) && 
        (dista > distance_buffer) && (distb > distance_buffer)){inside = true;}
    //y
    dista = fabs(maxy - position.Y());
    distb = fabs(position.Y() - miny);
    if (inside && (position.Y() > miny) && (position.Y() < maxy) &&
        (dista > distance_buffer) && (distb > distance_buffer)){inside = true;}
    else{inside = false;}
    //z
    dista = fabs(maxz - position.Z());
    distb = fabs(position.Z() - minz);
    if (inside && (position.Z() > minz) && (position.Z() < maxz) &&
        (dista > distance_buffer) && (distb > distance_buffer)){inside = true;}
    else{inside = false;}
  }//close if statement
  return inside;
}//Close the function

//Function return a map of number of hits per plane per track 
std::map<int,std::map<geo::PlaneID,int>> ana::TrackValidation::NumberofPlaneHitsPerTrack(
                                                               detinfo::DetectorClocksData const& clockData,
                                                               const std::vector<art::Ptr<recob::Hit> >& hits){
  //###### Declaration and Initiate variables #######
  std::map<int, std::map<geo::PlaneID, int> > HitNum;
   //Loop over the hits and find the IDE
   for(auto const& hit : hits){

     geo::WireID wireid = hit->WireID();
     geo::PlaneID  PlaneID = wireid.planeID();

     std::vector<sim::TrackIDE> trackIDEs = bt_serv->HitToTrackIDEs(clockData, hit);

     std::map<int,float> hitEnergies;

     //Loop over the IDEs associated to the hit and add up energies
     for(unsigned int idIt = 0; idIt < trackIDEs.size(); ++idIt){
       hitEnergies[std::abs(trackIDEs.at(idIt).trackID)] += trackIDEs.at(idIt).energy;
     }//close the loop over trackIDs

     //Find which track deposited the most energy.
     int   likelytrack = -9999;
     float MaxEnergy   = -9999;
     for(auto const& trkmapiter : hitEnergies){
       if(trkmapiter.second > MaxEnergy){
         MaxEnergy = trkmapiter.second;
         likelytrack = trkmapiter.first;
       }//Close if statement
     }//Close the loop over hitEnergies map

     ++HitNum[likelytrack][PlaneID];
   
   }//Close the loop over the hits
   return HitNum;
}//Close the function 

//Function for finding the track Best plane
int ana::TrackValidation::FindingTrackBeastPlane(detinfo::DetectorClocksData const& clockData,
                                                 const std::vector<art::Ptr<recob::Hit> >& trkhits){
  //###### Declaration and Initiate variables #######
  unsigned int BestPlaneId = 999;
  unsigned int bestPlaneNumHits = 0;
  float BestPlaneEnergy = -999;
  std::map<int,float> PlaneEnergies; //A map of planeId as a key and energy deposited in each plane as value
  std::map<geo::PlaneID::PlaneID_t, std::vector<art::Ptr<recob::Hit> > > planeHits;

  for(auto const& hit : trkhits){
    geo::WireID wireid = hit->WireID();
    unsigned int PlaneID = wireid.Plane;

    planeHits[PlaneID].push_back(hit);

    std::vector<sim::TrackIDE> trackIDEs = bt_serv->HitToTrackIDEs(clockData, hit);
    for(unsigned int idIt = 0; idIt < trackIDEs.size(); ++idIt){
      PlaneEnergies[PlaneID] += trackIDEs.at(idIt).energy;
    }//close the loop over trackIDEs
  }//cose the loop over the hits
  
  //Get the best plane for this track
  if(fUseBestPlane == 3){
    for(auto const& [plane, hits] : planeHits){
      unsigned int planeNumHits = hits.size();
      if(planeNumHits > bestPlaneNumHits){
        BestPlaneId      = plane;
        bestPlaneNumHits = planeNumHits;
      }
    }//close the loop over the planeHits map
  }else if(fUseBestPlane == 4){
    for(auto const& [planeIds, PEnergy] : PlaneEnergies){
      if(PEnergy > BestPlaneEnergy){
        BestPlaneId      = planeIds;
        BestPlaneEnergy  = PEnergy;
      }//close the if statement
    }//close the loop over the map
    //std::cout<< "BPEnergy = " <<BestPlaneId<<" U: "<<PlaneEnergies[0]<< " V: " <<PlaneEnergies[1]<<" Y: "<<PlaneEnergies[2]<<std::endl;
  }else{
    BestPlaneId = fUseBestPlane;
  }

  return BestPlaneId;
}//Close the function for finding the track best plane
//

//Function to finds the most probable track id associated with the set of hits from there true energy depositions.
//The pair returns the energy as well.
std::pair<int,double>ana::TrackValidation::pairTrackTrackIdToEnergyDep(detinfo::DetectorClocksData const& clockData,
                                                                    std::map<int,std::vector<int>> &TracksMothers,
                                                                    const std::vector<art::Ptr<recob::Hit> >& hits,
                                                                    int planeid){
  //###### Declaration and Initiate variables #######
  std::map<int, float> mapTrackIdToTrkHitEnergy; //Map of trackid and the energy deposited per the trackid for 3 planes
  std::map<int, float> mapTrackIdToEDep;//Map of trackid and the energy deposited per the trackid for best track planes
  //Looping over track hits
  for(auto const& hit : hits){
    geo::WireID wireid = hit->WireID();
    int trkPlaneId = wireid.Plane;
    //
    std::vector<sim::TrackIDE> TrkTrackIds = bt_serv->HitToTrackIDEs(clockData, hit);
    for(auto const& trkID : TrkTrackIds){
      mapTrackIdToTrkHitEnergy[std::abs(trkID.trackID)] += trkID.energy;
      if(trkPlaneId == planeid){
        mapTrackIdToEDep[std::abs(trkID.trackID)] += trkID.energy;
      }//Close if statement
    }//close the loop over trkids
  }//close the loop over hits

  //Find the total energy deposited by each track mother
  std::map<int, float> mapMotherIdToEnergyDep;//three planes
  std::map<int, float> mapMotherIdToEDep;    //Best track plane
  for(auto const& trkMother : TracksMothers){
    for(auto const& trkDaughter : trkMother.second){
      mapMotherIdToEnergyDep[trkMother.first] += mapTrackIdToTrkHitEnergy[trkDaughter];
      mapMotherIdToEDep[trkMother.first] += mapTrackIdToEDep[trkDaughter];
    }//close the loop over track daughter
  }//close the looping over track mothers

  {
    float maxEnergy = -99999;
    int TracktrkId  = -99999;
    for(auto const& mapmotherE : mapMotherIdToEnergyDep){
      float MEnergy = mapmotherE.second;
      int MotherTrkID = mapmotherE.first;
      if(MEnergy > maxEnergy){
        maxEnergy = MEnergy;
        TracktrkId = MotherTrkID;
      }//close if statement
    }//Close the loop over mapMotherIdToEDep
    if(maxEnergy == 0){
      return std::make_pair(-99999,-99999);
      std::cout <<"There is no energy in this track id \n";
    }
    return std::make_pair(TracktrkId, mapMotherIdToEDep[TracktrkId]);
  }//Close the block
}//close the function of TrueParticleIDFromTrueChain 

//Function return the track id of particle that has the most hits in a track.
int ana::TrackValidation::TrueParticleIDFromTotalRecoHits(detinfo::DetectorClocksData const& clockData,
                                                          const std::vector<art::Ptr<recob::Hit> >& hits){
  // Make a map of the tracks which are associated with this object and 
  // the number of hits they are the primary contributor to
  std::map<int,int> trackMap;
  for(auto const& hit : hits){
    int trackID = FunTrueParticleID(clockData, hit);
    trackMap[trackID]++;
  }//Close the loop over track hits
  
  // Pick the track which is the primary contributor to the most hits as the 'true track'
  int objectTrack = -99999;
  int highestCount = -1;
  int NHighestCounts = 0;
  for(auto const& mapIter : trackMap){
    if(mapIter.second > highestCount){
      highestCount = mapIter.second;
      objectTrack = mapIter.first;
      NHighestCounts = 1;
    }//close if statement
    else if(mapIter.second == highestCount){
      NHighestCounts++;
    }//
  }//
  if(NHighestCounts > 1){
    std::cout<<"TrueParticleIDFromTotalRecoHits - There are "<<NHighestCounts<<
      " particles which tie for highest number of contributing hits ("<< highestCount<<
      " hits).  Using TrueParticleIDFromTotalTrueEnergy instead."<<std::endl;
    objectTrack = TrueParticleIDFromTotalTrueEnergy(clockData, hits);
  }//Close the if statement
  return objectTrack;
}//close function return the track id of particle that has most hits in a track.

//Function return the track id of particle that deposite the most energy in the track
int ana::TrackValidation::TrueParticleIDFromTotalTrueEnergy(detinfo::DetectorClocksData const& clockData,
                                                            const std::vector<art::Ptr<recob::Hit> >& hits){
  //###### Declaration and Initiate variables #######
  std::map<int,double> trackIDToEDepMap;
  for(auto const& hit : hits){
    std::vector<sim::TrackIDE> trackIDs = bt_serv->HitToTrackIDEs(clockData, hit);
    for(unsigned int idIt = 0; idIt < trackIDs.size(); ++idIt){
      int unsigned_id = std::abs(trackIDs[idIt].trackID);
       trackIDToEDepMap[unsigned_id] += trackIDs[idIt].energy;
    }//close the loop over the trackids
  }//close the loop over the hits

  //Loop over the map and find the track which contributes the highest energy to the hit vector
  double max_energy = -1;
  int ParticleTrkId = -99999;
  for(auto const& mapIter : trackIDToEDepMap){
    double energy = mapIter.second;
    int    trkId  = mapIter.first;
    if(energy > max_energy){
      max_energy = energy;
      ParticleTrkId = trkId;
    }//close if statement
  }//close the loop over the map trackIDToEDepMap
  return ParticleTrkId;
}//close the function TrueParticleIDFromTotalTrueEnergy

//Function returns a float value of total energy deposited from a vector of hits for a specific plane
float ana::TrackValidation::TotalEnergyDepinHitsFun(detinfo::DetectorClocksData const& clockData,
                                                    const std::vector<art::Ptr<recob::Hit> >& hits,
                                                    int Plane){
  //###### Declaration and Initiate variables #######
  float DepEnergy = 0;
  //Loop over the hits and find the IDEs
  for(auto const& hit : hits){
    
    //Get the plane ID
    geo::WireID wireid = hit->WireID();
    int PlaneID = wireid.Plane;
    if(PlaneID != Plane){continue;}
    
    //Get the trackIds from the hit.
    std::vector<sim::TrackIDE> trackIDs = bt_serv->HitToTrackIDEs(clockData, hit);
    for(unsigned int idIt = 0; idIt < trackIDs.size(); ++idIt){
      //Find the true total energy deposited in a set of hits.
      DepEnergy += trackIDs.at(idIt).energy;
    }//Close the loop over trackIDs
  }//close the loop over hits

  return DepEnergy;

}//Close the function
//

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
