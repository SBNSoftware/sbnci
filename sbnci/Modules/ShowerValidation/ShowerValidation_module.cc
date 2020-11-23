////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Class:       ShowerValidation                                                                          //
// Modlue Type: Analyser                                                                                  //
// File         ShowerValidation_module.cc                                                                //
// Author:      Dominic Barker dominic.barker@sheffield.ac.uk                                             //
//              Edward Tyley e.tyley@sheffield.ac.uk                                                      //
//                                                                                                        //
// Usage:       The Validation module takes an array of shower modules and perfoms truth matching on      //
//              reconstructed EM showers. It calculates the truth information from geant information      //
//              such as the diretion and starting position and energy deposited. Metrics are then defined //
//              to compare the efficiency of the shower reconstruction modules. These are presented       //
//              as histograms and in terms of energy. Energy plots look at the metrics for specific       //
//              energies between +- fEnergyWidth of the values described in the fcl table                 //
//                                                                                                        //
// Updates:     31.10.2018  Clustering Validation Added                                                   //
//              13.11.2018  Hit Validation, 2D Histograms and refactoring of the code                     //
//              22.01.2019  TTree output added                                                            //
//              18.03.2019  Added PFParticle Validation                                                   //
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
#include "lardataobj/Simulation/SimChannel.h"
#include "lardataobj/Simulation/AuxDetSimChannel.h"
#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "MCRecoUtils/RecoUtils.h"
#include "MCRecoUtils/ShowerUtils.h"

//Root Includes
#include "TMath.h"
#include "TTree.h"
#include "TSystem.h"

//C++ Includes
#include <vector>
#include <iostream>
#include <fstream>

namespace ana {
  class ShowerValidation;
}

class ana::ShowerValidation : public art::EDAnalyzer {
  public:

    ShowerValidation(const fhicl::ParameterSet& pset);

    void analyze(const art::Event& evt);
    void endJob();
    void beginJob();

    void initTree(TTree* Tree, std::string branchName, float& Metric);
    void initTree(TTree* Tree, std::string branchName, std::vector<float>& Metric);
    void initTree(TTree* Tree, std::string branchName, std::map<std::string,std::vector<float> >& Metric, std::vector<std::string> fShowerModuleLabels);
    void initClusterTree(TTree* Tree,  std::string branchName,  std::map<std::string,std::vector<std::vector<std::vector<float> > > >& Metric, std::vector<std::string> fShowerModuleLabels);

    void ClusterValidation(std::vector<art::Ptr<recob::Cluster> >& clusters,
        const art::Event& evt,
        const detinfo::DetectorClocksData& clockData,
        art::Handle<std::vector<recob::Cluster> >& clusterHandle,
        std::map<int,std::vector<int> >& ShowerMotherTrackIDs,
        std::map<int,float>& MCTrack_Energy_map,
        std::map<art::ProductID,std::map<int,std::map<geo::PlaneID, int> > > & MCTrack_hit_map,
        int& TrueShowerID,
        float& trueShowerEnergy,
        const std::string & fShowerModuleLabel);

    void PFPValidation(std::vector<art::Ptr<recob::Cluster> >& clusters,
        art::Ptr<recob::PFParticle>  pfp,
        const art::Event& evt,
        const detinfo::DetectorClocksData& clockData,
        art::Handle<std::vector<recob::Cluster> >& clusterHandle,
        std::map<int,std::vector<int> >& ShowerMotherTrackIDs,
        std::map<int,float>& MCTrack_Energy_map,
        std::map<art::ProductID,std::map<int,std::map<geo::PlaneID, int> > > & MCTrack_hit_map,
        const std::string & fShowerModuleLabel);


  private:

    //fcl parameters
    std::string fGenieGenModuleLabel;
    std::string fLArGeantModuleLabel;
    std::string fHitsModuleLabel;
    std::string fTrackModuleLabel;
    std::string fPFParticleLabel;

    bool  fUseBiggestShower;
    bool  fFillOnlyClosestShower;
    bool  fRemoveNonContainedParticles;
    bool  fPFPValidation;
    bool  fUseNeutrinoShowers;
    bool  fClusterValidation;
    bool  fHitValidation;
    bool  fMatchShowersInTruth;
    int   fVerbose;
    int   fMinHitSize;
    float fSimEnergyCut;
    float fDensityCut;
    float fMaxSimEnergy;
    float fMinRecoEnergy;
    float fMatchedCut;

    std::vector<std::string> fShowerModuleLabels;
    std::vector<std::string> fHitModuleLabels;

    std::map<std::string,std::vector<float> > sDirX_TreeVal;
    std::map<std::string,std::vector<float> > sDirY_TreeVal;
    std::map<std::string,std::vector<float> > sDirZ_TreeVal;
    std::map<std::string,std::vector<float> > sTrueDirX_TreeVal;
    std::map<std::string,std::vector<float> > sTrueDirY_TreeVal;
    std::map<std::string,std::vector<float> > sTrueDirZ_TreeVal;
    std::map<std::string,std::vector<float> > sDirDiff_TreeVal;
    std::map<std::string,std::vector<float> > sStartX_TreeVal;
    std::map<std::string,std::vector<float> > sStartY_TreeVal;
    std::map<std::string,std::vector<float> > sStartZ_TreeVal;
    std::map<std::string,std::vector<float> > sStartDist_TreeVal;
    std::map<std::string,std::vector<float> > sLength_TreeVal;
    std::map<std::string,std::vector<float> > sOpeningAngle_TreeVal;
    std::map<std::string,std::vector<float> > sdEdx_TreeVal;
    std::map<std::string,std::vector<float> > sAllHitsPurity_TreeVal;
    std::map<std::string,std::vector<float> > sAllHitsComp_TreeVal;
    std::map<std::string,std::vector<float> > sEnergyPurityBestPlane_TreeVal;
    std::map<std::string,std::vector<float> > sEnergyCompBestPlane_TreeVal;
    std::map<std::string,std::vector<float> > sEnergy_TreeVal;
    std::map<std::string,std::vector<float> > sNumHits_TreeVal;
    std::map<std::string,std::vector<float> > sEnergyRat_TreeVal;
    std::map<std::string,std::vector<float> > sEnergyDiff_TreeVal;
    std::map<std::string,std::vector<float> > sEnergyDiffTrue_TreeVal;
    std::map<std::string,std::vector<float> > sTrueEnergy_TreeVal;
    std::map<std::string,std::vector<float> > sBestPlane_TreeVal;
    std::map<std::string,std::vector<float> > sGeoProjectionMatched_TreeVal;
    std::map<std::string,std::vector<float> > sTrackEnergyComp_TreeVal;
    std::map<std::string,std::vector<float> > sTrackEnergyPurity_TreeVal;
    std::map<std::string,std::vector<float> > sTrackHitComp_TreeVal;
    std::map<std::string,std::vector<float> > sTrackHitPurity_TreeVal;
    std::map<std::string,std::vector<float> > sTrackHitCompPurity_TreeVal;
    std::map<std::string,std::vector<float> > sTrueTrackLengh_TreeVal;
    std::map<std::string,std::vector<float> > sTrueTrackSpread_TreeVal;
    std::map<std::string,std::vector<float> > sPDG_TreeVal;
    std::map<std::string,std::vector<float> > sMother_TreeVal;
    std::map<std::string,std::vector<float> > sTrackPDG_TreeVal;

    std::map<std::string,std::vector<float> > pfpNeutrinos_TreeVal;
    std::map<std::string,std::vector<float> > pfpAllTracks_TreeVal;
    std::map<std::string,std::vector<float> > pfpAllShowers_TreeVal;
    std::map<std::string,std::vector<float> > pfpPrimaryTracks_TreeVal;
    std::map<std::string,std::vector<float> > pfpPrimaryShowers_TreeVal;
    std::map<std::string,std::vector<float> > pfpVertexDistX_TreeVal;
    std::map<std::string,std::vector<float> > pfpVertexDistY_TreeVal;
    std::map<std::string,std::vector<float> > pfpVertexDistZ_TreeVal;
    std::map<std::string,std::vector<float> > pfpVertexDistMag_TreeVal;
    std::map<std::string,std::vector<float> > pfpProjectionMatched_TreeVal;
    // include number of hits in pfp
    std::map<std::string,std::vector<float> > pfpHitsComp_TreeVal;
    std::map<std::string,std::vector<float> > pfpEnergyComp_TreeVal;
    std::map<std::string,std::vector<float> > pfpHitsPurity_TreeVal;
    std::map<std::string,std::vector<float> > pfpEnergyPurity_TreeVal;
    std::map<std::string,std::vector<float> > pfpTrackProjectionMatched_TreeVal;
    std::map<std::string,std::vector<float> > pfpTrackHitsComp_TreeVal;
    std::map<std::string,std::vector<float> > pfpTrackEnergyComp_TreeVal;
    std::map<std::string,std::vector<float> > pfpTrackHitsPurity_TreeVal;
    std::map<std::string,std::vector<float> > pfpTrackEnergyPurity_TreeVal;
    std::map<std::string,std::vector<float> > pfpShowerProjectionMatched_TreeVal;
    std::map<std::string,std::vector<float> > pfpShowerHitsComp_TreeVal;
    std::map<std::string,std::vector<float> > pfpShowerEnergyComp_TreeVal;
    std::map<std::string,std::vector<float> > pfpShowerHitsPurity_TreeVal;
    std::map<std::string,std::vector<float> > pfpShowerEnergyPurity_TreeVal;

    std::map<std::string,std::vector<float> > eSegmentation_TreeVal;
    std::map<std::string,std::vector<float> > eNumRecoTracks_TreeVal;
    std::map<std::string,std::vector<float> > eNumRecoShowers_TreeVal;

    float eNumTrueShowers_TreeVal;
    float eNumTrueShowersviaECut_TreeVal;
    float eNumTrueShowersviaDCut_TreeVal;
    std::vector<float> eTrueHitNum_TreeVal;
    std::vector<float> eTrueShowerE_TreeVal;
    std::vector<float> eTrueShowerEviaECut_TreeVal;
    std::vector<float> eTrueShowerEviaDCut_TreeVal;

    // std::map<std::string,std::vector<std::vector<std::vector<float> > > > cProjectionMatchedEnergy_TreeVal;
    std::map<std::string,std::vector<std::vector<std::vector<float> > > > cEnergyComp_TreeVal;
    std::map<std::string,std::vector<std::vector<std::vector<float> > > > cEnergyPurity_TreeVal;
    std::map<std::string,std::vector<std::vector<std::vector<float> > > > cHitsComp_TreeVal;
    std::map<std::string,std::vector<std::vector<std::vector<float> > > > cHitsPurity_TreeVal;

    std::map<std::string,std::vector<std::vector<float> > > hEnergyComp_TreeVal;
    std::map<std::string,std::vector<std::string> > sStartEndProcess_TreeVal;


    float EventRun_TreeVal;
    float EventSubrun_TreeVal;
    float EventNumber_TreeVal;

    //TTree
    TTree* Tree;

    //Service handlesacktracker
    art::ServiceHandle<cheat::BackTrackerService> backtracker;
    art::ServiceHandle<cheat::ParticleInventoryService> particleInventory;
    art::ServiceHandle<geo::Geometry> geom;
    art::ServiceHandle<art::TFileService> tfs;

    int numevents;
    int containmentCutNum;
    int numshowers;
    int numshowerspassTPC;
    int numshowerspassdensity;
    int numshoowerspassenergy;
    int numrecoshowers;
    int numrecoshowersana;

};

ana::ShowerValidation::ShowerValidation(const fhicl::ParameterSet& pset) : EDAnalyzer(pset){

  fGenieGenModuleLabel         = pset.get<std::string>("GenieGenModuleLabel");
  fLArGeantModuleLabel         = pset.get<std::string>("LArGeantModuleLabel");
  fHitsModuleLabel             = pset.get<std::string>("HitsModuleLabel");
  fTrackModuleLabel            = pset.get<std::string>("TrackModuleLabel");
  fPFParticleLabel             = pset.get<std::string>("PFParticleLabel");
  fShowerModuleLabels          = pset.get<std::vector<std::string> >("ShowerModuleLabels");
  fHitModuleLabels             = pset.get<std::vector<std::string> >("HitModuleLabels");
  fUseBiggestShower            = pset.get<bool>("UseBiggestShower");
  fFillOnlyClosestShower       = pset.get<bool>("FillOnlyClosestShower");
  fRemoveNonContainedParticles = pset.get<bool>("RemoveNonContainedParticles");
  fPFPValidation               = pset.get<bool>("PFPValidation");
  fClusterValidation           = pset.get<bool>("ClusterValidation");
  fHitValidation               = pset.get<bool>("HitValidation");
  fUseNeutrinoShowers          = pset.get<bool>("UseNeutrinoShowers");
  fMatchShowersInTruth         = pset.get<bool>("MatchShowersInTruth");
  fVerbose                     = pset.get<int>("Verbose");
  fMinHitSize                  = pset.get<int>("MinHitSize");
  fSimEnergyCut                = pset.get<float>("SimEnergyCut");
  fDensityCut                  = pset.get<float>("DensityCut");
  fMaxSimEnergy                = pset.get<float>("MaxSimEnergy");
  fMinRecoEnergy               = pset.get<float>("MinRecoEnergy");
  fMatchedCut                  = pset.get<float>("MatchedCut");

  if (fFillOnlyClosestShower && fUseBiggestShower) {
    throw cet::exception("ShowerValidation") << "Both fFillOnlyClosestShower and fUseBiggestShower are set. Please choose one of them" << std::endl;
  }
}

void ana::ShowerValidation::initTree(TTree* Tree, std::string branchName, std::vector<float>& Metric){
  const char* branchChar   = branchName.c_str();
  Tree->Branch(branchChar,"std::vector<float>", &Metric, 32000, 0);
}

void ana::ShowerValidation::initTree(TTree* Tree, std::string branchName, float& Metric){
  const char* branchChar   = branchName.c_str();
  Tree->Branch(branchChar, &Metric, 32000, 0);
}

void ana::ShowerValidation::initTree(TTree* Tree, std::string branchName, std::map<std::string,std::vector<float> >& Metric,   std::vector<std::string> fShowerModuleLabels){
  for (auto const& fShowerModuleLabel: fShowerModuleLabels){
    std::string branchString = branchName + "_" + fShowerModuleLabel;
    const char* branchChar   = branchString.c_str();
    Tree->Branch(branchChar,"std::vector<float>", &Metric[fShowerModuleLabel], 32000, 0);
  }
}

void ana::ShowerValidation::initClusterTree(TTree* Tree, std::string branchName,  std::map<std::string,std::vector<std::vector<std::vector<float> > > >& Metric, std::vector<std::string> fShowerModuleLabels){
  for (auto const& fShowerModuleLabel: fShowerModuleLabels){
    unsigned int numPlanes = geom->Nplanes();
    Metric[fShowerModuleLabel].resize(numPlanes);
    for (unsigned int plane=0; plane<numPlanes; plane++) {
      std::string branchString = "Plane_" + std::to_string(plane) + "_" + branchName + "_" + fShowerModuleLabel;
      const char* branchChar   = branchString.c_str();
      Tree->Branch(branchChar,"std::vector<std::vector<float> >", &Metric[fShowerModuleLabel].at(plane), 32000, 0);
    }
  }
}


void ana::ShowerValidation::beginJob() {

  numevents = 0;
  containmentCutNum = 0;
  numshowers = 0;
  numshowerspassTPC = 0;
  numshowerspassdensity = 0;
  numshoowerspassenergy = 0;
  numrecoshowers = 0;
  numrecoshowersana = 0;

  Tree = tfs->make<TTree>("MetricTree", "Tree Holding all metric information");
  gInterpreter->GenerateDictionary("vector<vector<float> > ","vector");
  Tree->Branch("EventRun", &EventRun_TreeVal, 32000, 0);
  Tree->Branch("EventSubrun", &EventSubrun_TreeVal, 32000, 0);
  Tree->Branch("EventNumber", &EventNumber_TreeVal, 32000, 0);

  initTree(Tree,"sDirX",sDirX_TreeVal,fShowerModuleLabels);
  initTree(Tree,"sDirY",sDirY_TreeVal,fShowerModuleLabels);
  initTree(Tree,"sDirZ",sDirZ_TreeVal,fShowerModuleLabels);
  initTree(Tree,"sTrueDirX",sTrueDirX_TreeVal,fShowerModuleLabels);
  initTree(Tree,"sTrueDirY",sTrueDirY_TreeVal,fShowerModuleLabels);
  initTree(Tree,"sTrueDirZ",sTrueDirZ_TreeVal,fShowerModuleLabels);
  initTree(Tree,"sDirDiff",sDirDiff_TreeVal,fShowerModuleLabels);
  initTree(Tree,"sStartX",sStartX_TreeVal,fShowerModuleLabels);
  initTree(Tree,"sStartY",sStartY_TreeVal,fShowerModuleLabels);
  initTree(Tree,"sStartZ",sStartZ_TreeVal,fShowerModuleLabels);
  initTree(Tree,"sStartDist",sStartDist_TreeVal,fShowerModuleLabels);
  initTree(Tree,"sLength",sLength_TreeVal,fShowerModuleLabels);
  initTree(Tree,"sOpeningAngle",sOpeningAngle_TreeVal,fShowerModuleLabels);
  initTree(Tree,"sdEdx",sdEdx_TreeVal,fShowerModuleLabels);
  initTree(Tree,"sAllHitsPurity",sAllHitsPurity_TreeVal,fShowerModuleLabels);
  initTree(Tree,"sAllHitsComp",sAllHitsComp_TreeVal,fShowerModuleLabels);
  initTree(Tree,"sEnergyPurityBestPlane",sEnergyPurityBestPlane_TreeVal,fShowerModuleLabels);
  initTree(Tree,"sEnergyCompBestPlane",sEnergyCompBestPlane_TreeVal,fShowerModuleLabels);
  initTree(Tree,"sEnergy",sEnergy_TreeVal,fShowerModuleLabels);
  initTree(Tree,"sNumHits",sNumHits_TreeVal,fShowerModuleLabels);
  initTree(Tree,"sEnergyRat",sEnergyRat_TreeVal,fShowerModuleLabels);
  initTree(Tree,"sEnergyDiff",sEnergyDiff_TreeVal,fShowerModuleLabels);
  initTree(Tree,"sEnergyDiffTrue",sEnergyDiffTrue_TreeVal,fShowerModuleLabels);
  initTree(Tree,"sTrueEnergy",sTrueEnergy_TreeVal,fShowerModuleLabels);
  initTree(Tree,"sBestPlane",sBestPlane_TreeVal,fShowerModuleLabels);
  initTree(Tree,"sGeoProjectionMatched",sGeoProjectionMatched_TreeVal,fShowerModuleLabels);
  initTree(Tree,"sTrackEnergyComp",sTrackEnergyComp_TreeVal,fShowerModuleLabels);
  initTree(Tree,"sTrackEnergyPurity",sTrackEnergyPurity_TreeVal,fShowerModuleLabels);
  initTree(Tree,"sTrackHitComp",sTrackHitComp_TreeVal,fShowerModuleLabels);
  initTree(Tree,"sTrackHitPurity",sTrackHitPurity_TreeVal,fShowerModuleLabels);
  initTree(Tree,"sTrackHitCompPurity",sTrackHitCompPurity_TreeVal,fShowerModuleLabels);
  initTree(Tree,"sTrueTrackLengh",sTrueTrackLengh_TreeVal,fShowerModuleLabels);
  initTree(Tree,"sTrueTrackSpread",sTrueTrackSpread_TreeVal,fShowerModuleLabels);
  initTree(Tree,"sPDG",sPDG_TreeVal,fShowerModuleLabels);
  initTree(Tree,"sPDGTrack",sTrackPDG_TreeVal,fShowerModuleLabels);
  initTree(Tree,"sMother",sMother_TreeVal,fShowerModuleLabels);

  initTree(Tree,"pfpNeutrinos",pfpNeutrinos_TreeVal,fShowerModuleLabels);
  initTree(Tree,"pfpAllTracks",pfpAllTracks_TreeVal,fShowerModuleLabels);
  initTree(Tree,"pfpAllShowers",pfpAllShowers_TreeVal,fShowerModuleLabels);
  initTree(Tree,"pfpTracks",pfpPrimaryTracks_TreeVal,fShowerModuleLabels);
  initTree(Tree,"pfpShowers",pfpPrimaryShowers_TreeVal,fShowerModuleLabels);
  initTree(Tree,"pfpVertexDistX",pfpVertexDistX_TreeVal,fShowerModuleLabels);
  initTree(Tree,"pfpVertexDistY",pfpVertexDistY_TreeVal,fShowerModuleLabels);
  initTree(Tree,"pfpVertexDistZ",pfpVertexDistZ_TreeVal,fShowerModuleLabels);
  initTree(Tree,"pfpVertexDistMag",pfpVertexDistMag_TreeVal,fShowerModuleLabels);
  if (fPFPValidation){
    initTree(Tree,"pfpProjectionMatched",pfpProjectionMatched_TreeVal,fShowerModuleLabels);
    initTree(Tree,"pfpHitsComp",pfpHitsComp_TreeVal,fShowerModuleLabels);
    initTree(Tree,"pfpEnergyComp",pfpEnergyComp_TreeVal,fShowerModuleLabels);
    initTree(Tree,"pfpHitsPurity",pfpHitsPurity_TreeVal,fShowerModuleLabels);
    initTree(Tree,"pfpEnergyPurity",pfpEnergyPurity_TreeVal,fShowerModuleLabels);
    initTree(Tree,"pfpTrackProjectionMatched",pfpTrackProjectionMatched_TreeVal,fShowerModuleLabels);
    initTree(Tree,"pfpTrackHitsComp",pfpTrackHitsComp_TreeVal,fShowerModuleLabels);
    initTree(Tree,"pfpTrackEnergyComp",pfpTrackEnergyComp_TreeVal,fShowerModuleLabels);
    initTree(Tree,"pfpTrackHitsPurity",pfpTrackHitsPurity_TreeVal,fShowerModuleLabels);
    initTree(Tree,"pfpTrackEnergyPurity",pfpTrackEnergyPurity_TreeVal,fShowerModuleLabels);
    initTree(Tree,"pfpShowerProjectionMatched",pfpShowerProjectionMatched_TreeVal,fShowerModuleLabels);
    initTree(Tree,"pfpShowerHitsComp",pfpShowerHitsComp_TreeVal,fShowerModuleLabels);
    initTree(Tree,"pfpShowerEnergyComp",pfpShowerEnergyComp_TreeVal,fShowerModuleLabels);
    initTree(Tree,"pfpShowerHitsPurity",pfpShowerHitsPurity_TreeVal,fShowerModuleLabels);
    initTree(Tree,"pfpShowerEnergyPurity",pfpShowerEnergyPurity_TreeVal,fShowerModuleLabels);
  }
  initTree(Tree,"eSegmentation",eSegmentation_TreeVal,fShowerModuleLabels);
  initTree(Tree,"eNumRecoTracks",eNumRecoTracks_TreeVal,fShowerModuleLabels);
  initTree(Tree,"eNumRecoShowers",eNumRecoShowers_TreeVal,fShowerModuleLabels);

  initTree(Tree,"eNumTrueShowers",eNumTrueShowers_TreeVal);
  initTree(Tree,"eNumTrueShowersviaECut",eNumTrueShowersviaECut_TreeVal);
  initTree(Tree,"eNumTrueShowersviaDCut",eNumTrueShowersviaDCut_TreeVal);
  initTree(Tree,"eTrueHitNum",eTrueHitNum_TreeVal);
  initTree(Tree,"eTrueShowerE",eTrueShowerE_TreeVal);
  initTree(Tree,"eTrueShowerEviaECut",eTrueShowerEviaECut_TreeVal);
  initTree(Tree,"eTrueShowerEviaDCut",eTrueShowerEviaDCut_TreeVal);

  // initClusterTree(Tree, "cProjectionMatchedEnergy", cProjectionMatchedEnergy_TreeVal, fShowerModuleLabels);
  if (fClusterValidation){
    initClusterTree(Tree, "cEnergyComp", cEnergyComp_TreeVal, fShowerModuleLabels);
    initClusterTree(Tree, "cEnergyPurity", cEnergyPurity_TreeVal, fShowerModuleLabels);
    initClusterTree(Tree, "cHitsComp", cHitsComp_TreeVal, fShowerModuleLabels);
    initClusterTree(Tree, "cHitsPurity", cHitsPurity_TreeVal, fShowerModuleLabels);
  }

  for (auto const& fHitModuleLabel: fHitModuleLabels){
    std::cout << geom->Nplanes() << std::endl;
    hEnergyComp_TreeVal[fHitModuleLabel].resize(geom->Nplanes());

    std::string hitString = "hEnergyComp_" + fHitModuleLabel;
    const char* hitChar   = hitString.c_str();

    Tree->Branch(hitChar,"std::vector<std::vector<float> > >", &hEnergyComp_TreeVal[fHitModuleLabel], 32000, 0);
  }


  for (auto const& fShowerModuleLabel: fShowerModuleLabels){

    std::string processString = "sStartEndProcess_" + fShowerModuleLabel;
    const char* processChar   = processString.c_str();

    Tree->Branch(processChar,"<std::string,std::vector<std::string>", &sStartEndProcess_TreeVal[fShowerModuleLabel], 32000, 0);
  }

  Tree->Print();
}


void ana::ShowerValidation::analyze(const art::Event& evt) {

  EventRun_TreeVal = evt.run();
  EventSubrun_TreeVal = evt.subRun();
  EventNumber_TreeVal = evt.event();

  // std::cout << "Analysing Event: " << EventNumber_TreeVal << EventSubrun_TreeVal << EventRun_TreeVal << std::endl;

  ++numevents;

  //Getting  MC truth information
  art::Handle< std::vector<simb::MCTruth> > mctruthListHandle;
  std::vector<art::Ptr<simb::MCTruth> > mclist;
  if(evt.getByLabel(fGenieGenModuleLabel,mctruthListHandle))
  {art::fill_ptr_vector(mclist, mctruthListHandle);}

  //Getting the SimWire Information
  //Get the SimChannels so that we can find the IDEs deposited on them.
  art::Handle<std::vector<sim::SimChannel> > simChannelHandle;
  std::vector<art::Ptr<sim::SimChannel> > simchannels;
  if(evt.getByLabel(fLArGeantModuleLabel,simChannelHandle))
  {art::fill_ptr_vector(simchannels, simChannelHandle);}

  //Getting the Hit Information
  art::Handle<std::vector<recob::Hit> > hitListHandle;
  std::vector<art::Ptr<recob::Hit> > hits;
  if(evt.getByLabel(fHitsModuleLabel,hitListHandle))
  {art::fill_ptr_vector(hits, hitListHandle);}

  //Get the track Information (hopfully you have pandora track)
  art::Handle<std::vector<recob::Track> > trackListHandle;
  std::vector<art::Ptr<recob::Track> > tracks;
  if(evt.getByLabel(fTrackModuleLabel,trackListHandle))
  {art::fill_ptr_vector(tracks, trackListHandle);}

  //I think that doing getManyByType kind of initalises the handles giving every particle product id. Doing this allows us to find handles for the individal hits later.

  //Get all the hits
  std::vector<art::Handle<std::vector<recob::Hit> > > hitHandles;
  evt.getManyByType(hitHandles);

  //Get all the clusters
  std::vector<art::Handle<std::vector<recob::Cluster> > > clusterHandles;
  evt.getManyByType(clusterHandles);

  //Get all the pfparticles
  std::vector<art::Handle<std::vector<recob::PFParticle> > > pfpHandles;
  evt.getManyByType(pfpHandles);


  //###############################################
  //### Get the Truth information for the event ###
  //###############################################

  //List the particles in the event
  const sim::ParticleList& particles = particleInventory->ParticleList();

  //Loop over the particles
  std::map<int,const simb::MCParticle*> trueParticles;
  std::map<int,bool> mcparticlescontained;
  std::map<int,float> trueParticleEnergy;
  std::map<int,std::vector<int> > showerMothers; //Mothers are the key Daughters are in the vector.
  std::map<int,float> MCTrack_Energy_map; // Track ID to energy
  std::map<art::ProductID,std::map<int,std::map<geo::PlaneID, int> > > MCTrack_hit_map;

  int numTrueShowers =0;

  float trueShowerEnergy=-99999999;

  auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(evt);
  if(fMatchShowersInTruth){

    //Get the MC Energy deposited for each MC track.
    MCTrack_Energy_map = RecoUtils::TrueEnergyDepositedFromMCTracks(simchannels);

    //Make a map of Track id and pdgcode
    for (const auto& particleIt: particles) {
      const simb::MCParticle* particle = particleIt.second;
      trueParticles[particle->TrackId()] = particle;

      if(fVerbose > 2){std::cout << "True Particle with track ID: " << particle->TrackId() << " Has code of: " << particle->PdgCode() << " and Energy of: " << particle->E() << " With Mother: " << particle->Mother() << " Proccess: " << particle->Process() << " End Process: "  << particle->EndProcess() << std::endl;}

    }

    //Get the shower mothers. The key of the map is the mother and the vector is the daughter chain.
    showerMothers = ShowerUtils::GetShowerMothersCandidates(trueParticles);

    eNumTrueShowers_TreeVal = showerMothers.size();
    for (auto const& showerMother: showerMothers){
      const simb::MCParticle* showerMotherParticle = trueParticles.at(showerMother.first);
      eTrueShowerE_TreeVal.push_back(showerMotherParticle->E());
    }
    if(fVerbose > 1)
      std::cout << "Inital shower candidates: " << showerMothers.size() << std::endl;

    if(fRemoveNonContainedParticles)
      ShowerUtils::RemoveNoneContainedParticles(showerMothers,trueParticles,MCTrack_Energy_map);
    if(fVerbose > 1)
      std::cout<< "Shower candidates after containment cut are: " << showerMothers.size() <<std::endl;

    //Cut on the shower mothers with energy.
    ShowerUtils::CutShowerMothersByE(showerMothers, trueParticles, fSimEnergyCut);
    eNumTrueShowersviaECut_TreeVal = showerMothers.size();
    for (auto const& showerMother: showerMothers){
      const simb::MCParticle* showerMotherParticle = trueParticles.at(showerMother.first);
      eTrueShowerEviaECut_TreeVal.push_back(showerMotherParticle->E());
    }

    if(fVerbose > 1)
      std::cout<< "shower candidates after E cut are: " << showerMothers.size() <<std::endl;


    //Cut on shower mothers with density
    ShowerUtils::CutShowerMothersByDensity(clockData, showerMothers, trueParticles, hits, fDensityCut);
    eNumTrueShowersviaDCut_TreeVal = showerMothers.size();
    for (auto const& showerMother: showerMothers){
      const simb::MCParticle* showerMotherParticle = trueParticles.at(showerMother.first);
      eTrueShowerEviaDCut_TreeVal.push_back(showerMotherParticle->E());
    }
    if(fVerbose > 1)
      std::cout<< "shower candidates after density cut are: " << showerMothers.size() <<std::endl;


    //If there are no true showers we can't validate
    if(showerMothers.size() == 0){
      if(fVerbose > 0){std::cout << "No shower mothers finishing" << std::endl;}
      return;
    }

    if(fVerbose > 2){
      for (auto const& showerMother: showerMothers){
        const simb::MCParticle *motherParticle = trueParticles[showerMother.first]; //test
        std::cout << " A Mother has track id: " << showerMother.first << " with Energy: "<< motherParticle->E() << std::endl;
        for (auto const& daughterID : showerMother.second){
          std::cout << "  Daughter of id: " << daughterID << std::endl;
        }
      }
    }

    numTrueShowers = showerMothers.size();
    numshowers += numTrueShowers; //Showers in all events

    //Get the number of hits associated wit each Track. This is done for every hits handle in the event.
    for(auto& handle : hitHandles){

      if(!handle.isValid()){
        mf::LogError("ShowerValidation") << "Bad hit handle from the all the hit handles" << std::endl;
        continue;
      }

      //RawHitFinder is a bit silly at the moment so ignore it for the time being - REMOVE
      if(handle.provenance()->moduleLabel() == "fasthit"){continue;}

      //Getting the Hit Information
      art::Handle<std::vector<recob::Hit> > hitHandle;
      std::vector<art::Ptr<recob::Hit> > hits_fromhandle;
      if(evt.getByLabel(handle.provenance()->moduleLabel(),hitHandle))
      {art::fill_ptr_vector(hits_fromhandle, hitHandle);}


      //Get a map of the number of hits per plane each track has deposited.
      MCTrack_hit_map[handle.id()] = RecoUtils::NumberofPlaneHitsPerTrack(clockData, hits_fromhandle);
    }
  }

  //######################
  //### Hit Validation ###
  //######################

  if(fHitValidation){

    //Loop over the hit handles
    for (auto const& fHitModuleLabel: fHitModuleLabels) {

      //Getting the Hit Information
      art::Handle<std::vector<recob::Shower> > hitValHandle;
      std::vector<art::Ptr<recob::Shower> > hits_fromValhandle;
      if(evt.getByLabel(fHitModuleLabel,hitValHandle))
      {art::fill_ptr_vector(hits_fromValhandle,hitValHandle);}

      //Calculate the total energy deposited
      float TotalEnergyDeposited = 0;
      for(auto const& track_iter: MCTrack_Energy_map){
        TotalEnergyDeposited += track_iter.second;
      }

      //Calculate how much energy was deposited in the hits
      for(geo::PlaneID plane_id: geom->IteratePlaneIDs()){

        float TotalEnergyDepinHits = RecoUtils::TotalEnergyDepinHits(clockData, hits,plane_id.Plane);
        float hitCompleteness = TotalEnergyDepinHits!=0 ? TotalEnergyDepinHits/TotalEnergyDeposited : -99999;
        hEnergyComp_TreeVal[fHitModuleLabel][plane_id.Plane].push_back(hitCompleteness);
      }
    }
    if(fVerbose > 1){std::cout << "Hit Validation Complete" << std::endl;}
  }

  //##############################################
  //Get the reconstructed info and Match with MC #
  //##############################################

  for (auto const& fShowerModuleLabel: fShowerModuleLabels){

    if(showerMothers.size() == 0 && fMatchShowersInTruth){
      if(fVerbose > 0){std::cout << "No Mothers to Match to the showers" << std::endl;}
      continue;
    }

    //Getting the Shower Information
    art::Handle<std::vector<recob::Shower> > showerListHandle;
    std::vector<art::Ptr<recob::Shower> > recoShowers;

    if(evt.getByLabel(fShowerModuleLabel,showerListHandle))
    {
      art::fill_ptr_vector(recoShowers,showerListHandle);
    }

    if(recoShowers.size() == 0){
      if(fVerbose > 0){std::cout << "No Shower in the Event" << std::endl;}
      continue;
    }

    art::Handle<std::vector<recob::PFParticle> > pfpListHandle;
    std::vector<art::Ptr<recob::PFParticle> > pfps;
    if(evt.getByLabel(fPFParticleLabel,pfpListHandle))
    {art::fill_ptr_vector(pfps,pfpListHandle);}

    art::Handle<std::vector<recob::Track> > showertrackHandle;
    std::vector<art::Ptr<recob::Track> > showertrack;
    if(evt.getByLabel(fShowerModuleLabel,showertrackHandle))
    {art::fill_ptr_vector(showertrack,showertrackHandle);}

    //Getting the Shower Information
    //Association between Showers and 2d Hits
    art::FindManyP<recob::Hit> fmh(showerListHandle, evt, fShowerModuleLabel);

    //Association between Showers and clusters
    art::FindManyP<recob::Cluster> fmch(showerListHandle, evt, fShowerModuleLabel);

    //Association between pfParticle and clusters
    art::FindManyP<recob::Cluster> fmpfc(pfpListHandle, evt, fPFParticleLabel);

    //Association between pfParticle and vertex
    art::FindManyP<recob::Vertex> fmpfv(pfpListHandle, evt, fPFParticleLabel);

    //Association between Showers and pfParticle
    art::FindManyP<recob::PFParticle> fmpf(showerListHandle, evt, fShowerModuleLabel);

    //Association between spacepoints and showers
    art::FindManyP<recob::SpacePoint> fmsp(showerListHandle, evt, fShowerModuleLabel);

    //######################
    //### PFP Validation ###
    //######################

    std::vector< art::Ptr<recob::Vertex> > neutrinoVertices;
    std::map<int, art::Ptr<recob::PFParticle> > pfpsMap;
    for (auto const& pfp: pfps){
      pfpsMap[pfp->Self()] = pfp;
    }


    if(fVerbose > 1) { std::cout<<"Doing PFP Validation"<<std::endl; };

    std::vector<int> pfpPrimaries;
    //Create a map between PFParticles and their IDs
    int pfpAllTrackCounter  = 0;
    int pfpAllShowerCounter = 0;
    // Only the primary pfps (Daughters of the pfp netrinos)
    int pfpTrackCounter     = 0;
    int pfpShowerCounter    = 0;
    int pfpNeutrinoCounter  = 0;

    if(fVerbose > 1) { std::cout<<"Number of PFP: "<<pfps.size()<<std::endl; };

    for (auto const& pfp: pfps) {
      if(fVerbose > 2){
        std::cout<<"PFParticle: "<<pfp->Self()<<" with: PDG: "<<pfp->PdgCode()
          <<" Parent: "<<pfp->Parent()<<std::endl;
      }
      if (pfp->PdgCode()==11) ++pfpAllShowerCounter;
      if (pfp->PdgCode()==13) ++pfpAllTrackCounter;
      if ((pfp->PdgCode()==12) ||(pfp->PdgCode()==14)){ //Find Neutrino and primary daughters
        // Gives a vector of Duaghter particle ID's, get Daughters using PFParticlesMap
        const std::vector<size_t> Daughters = pfp->Daughters();
        for (auto const& daughterIter: Daughters){
          const art::Ptr<recob::PFParticle>& daughter = pfpsMap[daughterIter];
          if (fPFPValidation){
            if(fmpfc.isValid() && fmpfc.size()!=0){
              art::Handle<std::vector<recob::Cluster > > clusterHandle;
              evt.get(fmpfc.at(0).front().id(),clusterHandle);
              if(clusterHandle.isValid()){
                std::vector<art::Ptr<recob::Cluster> > pfpClusters = fmpfc.at(daughter.key());
                ana::ShowerValidation::PFPValidation(pfpClusters,daughter,evt,clockData, clusterHandle,showerMothers,MCTrack_Energy_map,MCTrack_hit_map,fShowerModuleLabel);
              }
            }
          }
          // Split into shower like, 11, and track like, 13
          if (daughter->PdgCode() == 11) {
            ++pfpShowerCounter;
          } else if (daughter->PdgCode() == 13) {
            ++pfpTrackCounter;
          } else {
            std::cout<<"Something has gone horribly wrong, PFP PDG != 11||13"<<std::endl;
            return;
          }
          pfpPrimaries.push_back(daughter->Self());
        } // end loop over neutrino daughters
        // Get the PFP neutrino vertex
        if(fmpfv.isValid() && fmpfv.size()!=0 && fmpfv.at(0).size()!=0){
          art::Handle<std::vector<recob::Vertex > > vertexHandle;
          evt.get(fmpfv.at(0).front().id(),vertexHandle);

          if(vertexHandle.isValid()) {
            std::vector< art::Ptr<recob::Vertex> > pfpVertexVector = fmpfv.at(pfp.key());
            art::Ptr<recob::Vertex> pfpVertex = pfpVertexVector.at(0);
            neutrinoVertices.push_back(pfpVertex);
          }
        }
        ++pfpNeutrinoCounter;
      }
    }

    if(fVerbose > 1) { std::cout<<"Primary Tracks: "<<pfpTrackCounter<<" and Primary Showers: "<<pfpShowerCounter<<std::endl; };

    pfpNeutrinos_TreeVal[fShowerModuleLabel].push_back(pfpNeutrinoCounter);
    pfpAllTracks_TreeVal[fShowerModuleLabel].push_back(pfpAllTrackCounter);
    pfpAllShowers_TreeVal[fShowerModuleLabel].push_back(pfpAllShowerCounter);
    pfpPrimaryTracks_TreeVal[fShowerModuleLabel].push_back(pfpTrackCounter);
    pfpPrimaryShowers_TreeVal[fShowerModuleLabel].push_back(pfpShowerCounter);
    if(fVerbose > 1) { std::cout<<"Number of PFP Neutrinos:"<<pfpNeutrinoCounter
      <<" and vertex: "<<neutrinoVertices.size()<<std::endl; };

    //TODO: Make this more general than the vertex case
    for (auto const& pfpVertex: neutrinoVertices) {

      // get the pfp neutrino vertex
      double pfpvtx[3];
      pfpVertex->XYZ(pfpvtx);

      // TODO: here i just take the first particle becuase i'm lazy, needs fixing
      std::map<int,const simb::MCParticle*>::iterator trueInitialParticleIt = trueParticles.begin();
      const simb::MCParticle* initialParticle = (*trueInitialParticleIt).second;

      const TLorentzVector PositionTrajP = initialParticle->Position();
      double truevtx[3] = {PositionTrajP.X() ,PositionTrajP.Y(), PositionTrajP.Z()};

      double PFP_Start_Diff_X = pfpvtx[0] - truevtx[0];
      double PFP_Start_Diff_Y = pfpvtx[1] - truevtx[1];
      double PFP_Start_Diff_Z = pfpvtx[2] - truevtx[2];
      double PFP_Start_Diff = TMath::Sqrt(TMath::Power((pfpvtx[0] - truevtx[0]),2) + TMath::Power((pfpvtx[1] - truevtx[1]),2) + TMath::Power((pfpvtx[2] - truevtx[2]),2));

      pfpVertexDistX_TreeVal[fShowerModuleLabel].push_back(PFP_Start_Diff_X);
      pfpVertexDistY_TreeVal[fShowerModuleLabel].push_back(PFP_Start_Diff_Y);
      pfpVertexDistZ_TreeVal[fShowerModuleLabel].push_back(PFP_Start_Diff_Z);
      pfpVertexDistMag_TreeVal[fShowerModuleLabel].push_back(PFP_Start_Diff);

    }
    // } // End pfp validation

    if(fmh.size() == 0 || fmh.at(0).size() == 0){
      std::cout << " No hits in a recob shower. Association is made incorrectly. Bailing" << std::endl;
      continue;
    }

    //Get the ID of the shower hit module
    art::ProductID showerhit_productid = fmh.at(0).front().id();

    unsigned int max_hitnum=0;
    art::Ptr<recob::Shower> biggestShower;

    // Find the bigest shower (by number of hits)
    for (auto const& showerIter: recoShowers) {
      ++numrecoshowers;
      //Get the hits vector from the shower
      std::vector< art::Ptr<recob::Hit> > showerhits = fmh.at(showerIter.key());
      if(showerhits.size() > max_hitnum){
        max_hitnum = showerhits.size();
        biggestShower = showerIter;
      }
    }

    //Loop over the showers in the event
    for (auto const& shower: recoShowers) {
      if(fUseBiggestShower == true){
        if(shower!= biggestShower){
          if(fVerbose > 0){
            std::cout << "Removing shower as it is not the biggest shower" << std::endl;
          }
          continue;
        }
      }

      //Get the information for the shower
      TVector3  ShowerDirection                  = shower->Direction();
      TVector3  ShowerStart                      = shower->ShowerStart();//cm
      double ShowerLength                        = shower->Length();//cm
      double ShowerOpeningAngle                  = shower->OpenAngle();//cm
      std::vector< double >  ShowerEnergyPlanes  = shower->Energy();//MeV
      std::vector< double >  ShowerdEdX_vec      = shower->dEdx();//MeV/cm
      double max_energy = -9999;
      //Check the energy is above the limit

      if (ShowerEnergyPlanes.size()==0){
        if(fVerbose > 0){
          std::cout << "Shower energy not set. max_energy=-9999" << std::endl;
        }
      } else {
        max_energy = *max_element(std::begin(ShowerEnergyPlanes), std::end(ShowerEnergyPlanes));
      }
      if(max_energy < fMinRecoEnergy){
        if(fVerbose > 0){
          std::cout << "Removing shower as it is below the minimum reco energy cut" << std::endl;
        }
        continue;
      }

      //Check if user wants neutrino showers.
      if(fmpf.isValid() && fUseNeutrinoShowers){
        if(fmpf.at(shower.key()).size() > 0){
          if(pfpsMap.find(fmpf.at(shower.key()).at(0)->Parent()) != pfpsMap.end()){
            if(pfpsMap[fmpf.at(shower.key()).at(0)->Parent()]->PdgCode() != 12 && pfpsMap[fmpf.at(shower.key()).at(0)->Parent()]->PdgCode() != 14){
              if(fVerbose > 0){
                std::cout << "Removing shower as it is not daughter of the neutrino" << std::endl;
              }
              continue;
            }
          }
        }
      }

      //#########################
      //### Shower Validation ###
      //#########################

      //Get the hits vector from the shower
      std::vector< art::Ptr<recob::Hit> > showerhits = fmh.at(shower.key());
      int NumberofHitsinRecoShower = showerhits.size();
      if(fVerbose > 1){std::cout << "shower hits size: " << NumberofHitsinRecoShower << std::endl;}
      if(NumberofHitsinRecoShower < fMinHitSize) {
        if(fVerbose > 0){
          std::cout << "Removing shower as it doesn't have enough hits" << std::endl;
        }
        continue;
      }

      int ShowerBest_Plane = shower->best_plane();
      // Set the collection plane as default
      if(std::isnan(ShowerBest_Plane) || ShowerBest_Plane<0)
        ShowerBest_Plane = 2;

      //################################################
      //### Hit Completeness and purity Calculations ###
      //################################################

      //Function from RecoUtils, finds the most probable track ID associated with the set of hits from there true energy depositons. The pair returns the energy as well.
      std::pair<int,double> ShowerTrackInfo = ShowerUtils::TrueParticleIDFromTrueChain(clockData, showerMothers,showerhits,ShowerBest_Plane);

      int ShowerTrackID = ShowerTrackInfo.first;
      double TrueEnergyDepWithinShower_FromTrueShower = ShowerTrackInfo.second;
      double TrueEnergyDep_FromShower = 0;
      double TrueEnergyDep_WithinRecoShower = 0;
      int TrueHitDep_FromTrueShower = 0;
      int TrueHitsDep_WithinRecoShower = 0;

      //Check to see if the shower was correctly matched.
      if(TrueEnergyDepWithinShower_FromTrueShower == -99999 && fMatchShowersInTruth){
        if(fVerbose > 0){
          std::cout << "Reco Shower not matched to true shower. Think of reducing the energy threshold" << std::endl;
        }
        continue;
      }

      if(fMatchShowersInTruth){
        //Calculate the true Energy deposited By Shower
        for (auto const& daughter: showerMothers[ShowerTrackID]){
          TrueEnergyDep_FromShower += MCTrack_Energy_map[daughter];
          for(geo::PlaneID plane_id: geom->IteratePlaneIDs()){
            TrueHitDep_FromTrueShower +=  MCTrack_hit_map[showerhit_productid][daughter][plane_id];
          }
        }

        //Get the number of hits in the reco shower from the true shower.
        std::map<int,std::map<geo::PlaneID,int> >  MCTrack_showerhit_map = RecoUtils::NumberofPlaneHitsPerTrack(clockData, showerhits);
        for(auto const& planehit_map : MCTrack_showerhit_map){
          if(std::find(showerMothers[ShowerTrackID].begin(), showerMothers[ShowerTrackID].end(),planehit_map.first) == showerMothers[ShowerTrackID].end()){continue;}
          for(auto const& hit_num : planehit_map.second){
            TrueHitsDep_WithinRecoShower += hit_num.second;
          }
        }
      } // End if fMatchShowersInTruth

      //Loop over the hits and find the IDEs
      for (auto const& hitIt: showerhits){
        //Get the plane ID
        geo::WireID wireid = hitIt->WireID();
        int PlaneID = wireid.Plane;

        if(PlaneID != ShowerBest_Plane)
          continue;

        //Split the Hit into its IDE for each track it associates with.
        std::vector<sim::TrackIDE> trackIDEs = backtracker->HitToTrackIDEs(clockData, hitIt);
        for (auto const& trackIDE: trackIDEs){

          //Find the true total energy deposited in a set of hits.
          TrueEnergyDep_WithinRecoShower += trackIDE.energy;
        }
      }//Hit Loop

      double hitcompleteness = -99999;
      double hitpurity = -99999;
      double energycompleteness = -99999;
      double energypurity = -99999;

      if(TrueHitDep_FromTrueShower != 0)
        hitcompleteness = ((double) TrueHitsDep_WithinRecoShower)/(double) TrueHitDep_FromTrueShower;

      if(NumberofHitsinRecoShower != 0)
        hitpurity = (double) TrueHitsDep_WithinRecoShower/(double) NumberofHitsinRecoShower;

      if(TrueEnergyDep_FromShower != 0)
        energycompleteness = (TrueEnergyDepWithinShower_FromTrueShower)/TrueEnergyDep_FromShower;

      if(TrueEnergyDep_WithinRecoShower != 0)
        energypurity = TrueEnergyDepWithinShower_FromTrueShower/TrueEnergyDep_WithinRecoShower;

      sAllHitsComp_TreeVal[fShowerModuleLabel].push_back(hitcompleteness);
      sAllHitsPurity_TreeVal[fShowerModuleLabel].push_back(hitpurity);
      sEnergyCompBestPlane_TreeVal[fShowerModuleLabel].push_back(energycompleteness);
      sEnergyPurityBestPlane_TreeVal[fShowerModuleLabel].push_back(energypurity);


      //################################
      //### Initial track validation ###
      //################################

      bool doInitialTrackValidation = showertrackHandle.isValid();
      art::Ptr<recob::Track> initialTrack;

      double TrueTrackLength  = -99999;
      double initialTrackEnergyComp = -99999;
      double initialTrackEnergyPurity = -99999;
      double initialTrackHitComp = -99999;
      double initialTrackHitPurity = -99999;
      double initialTrackHitCompPurity = -999;

      if (doInitialTrackValidation) {
        art::FindManyP<recob::Track> fmit(showerListHandle, evt, fShowerModuleLabel);

        if(fmit.isValid()) {
          try{
            std::vector< art::Ptr<recob::Track> > initialTracks = fmit.at(shower.key());
            if(initialTracks.size() != 1){
              mf::LogError("ShowerValidation") << "Number of initial tracks: "<<initialTracks.size()
                <<". Skipping." << std::endl;
              doInitialTrackValidation = false;
            } else {
              initialTrack = initialTracks[0];
            }
          } catch (...) {
            doInitialTrackValidation = false;
          }
        } else {
          doInitialTrackValidation = false;
        }
      }

      if(doInitialTrackValidation) {

        //See if the initial track hit are created.
        int initialtrack_hits_from_mother      = 0;
        int true_num_hits_initialTrack         = 0;
        double initialtrack_energy             = 0;
        double true_energy_initialTrack        = 0;
        double initialtrack_energy_from_mother = 0;

        //Association between shower and initial track hits
        art::FindManyP<recob::Hit> fmith(showertrackHandle, evt, fShowerModuleLabel);

        //Get the hits.
        std::vector<art::Ptr<recob::Hit> > initialtrackhits  = fmith.at(initialTrack.key());

        //Cannot do truth matching if its turned off.
        if(fMatchShowersInTruth){
          //Where do the hits come from
          std::pair<int, double>  showerInitialTrackInfo = ShowerUtils::TrueParticleIDFromTrueChain(
              clockData, showerMothers, initialtrackhits, ShowerBest_Plane);

          if (showerInitialTrackInfo.first != ShowerTrackInfo.first){
            continue;
          }

          initialtrack_energy_from_mother = showerInitialTrackInfo.second;
          initialtrack_energy = RecoUtils::TotalEnergyDepinHits(clockData, initialtrackhits, ShowerBest_Plane);
          true_energy_initialTrack = MCTrack_Energy_map.at(ShowerTrackID);

          if(trueParticles[ShowerTrackID]->PdgCode() == 11){
            initialtrack_hits_from_mother = RecoUtils::NumberofPrimaryHitsFromTrack(clockData, ShowerTrackID,initialtrackhits);
            true_num_hits_initialTrack = RecoUtils::NumberofPrimaryHitsFromTrack(clockData, ShowerTrackID,showerhits);
          } else {
            //Get the ee+ pair daughters from the photon.
            std::vector<int> Electrons;
            for(auto const& daughter: showerMothers[ShowerTrackID]){
              if(trueParticles[daughter]->Mother() == ShowerTrackID && TMath::Abs(trueParticles[daughter]->PdgCode()) == 11){
                Electrons.push_back(daughter);
              }
            }
            initialtrack_hits_from_mother = RecoUtils::NumberofPrimaryHitsWithAllTracks(clockData, Electrons, initialtrackhits);
            true_num_hits_initialTrack    = RecoUtils::NumberofPrimaryHitsWithAllTracks(clockData, Electrons, showerhits);
          }
        }

        //Calculate metrics
        if(true_energy_initialTrack != 0)
          initialTrackEnergyComp = initialtrack_energy_from_mother / true_energy_initialTrack;
        if(initialtrack_energy != 0)
          initialTrackEnergyPurity = initialtrack_energy_from_mother / initialtrack_energy;

        if(true_num_hits_initialTrack != 0)
          initialTrackHitComp = (double) initialtrack_hits_from_mother/(double) true_num_hits_initialTrack;
        if(initialtrackhits.size() != 0)
          initialTrackHitPurity = (double) initialtrack_hits_from_mother/(double) initialtrackhits.size();
        if (initialTrackHitComp>0 && initialTrackHitPurity>0)
          initialTrackHitCompPurity = initialTrackHitComp*initialTrackHitPurity;
      }

      sTrackEnergyComp_TreeVal[fShowerModuleLabel].push_back(initialTrackEnergyComp);
      sTrackEnergyPurity_TreeVal[fShowerModuleLabel].push_back(initialTrackEnergyPurity);
      sTrackHitComp_TreeVal[fShowerModuleLabel].push_back(initialTrackHitComp);
      sTrackHitPurity_TreeVal[fShowerModuleLabel].push_back(initialTrackHitPurity);
      sTrackHitCompPurity_TreeVal[fShowerModuleLabel].push_back(initialTrackHitCompPurity);


      //### TODO

      TVector3 TrueShowerDirection = {-99999,-99999,-99999};
      TVector3 TrueShowerStart     = {-99999,-99999, -99999};
      if(fMatchShowersInTruth){

        TLorentzVector PositionTrajStart = {-99999,-99999,-99999, -99999};
        TLorentzVector PositionTrajEnd   = {-99999,-99999,-99999, -99999};

        const simb::MCParticle* MCShowerParticle = trueParticles.at(ShowerTrackID);

        int shower_pdg = MCShowerParticle->PdgCode();

        //Check if another particle actually gave more energy.
        int MostHitsTrackID = TMath::Abs(RecoUtils::TrueParticleIDFromTotalRecoHits(clockData, showerhits));
        int most_hit_pdg = -99999;
        if(trueParticles.find(MostHitsTrackID) != trueParticles.end()){
          most_hit_pdg = trueParticles.at(MostHitsTrackID)->PdgCode();
        }
        sTrackPDG_TreeVal[fShowerModuleLabel].push_back(most_hit_pdg);

        sPDG_TreeVal[fShowerModuleLabel].push_back(shower_pdg);
        sMother_TreeVal[fShowerModuleLabel].push_back(MCShowerParticle->Mother());

        //Find the Energy of the particle:
        float Energy = MCShowerParticle->E();
        sTrueEnergy_TreeVal[fShowerModuleLabel].push_back(Energy*1000);

        sStartEndProcess_TreeVal[fShowerModuleLabel].push_back(MCShowerParticle->EndProcess());

        //Get the number of Traj points to loop over
        unsigned int TrajPoints = MCShowerParticle->NumberTrajectoryPoints();

        // Select first traj point where the photon loses energy, last be default
        PositionTrajStart =  MCShowerParticle->Position(0);
        PositionTrajEnd   =  MCShowerParticle->Position(TrajPoints-1);

        // For phoyons say the shower started when the photon loses 10% of the energy
        if(MCShowerParticle->PdgCode() == 22){
          for (unsigned int trajPoint=0; trajPoint<TrajPoints;trajPoint++){
            if (MCShowerParticle->E(trajPoint) < MCShowerParticle->E()){
              PositionTrajStart = MCShowerParticle->Position(trajPoint);
              break;
            }
          }
        }

        //Get the track length of the initial track stub.
        double TrackStubLength = (PositionTrajStart-PositionTrajEnd).Vect().Mag();

        //How straight are the IDE points. First Define the line.
        TVector3 track_line = (PositionTrajStart.Vect()-PositionTrajEnd.Vect()).Unit();

        double TrackStubResidual = 0;
        double TrackStubResidual_iter = 0;
        //Get the energy depositions
        std::vector< const sim::IDE * > IDEs = backtracker->TrackIdToSimIDEs_Ps(ShowerTrackID);
        for(auto const& ide: IDEs){
          //Calculate the perpendicualr distance from the straight line
          if(ide->trackID == ShowerTrackID){
            TVector3 pos = {ide->x, ide->y,ide->z};
            TVector3 rel_pos = pos - PositionTrajStart.Vect();
            double   len  = rel_pos.Dot(track_line);
            if(len>TrackStubLength){continue;}
            double   perp = (rel_pos - len*track_line).Mag();
            TrackStubResidual += perp;
            ++TrackStubResidual_iter;
          }
        }

        sTrueTrackLengh_TreeVal[fShowerModuleLabel].push_back(TrackStubLength);
        sTrueTrackSpread_TreeVal[fShowerModuleLabel].push_back(TrackStubResidual/TrackStubResidual_iter);

        //The three vector for track length is the shower direction
        TrueShowerDirection = MCShowerParticle->Momentum().Vect();
        TrueShowerStart = PositionTrajStart.Vect();

        //Initial track lentgh of the shower. NEEDS UPDATING TODO
        TrueTrackLength = TrueShowerDirection.Mag();
      }


      //###########################
      //### Projection Matching ###
      //###########################

      float geoprojectionmatched_score = -99999;
      if(fmsp.isValid()){
        //Get the spacepoints
        std::vector<art::Ptr<recob::SpacePoint> > sps = fmsp.at(shower.key());

        if(sps.size() > 0){
          art::Handle<std::vector<recob::SpacePoint> > spHandle;
          evt.get(sps.front().id(),spHandle);

          if(spHandle.isValid()){
            art::FindManyP<recob::Hit> fmsph(spHandle, evt, spHandle.provenance()->moduleLabel());

            float geomatched = 0;
            for(auto const& sp: sps){

              //Get the the hit
              std::vector< art::Ptr<recob::Hit> > hit = fmsph.at(sp.key());
              if(hit.size() < 1)
                continue;

              //Get the true position
              std::vector<double> hitXYZ = {-999,-999,-999};
              try{
                hitXYZ = backtracker->HitToXYZ(clockData, hit[0]);
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
              if((trueposition-sp_postiion).Mag() < fMatchedCut)
                ++geomatched;
            }
            geoprojectionmatched_score = geomatched/(float) sps.size();
          }
        }
      }
      // std::cout<< "finished to spacepoint" << std::endl;

      //#############################
      //### Calculate the metrics ###
      //#############################

      //Get the angles between the direction
      float ShowerDirection_Xdiff = -99999;
      float ShowerDirection_Ydiff = -99999;
      float ShowerDirection_Zdiff = -99999;
      float ShowerDirection_XTrue = -99999;
      float ShowerDirection_YTrue = -99999;
      float ShowerDirection_ZTrue = -99999;
      float ShowerDirection_diff  = -99999;
      double Start_diff           = -99999;
      double Start_diff_X         = -99999;
      double Start_diff_Y         = -99999;
      double Start_diff_Z         = -99999;
      double sdEdx                = -99999;
      double sEnergy              = -99999;
      double sEnergyRat           = -99999;
      double sEnergyDiff          = -99999;

      if(fMatchShowersInTruth) {
        if(TMath::Sqrt((TrueShowerDirection.Y()*TrueShowerDirection.Y() + TrueShowerDirection.Z()*TrueShowerDirection.Z()))*TMath::Sqrt((ShowerDirection.Y()*ShowerDirection.Y() + ShowerDirection.Z()*ShowerDirection.Z())) !=0){
          ShowerDirection_Xdiff = (TrueShowerDirection.Y()*ShowerDirection.Y() + TrueShowerDirection.Z()*ShowerDirection.Z())/(TMath::Sqrt((TrueShowerDirection.Y()*TrueShowerDirection.Y() + TrueShowerDirection.Z()*TrueShowerDirection.Z()))*TMath::Sqrt((ShowerDirection.Y()*ShowerDirection.Y() + ShowerDirection.Z()*ShowerDirection.Z())));
          ShowerDirection_XTrue = TrueShowerDirection.X();
        }

        if(TMath::Sqrt((TrueShowerDirection.X()*TrueShowerDirection.X() + TrueShowerDirection.Z()*TrueShowerDirection.Z()))*TMath::Sqrt((ShowerDirection.X()*ShowerDirection.X() + ShowerDirection.Z()*ShowerDirection.Z())) !=0){
          ShowerDirection_Ydiff = (TrueShowerDirection.X()*ShowerDirection.X() + TrueShowerDirection.Z()*ShowerDirection.Z())/(TMath::Sqrt((TrueShowerDirection.X()*TrueShowerDirection.X() + TrueShowerDirection.Z()*TrueShowerDirection.Z()))*TMath::Sqrt((ShowerDirection.X()*ShowerDirection.X() + ShowerDirection.Z()*ShowerDirection.Z())));
          ShowerDirection_YTrue = TrueShowerDirection.Y();
        }

        if(TMath::Sqrt((TrueShowerDirection.Y()*TrueShowerDirection.Y() + TrueShowerDirection.X()*TrueShowerDirection.X()))*TMath::Sqrt((ShowerDirection.Y()*ShowerDirection.Y() + ShowerDirection.X()*ShowerDirection.X())) !=0){
          ShowerDirection_Zdiff = (TrueShowerDirection.Y()*ShowerDirection.Y() + TrueShowerDirection.X()*ShowerDirection.X())/(TMath::Sqrt((TrueShowerDirection.Y()*TrueShowerDirection.Y() + TrueShowerDirection.X()*TrueShowerDirection.X()))*TMath::Sqrt((ShowerDirection.Y()*ShowerDirection.Y() + ShowerDirection.X()*ShowerDirection.X())));
          ShowerDirection_ZTrue = TrueShowerDirection.Z();
        }

        if(TrueShowerDirection.Mag() != 0 || ShowerDirection.Mag() !=0){
          ShowerDirection_diff  = TrueShowerDirection.Dot(ShowerDirection)/(TrueShowerDirection.Mag()*ShowerDirection.Mag());
        }

        //Get the Error in the position. intialised as 0,0,0 this is a problem here.
        Start_diff = (TrueShowerStart - ShowerStart).Mag();
        Start_diff_X = (TrueShowerStart.X() - ShowerStart.X());
        Start_diff_Y = (TrueShowerStart.Y() - ShowerStart.Y());
        Start_diff_Z = (TrueShowerStart.Z() - ShowerStart.Z());
      }

      if (ShowerdEdX_vec.size()!=0)
        sdEdx = (ShowerdEdX_vec.at(ShowerBest_Plane));

      if (ShowerEnergyPlanes.size()!=0)
        sEnergy = (ShowerEnergyPlanes.at(ShowerBest_Plane));

      if(TrueEnergyDep_FromShower!=0 && sEnergy>0){
        sEnergyRat = sEnergy / TrueEnergyDep_FromShower;
        sEnergyDiff = (sEnergy - TrueEnergyDep_FromShower)/TrueEnergyDep_FromShower;
      }

      //Fill the histograms.
      sDirX_TreeVal[fShowerModuleLabel].push_back(ShowerDirection_Xdiff);
      sDirY_TreeVal[fShowerModuleLabel].push_back(ShowerDirection_Ydiff);
      sDirZ_TreeVal[fShowerModuleLabel].push_back(ShowerDirection_Zdiff);
      sTrueDirX_TreeVal[fShowerModuleLabel].push_back(ShowerDirection_XTrue);
      sTrueDirY_TreeVal[fShowerModuleLabel].push_back(ShowerDirection_YTrue);
      sTrueDirZ_TreeVal[fShowerModuleLabel].push_back(ShowerDirection_ZTrue);
      sDirDiff_TreeVal[fShowerModuleLabel].push_back(ShowerDirection_diff);

      sStartX_TreeVal[fShowerModuleLabel].push_back(Start_diff_X);
      sStartY_TreeVal[fShowerModuleLabel].push_back(Start_diff_Y);
      sStartZ_TreeVal[fShowerModuleLabel].push_back(Start_diff_Z);
      sStartDist_TreeVal[fShowerModuleLabel].push_back(Start_diff);

      sBestPlane_TreeVal[fShowerModuleLabel].push_back(ShowerBest_Plane);
      sLength_TreeVal[fShowerModuleLabel].push_back(ShowerLength);
      sOpeningAngle_TreeVal[fShowerModuleLabel].push_back(ShowerOpeningAngle);
      sNumHits_TreeVal[fShowerModuleLabel].push_back(NumberofHitsinRecoShower);
      sGeoProjectionMatched_TreeVal[fShowerModuleLabel].push_back(geoprojectionmatched_score);

      sEnergy_TreeVal[fShowerModuleLabel].push_back(sEnergy);
      sEnergyRat_TreeVal[fShowerModuleLabel].push_back(sEnergyRat);
      sEnergyDiff_TreeVal[fShowerModuleLabel].push_back(sEnergyDiff);

      sdEdx_TreeVal[fShowerModuleLabel].push_back(sdEdx);

      if(fVerbose > 1){
        std::cout << "#################################################" << std::endl;
        std::cout << "Global Event Information" << std::endl;
        std::cout << "#################################################" << std::endl;
        std::cout << "Shower Label: " <<  fShowerModuleLabel << std::endl;
        std::cout << "Number of Reco showers: " << recoShowers.size() << std::endl;
        std::cout << "Energy Simulated: " << trueShowerEnergy << std::endl;
        std::cout << "Number of True Showers Before Cuts: " <<  numTrueShowers << std::endl;
        std::cout << "Number of True Showers in the Event: " << showerMothers.size() << std::endl;
        std::cout << "#################################################" << std::endl;
        std::cout << "Reconstructed/Associated Truth Information" << std::endl;
        std::cout << "#################################################" << std::endl;
        std::cout << "Hit Size: " << NumberofHitsinRecoShower << std::endl;
        std::cout << "ShowerBest_Plane: " << ShowerBest_Plane << std::endl;
        std::cout << "True Start: " << TrueShowerStart.X() << " Shower Start: " << ShowerStart.X() << std::endl;
        std::cout << "X Poisition: " <<  ShowerStart.X() << "Y Position " << ShowerStart.Y() << " Z Poistion: " << ShowerStart.Z() << std::endl;
        std::cout << "TrueTrackLength: " << TrueTrackLength << " and ShowerLength: " << ShowerLength << std::endl;
        std::cout << "ShowerOpeningAngle: " << ShowerOpeningAngle << std::endl;
        std::cout << "Best Plane " << ShowerBest_Plane << std::endl;
        std::cout << "Best Plane Reco Shower Energy " << sEnergy << std::endl;
        std::cout << "Best Plane Reco Shower Energy " << sdEdx << std::endl;
        std::cout << "True Energy Deposited from true shower within Shower: " << TrueEnergyDepWithinShower_FromTrueShower << std::endl;
        std::cout << "True Energy Deposited From true associated Shower: " << TrueEnergyDep_FromShower << std::endl;
        std::cout << "True Energy Deposited by hits in shower: " <<  TrueEnergyDep_WithinRecoShower << std::endl;
        std::cout << "Energy Purity: " << energypurity << " Energy completeness: " << energycompleteness << std::endl;
        std::cout << "Hit Purity: " << hitpurity << "Hit Completeness: " << hitcompleteness << std::endl;
        std::cout << "#################################################" <<std::endl;
      }

      if(fVerbose > 0){std::cout << "Shower Validation Complete" << std::endl;}


      //##########################
      //### Cluster Validation ###
      //##########################

      if(fClusterValidation){

        //Get the clusters associated to the shower.Needs Work.
        if(fmch.isValid()){
          art::Handle<std::vector<recob::Cluster > > clusterHandle;
          evt.get(fmch.at(shower.key()).front().id(),clusterHandle);
          if(clusterHandle.isValid()){
            std::vector<art::Ptr<recob::Cluster> > showerclusters = fmch.at(shower.key());
            ana::ShowerValidation::ClusterValidation(showerclusters,evt,clockData,clusterHandle,showerMothers,MCTrack_Energy_map,MCTrack_hit_map,ShowerTrackID,trueShowerEnergy,fShowerModuleLabel);
          }
          else{
            mf::LogError("ShowerValidation") << "Cluster handle is stale. No clustering validation done" << std::endl;
          }
        }
        else if(fmpf.isValid()){
          //Find the Clusters associated to PF particle.
          art::Handle<std::vector<recob::PFParticle> > pfpHandle;
          evt.get(fmpf.at(shower.key()).front().id(),pfpHandle);
          if(pfpHandle.isValid()){
            art::FindManyP<recob::Cluster> fmcpf(pfpHandle, evt, pfpHandle.provenance()->moduleLabel());
            if(fmcpf.isValid()){
              art::Handle<std::vector<recob::Cluster > > clusterHandle;
              evt.get(fmcpf.at(0).front().id(),clusterHandle);
              if(clusterHandle.isValid()){
                std::vector< art::Ptr<recob::Cluster> > showerclusters = fmcpf.at(shower.key());
                ana::ShowerValidation::ClusterValidation(showerclusters,evt,clockData,clusterHandle,showerMothers,MCTrack_Energy_map,MCTrack_hit_map,ShowerTrackID,trueShowerEnergy,fShowerModuleLabel);
              }
              else{
                mf::LogError("ShowerValidation") << "Cluster handle is stale. No clustering validation done" << std::endl;
              }
            }
            else{
              mf::LogError("ShowerValidation") << "No Assosoication between pf particles and clusters was found for shower made from a pf particle. No clustering validation was done." << std::endl;
            }
          }
          else{
            mf::LogError("ShowerValidation") << "pf particle handle is stale" << std::endl;
          }
        }
        else{
          mf::LogError("ShowerValidation") << "No cluster or pandora association" << std::endl;
        }

        if(fVerbose > 1){std::cout << "Cluster Validation Complete" << std::endl;}

      }

      ++numrecoshowersana;


    }//Shower Loop

    //#########################
    //###   Event Metric    ###
    //#########################

    //Calculate the True Hit number
    for (auto const& showermother: showerMothers){
      int TrueHitDep_FromTrueShowers = 0;
      for (auto const& track_id: showermother.second){
        for(geo::PlaneID plane_id: geom->IteratePlaneIDs()){
          TrueHitDep_FromTrueShowers +=  MCTrack_hit_map[showerhit_productid][track_id][plane_id];
        }
      }
      eTrueHitNum_TreeVal.push_back(TrueHitDep_FromTrueShowers);
    }

    //This is rather pandora track specific for tuning pandora so lets add it if the shower module used pandora and the track module Pandora track is used. Very horrid and sorry for that.
    eNumRecoTracks_TreeVal[fShowerModuleLabel].push_back(tracks.size());

    //Whats the segementyness of the event.
    eSegmentation_TreeVal[fShowerModuleLabel].push_back((float)recoShowers.size()/(float) numTrueShowers);
    eNumRecoShowers_TreeVal[fShowerModuleLabel].push_back(recoShowers.size());

    if (fVerbose>2)
      std::cout<<"finished shower module: "<<fShowerModuleLabel<<std::endl;
}//Shower Module labels

//Fill the tree
Tree->Fill();

for (auto const& fShowerModuleLabel: fShowerModuleLabels) {
  sDirX_TreeVal[fShowerModuleLabel].clear();
  sDirY_TreeVal[fShowerModuleLabel].clear();
  sDirZ_TreeVal[fShowerModuleLabel].clear();
  sTrueDirX_TreeVal[fShowerModuleLabel].clear();
  sTrueDirY_TreeVal[fShowerModuleLabel].clear();
  sTrueDirZ_TreeVal[fShowerModuleLabel].clear();
  sDirDiff_TreeVal[fShowerModuleLabel].clear();
  sStartX_TreeVal[fShowerModuleLabel].clear();
  sStartY_TreeVal[fShowerModuleLabel].clear();
  sStartZ_TreeVal[fShowerModuleLabel].clear();
  sStartDist_TreeVal[fShowerModuleLabel].clear();
  sLength_TreeVal[fShowerModuleLabel].clear();
  sOpeningAngle_TreeVal[fShowerModuleLabel].clear();
  sEnergyRat_TreeVal[fShowerModuleLabel].clear();
  sdEdx_TreeVal[fShowerModuleLabel].clear();
  sEnergyCompBestPlane_TreeVal[fShowerModuleLabel].clear();
  sAllHitsPurity_TreeVal[fShowerModuleLabel].clear();
  sAllHitsComp_TreeVal[fShowerModuleLabel].clear();
  sEnergyPurityBestPlane_TreeVal[fShowerModuleLabel].clear();
  sEnergy_TreeVal[fShowerModuleLabel].clear();
  sNumHits_TreeVal[fShowerModuleLabel].clear();
  sEnergyDiff_TreeVal[fShowerModuleLabel].clear();
  sEnergyDiffTrue_TreeVal[fShowerModuleLabel].clear();
  sTrueEnergy_TreeVal[fShowerModuleLabel].clear();
  sBestPlane_TreeVal[fShowerModuleLabel].clear();
  sGeoProjectionMatched_TreeVal[fShowerModuleLabel].clear();
  sTrackEnergyComp_TreeVal[fShowerModuleLabel].clear();
  sTrackEnergyPurity_TreeVal[fShowerModuleLabel].clear();
  sTrackHitComp_TreeVal[fShowerModuleLabel].clear();
  sTrackHitPurity_TreeVal[fShowerModuleLabel].clear();
  sTrackHitCompPurity_TreeVal[fShowerModuleLabel].clear();
  sTrueTrackLengh_TreeVal[fShowerModuleLabel].clear();
  sTrueTrackSpread_TreeVal[fShowerModuleLabel].clear();
  sPDG_TreeVal[fShowerModuleLabel].clear();
  sTrackPDG_TreeVal[fShowerModuleLabel].clear();
  sMother_TreeVal[fShowerModuleLabel].clear();

  pfpNeutrinos_TreeVal[fShowerModuleLabel].clear();
  pfpPrimaryTracks_TreeVal[fShowerModuleLabel].clear();
  pfpPrimaryShowers_TreeVal[fShowerModuleLabel].clear();
  pfpAllTracks_TreeVal[fShowerModuleLabel].clear();
  pfpAllShowers_TreeVal[fShowerModuleLabel].clear();
  pfpVertexDistX_TreeVal[fShowerModuleLabel].clear();
  pfpVertexDistY_TreeVal[fShowerModuleLabel].clear();
  pfpVertexDistZ_TreeVal[fShowerModuleLabel].clear();
  pfpVertexDistMag_TreeVal[fShowerModuleLabel].clear();
  pfpProjectionMatched_TreeVal[fShowerModuleLabel].clear();
  pfpHitsComp_TreeVal[fShowerModuleLabel].clear();
  pfpEnergyComp_TreeVal[fShowerModuleLabel].clear();
  pfpHitsPurity_TreeVal[fShowerModuleLabel].clear();
  pfpEnergyPurity_TreeVal[fShowerModuleLabel].clear();
  pfpTrackProjectionMatched_TreeVal[fShowerModuleLabel].clear();
  pfpTrackHitsComp_TreeVal[fShowerModuleLabel].clear();
  pfpTrackEnergyComp_TreeVal[fShowerModuleLabel].clear();
  pfpTrackHitsPurity_TreeVal[fShowerModuleLabel].clear();
  pfpTrackEnergyPurity_TreeVal[fShowerModuleLabel].clear();
  pfpShowerProjectionMatched_TreeVal[fShowerModuleLabel].clear();
  pfpShowerHitsComp_TreeVal[fShowerModuleLabel].clear();
  pfpShowerEnergyComp_TreeVal[fShowerModuleLabel].clear();
  pfpShowerHitsPurity_TreeVal[fShowerModuleLabel].clear();
  pfpShowerEnergyPurity_TreeVal[fShowerModuleLabel].clear();

  eSegmentation_TreeVal[fShowerModuleLabel].clear();
  eNumRecoTracks_TreeVal[fShowerModuleLabel].clear();
  eNumRecoShowers_TreeVal[fShowerModuleLabel].clear();

  sStartEndProcess_TreeVal[fShowerModuleLabel].clear();

  if (fClusterValidation){
    for (unsigned int plane=0; plane<geom->Nplanes(); plane++) {
      // cProjectionMatchedEnergy_TreeVal[fShowerModuleLabel].clear();
      cEnergyComp_TreeVal[fShowerModuleLabel][plane].clear();
      cEnergyPurity_TreeVal[fShowerModuleLabel][plane].clear();
      cHitsComp_TreeVal[fShowerModuleLabel][plane].clear();
      cHitsPurity_TreeVal[fShowerModuleLabel][plane].clear();
    }
  }
}
eNumTrueShowers_TreeVal = -9999;
eNumTrueShowersviaECut_TreeVal = -9999;
eNumTrueShowersviaDCut_TreeVal = -9999;
eTrueHitNum_TreeVal.clear();
eTrueShowerE_TreeVal.clear();
eTrueShowerEviaECut_TreeVal.clear();
eTrueShowerEviaDCut_TreeVal.clear();

for(unsigned int hitlab_it=0; hitlab_it<fHitModuleLabels.size(); ++hitlab_it){
  for(unsigned int plane_it=0; plane_it<geom->Nplanes(); ++plane_it){
    hEnergyComp_TreeVal[fHitModuleLabels[hitlab_it]][plane_it].clear();
  }
}

return;
}


void ana::ShowerValidation::ClusterValidation(std::vector< art::Ptr<recob::Cluster> >& clusters,
    const art::Event& evt,
    const detinfo::DetectorClocksData& clockData,
    art::Handle<std::vector<recob::Cluster> >& clusterHandle,
    std::map<int,std::vector<int> >& ShowerMotherTrackIDs,
    std::map<int,float>& MCTrack_Energy_map,
    std::map<art::ProductID,std::map<int,std::map<geo::PlaneID, int> > >& MCTrack_hit_map,
    int& TrueShowerID,
    float& trueShowerEnergy,
    const std::string& fShowerModuleLabel){


  //Initialise Trees
  unsigned int numPlanes = (int)geom->Nplanes();
  std::vector<std::vector<float> > cProjectionMatchedEnergy_TreeVec(numPlanes);
  std::vector<std::vector<float> > cEnergyComp_TreeVec(numPlanes);
  std::vector<std::vector<float> > cEnergyPurity_TreeVec(numPlanes);
  std::vector<std::vector<float> > cHitsComp_TreeVec(numPlanes);
  std::vector<std::vector<float> > cHitsPurity_TreeVec(numPlanes);

  //Get the associated hits
  art::FindManyP<recob::Hit> fmhc(clusterHandle, evt, clusterHandle.provenance()->moduleLabel());

  //Holder for cluster its
  art::Handle<std::vector<recob::Hit > > hitHandle;

  // std::cout<<"Clusters size: "<<clusters.size()<<" and "<<fmhc.at(clusters.at(0).key()).size()<<" and "<<clusterHandle.provenance()->moduleLabel()<<std::endl;

  if (fmhc.at(clusters.at(0).key()).size()==0) {
    mf::LogError("ShowerValidation") << "Cluster has no hits. Trying next cluster."
      << std::endl;
    return;
  };

  //Get the Hits Handle used for this cluster type WARNING
  evt.get(fmhc.at(clusters.at(0).key()).front().id(),hitHandle);

  if(!hitHandle.isValid()){
    mf::LogError("ShowerValidation") << "Hits handle is stale. No clustering validation done" << std::endl;
    return;
  }

  //Get the hits vector from the shower
  for(auto const& cluster : clusters){

    std::vector< art::Ptr<recob::Hit> > clusterhits = fmhc.at(cluster.key());
    //Function from RecoUtils, finds the most probable track ID associated with the set of hits from there true energy depositons. The pair returns the energy in the hits as well.
    std::pair<int,double> ShowerTrackInfo = ShowerUtils::TrueParticleIDFromTrueChain(clockData, ShowerMotherTrackIDs,clusterhits, cluster->Plane().Plane);

    //Make sure the cluster has been matched.
    if(ShowerTrackInfo.second == -99999){
      if (fVerbose > 0){
        std::cout << "Reco cluster not matched to a True shower" << std::endl;
      }
      continue;
    }

    float TotalTrueEnergy = 0;
    float signalhits      = 0;
    float totalhits       = 0;

    std::map<int, std::map<geo::PlaneID, int> > hitPlaneMap = RecoUtils::NumberofPlaneHitsPerTrack(clockData, clusterhits);

    for (auto const& daughterID: ShowerMotherTrackIDs[ShowerTrackInfo.first]) {

      //Calculate the true Energy deposited By Shower
      TotalTrueEnergy += MCTrack_Energy_map[daughterID];
      // Get the map of planeIDs to number of hits

      auto hitPlaneMapIter = hitPlaneMap.find(daughterID);

      if (hitPlaneMapIter==hitPlaneMap.end()){
        if (fVerbose>2){
          std::cout<<"daughter: "<<daughterID<<" not found in cluster hit map"<<std::endl;
        }
        continue;
      }
      for (const auto& hitPlaneIt: hitPlaneMapIter->second){
        //std::cout<<"Plane ID: "<<(*hitPlaneMapit).first<<std::endl;
        //std::cout<<"SignalHits: "<<(*hitPlaneMapit).second<<std::endl;
        //std::cout<<"TotalHits: "<<MCTrack_hit_map[hitHandle.id()][daughterID][(*hitPlaneMapit).first]<<std::endl;

        //Count how many hits are from the true shower.
        signalhits += hitPlaneIt.second;
        //Count how many hits are missed in the plane from the true track id.
        totalhits += MCTrack_hit_map[hitHandle.id()][daughterID][hitPlaneIt.first];
      }
    }

    int projection_match = -99999;

    //Have we matched the 2D cluster to the correct shower correctly. In terms of Energy depositions:
    if(ShowerTrackInfo.first == TrueShowerID){projection_match = 1;}
    else{projection_match = 0;}

    //Calculate the purity and completeness metrics
    float completeness_hits   = -99999;
    float purity_hits         = -99999;
    float completeness_energy = -99999;
    float purity_energy       = -99999;

    float TotalEnergyDepinHits = RecoUtils::TotalEnergyDepinHits(clockData, clusterhits,cluster->Plane().Plane);

    cProjectionMatchedEnergy_TreeVec[cluster->Plane().Plane].push_back(projection_match);

    if(totalhits != 0)
      completeness_hits = (signalhits)/totalhits;

    if(clusterhits.size() != 0)
      purity_hits = signalhits/clusterhits.size();

    if(TotalTrueEnergy != 0)
      completeness_energy = (ShowerTrackInfo.second)/TotalTrueEnergy;

    if(TotalEnergyDepinHits != 0)
      purity_energy = ShowerTrackInfo.second/TotalEnergyDepinHits;

    cHitsComp_TreeVec[cluster->Plane().Plane].push_back(completeness_hits);
    cHitsPurity_TreeVec[cluster->Plane().Plane].push_back(purity_hits);
    cEnergyComp_TreeVec[cluster->Plane().Plane].push_back(completeness_energy);
    cEnergyPurity_TreeVec[cluster->Plane().Plane].push_back(purity_energy);

    if(fVerbose>1){
      std::cout << "#################################################"    << std::endl;
      std::cout << "               Cluster Metrics                   "    << std::endl;
      std::cout << "#################################################"    << std::endl;
      std::cout << "Projection matched:          " << projection_match    << std::endl;
      std::cout << "Cluster hit completeness:    " << completeness_hits   << std::endl;
      std::cout << "Cluster hit purity:          " << purity_hits         << std::endl;
      std::cout << "Cluster energy completeness: " << completeness_energy << std::endl;
      std::cout << "Cluster energy purity:       " << purity_energy       << std::endl;
      std::cout << "#################################################"    << std::endl;
    }
  }//Cluster Loop

  for (unsigned int plane=0; plane<numPlanes; plane++) {
    // cProjectionMatchedEnergy_TreeVal[fShowerModuleLabel][plane].push_back(cProjectionMatchedEnergy_TreeVec.at(plane));
    cEnergyComp_TreeVal[fShowerModuleLabel][plane].push_back(cEnergyComp_TreeVec.at(plane));
    cEnergyPurity_TreeVal[fShowerModuleLabel][plane].push_back(cEnergyPurity_TreeVec.at(plane));
    cHitsComp_TreeVal[fShowerModuleLabel][plane].push_back(cHitsComp_TreeVec.at(plane));
    cHitsPurity_TreeVal[fShowerModuleLabel][plane].push_back(cHitsPurity_TreeVec.at(plane));
  }
  return;
}


void ana::ShowerValidation::PFPValidation(std::vector<art::Ptr<recob::Cluster> >& clusters,
    art::Ptr<recob::PFParticle>  pfp,
    const art::Event& evt,
    const detinfo::DetectorClocksData& clockData,
    art::Handle<std::vector<recob::Cluster> >& clusterHandle,
    std::map<int,std::vector<int> >& ShowerMotherTrackIDs,
    std::map<int,float>& MCTrack_Energy_map,
    std::map<art::ProductID,std::map<int,std::map<geo::PlaneID, int> > > & MCTrack_hit_map,
    const std::string & fShowerModuleLabel){


  //Get the associated hits
  art::FindManyP<recob::Hit> fmhc(clusterHandle, evt, clusterHandle.provenance()->moduleLabel());

  //Holder for cluster hits

  //Get the Hits Handle used for this cluster type
  art::Handle<std::vector<recob::Hit > > hitHandle;
  evt.get(fmhc.at(clusters.at(0).key()).front().id(),hitHandle);

  if(!hitHandle.isValid()){
    mf::LogError("ShowerValidation") << "Hits handle is stale. No pfp validation done" << std::endl;
    return;
  }

  float TotalTrueEnergy      = 0;
  float TotalEnergyDepinHits = 0;
  float signalhits           = 0;
  float totalhits            = 0;
  float pfphits              = 0;
  float pfpenergy            = 0;

  //Calculate the purity and completeness metrics
  float completeness_hits   = -99999;
  float purity_hits         = -99999;
  float completeness_energy = -99999;
  float purity_energy       = -99999;

  int projectionMatched = 1;
  std::vector<int> projected_IDs;

  //Get the hits vector from the shower
  for(auto const& cluster : clusters){
    std::vector< art::Ptr<recob::Hit> > clusterhits = fmhc.at(cluster.key());

    std::map<int, std::map<geo::PlaneID, int> > hitPlaneMap = RecoUtils::NumberofPlaneHitsPerTrack(clockData, clusterhits);


    TotalEnergyDepinHits += RecoUtils::TotalEnergyDepinHits(clockData, clusterhits,cluster->Plane().Plane);
    pfphits      += clusterhits.size();

    std::pair<int,double> ShowerTrackInfo;
    std::vector<int>  daughters;

    // if shower-like use showerutils, else if track-like use recoutils
    int pdg = abs(pfp->PdgCode()); // Track or shower
    if (pdg==11 || pdg==22) {
      //Function from RecoUtils, finds the most probable track ID associated with the set of hits from there true energy depositons. The pair returns the energy in the hits as well.
      ShowerTrackInfo = ShowerUtils::TrueParticleIDFromTrueChain(clockData, ShowerMotherTrackIDs,clusterhits, cluster->Plane().Plane);


      //Make sure the cluster has been matched.
      if(ShowerTrackInfo.second == -99999){
        if (fVerbose>0){
          std::cout << "Reco cluster not matched to a True shower" << std::endl;
        }
        continue;
      }

      daughters = ShowerMotherTrackIDs[ShowerTrackInfo.first];

    } else {
      int PFPTrackInfo = RecoUtils::TrueParticleIDFromTotalRecoHits(clockData, clusterhits);

      if(PFPTrackInfo == -99999){
        if (fVerbose>0){
          std::cout << "Reco cluster not matched to a True shower" << std::endl;
        }
        continue;
      }
      double energy = RecoUtils::TotalEnergyDepinHitsFromTrack(clockData, clusterhits, PFPTrackInfo, cluster->Plane().Plane);//check with dom
      daughters = {PFPTrackInfo};

      ShowerTrackInfo = std::pair<int, double>(PFPTrackInfo,energy);
    }

    for (auto const& daughter: daughters) {
      TotalTrueEnergy += MCTrack_Energy_map[daughter];

      auto hitPlaneMapIter = hitPlaneMap.find(daughter);
      if (hitPlaneMapIter==hitPlaneMap.end()){
        if (fVerbose>2){
          std::cout<<"daughter: "<<daughter<<" not found in cluster hit map"<<std::endl;
        }
        continue;
      }
      for (const auto& hitPlaneIt: hitPlaneMapIter->second){
        //Count how many hits are from the true shower.
        signalhits += hitPlaneIt.second;
        //Count how many hits are missed in the plane from the true track id.
        totalhits += MCTrack_hit_map[hitHandle.id()][daughter][hitPlaneIt.first];
      }
      projected_IDs.push_back(ShowerTrackInfo.second);
    }
  }

  //Have we matched the 2D cluster to the correct shower correctly. In terms of Energy depositions:
  //if(ShowerTrackInfo.first == TrueShowerID){projection_match = 1;}
  //else{projection_match = 0;}

  for (auto ID : projected_IDs){
    if (ID != projected_IDs.at(0)) {projectionMatched=0;}
  }

  pfpProjectionMatched_TreeVal[fShowerModuleLabel].push_back(projectionMatched);
  if (pfp->PdgCode() == 11) { pfpShowerProjectionMatched_TreeVal[fShowerModuleLabel].push_back(projectionMatched);
  } else if (pfp->PdgCode() == 13) { pfpTrackProjectionMatched_TreeVal[fShowerModuleLabel].push_back(projectionMatched);
  }

  //cProjectionMatchedEnergy_TreeVec[cluster->Plane().Plane].push_back(projection_match);
  if(totalhits != 0){
    completeness_hits = (signalhits)/totalhits;
  }

  if(pfphits != 0){
    purity_hits = signalhits/pfphits;
  }

  if(TotalTrueEnergy != 0){
    completeness_energy = pfpenergy/TotalTrueEnergy;
  }

  if(TotalEnergyDepinHits != 0){
    purity_energy = pfpenergy/TotalEnergyDepinHits;
  }

  pfpHitsComp_TreeVal[fShowerModuleLabel].push_back(completeness_hits);
  pfpHitsPurity_TreeVal[fShowerModuleLabel].push_back(purity_hits);
  pfpEnergyComp_TreeVal[fShowerModuleLabel].push_back(completeness_energy);
  pfpEnergyPurity_TreeVal[fShowerModuleLabel].push_back(purity_energy);
  if (pfp->PdgCode() == 11) {
    pfpShowerHitsComp_TreeVal[fShowerModuleLabel].push_back(completeness_hits);
    pfpShowerHitsPurity_TreeVal[fShowerModuleLabel].push_back(purity_hits);
    pfpShowerEnergyComp_TreeVal[fShowerModuleLabel].push_back(completeness_energy);
    pfpShowerEnergyPurity_TreeVal[fShowerModuleLabel].push_back(purity_energy);
  } else if (pfp->PdgCode() == 13) {
    pfpTrackHitsComp_TreeVal[fShowerModuleLabel].push_back(completeness_hits);
    pfpTrackHitsPurity_TreeVal[fShowerModuleLabel].push_back(purity_hits);
    pfpEnergyComp_TreeVal[fShowerModuleLabel].push_back(completeness_energy);
    pfpEnergyPurity_TreeVal[fShowerModuleLabel].push_back(purity_energy);
  }


  if(fVerbose>1){
    std::cout << "#################################################"    << std::endl;
    std::cout << "                 PFP Metrics                     "    << std::endl;
    std::cout << "#################################################"    << std::endl;
    std::cout << "Projection matched:      " << projectionMatched   << std::endl;
    std::cout << "PFP hit completeness:    " << completeness_hits       << std::endl;
    std::cout << "PFP hit purity:          " << purity_hits             << std::endl;
    std::cout << "PFP energy completeness: " << completeness_energy     << std::endl;
    std::cout << "PFP energy purity:       " << purity_energy           << std::endl;
    std::cout << "#################################################"    << std::endl;
  }

  // need to add position difference

  return;
}

void ana::ShowerValidation::endJob() {

  std::cout << "Number of events ran over: " <<  numevents  << std::endl;
  std::cout << "Number of initial MC Showers: " <<  numshowers  << std::endl;
  std::cout << "Number of reco showers: " << numrecoshowers  << std::endl;
  std::cout << "Number of reco showers analysed: " << numrecoshowersana  << std::endl;

  if(numrecoshowersana == 0){
    std::ofstream myfile;
    myfile.open ("histcomp.txt");
    myfile << "0";
    myfile.close();

  }
}


DEFINE_ART_MODULE(ana::ShowerValidation)
