///////////////////////////////////////////////////////////////////////////
// Class:        RecoEff                                                 //
// Module Type:  Analyser                                                //
// File:         RecoEff_module.cc                                       //
// Author:       Henry Lay h.lay@lancaster.ac.uk                         //
//                                                                       //
// Purpose:      The module saves truth information regarding the        //
//               primary products of neutrino interactions. It also      //
//               saves their status in reconstruction. This allows       //
//               for plotting of a series of plots of different          //
//               reconstruction metrics.                                 //
///////////////////////////////////////////////////////////////////////////

//Standard
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

//Sim Base
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"

//Reco Base
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/Hit.h"

//Tools
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/FindOneP.h"

//LArSoft
#include "larsim/Utils/TruthMatchUtils.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"

//Root
#include "art_root_io/TFileService.h"
#include "TTree.h"

//SBNDCODE
#include "sbndcode/Geometry/GeometryWrappers/TPCGeoAlg.h"

constexpr int def_int     = -999;
constexpr float def_float = -999.0f;

class RecoEff;

class RecoEff : public art::EDAnalyzer {
public:
  explicit RecoEff(fhicl::ParameterSet const& pset);

  RecoEff(RecoEff const&) = delete;
  RecoEff(RecoEff&&) = delete;
  RecoEff& operator=(RecoEff const&) = delete;
  RecoEff& operator=(RecoEff&&) = delete;

  void analyze(art::Event const &e) override;

private:

  void ClearMaps();
  void ResetData();
  void SetupHitsMap(art::Event const &e);
  void ReconstructionProcessor(art::Event const &e);
  void TruthProcessor(art::Event const &e);

  std::vector<art::Ptr<recob::PFParticle> > GetPrimaryPFPs(art::Event const &e);
  art::Ptr<recob::PFParticle> GetPFP(art::Event const &e, long unsigned int const &id);

  float Purity(std::vector< art::Ptr<recob::Hit> > const &objectHits, int const &trackID);
  float Completeness(std::vector< art::Ptr<recob::Hit> > const &objectHits, int const &trackID);
  float TrackLength(art::Ptr<recob::Track> const &track);
  float TrueTrackLength(art::Ptr<simb::MCParticle> const &particle);


  detinfo::DetectorClocksData clockData = art::ServiceHandle<detinfo::DetectorClocksService>()->DataForJob();

  sbnd::TPCGeoAlg fTPCGeo;
  
  std::string fNuGenModuleLabel, fLArGeantModuleLabel, fPFParticleModuleLabel,
    fTrackModuleLabel, fShowerModuleLabel, fHitsModuleLabel;

  float fXEdgeCut, fYEdgeCut, fZFrontCut, fZBackCut,
    fXEdgeCutShowers, fYEdgeCutShowers, fZFrontCutShowers, fZBackCutShowers,
    fCathodeXCut;

  float fXOutEdge, fYEdge, fZFrontEdge, fZBackEdge,
    fXOutEdgeShowers, fYEdgeShowers, fZFrontEdgeShowers, fZBackEdgeShowers,
    fXCathodeEdge;

  TTree *fParticleTree;

  int mc_trackID, mc_PDG, reco_nTracks, reco_nShowers;
  float mc_x0, mc_y0, mc_z0, mc_xEnd, mc_yEnd, mc_zEnd, mc_pX0,
    mc_pY0, mc_pZ0, mc_energy0, mc_momentum, mc_pXEnd, mc_pYEnd, mc_pZEnd,
    mc_energyEnd, mc_mass, mc_theta_xy, mc_theta_yz, mc_theta_xz, mc_length, 
    reco_shower_purity, reco_shower_completeness, reco_track_purity, 
    reco_track_completeness, reco_track_length, reco_shower_dEdx;
  bool reco_isReconstructed;

  std::map<int,int> n_showers_map, n_tracks_map;
  std::map<int,float> shower_comp_map, shower_pur_map, track_comp_map,
    track_pur_map, track_length_map, shower_dEdx_map;
  std::map<int,int> hits_map;

};


RecoEff::RecoEff(fhicl::ParameterSet const &pset)
  : EDAnalyzer{pset},
  fNuGenModuleLabel (pset.get<std::string>("NuGenModuleLabel")),
  fLArGeantModuleLabel (pset.get<std::string>("LArGeantModuleLabel")),
  fPFParticleModuleLabel (pset.get<std::string>("PFParticleModuleLabel")),
  fTrackModuleLabel (pset.get<std::string>("TrackModuleLabel")),
  fShowerModuleLabel (pset.get<std::string>("ShowerModuleLabel")),
  fHitsModuleLabel (pset.get<std::string>("HitsModuleLabel")),
  fXEdgeCut (pset.get<float>("XEdgeCut",15.f)),
  fYEdgeCut (pset.get<float>("YEdgeCut",15.f)),
  fZFrontCut (pset.get<float>("ZFrontCut",30.f)),
  fZBackCut (pset.get<float>("ZBackCut",65.f)),
  fXEdgeCutShowers (pset.get<float>("XEdgeCutShowers",25.f)),
  fYEdgeCutShowers (pset.get<float>("YEdgeCutShowers",25.f)),
  fZFrontCutShowers (pset.get<float>("ZFrontCutShowers",30.f)),
  fZBackCutShowers (pset.get<float>("ZBackCutShowers",50.f)),
  fCathodeXCut (pset.get<float>("CathodeXCut",1.5f))

  {
    art::ServiceHandle<art::TFileService> tfs;
    fParticleTree = tfs->make<TTree>("ParticleTree","Particle data TTree");

    fParticleTree->Branch("mc_trackID",&mc_trackID);
    fParticleTree->Branch("mc_PDG",&mc_PDG);
    fParticleTree->Branch("mc_x0",&mc_x0);
    fParticleTree->Branch("mc_y0",&mc_y0);
    fParticleTree->Branch("mc_z0",&mc_z0);
    fParticleTree->Branch("mc_xEnd",&mc_xEnd);
    fParticleTree->Branch("mc_yEnd",&mc_yEnd);
    fParticleTree->Branch("mc_zEnd",&mc_zEnd);
    fParticleTree->Branch("mc_pX0",&mc_pX0);
    fParticleTree->Branch("mc_pY0",&mc_pY0);
    fParticleTree->Branch("mc_pZ0",&mc_pZ0);
    fParticleTree->Branch("mc_energy0",&mc_energy0);
    fParticleTree->Branch("mc_momentum",&mc_momentum);
    fParticleTree->Branch("mc_pXEnd",&mc_pXEnd);
    fParticleTree->Branch("mc_pYEnd",&mc_pYEnd);
    fParticleTree->Branch("mc_pZEnd",&mc_pZEnd);
    fParticleTree->Branch("mc_energyEnd",&mc_energyEnd);
    fParticleTree->Branch("mc_mass",&mc_mass);
    fParticleTree->Branch("mc_theta_xy",&mc_theta_xy);
    fParticleTree->Branch("mc_theta_yz",&mc_theta_yz);
    fParticleTree->Branch("mc_theta_xz",&mc_theta_xz);
    fParticleTree->Branch("mc_length",&mc_length);
    fParticleTree->Branch("reco_isReconstructed",&reco_isReconstructed);
    fParticleTree->Branch("reco_nTracks",&reco_nTracks);
    fParticleTree->Branch("reco_nShowers",&reco_nShowers);
    fParticleTree->Branch("reco_track_purity",&reco_track_purity);
    fParticleTree->Branch("reco_track_completeness",&reco_track_completeness);
    fParticleTree->Branch("reco_shower_purity",&reco_shower_purity);
    fParticleTree->Branch("reco_shower_completeness",&reco_shower_completeness);
    fParticleTree->Branch("reco_track_length",&reco_track_length);
    fParticleTree->Branch("reco_shower_dEdx",&reco_shower_dEdx);


    fXOutEdge   = fTPCGeo.MaxX() - fXEdgeCut;
    fYEdge      = fTPCGeo.MaxY() - fYEdgeCut;
    fZFrontEdge = fTPCGeo.MinZ() + fZFrontCut;
    fZBackEdge  = fTPCGeo.MaxZ() - fZBackCut;

    fXOutEdgeShowers   = fTPCGeo.MaxX() - fXEdgeCutShowers;
    fYEdgeShowers      = fTPCGeo.MaxY() - fYEdgeCutShowers;
    fZFrontEdgeShowers = fTPCGeo.MinZ() + fZFrontCutShowers;
    fZBackEdgeShowers  = fTPCGeo.MaxZ() - fZBackCutShowers;

    fXCathodeEdge = 0 + fCathodeXCut;

  }

void RecoEff::ClearMaps()
{
  n_tracks_map.clear(); track_comp_map.clear(); track_pur_map.clear();
  n_showers_map.clear(); shower_comp_map.clear(); shower_pur_map.clear();
  track_length_map.clear(); shower_dEdx_map.clear();
  hits_map.clear();
}

void RecoEff::ResetData()
{
  mc_trackID = def_int; mc_PDG = def_int; 

  mc_x0 = def_float; mc_y0 = def_float; mc_z0 = def_float; mc_xEnd = def_float; mc_yEnd = def_float; 
  mc_zEnd = def_float; mc_pX0 = def_float; mc_pY0 = def_float; mc_pZ0 = def_float; mc_energy0 = def_float; 
  mc_momentum = def_float; mc_pXEnd = def_float; mc_pYEnd = def_float; mc_pZEnd = def_float; 
  mc_energyEnd = def_float; mc_mass = def_float; mc_theta_xy = def_float; mc_theta_yz = def_float;
  mc_theta_xz = def_float; mc_length = def_float; 

  reco_nTracks = def_int; reco_nShowers = def_int; 
  
  reco_shower_purity = def_float; reco_shower_completeness = def_float; reco_track_purity = def_float; 
  reco_track_completeness = def_float; reco_track_length = def_float; reco_shower_dEdx = def_float; 

  reco_isReconstructed = false;
}

void RecoEff::SetupHitsMap(art::Event const &e)
{
  art::Handle<std::vector<recob::Hit> > handleHits;
  e.getByLabel(fHitsModuleLabel,handleHits);

  for(unsigned hit_i = 0; hit_i < handleHits->size(); ++hit_i) {
    const art::Ptr<recob::Hit> hit(handleHits,hit_i);
    hits_map[TruthMatchUtils::TrueParticleID(clockData,hit,true)]++;
  }
}

void RecoEff::ReconstructionProcessor(art::Event const &e)
{
  art::Handle<std::vector<recob::PFParticle> > handlePFPs;
  art::Handle<std::vector<recob::Track> > handleTracks;
  art::Handle<std::vector<recob::Shower> > handleShowers;
  e.getByLabel(fPFParticleModuleLabel,handlePFPs);
  e.getByLabel(fTrackModuleLabel,handleTracks);
  e.getByLabel(fShowerModuleLabel,handleShowers);

  art::FindManyP<recob::Track> pfpTrackAssn(handlePFPs,e,fTrackModuleLabel);
  art::FindManyP<recob::Shower> pfpShowerAssn(handlePFPs,e,fShowerModuleLabel);
  art::FindManyP<recob::Hit> trackHitAssn(handleTracks,e,fTrackModuleLabel);
  art::FindManyP<recob::Hit> showerHitAssn(handleShowers,e,fShowerModuleLabel);

  const std::vector<art::Ptr<recob::PFParticle> > primaries = GetPrimaryPFPs(e);

  for(auto primary : primaries){
    const std::vector<long unsigned int> daughterIDs = primary->Daughters();

    for(auto id: daughterIDs){
      const art::Ptr<recob::PFParticle> daughter = GetPFP(e,id);
      if(daughter.isNull()) continue;

      if(daughter->PdgCode() == 13) {
	const std::vector<art::Ptr<recob::Track> > tracksVec = pfpTrackAssn.at(daughter.key());
	if(tracksVec.size() != 1) continue;
	const art::Ptr<recob::Track> track = tracksVec[0];

	float x = track->Start().X(), y = track->Start().Y(), z = track->Start().Z();

	if(TMath::Abs(x) > fXOutEdge || TMath::Abs(x) < fXCathodeEdge ||
	   TMath::Abs(y) > fYEdge || z < fZFrontEdge || z > fZBackEdge) continue;

	std::vector<art::Ptr<recob::Hit> > trackHits = trackHitAssn.at(track.key());
	int trackID = TruthMatchUtils::TrueParticleIDFromTotalRecoHits(clockData,trackHits,true);
	float comp = Completeness(trackHits,trackID);
	float pur = Purity(trackHits,trackID);
	float length = TrackLength(track);

	if(n_tracks_map[trackID] == 0){
	  track_comp_map[trackID] = comp;
	  track_pur_map[trackID] = pur;
	  track_length_map[trackID] = length;
	}
	else if(comp > track_comp_map[trackID]) {
	  track_comp_map[trackID] = comp;
	  track_pur_map[trackID] = pur;
	  track_length_map[trackID] = length;
	}

	n_tracks_map[trackID]++;
      }
      else if(daughter->PdgCode() == 11) {
	const std::vector<art::Ptr<recob::Shower> > showersVec = pfpShowerAssn.at(daughter.key());
	if(showersVec.size() != 1) continue;
	const art::Ptr<recob::Shower> shower = showersVec[0];

	float x = shower->ShowerStart().X(), y = shower->ShowerStart().Y(), z = shower->ShowerStart().Z();

	if(TMath::Abs(x) > fXOutEdgeShowers || TMath::Abs(x) < fXCathodeEdge ||
	   TMath::Abs(y) > fYEdgeShowers || z < fZFrontEdgeShowers || z > fZBackEdgeShowers) continue;

	std::vector<art::Ptr<recob::Hit> > showerHits = showerHitAssn.at(shower.key());
	int trackID = TruthMatchUtils::TrueParticleIDFromTotalRecoHits(clockData,showerHits,true);
	float comp = Completeness(showerHits,trackID);
	float pur = Purity(showerHits,trackID);

	std::vector<double> dEdxVec = shower->dEdx();
	int best_plane = shower->best_plane();
	float dEdx = def_float;
	if(dEdxVec.size() != 0) dEdx = dEdxVec[best_plane];

	if(n_showers_map[trackID] == 0){
	  shower_comp_map[trackID] = comp;
	  shower_pur_map[trackID] = pur;
	  shower_dEdx_map[trackID] = dEdx;
	}
	else if(comp > shower_comp_map[trackID]) {
	  shower_comp_map[trackID] = comp;
	  shower_pur_map[trackID] = pur;
	  shower_dEdx_map[trackID] = dEdx;
	}

	n_showers_map[trackID]++;
      }
    }
  }
}

void RecoEff::TruthProcessor(art::Event const &e)
{
  art::Handle<std::vector<simb::MCTruth> > handleNeutrinos;
  e.getByLabel(fNuGenModuleLabel,handleNeutrinos);

  art::FindManyP<simb::MCParticle> nuParticleAssn(handleNeutrinos,e,fLArGeantModuleLabel);

  for(unsigned int nu_i = 0; nu_i < handleNeutrinos->size(); ++nu_i){
    const art::Ptr<simb::MCTruth> truthNeutrino(handleNeutrinos,nu_i);
    if(truthNeutrino.isNull()) continue;
    std::vector<art::Ptr<simb::MCParticle> > particles = nuParticleAssn.at(truthNeutrino.key());

    for(auto particle : particles){
      if(particle->Mother() != 0 || particle->StatusCode() != 1) continue;

      ResetData();

      mc_trackID = particle->TrackId();
      mc_PDG = particle->PdgCode();
      mc_x0 = particle->Vx();
      mc_y0 = particle->Vy();
      mc_z0 = particle->Vz();

      if(std::abs(mc_PDG) == 11 || std::abs(mc_PDG) == 22) {
	if(TMath::Abs(mc_x0) > fXOutEdgeShowers || TMath::Abs(mc_x0) < fXCathodeEdge ||
	   TMath::Abs(mc_y0) > fYEdgeShowers || mc_z0 < fZFrontEdgeShowers || mc_z0 > fZBackEdgeShowers) continue;
      }
      else if(TMath::Abs(mc_x0) > fXOutEdge || TMath::Abs(mc_x0) < fXCathodeEdge ||
	      TMath::Abs(mc_y0) > fYEdge || mc_z0 < fZFrontEdge || mc_z0 > fZBackEdge) continue;

      mc_xEnd = particle->EndX();
      mc_yEnd = particle->EndY();
      mc_zEnd = particle->EndZ();
      mc_pX0 = particle->Px();
      mc_pY0 = particle->Py();
      mc_pZ0 = particle->Pz();
      mc_energy0 = particle->E();
      mc_momentum = particle->P();
      mc_pXEnd = particle->EndPx();
      mc_pYEnd = particle->EndPy();
      mc_pZEnd = particle->EndPz();
      mc_energyEnd = particle->EndE();
      mc_mass = particle->Mass();
      mc_theta_xy = TMath::RadToDeg() * TMath::ATan(mc_pX0/mc_pY0);
      mc_theta_yz = TMath::RadToDeg() * TMath::ATan(mc_pY0/mc_pZ0);
      mc_theta_xz = TMath::RadToDeg() * TMath::ATan(mc_pX0/mc_pZ0);
      mc_length = TrueTrackLength(particle);

      if(n_tracks_map[mc_trackID] > 0 || n_showers_map[mc_trackID] > 0) reco_isReconstructed = true;
      reco_nTracks = n_tracks_map[mc_trackID];
      reco_nShowers = n_showers_map[mc_trackID];
      reco_shower_purity = shower_pur_map[mc_trackID];
      reco_shower_completeness = shower_comp_map[mc_trackID];
      reco_track_purity = track_pur_map[mc_trackID];
      reco_track_completeness = track_comp_map[mc_trackID];
      reco_track_length = track_length_map[mc_trackID];
      reco_shower_dEdx = shower_dEdx_map[mc_trackID];

      fParticleTree->Fill();
    }
  }
}

std::vector<art::Ptr<recob::PFParticle> > RecoEff::GetPrimaryPFPs(art::Event const &e)
{
  art::Handle<std::vector<recob::PFParticle> > handlePFPs;
  e.getByLabel(fPFParticleModuleLabel,handlePFPs);

  std::vector<art::Ptr<recob::PFParticle> > primaries;

  for(unsigned int pfp_i = 0; pfp_i < handlePFPs->size(); ++pfp_i) {
    const art::Ptr<recob::PFParticle> pfp(handlePFPs,pfp_i);
    if(pfp->IsPrimary()){
      if(std::abs(pfp->PdgCode()) == 12 || std::abs(pfp->PdgCode()) == 14) {
        primaries.push_back(pfp);
      }
    }
  }
  return primaries;
}

art::Ptr<recob::PFParticle> RecoEff::GetPFP(art::Event const &e, long unsigned int const &id)
{
  art::Handle<std::vector<recob::PFParticle> > handlePFPs;
  e.getByLabel(fPFParticleModuleLabel,handlePFPs);

  const art::Ptr<recob::PFParticle> nullReturn;

  for(unsigned int pfp_i = 0; pfp_i < handlePFPs->size(); ++pfp_i) {
    const art::Ptr<recob::PFParticle> pfp(handlePFPs,pfp_i);
    if(pfp->Self() == id) return pfp;
  }
  return nullReturn;
}

float RecoEff::Purity(std::vector< art::Ptr<recob::Hit> > const &objectHits, int const &trackID)
{
  std::map<int,int> objectHitsMap;

  for(unsigned int i = 0; i < objectHits.size(); ++i) {
    objectHitsMap[TruthMatchUtils::TrueParticleID(clockData,objectHits[i],true)]++;
  }
  return (objectHits.size() == 0) ? def_float : objectHitsMap[trackID]/static_cast<float>(objectHits.size());
}

float RecoEff::Completeness(std::vector< art::Ptr<recob::Hit> > const &objectHits, int const &trackID)
{
  std::map<int,int> objectHitsMap;

  for(unsigned int i = 0; i < objectHits.size(); ++i) {
    objectHitsMap[TruthMatchUtils::TrueParticleID(clockData,objectHits[i],true)]++;
  }
  return (hits_map[trackID] == 0) ? def_float : objectHitsMap[trackID]/static_cast<float>(hits_map[trackID]);
}

float RecoEff::TrackLength(art::Ptr<recob::Track> const &track)
{
  float length = 0;
  int nTrajPoints = track->NumberTrajectoryPoints();

  if(nTrajPoints < 2) return length;

  for(int point = 1; point < nTrajPoints; ++point) {
    TVector3 l = track->LocationAtPoint<TVector3>(point);
    if(!track->HasValidPoint(point)) break;

    TVector3 diff = track->LocationAtPoint<TVector3>(point) - track->LocationAtPoint<TVector3>(point-1);
    length += diff.Mag2();
  }
  return length;
}

float RecoEff::TrueTrackLength(art::Ptr<simb::MCParticle> const &particle)
{
  float length = 0;
  unsigned int nTrajPoints = particle->NumberTrajectoryPoints();

  if(nTrajPoints < 2) return length;

  for(unsigned int point = 1; point < nTrajPoints; ++point) {
    TVector3 l = particle->Position(point).Vect();
    if(l.X() > 200 || l.X() < -200 || l.Y() > 200 || l.Y() < -200 ||
       l.Z() > 500 || l.Z() < 0) break;

    TVector3 diff = particle->Position(point).Vect() - particle->Position(point-1).Vect();
    length += diff.Mag2();
  }
  return length;
}


void RecoEff::analyze(art::Event const &e)
{
  clockData = art::ServiceHandle<detinfo::DetectorClocksService>()->DataFor(e);

  ClearMaps();
  SetupHitsMap(e);
  ReconstructionProcessor(e);
  TruthProcessor(e);
}

DEFINE_ART_MODULE(RecoEff)
