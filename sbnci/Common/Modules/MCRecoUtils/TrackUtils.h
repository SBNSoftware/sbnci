#ifndef TRACKUTILS_H_SEEN
#define TRACKUTILS_H_SEEN


/////////////////////////////////////////////////////
// TrackUtils.h 
//
// A few reco utilities like truth matching
// Functions for track specific.
// 
////////////////////////////////////////////////////


//Framework
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Handle.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/PtrVector.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "canvas/Persistency/Common/FindManyP.h"

//LArsoft
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "lardataobj/RecoBase/Hit.h"
#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/PlaneGeo.h"
#include "larcorealg/Geometry/WireGeo.h"
#include "sbnci/Common/Modules/MCRecoUtils/RecoUtils.h"
// c++
#include <vector>
#include <map>
#include <utility>
//
// ROOT
#include "TTree.h"


namespace TrackUtils{
  //Function retunr a vector of event tracks' trackids
  std::vector<int> EventMCTrkIdsFun(art::Event const& evt,
                                    std::map<int, const simb::MCParticle*> TrueParticlesMap,
                                    std::vector<unsigned int> ParTrkPdgCode);
  
  //Function return a vector of true hits of a reco particle
  std::vector<art::Ptr<recob::Hit> > TrueHitsOfParticle(detinfo::DetectorClocksData const& clockData,
                                                        int TrackTrkIds,
                                                        std::vector<art::Ptr<recob::Hit> > RecoParticleHits);
 
  //Function remove non contained particles
  void RemoveNoneContainedParticles(std::vector<int>& EventMCTracks,
                                    std::map<int,const simb::MCParticle*>& trueParticlesMap,
                                    std::map<int,float>& MCTrackEnergyMap);
  //
  //Function cut on the Track mothers with energy
  void CutTrkMothersByE(std::vector<int>& EventMCTracks,
                        std::map<int,const simb::MCParticle*>& trueParticlesMap,
                        float& EnergyCut);
  //

  //Function for cut on track via length 
  void CalculateTrackLengthCut(std::vector<int>&  EventMCTracks,
                               std::map<int,const simb::MCParticle*>& trueParticlesMap,
                               float& LengthCut);
  //Function for finding the track Best plane
  int FindingTrackBeastPlane(detinfo::DetectorClocksData const& clockData,
                             const std::vector<art::Ptr<recob::Hit> >& trkhits,
                             unsigned int fUseBestPlane);
  //Function to finds the most probable track id associated with the set of hits from 
  //there true energy depositions. The pair returns the energy as well.
  std::pair<int,double> pairTrackTrackIdToEnergyDep(detinfo::DetectorClocksData const& clockData,
                                                    std::vector<int>& EventMCTracks,
                                                    const std::vector<art::Ptr<recob::Hit> >& hits,
                                                    int planeid);
  
}

#endif
