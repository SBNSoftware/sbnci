#include "TrackUtils.h"


//Function retunr a vector of event tracks' trackid.
std::vector<int> TrackUtils::EventMCTrkIdsFun(art::Event const& evt,
                                              std::map<int, const simb::MCParticle*> TrueParticlesMap,
                                              std::vector<unsigned int> ParTrkPdgCode){

  //###### Declaration and Initiate variables #######
  std::vector<int> EventTracksIDsVec;
  //
  for(auto const& mapiter : TrueParticlesMap){
    //
    //Declaration and Initiate variables
    const simb::MCParticle *particle = mapiter.second;
    // looking at the PdgCode of the particles
    unsigned int particlePdg = std::abs(particle->PdgCode());
    // Make sure that we are looking at track particles and daughter
    if(std::find(ParTrkPdgCode.begin(), ParTrkPdgCode.end(), particlePdg) == ParTrkPdgCode.end()){continue;}
    //Add to the mother chain.
    EventTracksIDsVec.push_back(std::abs(particle->TrackId()));
  }//close the loop over true mc particles
  return EventTracksIDsVec;
}//close MotherTrackIdsTodauIdsFun()

//Function return a vector of true hits of a reco particle
std::vector<art::Ptr<recob::Hit> > TrackUtils::TrueHitsOfParticle(detinfo::DetectorClocksData const& clockData,
                                                                  int TrackTrkIds,
                                                                  std::vector<art::Ptr<recob::Hit> > RHits){

  //###### Declaration and Initiate variables #######
  std::vector<art::Ptr<recob::Hit> >  hitsOfTrkIds;
  for(auto const& hit : RHits){
    int TrkIdOfHit = RecoUtils::TrueParticleID(clockData, hit);
    if(TrkIdOfHit == TrackTrkIds){
        hitsOfTrkIds.push_back(hit);
    }//colse condition statement
  }//close looping over all reco hit 
  return hitsOfTrkIds;
}//Close the function that return a vector of ture hits of a reco particle

//Function remove non contained particles
void TrackUtils::RemoveNoneContainedParticles(std::vector<int>&  EventMCTracks,
                                              std::map<int,const simb::MCParticle*>& trueParticlesMap,
                                              std::map<int,float>& MCTrackEnergyMap){

  //###### Declaration and Initiate variables #######
  std::vector<int>::iterator lastMCtrack = EventMCTracks.end();
  for(std::vector<int>::iterator iter = EventMCTracks.begin(); iter != lastMCtrack; ++iter){
    //Copy the iterator value
    int mctrk = (*iter);
    //Loop over the daughter
    float deposited_energy = 0.f;
    float sim_energy       = (trueParticlesMap[mctrk]->E())*1000;
    deposited_energy = MCTrackEnergyMap[mctrk];
    //If the energy of track is over 90% seen on the wires in truth, It is contained.
    if((deposited_energy/sim_energy) < 0.9){
      EventMCTracks.erase(iter);
    }  
  }
  return;
}//Close the function that remove non contained particles

//Function cut on the Track mothers with energy
void TrackUtils::CutTrkMothersByE(std::vector<int>&  EventMCTracks,
                                  std::map<int,const simb::MCParticle*>& trueParticlesMap,
                                  float& EnergyCut){   
  //###### Declaration and Initiate variables #######
  std::vector<int>::iterator lastMCtrack = EventMCTracks.end();
  //Cut the true tracks and make sure they are a tracks.
  for(std::vector<int>::iterator itera = EventMCTracks.begin(); itera != lastMCtrack; ++itera){
    //Copy the iterator value
    int mctrk = (*itera);
    const simb::MCParticle* TrkMotherParticle = trueParticlesMap[mctrk];
    if(TrkMotherParticle->E() < EnergyCut){
      EventMCTracks.erase(itera);
    }
  }//close the loop over mothers
  return;
}//close the function

//Function for cut on track via length 
void TrackUtils::CalculateTrackLengthCut(std::vector<int>&  EventMCTracks,
                                         std::map<int,const simb::MCParticle*>& trueParticlesMap,
                                         float& LengthCut){
  //###### Declaration and Initiate variables #######
  std::vector<int>::iterator lastMCtrack = EventMCTracks.end();
  //Cut the true tracks and make sure they are a tracks.
  for(std::vector<int>::iterator iterat = EventMCTracks.begin(); iterat != lastMCtrack; ++iterat){
  //
    float TrackLength = 0.f;
    bool startpc = false;
    bool endtpc  = false;
    //Copy the iterator value
    int mctrks = (*iterat);
    //
    const simb::MCParticle* TrksMCParticle = trueParticlesMap[mctrks];
    //Get the number of Traj points to loop over
    unsigned int TrajPoints = TrksMCParticle->NumberTrajectoryPoints();
    // Select first traj point where the track loses energy, last be default
    TLorentzVector PositionTrajStart =  TrksMCParticle->Position(0);
    TLorentzVector PositionTrajEnd   =  TrksMCParticle->Position(TrajPoints-1);
    //Check if the start and end position of track is inside the tpc
    startpc = RecoUtils::IsInsideTPC((PositionTrajStart.Vect()), 0);
    endtpc  = RecoUtils::IsInsideTPC((PositionTrajEnd.Vect()), 0);
    //Get the track length.
    TrackLength = (PositionTrajStart-PositionTrajEnd).Vect().Mag();
    if((!startpc) || (!endtpc) || (TrackLength < LengthCut)){
    //
      if(EventMCTracks.size() == 1){EventMCTracks.clear();}
      else {EventMCTracks.erase(iterat);}
    }
  }//close the loop over tracks mothers
  return;
}//close the function

//Function for finding the track Best plane
int TrackUtils::FindingTrackBeastPlane(detinfo::DetectorClocksData const& clockData,
                                       const std::vector<art::Ptr<recob::Hit> >& trkhits,
                                       unsigned int fUseBestPlane){
  //###### Declaration and Initiate variables ######
  unsigned int BestPlaneId = 999;
  unsigned int bestPlaneNumHits = 0;
  float BestPlaneEnergy = -999;
  std::map<int,float> PlaneEnergies; //A map of planeId as a key and energy deposited in each plane as value
  std::map<geo::PlaneID::PlaneID_t, unsigned int> planeHits;
  art::ServiceHandle<cheat::BackTrackerService> bt_serv;
  for(auto const& hit : trkhits){
    geo::WireID wireid = hit->WireID();
    unsigned int PlaneID = wireid.Plane;
    //planeHits[PlaneID].push_back(hit);
    ++planeHits[PlaneID];
    std::vector<sim::TrackIDE> trackIDEs = bt_serv->HitToTrackIDEs(clockData, hit);
    for(unsigned int idIt = 0; idIt < trackIDEs.size(); ++idIt){
    
      PlaneEnergies[PlaneID] += trackIDEs.at(idIt).energy;
    }//close the loop over trackIDEs
  }//cose the loop over the hits
  //Get the best plane for this track
  if(fUseBestPlane == 3){
    for(auto const& [plane, Nohits] : planeHits){
      unsigned int planeNumHits = Nohits;
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
  }else{
    BestPlaneId = fUseBestPlane;
  }

  return BestPlaneId;
}//Close the function for finding the track best plane


//Function to finds the most probable track id associated with the set of hits from there true energy depositions.
//The pair returns the energy as well.
std::pair<int,double> TrackUtils::pairTrackTrackIdToEnergyDep(detinfo::DetectorClocksData const& clockData,
                                                              std::vector<int>& EventMCTracks,
                                                              const std::vector<art::Ptr<recob::Hit> >& hits,
                                                              int planeid){
  //###### Declaration and Initiate variables #######
  std::map<int, float> mapTrackIdToTrkHitEnergy; //Map of trackid and the energy deposited per the trackid for 3 planes
  std::map<int, float> mapTrackIdToEDep;//Map of trackid and the energy deposited per the trackid for best track planes
  art::ServiceHandle<cheat::BackTrackerService> bt_serv;
  //Looping over track hits
  for(auto const& hit : hits){
    geo::WireID wireid = hit->WireID();
    int trkPlaneId = wireid.Plane;
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
  for(auto const& mctrks : EventMCTracks){
    mapMotherIdToEnergyDep[mctrks] += mapTrackIdToTrkHitEnergy[mctrks];
    mapMotherIdToEDep[mctrks] += mapTrackIdToEDep[mctrks];
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


