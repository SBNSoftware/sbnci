/*   //         PlotdEdx_module.cc         \\
    //   Andrew Scarff (Uni of Sheffield)   \\
   //========================================\\
  // Module to plot the dEdx                  \\*/



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

//Larsoft Includes 

#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCFlux.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "lardataobj/Simulation/AuxDetSimChannel.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "lardataobj/AnalysisBase/ParticleID.h"
#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/RawData/BeamInfo.h"
#include "larcoreobj/SummaryData/POTSummary.h"
#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/EndPoint2D.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "lardataobj/AnalysisBase/CosmicTag.h"
#include "sbndcode/RecoUtils/RecoUtils.h"
#include "lardataobj/RawData/raw.h"

//C++ Includes
#include <vector>
#include <string>
#include "TH1.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TVirtualFFT.h"
#include "TLegend.h"

namespace ana {
  class PlotdEdx;
}

class ana::PlotdEdx : public art::EDAnalyzer {
public:

  PlotdEdx(const fhicl::ParameterSet& pset);

  //Larsoft specific funcitons 
  void beginJob();
  void analyze(const art::Event& evt);
  void endJob();

  //Class Functions

private:

  std::string fPlotdEdxModuleLabel;
  std::string fTrackModuleLabel;
  std::string fCaloModuleLabel;
  art::ServiceHandle<art::TFileService> tfs;
  int nevents=0;
  
  TH1F* dEdxPlotP0;
  TH1F* dEdxPlotP1;
  TH1F* dEdxPlotP2;

  std::map<int,TH1F*> plotMap;

};

ana::PlotdEdx::PlotdEdx(const fhicl::ParameterSet& pset) 
  : EDAnalyzer(pset)
{
  //Initialise set the fcl values 
  fPlotdEdxModuleLabel = pset.get<std::string>("PlotdEdxModuleLabel");
  fTrackModuleLabel = pset.get<std::string>("TrackModuleLabel");  
  fCaloModuleLabel = pset.get<std::string>("CaloModuleLabel");

}

void ana::PlotdEdx::beginJob() {
  std::cout << "Beginning job!" << std::endl;
  dEdxPlotP0 = tfs->make<TH1F>("dEdx-P0","",100,0.0,6);
  dEdxPlotP1 = tfs->make<TH1F>("dEdx-P1","",100,0.0,6);
  dEdxPlotP2 = tfs->make<TH1F>("dEdx-P2","",100,0.0,6);
  //dEdxPlotT1P0 = tfs->make<TH1F>("dE/dx T1 P0","",100,0.0,6);
  //dEdxPlotT1P1 = tfs->make<TH1F>("dE/dx T1 P1","",100,0.0,6);
  //dEdxPlotT1P2 = tfs->make<TH1F>("dE/dx T1 P2","",100,0.0,6);

  dEdxPlotP0->SetTitle("dE/dx - Horizontal Muons - Plane 0");
  dEdxPlotP0->GetXaxis()->SetTitle("dE/dx (MeV/cm)");
  dEdxPlotP0->GetYaxis()->SetTitle("Counts");
  dEdxPlotP1->SetTitle("dE/dx - Horizontal Muons - Plane 1");
  dEdxPlotP1->GetXaxis()->SetTitle("dE/dx (MeV/cm)");
  dEdxPlotP1->GetYaxis()->SetTitle("Counts");
  dEdxPlotP2->SetTitle("dE/dx - Horizontal Muons - Plane 2");
  dEdxPlotP2->GetXaxis()->SetTitle("dE/dx (MeV/cm)");
  dEdxPlotP2->GetYaxis()->SetTitle("Counts");

  plotMap[0]=dEdxPlotP0;
  plotMap[1]=dEdxPlotP1;
  plotMap[2]=dEdxPlotP2;
  

}

void ana::PlotdEdx::analyze(const art::Event& evt) {

  std::cout << "Analysing Event " << evt.event() << std::endl;
  nevents++;
  
  art::Handle< std::vector<recob::Track> > trackListHandle;
  std::vector< art::Ptr<recob::Track> > trackList;
  
// Filling the tracklist from the tracklistHandle                                                          
  if (evt.getByLabel(fTrackModuleLabel,trackListHandle))
    {art::fill_ptr_vector(trackList, trackListHandle);}

  art::FindManyP<recob::Hit> findManyHits(trackListHandle, evt, fTrackModuleLabel);
  art::FindManyP<anab::Calorimetry> findManyCalo(trackListHandle, evt, fCaloModuleLabel);


  for(auto const& track: trackList){ //loop over tracks
    
    std::vector< art::Ptr<anab::Calorimetry> > trackCaloInfo = findManyCalo.at(track.key());
   
    for(auto const &caloinfoplane: trackCaloInfo){ //loop over planes
      const std::vector<float> dEdx = caloinfoplane->dEdx();
      int planeID = caloinfoplane->PlaneID().Plane;
      
      for(auto const& dedx : dEdx){
	//fill histo 
	plotMap[planeID]->Fill(dedx);
      }
      
    }
  }
  std::cout << "done mate" << std::endl;
  return;
}


void ana::PlotdEdx::endJob() {
  



  std::cout << "***************************************************" << std::endl;
  std::cout << "Analysing finished." << std::endl;
  std::cout << "***************************************************" << std::endl;

}


DEFINE_ART_MODULE(ana::PlotdEdx)
