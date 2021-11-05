/*   //         NoiseTest_module.cc        \\
    //   Andrew Scarff (Uni of Sheffield)   \\
   //========================================\\
  // Module to find the FFT from noise events \\
 //  created using the MicroBooNE noise method \\
//   to test the implementation into sbndcode.  \\*/



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
//#include "art/Framework/Services/Optional/TFileService.h"
//#include "art/Framework/Services/Optional/TFileDirectory.h"
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
  class NoiseTest;
}

class ana::NoiseTest : public art::EDAnalyzer {
public:

  NoiseTest(const fhicl::ParameterSet& pset);

  //Larsoft specific funcitons 
  void analyze(const art::Event& evt);
  void beginJob();
  void endJob();

  //Class Functions

private:

  std::string fNoiseModuleLabel;
  art::ServiceHandle<art::TFileService> tfs;
  
  //  std::vector<float> rms;
  float rms[11276];
  float rmsT0P0[1986] = { 0 };
  float rmsT0P1[1986] = { 0 };
  float rmsT0P2[1666] = { 0 };
  float rmsT1P0[1986] = { 0 };
  float rmsT1P1[1986] = { 0 };
  float rmsT1P2[1666] = { 0 };
  //float rmsP0Wrong[1986] = { 0 };
  //float rmsP1Wrong[1986] = { 0 };
  //float rmsP2Wrong[1666] = { 0 };
  float wireT0P0[1986] = { 0 };
  float wireT0P1[1986] = { 0 };
  float wireT0P2[1666] = { 0 };
  float wireT1P0[1986] = { 0 };
  float wireT1P1[1986] = { 0 };
  float wireT1P2[1666] = { 0 };
  float wirelengthT0P0[1986] = { 0 };
  float wirelengthT0P1[1986] = { 0 };
  float wirelengthT0P2[1666] = { 0 };
  float wirelengthT1P0[1986] = { 0 };
  float wirelengthT1P1[1986] = { 0 };
  float wirelengthT1P2[1666] = { 0 };

  int nevents = 0;
  int t0p0iter, t0p1iter, t0p2iter;
  int t1p0iter, t1p1iter, t1p2iter;
  Int_t N;

  Double_t *fftoutputT0P0 = new Double_t[3000];
  Double_t *fftoutputT0P1 = new Double_t[3000];
  Double_t *fftoutputT0P2 = new Double_t[3000];
  Double_t *fftoutputT1P0 = new Double_t[3000];
  Double_t *fftoutputT1P1 = new Double_t[3000];
  Double_t *fftoutputT1P2 = new Double_t[3000];

  
  //TH1F* FFTHist;
  
  //fcl Parameters - usually have start with f
  /*  
  std::string fHitsModuleLabel;
  std::string fTrackModuleLabel;
  int         fParameter;
  TH1F*       hist = new TH1F("hist", "hist", 1000,0,50);
  TCanvas*    c1 = new TCanvas("c1","c1",800,1000);
  TFile*      myfile = new TFile("test.root","RECREATE");

  // This is just a class you can initialise anything
  int numtracks  = 0;
  */
};

ana::NoiseTest::NoiseTest(const fhicl::ParameterSet& pset) : EDAnalyzer(pset) {
  //Initialise set the fcl values 
  fNoiseModuleLabel = pset.get<std::string>("NoiseModuleLabel");
}

void ana::NoiseTest::beginJob() {
  std::cout << "Beginning job!" << std::endl;
}

void ana::NoiseTest::analyze(const art::Event& evt) {

  std::cout << "Analysing Event " << evt.event() << std::endl;
  nevents++;
  // Read in the digit List object(s).                                                                                                 
  art::Handle< std::vector<raw::RawDigit> > digitVecHandle;
  std::vector<art::Ptr<raw::RawDigit> > rawdigits; 

  if(evt.getByLabel(fNoiseModuleLabel, digitVecHandle)){
    art::fill_ptr_vector(rawdigits, digitVecHandle); 
  }
  std::vector<short> rawwaveform; 
  //  Int_t nWires = rawdigits.size();
  //std::cout << rawdigits.size() << std::endl;
  
  t0p0iter = 0;
  t0p1iter = 0;
  t0p2iter = 0;
  t1p0iter = 0;
  t1p1iter = 0;
  t1p2iter = 0;

  

  //Loop over all the wires (>11,000)
  for(size_t rdIter = 0; rdIter < rawdigits.size(); ++rdIter){
  
      art::Ptr<raw::RawDigit> digitVec = rawdigits[rdIter];
      rawwaveform.resize(digitVec->Samples());
      raw::Uncompress(digitVec->ADCs(), rawwaveform, digitVec->Compression());
      
      
      Int_t channel = digitVec->Channel(); //This = rdIter, keep as channel for ease of reading.

      ///Find wire length
      art::ServiceHandle<geo::Geometry> geom;
      std::vector<geo::WireID> wireIDs = geom->ChannelToWire(channel);
      geo::WireGeo const& wire = geom->Wire(wireIDs.front());      
      double wirelength = wire.Length();
      
      ///Get wire and plane numbers.
      unsigned int wireID = wireIDs.front().Wire;
      unsigned int planeID = wireIDs.front().Plane;
      unsigned int tpcID = wireIDs.front().TPC;
          
      N = rawwaveform.size();//3000 ticks // NUMBER OF ENTRIES IN WAVEFORM
      
      //std::string rawNameStr = "RawWaveform_" + std::to_string(evt.event());
      //const char* rawName = rawNameStr.c_str();
      
      std::string RMSNameStr = "RMSHist_" + std::to_string(evt.event());
      const char* RMSName = RMSNameStr.c_str();
      
      //TH1F* FFTHist = tfs->make<TH1F>(name,"",N/2,0.0,1.0);  // In MHz (Nyquist frequency is 1 MHz for 2 MHz digitization)                               
      //TH1F* RawWaveform = tfs->make<TH1F>(rawName,"",N,0.0,N); 
      TH1F* RMSHist = tfs->make<TH1F>(RMSName,"",N,0.0,N); 

      Double_t *fftinput = new Double_t[N];

      
      for(int j = 0; j<N; j++)
	{
	  fftinput[j] = rawwaveform[j];
	  //RawWaveform->SetBinContent(j+1,fftinput[j]);
	  RMSHist->Fill(rawwaveform[j]);	  
	}   
      float temprms =  RMSHist->GetRMS();
      
      TVirtualFFT *tempfft = TVirtualFFT::FFT(1,&N,"R2C");
      tempfft->SetPoints(fftinput);
      tempfft->Transform();
      

      
      if(planeID==0){
	if(tpcID==0){
      	  rmsT0P0[t0p0iter] += temprms;
	  wireT0P0[t0p0iter] = wireID;
	  wirelengthT0P0[t0p0iter] = wirelength;

	  ///Calculating the FFT.
	  Double_t re, im;
	  for(int k = 0; k < N; k++)
	    {
	      tempfft->GetPointComplex(k,re,im);
	      fftoutputT0P0[k] += sqrt(pow(re,2)+pow(im,2));
	      //	      std::cout<<"k: " << k << ", fftoutputT0P0[k]: " << fftoutputT0P0[k] << ", fftinput[k]: " << fftinput[k] << ", rawwaveform[k]: " << rawwaveform[k] << std::endl;
      
	    }
	  fftoutputT0P0[0] = 0.0;
      
	  t0p0iter++;
	}
	
	if(tpcID==1){
      	  rmsT1P0[t1p0iter] += temprms;
	  wireT1P0[t1p0iter] = wireID;
	  wirelengthT1P0[t1p0iter] = wirelength;

	  ///Calculating the FFT.
	  Double_t re, im;
	  for(int k = 0; k < N; k++)
	    {
	      tempfft->GetPointComplex(k,re,im);
	      fftoutputT1P0[k] += sqrt(pow(re,2)+pow(im,2));
	    }
	  fftoutputT1P0[0] = 0.0;
	  t1p0iter++;
	}
      }
	
      if(planeID==1){
	if(tpcID==0){
	  rmsT0P1[t0p1iter] += temprms;
	  wireT0P1[t0p1iter] = wireID;
	  wirelengthT0P1[t0p1iter] = wirelength;

	  ///Calculating the FFT.
	  Double_t re, im;
	  for(int k = 0; k < N; k++)
	    {
	      tempfft->GetPointComplex(k,re,im);
	      fftoutputT0P1[k] += sqrt(pow(re,2)+pow(im,2));
	    }
	  fftoutputT0P1[0] = 0.0;

	  t0p1iter++;
	}
	if(tpcID==1){
	  rmsT1P1[t1p1iter] += temprms;
	  wireT1P1[t1p1iter] = wireID;
	  wirelengthT1P1[t1p1iter] = wirelength;

	  ///Calculating the FFT.
	  Double_t re, im;
	  for(int k = 0; k < N; k++)
	    {
	      tempfft->GetPointComplex(k,re,im);
	      fftoutputT1P1[k] += sqrt(pow(re,2)+pow(im,2));
	    }
	  fftoutputT1P1[0] = 0.0;

	  t1p1iter++;
	}
      }
      if(planeID==2){
	if(tpcID==0){
	  rmsT0P2[t0p2iter] += temprms;
	  wireT0P2[t0p2iter] = wireID;
	  wirelengthT0P2[t0p2iter] = wirelength;

	  ///Calculating the FFT.
	  Double_t re, im;
	  for(int k = 0; k < N; k++)
	    {
	      tempfft->GetPointComplex(k,re,im);
	      fftoutputT0P2[k] += sqrt(pow(re,2)+pow(im,2));
	    }
	  fftoutputT0P2[0] = 0.0;

	  t0p2iter++;
	}
	if(tpcID==1){
	  rmsT1P2[t1p2iter] += temprms;
	  wireT1P2[t1p2iter] = wireID;
	  wirelengthT1P2[t1p2iter] = wirelength;

	  ///Calculating the FFT.
	  Double_t re, im;
	  for(int k = 0; k < N; k++)
	    {
	      tempfft->GetPointComplex(k,re,im);
	      fftoutputT1P2[k] += sqrt(pow(re,2)+pow(im,2));
	    }
	  fftoutputT1P2[0] = 0.0;

	  t1p2iter++;
	}
      }
      
      
      //delete RawWaveform;
      delete RMSHist;
      delete[] fftinput;    

  }//end wire loop
  
}


void ana::NoiseTest::endJob() {
  

  // RMS Plots

  float ADCtoElectron = (6241000/10696) * (4.7/14);; //New corrected ADCTicksPerPCAtLowestASICGainSetting  

  // Plane 0 - TPC 0

  int rmsT0P0Len = sizeof(rmsT0P0)/sizeof(*rmsT0P0);
  TH1F* RMSvsWireT0P0 = tfs->make<TH1F>("RMS vs Wire, TPC 0, Plane 0","",rmsT0P0Len,0.0,rmsT0P0Len);   
  TH1F* WirelengthvsWireT0P0 = tfs->make<TH1F>("Wirelength vs Wire, TPC 0, Plane 0","",rmsT0P0Len,0.0,rmsT0P0Len);   

  for(int i=0; i<rmsT0P0Len-1; ++i){
    rmsT0P0[i] /= nevents;
    rmsT0P0[i] *= ADCtoElectron;
    RMSvsWireT0P0->SetBinContent(wireT0P0[i],rmsT0P0[i]);
    WirelengthvsWireT0P0->SetBinContent(wireT0P0[i],wirelengthT0P0[i]);
  }
  
  RMSvsWireT0P0->GetXaxis()->SetTitle("Wire Number");
  RMSvsWireT0P0->GetYaxis()->SetTitle("ENC (Electrons)");
  RMSvsWireT0P0->SetTitle("RMS vs Wire, TPC 0, Plane 0");
  RMSvsWireT0P0->SetStats(0);
  RMSvsWireT0P0->SetLineColor(kRed);  

/*
  TLegend* legend = new TLegend(0.1,0.1,0.48,0.3);
  legend->AddEntry(RMSvsWireP0,"ENC calculated using correct conversion","l");
  legend->AddEntry(RMSvsWireP0Wrong,"ENC using current conversion","l");
  legend->Draw();
  */

  WirelengthvsWireT0P0->GetXaxis()->SetTitle("Wire Number");
  WirelengthvsWireT0P0->GetYaxis()->SetTitle("Wire Length (cm)");
  WirelengthvsWireT0P0->SetTitle("Wire Length vs Wire, TPC 0, Plane 0");
  WirelengthvsWireT0P0->SetStats(0);

  //Plane 0 - TPC 1

  int rmsT1P0Len = sizeof(rmsT1P0)/sizeof(*rmsT1P0);
  TH1F* RMSvsWireT1P0 = tfs->make<TH1F>("RMS vs Wire, TPC 1, Plane 0","",rmsT1P0Len,0.0,rmsT1P0Len);   
  TH1F* WirelengthvsWireT1P0 = tfs->make<TH1F>("Wirelength vs Wire, TPC 1, Plane 0","",rmsT1P0Len,0.0,rmsT1P0Len);   

  for(int i=0; i<rmsT1P0Len-1; ++i){
    rmsT1P0[i] /= nevents;
    rmsT1P0[i] *= ADCtoElectron;
    RMSvsWireT1P0->SetBinContent(wireT1P0[i],rmsT1P0[i]);
    WirelengthvsWireT1P0->SetBinContent(wireT1P0[i],wirelengthT1P0[i]);
  }
  
  RMSvsWireT1P0->GetXaxis()->SetTitle("Wire Number");
  RMSvsWireT1P0->GetYaxis()->SetTitle("ENC (Electrons)");
  RMSvsWireT1P0->SetTitle("RMS vs Wire, TPC 1, Plane 0");
  RMSvsWireT1P0->SetStats(0);
  RMSvsWireT1P0->SetLineColor(kRed);
  
  WirelengthvsWireT1P0->GetXaxis()->SetTitle("Wire Number");
  WirelengthvsWireT1P0->GetYaxis()->SetTitle("Wire Length (cm)");
  WirelengthvsWireT1P0->SetTitle("Wire Length vs Wire, TPC 1, Plane 0");
  WirelengthvsWireT1P0->SetStats(0);


  // Plane 1 - TPC 0

  int rmsT0P1Len = sizeof(rmsT0P1)/sizeof(*rmsT0P1);
  TH1F* RMSvsWireT0P1 = tfs->make<TH1F>("RMS vs Wire, TPC 0, Plane 1","",rmsT0P1Len,0.0,rmsT0P1Len);   
  TH1F* WirelengthvsWireT0P1 = tfs->make<TH1F>("Wirelength vs Wire, TPC 0, Plane 1","",rmsT0P1Len,0.0,rmsT0P1Len);   

  for(int i=0; i<rmsT0P1Len-1; ++i){
    rmsT0P1[i] /= nevents;
    rmsT0P1[i] *= ADCtoElectron;
    RMSvsWireT0P1->SetBinContent(wireT0P1[i],rmsT0P1[i]);
    WirelengthvsWireT0P1->SetBinContent(wireT0P1[i],wirelengthT0P1[i]);
  }
  RMSvsWireT0P1->GetXaxis()->SetTitle("Wire Number");
  RMSvsWireT0P1->GetYaxis()->SetTitle("ENC (Electrons)");
  RMSvsWireT0P1->SetTitle("RMS vs Wire, TPC 0, Plane 1");
  RMSvsWireT0P1->SetStats(0);
  RMSvsWireT0P1->SetLineColor(kRed);

  WirelengthvsWireT0P1->GetXaxis()->SetTitle("Wire Number");
  WirelengthvsWireT0P1->GetYaxis()->SetTitle("Wire Length (cm)");
  WirelengthvsWireT0P1->SetTitle("Wire Length vs Wire, TPC 0, Plane 1");
  WirelengthvsWireT0P1->SetStats(0);

  // Plane 1 - TPC 1

  int rmsT1P1Len = sizeof(rmsT1P1)/sizeof(*rmsT1P1);
  TH1F* RMSvsWireT1P1 = tfs->make<TH1F>("RMS vs Wire, TPC 1, Plane 1","",rmsT1P1Len,0.0,rmsT1P1Len);   
  TH1F* WirelengthvsWireT1P1 = tfs->make<TH1F>("Wirelength vs Wire, TPC 1, Plane 1","",rmsT1P1Len,0.0,rmsT1P1Len);   

  for(int i=0; i<rmsT1P1Len-1; ++i){
    rmsT1P1[i] /= nevents;
    rmsT1P1[i] *= ADCtoElectron;
    RMSvsWireT1P1->SetBinContent(wireT1P1[i],rmsT1P1[i]);
    WirelengthvsWireT1P1->SetBinContent(wireT1P1[i],wirelengthT1P1[i]);
  }

  RMSvsWireT1P1->GetXaxis()->SetTitle("Wire Number");
  RMSvsWireT1P1->GetYaxis()->SetTitle("ENC (Electrons)");
  RMSvsWireT1P1->SetTitle("RMS vs Wire, TPC 1, Plane 1");
  RMSvsWireT1P1->SetStats(0);
  RMSvsWireT1P1->SetLineColor(kRed);

  WirelengthvsWireT1P1->GetXaxis()->SetTitle("Wire Number");
  WirelengthvsWireT1P1->GetYaxis()->SetTitle("Wire Length (cm)");
  WirelengthvsWireT1P1->SetTitle("Wire Length vs Wire, TPC 1, Plane 1");
  WirelengthvsWireT1P1->SetStats(0);


  //Plane 2 - TPC 0

  int rmsT0P2Len = sizeof(rmsT0P2)/sizeof(*rmsT0P2);
  TH1F* RMSvsWireT0P2 = tfs->make<TH1F>("RMS vs Wire, TPC 0, Plane 2","",rmsT0P2Len,0.0,rmsT0P2Len); 
  TH1F* WirelengthvsWireT0P2 = tfs->make<TH1F>("Wirelength vs Wire, TPC 0, Plane 2","",rmsT0P2Len,0.0,rmsT0P2Len);   

  for(int i=0; i<rmsT0P2Len-1; ++i){
    rmsT0P2[i] /= nevents;
    rmsT0P2[i] *= ADCtoElectron;
    RMSvsWireT0P2->SetBinContent(wireT0P2[i],rmsT0P2[i]);
    WirelengthvsWireT0P2->SetBinContent(wireT0P2[i],wirelengthT0P2[i]);
  }

  RMSvsWireT0P2->GetXaxis()->SetTitle("Wire Number");
  RMSvsWireT0P2->GetYaxis()->SetTitle("ENC (Electrons)");
  RMSvsWireT0P2->SetTitle("RMS vs Wire, TPC 0, Plane 2");
  RMSvsWireT0P2->SetStats(0);
  RMSvsWireT0P2->SetLineColor(kRed);

  WirelengthvsWireT0P2->GetXaxis()->SetTitle("Wire Number");
  WirelengthvsWireT0P2->GetYaxis()->SetTitle("Wire Length (cm)");
  WirelengthvsWireT0P2->SetTitle("Wire Length vs Wire, TPC 0, Plane 2");
  WirelengthvsWireT0P2->SetStats(0);

  //Plane 2 - TPC 1

  int rmsT1P2Len = sizeof(rmsT1P2)/sizeof(*rmsT1P2);
  TH1F* RMSvsWireT1P2 = tfs->make<TH1F>("RMS vs Wire, TPC 1, Plane 2","",rmsT1P2Len,0.0,rmsT1P2Len); 
  TH1F* WirelengthvsWireT1P2 = tfs->make<TH1F>("Wirelength vs Wire, TPC 1, Plane 2","",rmsT1P2Len,0.0,rmsT1P2Len);   

  for(int i=0; i<rmsT1P2Len-1; ++i){
    rmsT1P2[i] /= nevents;
    rmsT1P2[i] *= ADCtoElectron;
    RMSvsWireT1P2->SetBinContent(wireT1P2[i],rmsT1P2[i]);
    WirelengthvsWireT1P2->SetBinContent(wireT1P2[i],wirelengthT1P2[i]);
  }

  RMSvsWireT1P2->GetXaxis()->SetTitle("Wire Number");
  RMSvsWireT1P2->GetYaxis()->SetTitle("ENC (Electrons)");
  RMSvsWireT1P2->SetTitle("RMS vs Wire, TPC 1, Plane 2");
  RMSvsWireT1P2->SetStats(0);
  RMSvsWireT1P2->SetLineColor(kRed);

  WirelengthvsWireT1P2->GetXaxis()->SetTitle("Wire Number");
  WirelengthvsWireT1P2->GetYaxis()->SetTitle("Wire Length (cm)");
  WirelengthvsWireT1P2->SetTitle("Wire Length vs Wire, TPC 1, Plane 2");
  WirelengthvsWireT1P2->SetStats(0);


  // FFT Calculation
  for(int k = 0; k < N; k++)
    {
      
      if(fftoutputT0P0[k]>1e+10 || fftoutputT0P0[k]< 0){fftoutputT0P0[k] = 0;}
      if(fftoutputT0P1[k]>1e+10 || fftoutputT0P1[k]< 0){fftoutputT0P1[k] = 0;}
      if(fftoutputT0P2[k]>1e+10 || fftoutputT0P2[k]< 0){fftoutputT0P2[k] = 0;}
      if(fftoutputT1P0[k]>1e+10 || fftoutputT1P0[k]< 0){fftoutputT1P0[k] = 0;}
      if(fftoutputT1P1[k]>1e+10 || fftoutputT1P1[k]< 0){fftoutputT1P1[k] = 0;}
      if(fftoutputT1P2[k]>1e+10 || fftoutputT1P2[k]< 0){fftoutputT1P2[k] = 0;}
      
      std::cout<<"BEFORE DIVISION - fftoutputT0P1[k]: " << fftoutputT0P1[k] << ", t0p1iter: " << t0p1iter << ", nevents: " << nevents << std::endl;

      fftoutputT0P0[k] /= t0p0iter*nevents;
      fftoutputT0P1[k] /= t0p1iter*nevents;
      fftoutputT0P2[k] /= t0p2iter*nevents;
      fftoutputT1P0[k] /= t1p0iter*nevents;
      fftoutputT1P1[k] /= t1p1iter*nevents;
      fftoutputT1P2[k] /= t1p2iter*nevents;
      
      std::cout<<"AFTER - fftoutputT0P0[k]: " << fftoutputT0P0[k] << ", t0p0iter: " << t0p0iter << ", nevents: " << nevents << std::endl;

      
    }
  
  TH1F* FFTHistT0P0 = tfs->make<TH1F>("FFT Hist, TPC 0, Plane 0","",N/2,0.0,1.0);  // In MHz (Nyquist frequency is 1 MHz for 2 MHz digitization)
  TH1F* FFTHistT0P1 = tfs->make<TH1F>("FFT Hist, TPC 0, Plane 1","",N/2,0.0,1.0);  // In MHz (Nyquist frequency is 1 MHz for 2 MHz digitization)
  TH1F* FFTHistT0P2 = tfs->make<TH1F>("FFT Hist, TPC 0, Plane 2","",N/2,0.0,1.0);  // In MHz (Nyquist frequency is 1 MHz for 2 MHz digitization)
  TH1F* FFTHistT1P0 = tfs->make<TH1F>("FFT Hist, TPC 1, Plane 0","",N/2,0.0,1.0);  // In MHz (Nyquist frequency is 1 MHz for 2 MHz digitization)
  TH1F* FFTHistT1P1 = tfs->make<TH1F>("FFT Hist, TPC 1, Plane 1","",N/2,0.0,1.0);  // In MHz (Nyquist frequency is 1 MHz for 2 MHz digitization)
  TH1F* FFTHistT1P2 = tfs->make<TH1F>("FFT Hist, TPC 1, Plane 2","",N/2,0.0,1.0);  // In MHz (Nyquist frequency is 1 MHz for 2 MHz digitization)
       
  ///Creating the FFT plot.
  
  for(int k = 0; k < N/2; k++)
    {
      FFTHistT0P0->SetBinContent(k+1,fftoutputT0P0[k]);
      FFTHistT0P1->SetBinContent(k+1,fftoutputT0P1[k]);
      FFTHistT0P2->SetBinContent(k+1,fftoutputT0P2[k]);
      FFTHistT1P0->SetBinContent(k+1,fftoutputT1P0[k]);
      FFTHistT1P1->SetBinContent(k+1,fftoutputT1P1[k]);
      FFTHistT1P2->SetBinContent(k+1,fftoutputT1P2[k]);
    }
  FFTHistT0P0->SetBinContent(1,0.0);
  FFTHistT0P1->SetBinContent(1,0.0);
  FFTHistT0P2->SetBinContent(1,0.0);
  FFTHistT1P0->SetBinContent(1,0.0);
  FFTHistT1P1->SetBinContent(1,0.0);
  FFTHistT1P2->SetBinContent(1,0.0);

  FFTHistT0P0->SetTitle("FFT, TPC 0, Plane 0");
  FFTHistT0P0->GetXaxis()->SetTitle("Frequency (MHz)");
  FFTHistT0P0->GetYaxis()->SetTitle("Magnitude (Arbitrary Units)");
  FFTHistT0P1->SetTitle("FFT, TPC 0, Plane 1");
  FFTHistT0P1->GetXaxis()->SetTitle("Frequency (MHz)");
  FFTHistT0P1->GetYaxis()->SetTitle("Magnitude (Arbitrary Units)");
  FFTHistT0P2->SetTitle("FFT, TPC 0, Plane 2");
  FFTHistT0P2->GetXaxis()->SetTitle("Frequency (MHz)");
  FFTHistT0P2->GetYaxis()->SetTitle("Magnitude (Arbitrary Units)");

  FFTHistT1P0->SetTitle("FFT, TPC 1, Plane 0");
  FFTHistT1P0->GetXaxis()->SetTitle("Frequency (MHz)");
  FFTHistT1P0->GetYaxis()->SetTitle("Magnitude (Arbitrary Units)");
  FFTHistT1P1->SetTitle("FFT, TPC 1, Plane 1");
  FFTHistT1P1->GetXaxis()->SetTitle("Frequency (MHz)");
  FFTHistT1P1->GetYaxis()->SetTitle("Magnitude (Arbitrary Units)");
  FFTHistT1P2->SetTitle("FFT, TPC 1, Plane 2");
  FFTHistT1P2->GetXaxis()->SetTitle("Frequency (MHz)");
  FFTHistT1P2->GetYaxis()->SetTitle("Magnitude (Arbitrary Units)");
  
  //delete[] fftoutputP0;
  //delete FFTHist;
  



  std::cout << "***************************************************" << std::endl;
  std::cout << "NoiseTest finished." << std::endl;
  std::cout << "***************************************************" << std::endl;

}


DEFINE_ART_MODULE(ana::NoiseTest)
