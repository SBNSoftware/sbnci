#include "TCanvas.h"
#include "TStyle.h"
#include "TH1.h"
#include "TH2.h"
#include "TGaxis.h"
#include "TRandom.h"
#include "TFile.h"
#include "TLegend.h"
#include "TCollection.h"
#include <iostream>
#include "THStack.h"
#include "TAxis.h"
#include "TGraphErrors.h"
#include <sstream>
#include <vector>
#include "TString.h"

void SetDetectorRange(float &minx, float &miny, float &minz,
                      float &maxx, float &maxy, float &maxz,
                      int &binsx, int &binsy, int &binsz,
                      bool isND)
{

  float min_nd = {-200.,-200.,0.};
  float max_nd = { 200., 200.,500.};
  float min_fd = {-400.,-200.,-1000.};
  float max_fd = { 400., 150., 1000.};
  int bins_nd  = { 40, 40, 50};
  int bins_fd  = { 80, 35, 100};

  float min = min_nd;
  float max = max_nd;
  int bins  = bins_nd;

  if(!isND){
    min  = min_fd;
    max  = max_fd;
    bins = bins_fd;
  }

  minx = min[0];
  miny = min[1];
  minz = min[2];
  maxx = max[0];
  maxy = max[1];
  maxz = max[2];
  binsx = bins[0];
  binsy = bins[1];
  binsz = bins[2];

}

void PlotHists(TFile *inFile, TFile *outFile, TString histName, TString histTitle, TString histXaxis, TString histYaxis, int nBins, int xMin, int xMax, TString metricToPlot, TString dirName);

void sbnciplot_trackvalidation(TString inputFile)
{
  
  TFile *inFile = TFile::Open(inputFile.Data());
  TFile *outFile = new TFile("ci_validation_histos.root", "RECREATE");

  TString histName, histTitle, histXaxis, metricToPlot;
  int xMin, xMax;
  TString histYaxis = "Counts";
  int nBins = 100;

  TString dirName = "pandoraTrack";

  outFile->cd();
  gDirectory->mkdir(dirName);
  gSystem->Exec(Form("mkdir -v %s", dirName.Data()));

  float minX, minY, minZ;
  float maxX, maxY, maxZ;
  int nBinsX, nBinsY, nBinsZ;
  SetDetectorRange(minX, minY, minZ, maxX, maxY, maxZ, nBinsX, nBinsY, nBinsZ, true);

  PlotHists(inFile, outFile, "track_lenght",       "Reco Track Lenght",           "Reco Track Lenght (cm)",            histYaxis, 70,    0, 700,      "track_lenght",       dirName);
  PlotHists(inFile, outFile, "track_startdirphi",  "Reco Track Start Direction",  "Reco Track Start Direction (#Phi)", histYaxis, 70, -3.5, 3.5,      "track_startdirphi",  dirName);
  PlotHists(inFile, outFile, "track_startx",       "Reco Track Start Position X", "Reco Track Start Position X (cm)",  histYaxis, nBinsX, minX, maxX, "track_startx",       dirName);
  PlotHists(inFile, outFile, "track_starty",       "Reco Track Start Position Y", "Reco Track Start Position Y (cm)",  histYaxis, nBinsY, minY, maxY, "track_starty",       dirName);
  PlotHists(inFile, outFile, "track_startz",       "Reco Track Start Position Z", "Reco Track Start Position Z (cm)",  histYaxis, nBinsZ, minZ, maxZ, "track_startz",       dirName);
  PlotHists(inFile, outFile, "track_nhits",        "Reco Track Number of Hits",   "Reco Track Number of Hits",         histYaxis, 500, 0, 5000,       "track_nhits",        dirName);
  PlotHists(inFile, outFile, "track_purity",       "Track Purity",                "Track Purity",                      histYaxis, 50, 0, 1.,          "track_purity",       dirName);
  PlotHists(inFile, outFile, "track_completeness", "Track Completeness",          "Track Completeness",                histYaxis, 50, 0, 1.,          "track_completeness", dirName);

  PlotHists(inFile, outFile, "cluster_nhits",        "Cluster Number of Hits", "Cluster Number of Hits", histYaxis, nBins, 0, 1500, "cluster_nhits", dirName);
  PlotHists(inFile, outFile, "cluster_purity",       "Cluster Purity",         "Cluster Purity",         histYaxis, nBins, 0, 1500, "cluster_purity", dirName);
  PlotHists(inFile, outFile, "cluster_completeness", "Cluster Completeness",   "Cluster Completeness",   histYaxis, nBins, 0, 1500, "cluster_completeness", dirName);
  
}


void PlotHists(TFile *inFile, TFile *outFile, TString histName, TString histTitle, TString histXaxis, TString histYaxis, int nBins, int xMin, int xMax, TString metricToPlot, TString dirName)
{

  TH1F* hist = new TH1F(histName, histTitle+";"+histXaxis+";"+histYaxis, nBins, xMin, xMax); 

  TTree* metricTree; 
  inFile->GetObject("trackvalidation", metricTree);
  std::cout<<"histName: " << histName << std::endl;
  metricTree->Draw(metricToPlot + " >> " + histName);
  gStyle->SetOptStat(0);

  outFile->cd();
  outFile->cd(dirName.Data());

  hist->Write("", TObject::kOverwrite);
  
  TCanvas *c1 = new TCanvas("c1","c1",800,1000);
  hist->Draw();
  c1->Print(Form("%s/%s.png", dirName.Data(), histName.Data()));
  
  delete c1;
  delete hist;
  delete metricTree;
  return;
}