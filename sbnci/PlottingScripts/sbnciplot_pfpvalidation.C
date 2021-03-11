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

void PlotHists(TFile *inFile, TFile *outFile, TString treeName, TString histName, TString histTitle, TString histXaxis, TString histYaxis, int nBins, int xMin, int xMax, TString metricToPlot, TString dirName);

void sbnciplot_pfpvalidation(TString inputFile)
{
  
  TFile *inFile = TFile::Open(inputFile.Data());
  TFile *outFile = new TFile("ci_validation_histos.root", "RECREATE");

  TString histName, histTitle, histXaxis, metricToPlot;
  int xMin, xMax;
  TString histYaxis = "Counts";
  int nBins = 100;

  TString dirName = "pfpvalidation";

  outFile->cd();
  gDirectory->mkdir(dirName);
  gSystem->Exec(Form("mkdir -v %s", dirName.Data()));


  PlotHists(inFile, outFile, "pfpvalidation/pfpTree", "pfpNumHits", "PFP nHits", "PFP nHits", histYaxis, nBins, 0, 2000, "pfpNumHits_pandora", dirName);
  PlotHists(inFile, outFile, "pfpvalidation/pfpTree", "pfpTrackScore", "PFP track score", "PFP track score", histYaxis, nBins, 0, 1, "pfpTrackScore_pandora", dirName);
  PlotHists(inFile, outFile, "pfpvalidation/pfpTree", "pfpHitPurity", "PFP hit purity", "PFP hit purity", histYaxis, nBins, 0, 1, "pfpHitPurity_pandora", dirName);
  PlotHists(inFile, outFile, "pfpvalidation/pfpTree", "pfpHitComp", "PFP hit completeness", "PFP hit completeness", histYaxis, nBins, 0, 1, "pfpHitComp_pandora", dirName);
  PlotHists(inFile, outFile, "pfpvalidation/pfpTree", "pfpEnergyPurity", "PFP energy purity", "PFP energy purity", histYaxis, nBins, 0, 1, "pfpEnergyPurity_pandora", dirName);
  PlotHists(inFile, outFile, "pfpvalidation/pfpTree", "pfpEnergyComp", "PFP energy completeness", "PFP energy completeness", histYaxis, nBins, 0, 1, "pfpEnergyComp_pandora", dirName);
  PlotHists(inFile, outFile, "pfpvalidation/pfpTree", "pfpHitSPRatio", "PFP Space point to hit ratio", "PFP Space point to hit ratio", histYaxis, nBins, 0, 1, "pfpHitSPRatio_pandora", dirName);



  PlotHists(inFile, outFile, "pfpvalidation/trueTree", "recoPFPs", "Reco PFPs per true particle", "Reco PFPs per true particle", histYaxis, 10, 0, 10, "recoPFPs_pandora", dirName);
  PlotHists(inFile, outFile, "pfpvalidation/trueTree", "recoPFPTracks", "Reco tracks per true particle", "Reco tracks per true particle", histYaxis, 10, 0, 10, "recoPFPTracks_pandora", dirName);
  PlotHists(inFile, outFile, "pfpvalidation/trueTree", "recoPFPShowers", "Reco showers per true particle", "Reco showers per true particle", histYaxis, 10, 0, 10, "recoPFPShowers_pandora", dirName);
  PlotHists(inFile, outFile, "pfpvalidation/trueTree", "hitPurity", "Hit purity per true particle", "Hit purity per true particle", histYaxis, nBins, 0, 1, "hitPurity_pandora", dirName);
  PlotHists(inFile, outFile, "pfpvalidation/trueTree", "hitComp", "Hit completeness per true particle", "Hit completeness per true particle", histYaxis, nBins, 0, 1, "hitComp_pandora", dirName);
  PlotHists(inFile, outFile, "pfpvalidation/trueTree", "energyPurity", "Energy purity per true particle", "Energy purity per true particle", histYaxis, nBins, 0, 1, "energyPurity_pandora", dirName);
  PlotHists(inFile, outFile, "pfpvalidation/trueTree", "energyComp", "Energy completeness per true particle", "Energy completeness per true particle", histYaxis, nBins, 0, 1, "energyComp_pandora", dirName);
  PlotHists(inFile, outFile, "pfpvalidation/trueTree", "hitSPRatio", "Space point to hit ratio per true particle", "Space point to hit ratio per true particle", histYaxis, nBins, 0, 1, "hitPurity_pandora", dirName);

  PlotHists(inFile, outFile, "pfpvalidation/eventTree", "numPFP", "Number of PFPs per event", "Number of PFPs per event", histYaxis, 20, 0, 20, "numPFP_pandora", dirName);
  PlotHists(inFile, outFile, "pfpvalidation/eventTree", "numPFPNu", "Number of neutrino PFPs per event", "Number of neutrino PFPs per event", histYaxis, 10, 0, 10, "numPFPNu_pandora", dirName);
  PlotHists(inFile, outFile, "pfpvalidation/eventTree", "numPFPShower", "Number of showers per event", "Number of showers per event", histYaxis, 10, 0, 10, "numPFPShower_pandora", dirName);
  PlotHists(inFile, outFile, "pfpvalidation/eventTree", "numPFPTrack", "Number of tracks per event", "Number of tracks per event", histYaxis, 10, 0, 10, "numPFPTrack_pandora", dirName);
}


void PlotHists(TFile *inFile, TFile *outFile, TString treeName, TString histName, TString histTitle, TString histXaxis, TString histYaxis, int nBins, int xMin, int xMax, TString metricToPlot, TString dirName)
{

  TH1F* hist = new TH1F(histName, histTitle+";"+histXaxis+";"+histYaxis, nBins, xMin, xMax); 

  TTree* metricTree; 
  inFile->GetObject(treeName, metricTree);
  std::cout<<"metricToPlot: " << metricToPlot << ". histName: " << histName << std::endl;
  metricTree->Draw(metricToPlot + " >> " + histName,metricToPlot + " != -999");
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