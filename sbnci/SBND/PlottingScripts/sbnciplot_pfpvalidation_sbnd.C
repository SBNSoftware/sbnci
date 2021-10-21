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

void PlotHists(TFile *inFile, TFile *outFile, TString treeName, TString histName, TString histTitle, TString histXaxis, TString histYaxis, int nBins, int xMin, int xMax, TString metricToPlot, TString condition, TString dirName);

void sbnciplot_pfpvalidation_sbnd(TString inputFile)
{
  
  TFile *inFile = TFile::Open(inputFile.Data());
  TFile *outFile = new TFile("ci_validation_histos.root", "RECREATE");

  TString histName, histTitle, histXaxis, metricToPlot;
  int xMin, xMax;
  TString histYaxis = "Counts";
  int nBins = 40;

  TString dirName = "pfpvalidation";

  outFile->cd();
  gDirectory->mkdir(dirName);
  gSystem->Exec(Form("mkdir -v %s", dirName.Data()));


  PlotHists(inFile, outFile, "pfpvalidation/pfpTree", "pfpNumHits", "PFP nHits", "PFP nHits", histYaxis, nBins, 0, 2000, "pfpNumHits_pandora", "", dirName);
  PlotHists(inFile, outFile, "pfpvalidation/pfpTree", "pfpTrackScore", "PFP track score", "PFP track score", histYaxis, nBins, 0, 1, "pfpTrackScore_pandora", "", dirName);
  PlotHists(inFile, outFile, "pfpvalidation/pfpTree", "pfpHitPurity", "PFP hit purity", "PFP hit purity", histYaxis, nBins, 0, 1, "pfpHitPurity_pandora", "", dirName);
  PlotHists(inFile, outFile, "pfpvalidation/pfpTree", "pfpHitComp", "PFP hit completeness", "PFP hit completeness", histYaxis, nBins, 0, 1, "pfpHitComp_pandora", "", dirName);
  PlotHists(inFile, outFile, "pfpvalidation/pfpTree", "pfpEnergyPurity", "PFP energy purity", "PFP energy purity", histYaxis, nBins, 0, 1, "pfpEnergyPurity_pandora", "", dirName);
  PlotHists(inFile, outFile, "pfpvalidation/pfpTree", "pfpEnergyComp", "PFP energy completeness", "PFP energy completeness", histYaxis, nBins, 0, 1, "pfpEnergyComp_pandora", "", dirName);
  PlotHists(inFile, outFile, "pfpvalidation/pfpTree", "pfpHitSPRatio", "PFP Space point to hit ratio", "PFP Space point to hit ratio", histYaxis, nBins, 0, 1, "pfpHitSPRatio_pandora", "", dirName);

  PlotHists(inFile, outFile, "pfpvalidation/trueTree", "muonRecoPFPs", "Reco PFPs per true muon", "Reco PFPs per true muon", histYaxis, 10, 0, 10, "recoPFPs_pandora", "truePdg==13", dirName);
  PlotHists(inFile, outFile, "pfpvalidation/trueTree", "muonRecoPFPTracks", "Reco tracks per true muon", "Reco tracks per true muon", histYaxis, 10, 0, 10, "recoPFPTracks_pandora", "truePdg==13", dirName);
  PlotHists(inFile, outFile, "pfpvalidation/trueTree", "muonRecoPFPShowers", "Reco showers per true muon", "Reco showers per true muon", histYaxis, 10, 0, 10, "recoPFPShowers_pandora", "truePdg==13", dirName);
  PlotHists(inFile, outFile, "pfpvalidation/trueTree", "muonHitPurity", "Hit purity per true muon", "Hit purity per true muon", histYaxis, nBins, 0, 1, "hitPurity_pandora", "truePdg==13", dirName);
  PlotHists(inFile, outFile, "pfpvalidation/trueTree", "muonHitComp", "Hit completeness per true muon", "Hit completeness per true muon", histYaxis, nBins, 0, 1, "hitComp_pandora", "truePdg==13", dirName);
  PlotHists(inFile, outFile, "pfpvalidation/trueTree", "muonEnergyPurity", "Energy purity per true muon", "Energy purity per true muon", histYaxis, nBins, 0, 1, "energyPurity_pandora", "truePdg==13", dirName);
  PlotHists(inFile, outFile, "pfpvalidation/trueTree", "muonEnergyComp", "Energy completeness per true muon", "Energy completeness per true muon", histYaxis, nBins, 0, 1, "energyComp_pandora", "truePdg==13", dirName);
  PlotHists(inFile, outFile, "pfpvalidation/trueTree", "muonHitSPRatio", "Space point to hit ratio per true muon", "Space point to hit ratio per true muon", histYaxis, nBins, 0, 1, "hitPurity_pandora", "truePdg==13", dirName);

  PlotHists(inFile, outFile, "pfpvalidation/trueTree", "protonRecoPFPs", "Reco PFPs per true proton", "Reco PFPs per true proton", histYaxis, 10, 0, 10, "recoPFPs_pandora", "truePdg==2212", dirName);
  PlotHists(inFile, outFile, "pfpvalidation/trueTree", "protonRecoPFPTracks", "Reco tracks per true proton", "Reco tracks per true proton", histYaxis, 10, 0, 10, "recoPFPTracks_pandora", "truePdg==2212", dirName);
  PlotHists(inFile, outFile, "pfpvalidation/trueTree", "protonRecoPFPShowers", "Reco showers per true proton", "Reco showers per true proton", histYaxis, 10, 0, 10, "recoPFPShowers_pandora", "truePdg==2212", dirName);
  PlotHists(inFile, outFile, "pfpvalidation/trueTree", "protonHitPurity", "Hit purity per true proton", "Hit purity per true proton", histYaxis, nBins, 0, 1, "hitPurity_pandora", "truePdg==2212", dirName);
  PlotHists(inFile, outFile, "pfpvalidation/trueTree", "protonHitComp", "Hit completeness per true proton", "Hit completeness per true proton", histYaxis, nBins, 0, 1, "hitComp_pandora", "truePdg==2212", dirName);
  PlotHists(inFile, outFile, "pfpvalidation/trueTree", "protonEnergyPurity", "Energy purity per true proton", "Energy purity per true proton", histYaxis, nBins, 0, 1, "energyPurity_pandora", "truePdg==2212", dirName);
  PlotHists(inFile, outFile, "pfpvalidation/trueTree", "protonEnergyComp", "Energy completeness per true proton", "Energy completeness per true proton", histYaxis, nBins, 0, 1, "energyComp_pandora", "truePdg==2212", dirName);
  PlotHists(inFile, outFile, "pfpvalidation/trueTree", "protonHitSPRatio", "Space point to hit ratio per true proton", "Space point to hit ratio per true proton", histYaxis, nBins, 0, 1, "hitPurity_pandora", "truePdg==2212", dirName);

  PlotHists(inFile, outFile, "pfpvalidation/trueTree", "electronRecoPFPs", "Reco PFPs per true electron", "Reco PFPs per true electron", histYaxis, 10, 0, 10, "recoPFPs_pandora", "truePdg==11", dirName);
  PlotHists(inFile, outFile, "pfpvalidation/trueTree", "electronRecoPFPTracks", "Reco tracks per true electron", "Reco tracks per true electron", histYaxis, 10, 0, 10, "recoPFPTracks_pandora", "truePdg==11", dirName);
  PlotHists(inFile, outFile, "pfpvalidation/trueTree", "electronRecoPFPShowers", "Reco showers per true electron", "Reco showers per true electron", histYaxis, 10, 0, 10, "recoPFPShowers_pandora", "truePdg==11", dirName);
  PlotHists(inFile, outFile, "pfpvalidation/trueTree", "electronHitPurity", "Hit purity per true electron", "Hit purity per true electron", histYaxis, nBins, 0, 1, "hitPurity_pandora", "truePdg==11", dirName);
  PlotHists(inFile, outFile, "pfpvalidation/trueTree", "electronHitComp", "Hit completeness per true electron", "Hit completeness per true electron", histYaxis, nBins, 0, 1, "hitComp_pandora", "truePdg==11", dirName);
  PlotHists(inFile, outFile, "pfpvalidation/trueTree", "electronEnergyPurity", "Energy purity per true electron", "Energy purity per true electron", histYaxis, nBins, 0, 1, "energyPurity_pandora", "truePdg==11", dirName);
  PlotHists(inFile, outFile, "pfpvalidation/trueTree", "electronEnergyComp", "Energy completeness per true electron", "Energy completeness per true electron", histYaxis, nBins, 0, 1, "energyComp_pandora", "truePdg==11", dirName);
  PlotHists(inFile, outFile, "pfpvalidation/trueTree", "electronHitSPRatio", "Space point to hit ratio per true electron", "Space point to hit ratio per true electron", histYaxis, nBins, 0, 1, "hitPurity_pandora", "truePdg==11", dirName);

  PlotHists(inFile, outFile, "pfpvalidation/eventTree", "numPFP", "Number of PFPs per event", "Number of PFPs per event", histYaxis, 20, 0, 20, "numPFP_pandora", "", dirName);
  PlotHists(inFile, outFile, "pfpvalidation/eventTree", "numPFPNu", "Number of neutrino PFPs per event", "Number of neutrino PFPs per event", histYaxis, 10, 0, 10, "numPFPNu_pandora", "", dirName);
  PlotHists(inFile, outFile, "pfpvalidation/eventTree", "numPFPShower", "Number of showers per event", "Number of showers per event", histYaxis, 10, 0, 10, "numPFPShower_pandora", "", dirName);
  PlotHists(inFile, outFile, "pfpvalidation/eventTree", "numPFPTrack", "Number of tracks per event", "Number of tracks per event", histYaxis, 10, 0, 10, "numPFPTrack_pandora", "", dirName);
}


void PlotHists(TFile *inFile, TFile *outFile, TString treeName, TString histName, TString histTitle, TString histXaxis, TString histYaxis, int nBins, int xMin, int xMax, TString metricToPlot, TString condition, TString dirName)
{

  TH1F* hist = new TH1F(histName, histTitle+";"+histXaxis+";"+histYaxis, nBins, xMin, xMax); 

  TTree* metricTree; 
  inFile->GetObject(treeName, metricTree);
  std::cout<<"metricToPlot: " << metricToPlot << ". histName: " << histName << std::endl;
  metricTree->Draw(metricToPlot + " >> " + histName, condition);
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
