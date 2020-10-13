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

void PlotHists(TFile *inFile, TString histName, TString histTitle, TString histXaxis, TString histYaxis, int nBins, int xMin, int xMax, TString metricToPlot);

void sbnciplot_showervalidation(TString inputFile)
{
  
  //TString inputFile = "root://fndca1.fnal.gov:1094/pnfs/fnal.gov/usr/sbnd/persistent/users/ascarff/ciplots/ShowerValidation_trees.root";
  TFile *inFile = TFile::Open(inputFile.Data());
  TString histName = "sTrueEnergy";
  TString histTitle = "True Shower Energy";
  TString histXaxis = "True Shower Energy";
  TString histYaxis = "Counts";
  int nBins = 100;
  int xMin = 0;
  int xMax = 2000;

  TString metricToPlot = "sTrueEnergy_tracs";

  TFile *outFile = new TFile("ShowerValidationPlots.root", "RECREATE");

  PlotHists(inFile,  histName, histTitle, histXaxis, histYaxis, nBins, xMin, xMax, metricToPlot);
  PlotHists(inFile,  "pfpShowerHitsPurity", "Shower Hit Purity", "Purity", histYaxis, nBins, 0, 1, "pfpShowerHitsPurity_tracs");
  
}


void PlotHists(TFile *inFile, TString histName, TString histTitle, TString histXaxis, TString histYaxis, int nBins, int xMin, int xMax, TString metricToPlot)
{
  TH1F* hist = new TH1F(histName, histTitle+";"+histXaxis+";"+histYaxis, nBins, xMin, xMax); 

  TTree* metricTree; 
  inFile->GetObject("showerValidation/MetricTree", metricTree);
  metricTree->Draw(metricToPlot + " >> " + histName);
  
  gStyle->SetOptStat(0);
  hist->Write();
  
  TCanvas *c1 = new TCanvas("c1","c1",800,1000);
  hist->Draw();
  c1->Print(histName+".png");
  
  delete c1;
  delete hist;
  return;
}
