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

void sbnciplot_showervalidation()
{
  int PlotHists(TFile *inFile, TString histName, TString histTitle, TString histXaxis, TString histYaxis, int nBins, int xMin, int xMax, TString metricToPlot);
  
  TString inputFile = "/sbnd/data/users/ascarff/ci/ShowerValidation_trees.root";
  TFile *inFile = new TFile(inputFile);

  TString histName = "sTrueEnergy";
  TString histTitle = "True Shower Energy";
  TString histXaxis = "True Shower Energy";
  TString histYaxis = "Counts";
  int nBins = 100;
  int xMin = 0;
  int xMax = 2000;

  TString metricToPlot = "sTrueEnergy_tracs";

  PlotHists(inFile,  histName, histTitle, histXaxis, histYaxis, nBins, xMin, xMax, metricToPlot);
  PlotHists(inFile,  "pfpShowerHitsPurity", "Shower Hit Purity", "Purity", histYaxis, nBins, 0, 1, "pfpShowerHitsPurity_tracs");


}


int PlotHists(TFile *inFile, TString histName, TString histTitle, TString histXaxis, TString histYaxis, int nBins, int xMin, int xMax, TString metricToPlot)
{
  TH1F* hist = new TH1F(histName, histTitle+";"+histXaxis+";"+histYaxis, nBins, xMin, xMax); 

  TTree* metricTree; 
  inFile->GetObject("showerValidation/MetricTree", metricTree);
  metricTree->Draw(metricToPlot + " >> " + histName);
  
  TFile *outFile = new TFile("ShowerValidation_plots.root","RECREATE");
  gStyle->SetOptStat(0);
  hist->Write();
  
  return(1);
}
