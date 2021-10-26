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

void PlotHists(TFile *inFile, TFile *outFile, TString histName, TString histTitle, TString histXaxis, TString histYaxis, int nBins, int xMin, int xMax, TString metricToPlot, TString dirName);

void sbnciplot_showervalidation(TString inputFile)
{
  
  TFile *inFile = TFile::Open(inputFile.Data());
  TFile *outFile = new TFile("ci_validation_histos.root", "RECREATE");

  TString histName, histTitle, histXaxis, metricToPlot;
  int xMin, xMax;
  TString histYaxis = "Counts";
  int nBins = 30;

  TString dirName = "pandoraShower";

  outFile->cd();
  gDirectory->mkdir(dirName);
  gSystem->Exec(Form("mkdir -v %s", dirName.Data()));


  PlotHists(inFile, outFile, "sTrueEnergy", "True Shower Energy", "True Shower Energy", histYaxis, nBins, 0, 2000, "sTrueEnergy_pandoraShower", dirName);
  PlotHists(inFile, outFile, "pfpShowerHitsPurity", "Shower Hit Purity", "Purity", histYaxis, nBins, 0, 1, "pfpShowerHitsPurity_pandoraShower", dirName);
  PlotHists(inFile, outFile, "pfpShowerHitsComp", "Shower Hit Completeness", "Completeness", histYaxis, nBins, 0, 1, "pfpShowerHitsComp_pandoraShower", dirName);
  PlotHists(inFile, outFile, "sDirDiff", "Shower Start Direction Difference", "X", histYaxis, nBins, -1, 1, "sDirDiff_pandoraShower", dirName);
  PlotHists(inFile, outFile, "sStartDist", "Shower Start Position Difference", "Y", histYaxis, nBins, 0, 60, "sStartDist_pandoraShower", dirName);
  PlotHists(inFile, outFile, "sdEdx", "Shower Start dE/dx", "dE/dx", histYaxis, nBins, 0, 20, "sdEdx_pandoraShower", dirName);
  PlotHists(inFile, outFile, "sEnergy", "Shower Start Energy", "Energy", histYaxis, nBins, 0, 1500, "sEnergy_pandoraShower", dirName);
  
}


void PlotHists(TFile *inFile, TFile *outFile, TString histName, TString histTitle, TString histXaxis, TString histYaxis, int nBins, int xMin, int xMax, TString metricToPlot, TString dirName)
{

  TH1F* hist = new TH1F(histName, histTitle+";"+histXaxis+";"+histYaxis, nBins, xMin, xMax); 

  TTree* metricTree; 
  inFile->GetObject("showervalidation/MetricTree", metricTree);
  std::cout<<"metricToPlot: " << metricToPlot << ". histName: " << histName << std::endl;
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
