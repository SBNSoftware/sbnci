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

void sbnciplot_tpcsimvalidation_icarus(TString inputFile)
{
  
  TFile *inFile = TFile::Open(inputFile.Data());
  TFile *outFile = new TFile("ci_validation_histos.root", "RECREATE");

  TString histName, histTitle, histXaxis, metricToPlot;
  int xMin, xMax;
  TString histYaxis = "Counts";
  int nBins = 50;

  TString dirName = "tpcsimvalidation";

  outFile->cd();
  gDirectory->mkdir(dirName);
  gSystem->Exec(Form("mkdir -v %s", dirName.Data()));


  PlotHists(inFile, outFile, "tpcsimvalidation/tpcsimTree", "NADC","", "NADC", histYaxis, nBins, 0, 5000, "NADC", "", dirName);
  PlotHists(inFile, outFile, "tpcsimvalidation/tpcsimTree", "Integrated ADC Counts","", "Integrated ADC Cts/Channel", histYaxis, nBins, 1, 5000, "ADC", "", dirName);
  PlotHists(inFile, outFile, "tpcsimvalidation/tpcsimTree", "Samples","", "Samples", histYaxis, nBins, 0, 5000, "Samples", "", dirName);
  PlotHists(inFile, outFile, "tpcsimvalidation/tpcsimTree", "Pedestal","", "Pedestal", histYaxis, nBins, 0, 5000, "Pedestal", "", dirName);
  PlotHists(inFile, outFile, "tpcsimvalidation/tpcsimTree", "Sigma","", "Sigma", histYaxis, 10, 0, 10, "Sigma", "", dirName);
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
