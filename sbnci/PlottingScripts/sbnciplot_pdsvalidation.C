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

void sbnciplot_pdsvalidation(TString inputFile)
{
  
  TFile *inFile = TFile::Open(inputFile.Data());
  TFile *outFile = new TFile("ci_validation_histos.root", "RECREATE");

  TString histName, histTitle, histXaxis, metricToPlot;
  int xMin, xMax;
  TString histYaxis = "Counts";
  int nBins = 50;

  TString dirName = "pdsvalidation";

  outFile->cd();
  gDirectory->mkdir(dirName);
  gSystem->Exec(Form("mkdir -v %s", dirName.Data()));


  PlotHists(inFile, outFile, "pdsvalidation/pdsTree", "TotalFlashPE","", "TotalFlashPE", histYaxis, nBins, 0, 2000, "TotalFlashPE", "", dirName);
  PlotHists(inFile, outFile, "pdsvalidation/pdsTree", "FlashTime","", "FlashTime", histYaxis, nBins, 0, 100, "FlashTime", "", dirName);
  PlotHists(inFile, outFile, "pdsvalidation/pdsTree", "FlashTimeWidth","", "FlashTimeWidth", histYaxis, nBins, 0, 100, "FlashTimeWidth", "", dirName);
  PlotHists(inFile, outFile, "pdsvalidation/pdsTree", "YCenter","", "YCenter", histYaxis, nBins, 0, 500, "YCenter", "", dirName);
  PlotHists(inFile, outFile, "pdsvalidation/pdsTree", "ZCenter","", "ZCenter", histYaxis, nBins, 0, 500, "ZCenter", "", dirName);
  PlotHists(inFile, outFile, "pdsvalidation/pdsTree", "YWidth","","YWidth", histYaxis, nBins, 0, 500, "YWidth", "", dirName);
  PlotHists(inFile, outFile, "pdsvalidation/pdsTree", "ZWidth","", "ZWidth", histYaxis, nBins, 0, 500, "ZWidth", "", dirName);
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
