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

void PlotHists(TFile *inFile, TFile *outFile, TString treeName, TString histName, TString histTitle, TString histXaxis, TString histYaxis, int nBins, double xMin, double xMax, TString metricToPlot, TString condition, TString dirName);

void sbnciplot_pdssimvalidation_sbnd(TString inputFile)
{

  TFile *inFile = TFile::Open(inputFile.Data());
  TFile *outFile = new TFile("ci_validation_histos.root", "RECREATE");

  TString histName, histTitle, histXaxis, metricToPlot;
  int xMin, xMax;
  TString histYaxis = "Counts";
  int nBins = 50;

  TString dirName = "pdssimvalidation";

  outFile->cd();
  gDirectory->mkdir(dirName);
  gSystem->Exec(Form("mkdir -v %s", dirName.Data()));

  int nBinsPE=25;
  int nBinsPEOpCh=60;
  PlotHists(inFile, outFile, "pdssimvalidation/pdsSimTree", "NPhotonsVec","", "#SimPhotons (per OpCh)", histYaxis, nBinsPEOpCh, 0, 10000, "NPhotonsVec", "", dirName);
  PlotHists(inFile, outFile, "pdssimvalidation/pdsSimTree", "NPhotons","", "#SimPhotons", histYaxis, nBinsPE, 0, 60000, "NPhotons", "", dirName);
  PlotHists(inFile, outFile, "pdssimvalidation/pdsSimTree", "NPhotonsPMTCoated","", "#SimPhotons_{CoatedPMTs}", histYaxis, nBinsPE, 0, 40000, "NPhotonsPMTCoated", "", dirName);
  PlotHists(inFile, outFile, "pdssimvalidation/pdsSimTree", "NPhotonsPMTUncoated","", "#SimPhotons_{UncoatedPMTs}", histYaxis, nBinsPE, 0, 5000, "NPhotonsPMTUncoated", "", dirName);
  PlotHists(inFile, outFile, "pdssimvalidation/pdsSimTree", "NPhotonsXARAPUCAVuv","", "#SimPhotons_{XARAPUCAVuv}", histYaxis, nBinsPE, 0, 5000, "NPhotonsXARAPUCAVuv", "", dirName);
  PlotHists(inFile, outFile, "pdssimvalidation/pdsSimTree", "NPhotonsXARAPUCAVis","", "#SimPhotons_{XARAPUCAVis}", histYaxis, nBinsPE, 0, 5000, "NPhotonsXARAPUCAVis", "", dirName);

}

void PlotHists(TFile *inFile, TFile *outFile, TString treeName, TString histName, TString histTitle, TString histXaxis, TString histYaxis, int nBins, double xMin, double xMax, TString metricToPlot, TString condition, TString dirName)
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
