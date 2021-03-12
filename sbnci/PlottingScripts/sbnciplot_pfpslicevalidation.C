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

void sbnciplot_pfpslicevalidation(TString inputFile)
{
  
  TFile *inFile = TFile::Open(inputFile.Data());
  TFile *outFile = new TFile("ci_validation_histos.root", "RECREATE");

  TString histName, histTitle, histXaxis, metricToPlot;
  int xMin, xMax;
  TString histYaxis = "Counts";
  int nBins = 100;

  TString dirName = "pfpslicevalidation";

  outFile->cd();
  gDirectory->mkdir(dirName);
  gSystem->Exec(Form("mkdir -v %s", dirName.Data()));


  PlotHists(inFile, outFile, "pfpslicevalidation/trueTree", "numSlices", "Slices", "Slices", histYaxis, 10, 0, 10, "numSlices_pandora", dirName);
  PlotHists(inFile, outFile, "pfpslicevalidation/trueTree", "numNeutrinos", "Neutrino slices", "Neutrino slices", histYaxis, 10, 0, 10, "numNeutrinos_pandora", dirName);
  PlotHists(inFile, outFile, "pfpslicevalidation/trueTree", "purity", "Slice purity", "Slice purity", histYaxis, nBins, 0, 1, "purity_pandora", dirName);
  PlotHists(inFile, outFile, "pfpslicevalidation/trueTree", "comp", "Slice completeness", "Slice completness", histYaxis, nBins, 0, 1, "comp_pandora", dirName);
  PlotHists(inFile, outFile, "pfpslicevalidation/trueTree", "score", "Slice ID score", "Slice ID score", histYaxis, nBins, 0, 1, "score_pandora", dirName);
  PlotHists(inFile, outFile, "pfpslicevalidation/trueTree", "recoPDG", "Slice PDG ID", "Slice PDG ID", histYaxis, 4, 11, 15, "recoPDG_pandora", dirName);
  PlotHists(inFile, outFile, "pfpslicevalidation/trueTree", "pfpVertexX", "Reconstructed vertex x", "Reconstructed vertex x (cm)", histYaxis, nBins, -200,200, "pfpVertexX_pandora", dirName);
  PlotHists(inFile, outFile, "pfpslicevalidation/trueTree", "pfpVertexY", "Reconstructed vertex y", "Reconstructed vertex y (cm)", histYaxis, nBins, -200,200, "pfpVertexY_pandora", dirName);
  PlotHists(inFile, outFile, "pfpslicevalidation/trueTree", "pfpVertexZ", "Reconstructed vertex z", "Reconstructed vertex z (cm)", histYaxis, nBins, 0, 500, "pfpVertexZ_pandora", dirName);
  PlotHists(inFile, outFile, "pfpslicevalidation/trueTree", "pfpVertexDistX", "(reco - true) vertex x", "(reco - true) vertex x (cm)", histYaxis, nBins, -10, 10, "pfpVertexDistX_pandora", dirName);
  PlotHists(inFile, outFile, "pfpslicevalidation/trueTree", "pfpVertexDistY", "(reco - true) vertex y", "(reco - true) vertex y (cm)", histYaxis, nBins, -10, 10, "pfpVertexDistY_pandora", dirName);
  PlotHists(inFile, outFile, "pfpslicevalidation/trueTree", "pfpVertexDistZ", "(reco - true) vertex z", "(reco - true) vertex z (cm)", histYaxis, nBins, -10, 10, "pfpVertexDistZ_pandora", dirName);

  PlotHists(inFile, outFile, "pfpslicevalidation/eventTree", "cosmicScores", "Cosmic slice ID scores", "Cosmic slice ID scores", histYaxis, nBins, 0, 1, "cosmicScores_pandora", dirName);
  PlotHists(inFile, outFile, "pfpslicevalidation/eventTree", "nuScores", "Neutrino slice ID scores", "Neutrino slice ID scores", histYaxis, nBins, 0, 1, "nuScores_pandora", dirName);

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
