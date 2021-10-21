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

void sbnciplot_crtvalidation_icarus(TString inputFile)
{
  
  TFile *inFile = TFile::Open(inputFile.Data());
  TFile *outFile = new TFile("ci_validation_histos.root", "RECREATE");

  TString histName, histTitle, histXaxis, metricToPlot;
  int xMin, xMax;
  TString histYaxis = "Counts";
  int nBins = 50;

  TString dirName = "crtvalidation";

  outFile->cd();
  gDirectory->mkdir(dirName);
  gSystem->Exec(Form("mkdir -v %s", dirName.Data()));


  PlotHists(inFile, outFile, "crtvalidation/crtTree", "HitPESHit","", "HitPESHit", histYaxis, nBins, 0, 1000, "HitPESHit", "", dirName);
  PlotHists(inFile, outFile, "crtvalidation/crtTree", "x_pos","", "x_pos", histYaxis, nBins, -500, 500, "x_pos", "", dirName);
  PlotHists(inFile, outFile, "crtvalidation/crtTree", "x_err","", "x_err", histYaxis, nBins, 0, 500, "x_err", "", dirName);
  PlotHists(inFile, outFile, "crtvalidation/crtTree", "y_pos","", "y_pos", histYaxis, nBins, -500, 500, "y_pos", "", dirName);
  PlotHists(inFile, outFile, "crtvalidation/crtTree", "y_err","", "y_err", histYaxis, nBins, 0, 500, "y_err", "", dirName);
  PlotHists(inFile, outFile, "crtvalidation/crtTree", "z_pos","", "z_pos", histYaxis, nBins, -1000, 1000, "z_pos", "", dirName);
  PlotHists(inFile, outFile, "crtvalidation/crtTree", "z_err","", "z_err", histYaxis, nBins, 0, 500, "z_err", "", dirName);

  PlotHists(inFile, outFile, "crtvalidation/crtTree", "TrackPESHit","", "TrackPESHit", histYaxis, nBins, 0, 1000, "HitPESHit", "", dirName);
  PlotHists(inFile, outFile, "crtvalidation/crtTree", "x1_pos","", "x2_pos", histYaxis, nBins, -500, 500, "x1_pos", "", dirName);
  PlotHists(inFile, outFile, "crtvalidation/crtTree", "x1_err","", "x2_err", histYaxis, nBins, 0, 500, "x1_err", "", dirName);
  PlotHists(inFile, outFile, "crtvalidation/crtTree", "y1_pos","", "y2_pos", histYaxis, nBins, -500, 500, "y1_pos", "", dirName);
  PlotHists(inFile, outFile, "crtvalidation/crtTree", "y1_err","", "y2_err", histYaxis, nBins, 0, 500, "y1_err", "", dirName);
  PlotHists(inFile, outFile, "crtvalidation/crtTree", "z1_pos","", "z2_pos", histYaxis, nBins, -1000, 1000, "z1_pos", "", dirName);
  PlotHists(inFile, outFile, "crtvalidation/crtTree", "z1_err","", "z2_err", histYaxis, nBins, 0, 500, "z1_err", "", dirName);
  PlotHists(inFile, outFile, "crtvalidation/crtTree", "x2_pos","", "x2_pos", histYaxis, nBins, -500, 500, "x2_pos", "", dirName);
  PlotHists(inFile, outFile, "crtvalidation/crtTree", "x2_err","", "x2_err", histYaxis, nBins, 0, 500, "x2_err", "", dirName);
  PlotHists(inFile, outFile, "crtvalidation/crtTree", "y2_pos","", "y2_pos", histYaxis, nBins, -500, 500, "y2_pos", "", dirName);
  PlotHists(inFile, outFile, "crtvalidation/crtTree", "y2_err","", "y2_err", histYaxis, nBins, 0, 500, "y2_err", "", dirName);
  PlotHists(inFile, outFile, "crtvalidation/crtTree", "z2_pos","", "z2_pos", histYaxis, nBins, -1000, 1000, "z2_pos", "", dirName);
  PlotHists(inFile, outFile, "crtvalidation/crtTree", "z2_err","", "z2_err", histYaxis, nBins, 0, 500, "z2_err", "", dirName);
  PlotHists(inFile, outFile, "crtvalidation/crtTree", "Length","", "Length", histYaxis, nBins, 0, 2000, "Length", "", dirName);
  PlotHists(inFile, outFile, "crtvalidation/crtTree", "Thetaxy","", "Thetaxy", histYaxis, nBins, -10, 10, "Thetaxy", "", dirName);
  PlotHists(inFile, outFile, "crtvalidation/crtTree", "Phizy","", "Phizy", histYaxis, nBins, -10, 10, "Phizy", "", dirName);
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
