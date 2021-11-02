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

void sbnciplot_pdsvalidation_sbnd(TString inputFile)
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

  int nBinsPE=25;
  int nBinsPEOpCh=60;
  PlotHists(inFile, outFile, "pdsvalidation/pdsTree", "NPhotonsVec","", "#SimPhotons (per OpCh)", histYaxis, nBinsPEOpCh, 0, 10000, "NPhotonsVec", "", dirName);
  PlotHists(inFile, outFile, "pdsvalidation/pdsTree", "NPhotons","", "#SimPhotons", histYaxis, nBinsPE, 0, 60000, "NPhotons", "", dirName);
  PlotHists(inFile, outFile, "pdsvalidation/pdsTree", "NPhotonsPMTCoated","", "#SimPhotons_{CoatedPMTs}", histYaxis, nBinsPE, 0, 40000, "NPhotonsPMTCoated", "", dirName);
  PlotHists(inFile, outFile, "pdsvalidation/pdsTree", "NPhotonsPMTUncoated","", "#SimPhotons_{UncoatedPMTs}", histYaxis, nBinsPE, 0, 5000, "NPhotonsPMTUncoated", "", dirName);
  PlotHists(inFile, outFile, "pdsvalidation/pdsTree", "NPhotonsXARAPUCAVuv","", "#SimPhotons_{XARAPUCAVuv}", histYaxis, nBinsPE, 0, 5000, "NPhotonsXARAPUCAVuv", "", dirName);
  PlotHists(inFile, outFile, "pdsvalidation/pdsTree", "NPhotonsXARAPUCAVis","", "#SimPhotons_{XARAPUCAVis}", histYaxis, nBinsPE, 0, 5000, "NPhotonsXARAPUCAVis", "", dirName);


  PlotHists(inFile, outFile, "pdsvalidation/pdsTree", "PMTNHit","", "PMTNHit", histYaxis, nBins, 0, 5000, "PMTNHit", "", dirName);
  PlotHists(inFile, outFile, "pdsvalidation/pdsTree", "PMTPeakTime","", "PMTPeakTime-NuG4Time [#mus]", histYaxis, 200, -10, 40, "PMTPeakTime", "", dirName);
  PlotHists(inFile, outFile, "pdsvalidation/pdsTree", "PMTWidth","", "PMTWidth [#mus] ", histYaxis, nBins, 0, 3, "PMTWidth", "", dirName);
  PlotHists(inFile, outFile, "pdsvalidation/pdsTree", "PMTAmplitude","", "PMTAmplitude [ADC]", histYaxis, nBins, 0, 3050, "PMTAmplitude", "", dirName); //XRange upper value set to PMT ADC saturation value
  PlotHists(inFile, outFile, "pdsvalidation/pdsTree", "PMTPE","", "PMTPE", histYaxis, nBinsPEOpCh, 0, 10000, "PMTPE", "", dirName);
  PlotHists(inFile, outFile, "pdsvalidation/pdsTree", "PMTPEResolutionVec","", "(PMTPE-PMTSimPhot)/PMTSimPhot (per OpCh)", histYaxis, nBins, -1, 1.5, "PMTPEResolutionVec", "", dirName);


  //PlotHists(inFile, outFile, "pdsvalidation/pdsTree", "ArapucaNHit","", "ArapucaNHit", histYaxis, nBins, 0, 500, "ArapucaNHit", "", dirName);
  //PlotHists(inFile, outFile, "pdsvalidation/pdsTree", "ArapucaPeakTime","", "ArapucaPeakTime", histYaxis, nBins, 0, 10, "ArapucaPeakTime", "", dirName);
  //PlotHists(inFile, outFile, "pdsvalidation/pdsTree", "ArapucaWidth","", "ArapucaWidth", histYaxis, nBins, 0, 5, "ArapucaWidth", "", dirName);
  //PlotHists(inFile, outFile, "pdsvalidation/pdsTree", "ArapucaAmplitude","", "ArapucaAmplitude", histYaxis, nBins, 0, 4000, "ArapucaAmplitude", "", dirName);
  //PlotHists(inFile, outFile, "pdsvalidation/pdsTree", "ArapucaPE","", "ArapucaPE", histYaxis, nBinsPEOpCh, 0, 500, "ArapucaPE", "", dirName);
  //PlotHists(inFile, outFile, "pdsvalidation/pdsTree", "ArapucaPEResolutionVec","", "(ArapucaPE-ArapucaSimPhot)/ArapucaSimPhot", histYaxis, nBins, -1, 1.5, "ArapucaPEResolutionVec", "", dirName);


  PlotHists(inFile, outFile, "pdsvalidation/pdsTree", "NFlash","", "NFlash", histYaxis, 10, 0, 10, "NFlash", "", dirName);
  PlotHists(inFile, outFile, "pdsvalidation/pdsTree", "TotalFlashPE","", "TotalFlashPE", histYaxis, nBinsPE, 0,50000, "TotalFlashPE", "", dirName);
  PlotHists(inFile, outFile, "pdsvalidation/pdsTree", "FlashTime","", "FlashTime [#mus]", histYaxis, 20, 0, 2, "FlashTime", "", dirName);
  PlotHists(inFile, outFile, "pdsvalidation/pdsTree", "FlashTimeResolution","", "FlashTime-NuG4Time [ns]", histYaxis, nBins, -0.05, 0.1, "FlashTimeResolution", "", dirName);
  PlotHists(inFile, outFile, "pdsvalidation/pdsTree", "YCenter","", "YCenter [cm]", histYaxis, 20, -200, 200, "YCenter", "", dirName);
  PlotHists(inFile, outFile, "pdsvalidation/pdsTree", "ZCenter","", "ZCenter [cm]", histYaxis, 25, 0, 500, "ZCenter", "", dirName);

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
