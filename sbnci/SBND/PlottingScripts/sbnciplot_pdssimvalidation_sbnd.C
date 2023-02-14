#include "TCanvas.h"
#include "TStyle.h"
#include "TH1.h"
#include "TH2.h"
#include "TGaxis.h"
#include "TRandom.h"
#include "TFile.h"
#include "TLegend.h"
#include "TCollection.h"
#include "THStack.h"
#include "TAxis.h"
#include "TGraphErrors.h"
#include "TString.h"

#include <iostream>
#include <map>
#include <sstream>
#include <vector>

void PlotHists(TFile *inFile, TFile *outFile,
               TString treeName, TString histName, TString histTitle,
               TString histXaxis, TString histYaxis,
               int nBins, double xMin, double xMax,
               TString metricToPlot, TString condition, TString dirName);
void PlotMapHists(TFile *inFile, TFile *outFile,
                  TString treeName, TString mapName, TString histTitle,
                  TString histXaxis, TString histYaxis,
                  int nBins, double xMin, double xMax,
                  TString metricToPlot, TString condition, TString dirName);


void sbnciplot_pdssimvalidation_sbnd(TString inputFile)
{
  TFile *inFile = TFile::Open(inputFile.Data());
  TFile *outFile = new TFile("ci_validation_histos.root", "RECREATE");

  TString histYaxis = "Counts";
  TString dirName = "pdssimvalidation";
  TString treeName = dirName + "/pdsSimTree";
  int nBinsPE = 25;
  int nBinsPEOpCh = 50;

  outFile->cd();
  gDirectory->mkdir(dirName);
  gSystem->Exec(Form("mkdir -v %s", dirName.Data()));

  PlotHists(inFile, outFile, treeName, "NPhotonsVec",         "", "#SimPhotons (per OpCh)",    histYaxis, nBinsPEOpCh, 0, 500, "NPhotonsVec", "", dirName);
  PlotHists(inFile, outFile, treeName, "NPhotons",            "", "#SimPhotons",               histYaxis, nBinsPE, 0, 60000, "NPhotons", "", dirName);
  PlotHists(inFile, outFile, treeName, "NPhotonsPMTCoated",   "", "#SimPhotons_{CoatedPMTs}",  histYaxis, nBinsPE, 0, 40000, "NPhotonsPMTCoated", "", dirName);
  PlotHists(inFile, outFile, treeName, "NPhotonsPMTUncoated", "", "#SimPhotons_{UncoatdPMTs}", histYaxis, nBinsPE, 0, 5000, "NPhotonsPMTUncoated", "", dirName);
  PlotHists(inFile, outFile, treeName, "NPhotonsXARAPUCAVuv", "", "#SimPhotons_{XARAPUCAVuv}", histYaxis, nBinsPE, 0, 5000, "NPhotonsXARAPUCAVuv", "", dirName);
  PlotHists(inFile, outFile, treeName, "NPhotonsXARAPUCAVis", "", "#SimPhotons_{XARAPUCAVis}", histYaxis, nBinsPE, 0, 5000, "NPhotonsXARAPUCAVis", "", dirName);
  PlotMapHists(inFile, outFile, treeName, "NTimePhotonsMap",  "", "T_{SimPhoton} - T_{0} [ns]", histYaxis, 40, 0, 6000, "NTimePhotonsHist", "", dirName);  

}

void PlotHists(TFile *inFile, TFile *outFile,
               TString treeName, TString histName, TString histTitle,
               TString histXaxis, TString histYaxis,
               int nBins, double xMin, double xMax,
               TString metricToPlot, TString condition, TString dirName)
{

  TH1F* hist = new TH1F(histName, histTitle+";"+histXaxis+";"+histYaxis, nBins, xMin, xMax);

  TTree* metricTree;
  inFile->GetObject(treeName, metricTree);
  std::cout << "metricToPlot: " << metricToPlot
            << ". histName: " << histName
            << ". condition: " << condition << std::endl;
  metricTree->Draw(metricToPlot + " >> " + histName, condition);
  gStyle->SetOptStat(0);

  outFile->cd();
  outFile->cd(dirName.Data());

  hist->Write("", TObject::kOverwrite);

  TCanvas *c1 = new TCanvas("c1"+histName,"c1",800,1000);
  hist->Draw();
  c1->Print(Form("%s/%s.png", dirName.Data(), histName.Data()));

  delete c1;
  delete hist;
  delete metricTree;
  return;
}

void PlotMapHists(TFile *inFile, TFile *outFile,
                  TString treeName, TString mapName, TString histTitle,
                  TString histXaxis, TString histYaxis,
                  int nBins, double xMin, double xMax,
                  TString metricToPlot, TString condition, TString dirName)
{
  double xMinFixed = 0.5;
  xMin = (xMin<xMinFixed) ? xMinFixed: xMin;
  nBins = (nBins<10) ? 10 : nBins;
  TString histName = metricToPlot;
  std::vector<double> bin_lims(nBins);
  bin_lims[0] = xMin, bin_lims[1] = 1, bin_lims[2] = 3, bin_lims[3] = 5;
  bin_lims[4] = 8, bin_lims[5] = 12,  bin_lims[6] = 17, bin_lims[7] = 23;

  int idx = 8;
  float logmin = log10(30);
  float logmax = log10(xMax);
  int index = 0;
  for(int i=idx; i<bin_lims.size() ; i++) {
    bin_lims[i] =  pow(10., (logmin + index*(logmax-logmin)/(bin_lims.size()-idx-1)) );
    // std::cout << "index:\t" << index << "\t";
    // std::cout << "i:\t" << i << "\t";
    // std::cout << "bin_lims[i]:\t" << bin_lims[i] << "\n";
    index++;
  }

  // TH1F* hist = new TH1F(histName, histTitle+";"+histXaxis+";"+histYaxis, nBins, xMin, xMax);
  TH1F* hist = new TH1F(histName, histTitle+";"+histXaxis+";"+histYaxis,
                        bin_lims.size()-1, bin_lims.data());

  std::cout << "metricToPlot: " << metricToPlot
            << ". histName: " << histName
            << ". condition: " << condition << std::endl;

  TTree* metricTree;
  inFile->GetObject(treeName, metricTree);
  std::map<int,int>* mapp = 0;
  metricTree->SetBranchAddress(mapName, &mapp);
  for (int iEntry = 0; metricTree->LoadTree(iEntry) >= 0; ++iEntry) {
    metricTree->GetEntry(iEntry);
    for(const auto& ntp : *mapp) {
      // Include under and overflow values
      if(ntp.first <= xMinFixed) hist->AddBinContent(hist->FindBin(xMin), ntp.second);
      else if(ntp.first > xMax) hist->AddBinContent(hist->FindBin(xMax), ntp.second);
      else hist->AddBinContent(hist->FindBin(ntp.first), ntp.second);
    }
  }

  gStyle->SetOptStat(0);

  outFile->cd();
  outFile->cd(dirName.Data());

  TCanvas *c1 = new TCanvas("c1"+histName,"c1",1200,1000);
  c1->SetLogx();

  hist->Write("", TObject::kOverwrite);

  hist->Draw("");
  c1->Print(Form("%s/%s.png", dirName.Data(), histName.Data()));

  delete c1;
  delete hist;
  delete metricTree;
  return;
}
