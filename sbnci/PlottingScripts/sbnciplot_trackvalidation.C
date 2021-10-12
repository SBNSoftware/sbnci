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
void Plot2DHists(TFile *inFile, TFile *outFile, TString histName, TString histTitle, TString histXaxis, TString histYaxis, int nBins, int xMin, int xMax, int yMin, int yMax, TString metricToPlotonX, TString metricToPlotonY, TString dirName);

//Global variables initialisation
std::vector<TString> trkparticles = {"All", "Muon", "Pion", "Proton", "Other"};
//
void sbnciplot_trackvalidation(TString inputFile){
  TFile *inFile = TFile::Open(inputFile.Data());
  TFile *outFile = new TFile("ci_validation_histos.root", "RECREATE");

  TString histName, histTitle, histXaxis, metricToPlot;
  int xMin, xMax;
  TString histYaxis = "Counts";
  int nBins = 100;

  TString dirName = "pandoraTrack";

  outFile->cd();
  gDirectory->mkdir(dirName);
  gSystem->Exec(Form("mkdir -v %s", dirName.Data()));
  
  PlotHists(inFile, outFile, "trkTrueEnergy", "True Track Energy", "True Track Energy", histYaxis, nBins, 0, 2000, "trkTrueEnergy_pandoraTrack", dirName);
  PlotHists(inFile, outFile, "pfpTrackHitsPurity", "Track Hit Purity", "Purity", histYaxis, nBins, 0, 1, "PFPTrackHitsPurity_pandoraTrack", dirName);
  PlotHists(inFile, outFile, "pfpTrackHitsComp", "Track Hit Completeness", "Completeness", histYaxis, nBins, 0, 1, "PFPTrackHitsCompl_pandoraTrack", dirName);
  PlotHists(inFile, outFile, "trkDirDiff", "Track Start Direction Difference", "X", histYaxis, nBins, -1, 1, "trkStratDirDiff_pandoraTrack", dirName);
  PlotHists(inFile, outFile, "trkEndDirDiff", "Track End Direction Difference", "X", histYaxis, nBins, -1, 1, "trkEndDirDiff_pandoraTrack", dirName);
  PlotHists(inFile, outFile, "trkStartDist", "Track Start Position Difference", "Y", histYaxis, nBins, 0, 60, "trkStartDist_pandoraTrack", dirName);
  //Track energy information
  PlotHists(inFile, outFile, "trkdEdxBP", "Track Reconstruced dE/dx for Best Plane", "dE/dx (MeV/cm)", histYaxis, nBins, 0, 20, "trkdEdx_bp_pandoraTrack", dirName);
  PlotHists(inFile, outFile, "trkdEdxU", "Track Reconstruced dE/dx for U Plane", "dE/dx (MeV/cm)", histYaxis, nBins, 0, 20, "trkdEdx_u_pandoraTrack", dirName);
  PlotHists(inFile, outFile, "trkdEdxV", "Track Reconstruced dE/dx for V Plane", "dE/dx (MeV/cm)", histYaxis, nBins, 0, 20, "trkdEdx_v_pandoraTrack", dirName);
  PlotHists(inFile, outFile, "trkdEdxY", "Track Reconstruced dE/dx for Y Plane", "dE/dx (MeV/cm)", histYaxis, nBins, 0, 20, "trkdEdx_y_pandoraTrack", dirName);
  //Track kinetic energy information 
  PlotHists(inFile, outFile, "trkKEBP", "Track Reconstructed Kinetic Energy for the Best Plane", "Kinetic Energy(MeV)", histYaxis, nBins, 0, 1500, "trkKEBP_pandoraTrack", dirName);
  PlotHists(inFile, outFile, "trkKEU", "Track Reconstructed Kinetic Energy for U Plane", "Kinetic Energy(MeV)", histYaxis, nBins, 0, 1500, "trkKEU_pandoraTrack", dirName);
  PlotHists(inFile, outFile, "trkKEV", "Track Reconstructed Kinetic Energy for V Plane", "Kinetic Energy(MeV)", histYaxis, nBins, 0, 1500, "trkKEV_pandoraTrack", dirName);
  PlotHists(inFile, outFile, "trkKEY", "Track Reconstructed Kinetic Energy for Y Plane", "Kinetic Energy(MeV)", histYaxis, nBins, 0, 1500, "trkKEY_pandoraTrack", dirName);
  PlotHists(inFile, outFile, "trkTrueKinEnergy", "Track True Kinetic Energy", "Energy(MeV)", histYaxis, nBins, 0, 1500, "trkTrueKinEnergy_pandoraTrack", dirName);
  //Track length infromation
  PlotHists(inFile, outFile, "trkTrueLength", "Track True Length", "Length(cm)", histYaxis, nBins, 0, 200, "trkTrueTrackLength_pandoraTrack", dirName);
  PlotHists(inFile, outFile, "trkRecoLength", "Track Reconstructed Length", "Length(cm)", histYaxis, nBins, 0, 200, "trkLength_pandoraTrack", dirName);
  //Track momentum information 
  //Momentum by range plot
  PlotHists(inFile, outFile, "trkRecoMomentumByRange", "Track Reconstructed Momentum by Range", "Momentum by Range (GeV)", histYaxis, nBins, 0, 2, "recotrk_range_p_pandoraTrack", dirName);
  //Momentum coulomb scattering; Best Momentum
  PlotHists(inFile, outFile, "trkRecoMomMCS", "Track Reconstructed Momentum MCS", "Momentum MCS (GeV)", histYaxis, nBins, 0, 2, "recotrk_bestMomentum_p_pandoraTrack", dirName);
  //
  //Track dEdx vector and residual range vector plots
  PlotHists(inFile, outFile, "trkRecodEdxVecBP", "Track Reconstructed dEdx vector", "dEdx (MeV/cm)", histYaxis, nBins, 0, 15, "trkdEdxVecBP_pandoraTrack", dirName);
  PlotHists(inFile, outFile, "trkRecoResidualRangeVecBP", "Track Reconstructed Residual Range", "Residual Range(cm)", histYaxis, nBins, 0, 50, "trkResRangeBP_pandoraTrack", dirName);
  //2D histograms
  //Track dEdx vs Residulal range for 3 planes and Best Plane
  Plot2DHists(inFile, outFile, "dEdx_vs_Residual_Range_BP", "dEdx vs Residual Range for Best Plane", "Residual Range(cm)","dEdx (MeV/cm)", nBins, 0, 50, 0, 30,"trkResRangeBP_pandoraTrack","trkdEdxVecBP_pandoraTrack",dirName);
  Plot2DHists(inFile, outFile, "dEdx_vs_Residual_Range_U", "dEdx vs Residual Range for U plane", "Residual Range(cm)","dEdx (MeV/cm)", nBins, 0, 50, 0, 30,"trkResRangeY_pandoraTrack","trkdEdxVecY_pandoraTrack",dirName);
  Plot2DHists(inFile, outFile, "dEdx_vs_Residual_Range_V", "dEdx vs Residual Range for V Plane", "Residual Range(cm)","dEdx (MeV/cm)", nBins, 0, 50, 0, 30,"trkResRangeV_pandoraTrack","trkdEdxVecV_pandoraTrack",dirName);
  Plot2DHists(inFile, outFile, "dEdx_vs_Residual_Range_Y", "dEdx vs Residual Range for Y Plane", "Residual Range(cm)","dEdx (MeV/cm)", nBins, 0, 50, 0, 30,"trkResRangeY_pandoraTrack","trkdEdxVecY_pandoraTrack",dirName);
 
}

void PlotHists(TFile *inFile, TFile *outFile, TString histName, TString histTitle, TString histXaxis, TString histYaxis, int nBins, int xMin, int xMax, TString metricToPlot, TString dirName){

  TTree* metricTree;
  inFile->GetObject("trackvalidation/MetricTree", metricTree);
  std::cout<<"metricToPlot: " << metricToPlot << ". histName: " << histName << std::endl;
  
  for(auto trkpar : trkparticles){
    TString BraHisName = trkpar + "_" + histName;
    TTree* BranchTree;
    TH1F* hist   = new TH1F(BraHisName, (histTitle + " for " + trkpar)+";"+histXaxis+";"+histYaxis, nBins, xMin, xMax);

    if(trkpar == "All"){
      BranchTree = metricTree->CloneTree();
      BranchTree->Draw(metricToPlot + " >> " + BraHisName);
    }
    if(trkpar == "Muon"){
      BranchTree = metricTree->CopyTree("trkPDG_pandoraTrack==13");
      BranchTree->Draw(metricToPlot + " >> " + BraHisName);
    }
    if(trkpar == "Pion"){
      BranchTree = metricTree->CopyTree("trkPDG_pandoraTrack==211");
      BranchTree->Draw(metricToPlot + " >> " + BraHisName);
    }
    if(trkpar == "Proton"){
      BranchTree = metricTree->CopyTree("trkPDG_pandoraTrack==2221");
      BranchTree->Draw(metricToPlot + " >> " + BraHisName);
    }
    if(trkpar == "Other"){
      BranchTree = metricTree->CopyTree("trkPDG_pandoraTrack!=13 && trkPDG_pandoraTrack!=211 && trkPDG_pandoraTrack!=2212");
      BranchTree->Draw(metricToPlot + " >> " + BraHisName);
    }

    gStyle->SetOptStat(0);
    
    outFile->cd();
    outFile->cd(dirName.Data());
    
    hist->Write("", TObject::kOverwrite);
    TCanvas *c1 = new TCanvas(trkpar,trkpar,800,1000);
    hist->Draw();
    
    c1->Print(Form("%s/%s.png", dirName.Data(), BraHisName.Data()));
    
    std::cout<<"############### End of the loop over track particles ###############\n";
    delete c1;
    delete hist;
    delete BranchTree;

  }
  delete metricTree;
  return;
}

void Plot2DHists(TFile *inFile, TFile *outFile, TString histName, TString histTitle, TString histXaxis, TString histYaxis, int nBins, int xMin, int xMax, int yMin, int yMax, TString metricToPlotonX, TString metricToPlotonY, TString dirName){
  
  TTree* metricTree;
  inFile->GetObject("trackvalidation/MetricTree", metricTree);
  std::cout<<"metricToPlot on X: " << metricToPlotonX << " metricToPlot on Y: " << metricToPlotonY <<". histName: " << histName << std::endl;

  for(auto trkpar : trkparticles){
    TString BraHisName = trkpar + "_" + histName;
    TTree* BranchTree;
    
    TH2F* hist2d = new TH2F(BraHisName, (histTitle + " for " + trkpar)+";"+histXaxis+";"+histYaxis, nBins, xMin, xMax, nBins, yMin, yMax);

    if(trkpar == "All"){
      BranchTree = metricTree->CloneTree();
      BranchTree->Draw(metricToPlotonY + ":" + metricToPlotonX + " >> " + BraHisName);
    }
    if(trkpar == "Muon"){
      BranchTree = metricTree->CopyTree("trkPDG_pandoraTrack==13");
      BranchTree->Draw(metricToPlotonY + ":" + metricToPlotonX + " >> " + BraHisName);
    }
    if(trkpar == "Pion"){
      BranchTree = metricTree->CopyTree("trkPDG_pandoraTrack==211");
      BranchTree->Draw(metricToPlotonY + ":" + metricToPlotonX + " >> " + BraHisName);
    }
    if(trkpar == "Proton"){
      BranchTree = metricTree->CopyTree("trkPDG_pandoraTrack==2221");
      BranchTree->Draw(metricToPlotonY + ":" + metricToPlotonX + " >> " + BraHisName);
    }
    if(trkpar == "Other"){
      BranchTree = metricTree->CopyTree("trkPDG_pandoraTrack!=13 && trkPDG_pandoraTrack!=211 && trkPDG_pandoraTrack!=2212");
      BranchTree->Draw(metricToPlotonY + ":" + metricToPlotonX + " >> " + BraHisName);
    }
    
    gStyle->SetOptStat(0);

    outFile->cd();
    outFile->cd(dirName.Data());
    hist2d->Write("", TObject::kOverwrite);
    TCanvas *c2 = new TCanvas("c2","c2",800,1000);
    hist2d->Draw("colz");

    c2->Print(Form("%s/%s.png", dirName.Data(), BraHisName.Data()));
    delete c2;
    delete hist2d;
    delete BranchTree;

  }

  delete metricTree;
  return;
}



