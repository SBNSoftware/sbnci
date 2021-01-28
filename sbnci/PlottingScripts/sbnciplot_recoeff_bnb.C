#include "TTree.h"
#include "TFile.h"
#include "TH1F.h"
#include "TString.h"

void Plot(TTree *tree, TFile* outputFile, TString directoryName, TString PDG, TString variable, 
	  int nBins, float xMin, float xMax, TString title, bool goodReco = false);

void sbnciplot_recoeff_bnb(TString inputFileName)
{
  TFile *inputFile = TFile::Open(inputFileName.Data());
  TTree *tree = (TTree*) inputFile->Get("recoEff/ParticleTree");

  TFile *outputFile = new TFile("ci_validation_histos.root","RECREATE");

  TString directoryName = "efficiencyPlots";

  outputFile->cd();
  gDirectory->mkdir(directoryName);
  gSystem->Exec(Form("mkdir -v %s",directoryName.Data()));

  //Muon Plots
  Plot(tree,outputFile,directoryName,"13","mc_theta_xz",30,-180,180,";True #theta_{xz} (#circ);Fraction;");
  Plot(tree,outputFile,directoryName,"13","mc_theta_yz",30,-180,180,";True #theta_{yz} (#circ);Fraction;");
  Plot(tree,outputFile,directoryName,"13","mc_theta_xy",30,-180,180,";True #theta_{xy} (#circ);Fraction;");
  Plot(tree,outputFile,directoryName,"13","mc_momentum",30,0,1.5,";True p (GeV/c);Fraction;");
  Plot(tree,outputFile,directoryName,"13","mc_energy0",50,0,2,";True E (GeV);Fraction;");
  Plot(tree,outputFile,directoryName,"13","mc_length",40,0,400,";True track length (cm);Fraction;");

  //Muon High Quality Plots
  Plot(tree,outputFile,directoryName,"13","mc_theta_xz",30,-180,180,";True #theta_{xz} (#circ);Fraction;",true);
  Plot(tree,outputFile,directoryName,"13","mc_theta_yz",30,-180,180,";True #theta_{yz} (#circ);Fraction;",true);
  Plot(tree,outputFile,directoryName,"13","mc_theta_xy",30,-180,180,";True #theta_{xy} (#circ);Fraction;",true);
  Plot(tree,outputFile,directoryName,"13","mc_momentum",30,0,1.5,";True p (GeV/c);Fraction;",true);
  Plot(tree,outputFile,directoryName,"13","mc_energy0",50,0,2,";True E (GeV);Fraction;",true);
  Plot(tree,outputFile,directoryName,"13","mc_length",40,0,400,";True track length (cm);Fraction;",true);


  // Proton Plots
  Plot(tree,outputFile,directoryName,"2212","mc_theta_xz",30,-180,180,";True #theta_{xz} (#circ);Fraction;");
  Plot(tree,outputFile,directoryName,"2212","mc_theta_yz",30,-180,180,";True #theta_{yz} (#circ);Fraction;");
  Plot(tree,outputFile,directoryName,"2212","mc_theta_xy",30,-180,180,";True #theta_{xy} (#circ);Fraction;");
  Plot(tree,outputFile,directoryName,"2212","mc_momentum",30,0,1.5,";True p (GeV/c);Fraction;");
  Plot(tree,outputFile,directoryName,"2212","mc_energy0",50,0,2,";True E (GeV);Fraction;");
  Plot(tree,outputFile,directoryName,"2212","mc_length",50,0,50,";True track length (cm);Fraction;");

  // Proton High Quality Plots
  Plot(tree,outputFile,directoryName,"2212","mc_theta_xz",30,-180,180,";True #theta_{xz} (#circ);Fraction;",true);
  Plot(tree,outputFile,directoryName,"2212","mc_theta_yz",30,-180,180,";True #theta_{yz} (#circ);Fraction;",true);
  Plot(tree,outputFile,directoryName,"2212","mc_theta_xy",30,-180,180,";True #theta_{xy} (#circ);Fraction;",true);
  Plot(tree,outputFile,directoryName,"2212","mc_momentum",30,0,1.5,";True p (GeV/c);Fraction;",true);
  Plot(tree,outputFile,directoryName,"2212","mc_energy0",50,0,2,";True E (GeV);Fraction;",true);
  Plot(tree,outputFile,directoryName,"2212","mc_length",50,0,50,";True track length (cm);Fraction;",true);


  //Pi Plus Plots
  Plot(tree,outputFile,directoryName,"211","mc_theta_xz",30,-180,180,";True #theta_{xz} (#circ);Fraction;");
  Plot(tree,outputFile,directoryName,"211","mc_theta_yz",30,-180,180,";True #theta_{yz} (#circ);Fraction;");
  Plot(tree,outputFile,directoryName,"211","mc_theta_xy",30,-180,180,";True #theta_{xy} (#circ);Fraction;");
  Plot(tree,outputFile,directoryName,"211","mc_momentum",30,0,1.5,";True p (GeV/c);Fraction;");
  Plot(tree,outputFile,directoryName,"211","mc_energy0",50,0,2,";True E (GeV);Fraction;");
  Plot(tree,outputFile,directoryName,"211","mc_length",40,0,100,";True track length (cm);Fraction;");

  //Pi Plus High Quality Plots
  Plot(tree,outputFile,directoryName,"211","mc_theta_xz",30,-180,180,";True #theta_{xz} (#circ);Fraction;",true);
  Plot(tree,outputFile,directoryName,"211","mc_theta_yz",30,-180,180,";True #theta_{yz} (#circ);Fraction;",true);
  Plot(tree,outputFile,directoryName,"211","mc_theta_xy",30,-180,180,";True #theta_{xy} (#circ);Fraction;",true);
  Plot(tree,outputFile,directoryName,"211","mc_momentum",30,0,1.5,";True p (GeV/c);Fraction;",true);
  Plot(tree,outputFile,directoryName,"211","mc_energy0",50,0,2,";True E (GeV);Fraction;",true);
  Plot(tree,outputFile,directoryName,"211","mc_length",40,0,100,";True track length (cm);Fraction;",true);


  //Pi Minus Plots
  Plot(tree,outputFile,directoryName,"-211","mc_theta_xz",30,-180,180,";True #theta_{xz} (#circ);Fraction;");
  Plot(tree,outputFile,directoryName,"-211","mc_theta_yz",30,-180,180,";True #theta_{yz} (#circ);Fraction;");
  Plot(tree,outputFile,directoryName,"-211","mc_theta_xy",30,-180,180,";True #theta_{xy} (#circ);Fraction;");
  Plot(tree,outputFile,directoryName,"-211","mc_momentum",30,0,1.5,";True p (GeV/c);Fraction;");
  Plot(tree,outputFile,directoryName,"-211","mc_energy0",50,0,2,";True E (GeV);Fraction;");
  Plot(tree,outputFile,directoryName,"-211","mc_length",40,0,100,";True track length (cm);Fraction;");

  //Pi Minus High Quality Plots
  Plot(tree,outputFile,directoryName,"-211","mc_theta_xz",30,-180,180,";True #theta_{xz} (#circ);Fraction;",true);
  Plot(tree,outputFile,directoryName,"-211","mc_theta_yz",30,-180,180,";True #theta_{yz} (#circ);Fraction;",true);
  Plot(tree,outputFile,directoryName,"-211","mc_theta_xy",30,-180,180,";True #theta_{xy} (#circ);Fraction;",true);
  Plot(tree,outputFile,directoryName,"-211","mc_momentum",30,0,1.5,";True p (GeV/c);Fraction;",true);
  Plot(tree,outputFile,directoryName,"-211","mc_energy0",50,0,2,";True E (GeV);Fraction;",true);
  Plot(tree,outputFile,directoryName,"-211","mc_length",40,0,100,";True track length (cm);Fraction;",true);



  outputFile->Close();

  delete inputFile, outputFile, tree;
}
  
void Plot(TTree* tree, TFile* outputFile, TString directoryName, TString PDG, TString variable, 
	  int nBins, float xMin, float xMax, TString title, bool goodReco = false)
{
  TString recoDef = "";
  TString histName = "";
  TString short_variable(variable);
  short_variable.Remove(0,3);

  if(PDG == "13"){
    histName += "muon";
    recoDef += "reco_nTracks > 0";
    if(goodReco) recoDef += " && reco_track_completeness > 0.9 && reco_track_purity > 0.9";
  }
  else if(PDG == "2212"){
    histName += "proton";
    recoDef += "reco_nTracks > 0";
    if(goodReco) recoDef += " && reco_track_completeness > 0.8 && reco_track_purity > 0.8";
  }
  else if(PDG == "211"){
    histName += "piplus";
    recoDef += "reco_nTracks > 0";
    if(goodReco) recoDef += " && reco_track_completeness > 0.8 && reco_track_purity > 0.8";
  }
  else if(PDG == "-211"){
    histName += "piminus";
    recoDef += "reco_nTracks > 0";
    if(goodReco) recoDef += " && reco_track_completeness > 0.8 && reco_track_purity > 0.8";
  }

  histName += "_";
  histName += short_variable;
  if(goodReco) histName += "_high_quality";

  TH1F *tHist = new TH1F("tHist","",nBins,xMin,xMax);
  TH1F *rHist = new TH1F(histName,"",nBins,xMin,xMax);

  tree->Draw(variable + ">>tHist","mc_PDG == " + PDG);
  tree->Draw(variable + ">>" + histName,"mc_PDG == " + PDG + " && " + recoDef);

  TEfficiency eff(*rHist,*tHist);
  eff.SetTitle(title);
  eff.SetName(histName);

  outputFile->cd();
  outputFile->cd(directoryName.Data());

  eff.Write(histName,TObject::kOverwrite);

  TCanvas *canvas = new TCanvas("canvas","canvas");
  eff.Draw();
  canvas->Print(Form("%s/%s.png",directoryName.Data(),histName.Data()));

  delete tHist, rHist, canvas;
}
