#include "TTree.h"
#include "TFile.h"
#include "TH1F.h"
#include "TString.h"

void EffPlot(TTree *tree, TFile* outputFile, TString directoryName, TString PDG, TString variable, 
	     int nBins, float xMin, float xMax, TString title, bool goodReco = false);

void Plot(TTree* tree, TFile* outputFile, TString directoryName, TString PDG, TString variable, 
          int nBins, float xMin, float xMax, TString title, TString variableName);

void sbnciplot_recoeff(TString inputFileName)
 {
  TFile *inputFile = TFile::Open(inputFileName.Data());
  TTree *tree = (TTree*) inputFile->Get("recoeff/ParticleTree");

  TFile *outputFile = new TFile("ci_validation_histos.root","RECREATE");
  gStyle->SetOptStat(0);
  gROOT->ForceStyle();


  outputFile->cd();
  gDirectory->mkdir("muon");
  gSystem->Exec(Form("mkdir -v %s","muon"));
  gDirectory->mkdir("electron");
  gSystem->Exec(Form("mkdir -v %s","electron"));
  gDirectory->mkdir("proton");
  gSystem->Exec(Form("mkdir -v %s","proton"));

  //Muon Plots
  EffPlot(tree,outputFile,"muon","13","mc_theta_xz",30,-90,90,";True #theta_{xz} (#circ);Fraction;");
  EffPlot(tree,outputFile,"muon","13","mc_theta_yz",30,-90,90,";True #theta_{yz} (#circ);Fraction;");
  EffPlot(tree,outputFile,"muon","13","mc_theta_xy",30,-90,90,";True #theta_{xy} (#circ);Fraction;");
  EffPlot(tree,outputFile,"muon","13","mc_momentum",30,0,1.5,";True p (GeV/c);Fraction;");
  EffPlot(tree,outputFile,"muon","13","mc_energy0",50,0,1.2,";True E (GeV);Fraction;");
  EffPlot(tree,outputFile,"muon","13","mc_length",40,0,400,";True track length (cm);Fraction;");
  Plot(tree,outputFile,"muon","13","reco_track_completeness",40,0,1,";Track completeness;Entries","completeness");
  Plot(tree,outputFile,"muon","13","reco_track_purity",40,0,1,";Track purity;Entries","purity");
  Plot(tree,outputFile,"muon","13","(reco_track_length - mc_length)/mc_length",40,-1,1,";#frac{reco - true}{true} track length (cm);Entries","length_metric");

  //Muon High Quality Plots
  EffPlot(tree,outputFile,"muon","13","mc_theta_xz",30,-90,90,";True #theta_{xz} (#circ);Fraction;",true);
  EffPlot(tree,outputFile,"muon","13","mc_theta_yz",30,-90,90,";True #theta_{yz} (#circ);Fraction;",true);
  EffPlot(tree,outputFile,"muon","13","mc_theta_xy",30,-90,90,";True #theta_{xy} (#circ);Fraction;",true);
  EffPlot(tree,outputFile,"muon","13","mc_momentum",30,0,1.5,";True p (GeV/c);Fraction;",true);
  EffPlot(tree,outputFile,"muon","13","mc_energy0",50,0,1.2,";True E (GeV);Fraction;",true);
  EffPlot(tree,outputFile,"muon","13","mc_length",40,0,400,";True track length (cm);Fraction;",true);


  //Electron Plots
  EffPlot(tree,outputFile,"electron","11","mc_theta_xz",30,-90,90,";True #theta_{xz} (#circ);Fraction;");
  EffPlot(tree,outputFile,"electron","11","mc_theta_yz",30,-90,90,";True #theta_{yz} (#circ);Fraction;");
  EffPlot(tree,outputFile,"electron","11","mc_theta_xy",30,-90,90,";True #theta_{xy} (#circ);Fraction;");
  EffPlot(tree,outputFile,"electron","11","mc_momentum",30,0,2,";True p (GeV/c);Fraction;");
  EffPlot(tree,outputFile,"electron","11","mc_energy0",30,0,2,";True E (GeV);Fraction;");
  Plot(tree,outputFile,"electron","11","reco_shower_completeness",40,0,1,";Shower completeness;Entries","completeness");
  Plot(tree,outputFile,"electron","11","reco_shower_purity",40,0,1,";Shower purity;Entries","purity");

  //Electron High Quality Plots
  EffPlot(tree,outputFile,"electron","11","mc_theta_xz",30,-90,90,";True #theta_{xz} (#circ);Fraction;",true);
  EffPlot(tree,outputFile,"electron","11","mc_theta_yz",30,-90,90,";True #theta_{yz} (#circ);Fraction;",true);
  EffPlot(tree,outputFile,"electron","11","mc_theta_xy",30,-90,90,";True #theta_{xy} (#circ);Fraction;",true);
  EffPlot(tree,outputFile,"electron","11","mc_momentum",30,0,2,";True p (GeV/c);Fraction;",true);
  EffPlot(tree,outputFile,"electron","11","mc_energy0",30,0,2,";True E (GeV);Fraction;",true);


  // Proton Plots
  EffPlot(tree,outputFile,"proton","2212","mc_theta_xz",30,-90,90,";True #theta_{xz} (#circ);Fraction;");
  EffPlot(tree,outputFile,"proton","2212","mc_theta_yz",30,-90,90,";True #theta_{yz} (#circ);Fraction;");
  EffPlot(tree,outputFile,"proton","2212","mc_theta_xy",30,-90,90,";True #theta_{xy} (#circ);Fraction;");
  EffPlot(tree,outputFile,"proton","2212","mc_momentum",30,0,1.2,";True p (GeV/c);Fraction;");
  EffPlot(tree,outputFile,"proton","2212","mc_energy0",50,0,1.5,";True E (GeV);Fraction;");
  EffPlot(tree,outputFile,"proton","2212","mc_length",50,0,40,";True track length (cm);Fraction;");
  Plot(tree,outputFile,"proton","2212","reco_track_completeness",40,0,1,";Track completeness;Entries","completeness");
  Plot(tree,outputFile,"proton","2212","reco_track_purity",40,0,1,";Track purity;Entries","purity");
  Plot(tree,outputFile,"proton","2212","(reco_track_length - mc_length)/mc_length",40,-1,1,";#frac{reco - true}{true} track length (cm);Entries","length_metric");

  // Proton High Quality Plots
  EffPlot(tree,outputFile,"proton","2212","mc_theta_xz",30,-90,90,";True #theta_{xz} (#circ);Fraction;",true);
  EffPlot(tree,outputFile,"proton","2212","mc_theta_yz",30,-90,90,";True #theta_{yz} (#circ);Fraction;",true);
  EffPlot(tree,outputFile,"proton","2212","mc_theta_xy",30,-90,90,";True #theta_{xy} (#circ);Fraction;",true);
  EffPlot(tree,outputFile,"proton","2212","mc_momentum",30,0,1.2,";True p (GeV/c);Fraction;",true);
  EffPlot(tree,outputFile,"proton","2212","mc_energy0",50,0,1.5,";True E (GeV);Fraction;",true);
  EffPlot(tree,outputFile,"proton","2212","mc_length",50,0,40,";True track length (cm);Fraction;",true);


  outputFile->Close();

  delete inputFile, outputFile, tree;
 }

 void EffPlot(TTree* tree, TFile* outputFile, TString directoryName, TString PDG, TString variable, 
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
   else if(PDG == "11"){
     histName += "electron";
     recoDef += "reco_nShowers > 0";
     if(goodReco) recoDef += " && reco_shower_completeness > 0.8 && reco_shower_purity > 0.8";
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

void Plot(TTree* tree, TFile* outputFile, TString directoryName, TString PDG, TString variable, 
	  int nBins, float xMin, float xMax, TString title, TString variableName)
{
  TString recoDef = "";
  TString histName = "";

  if(PDG == "13"){
    histName += "muon";
    recoDef += "reco_nTracks > 0";
  }
  else if(PDG == "11"){
    histName += "electron";
    recoDef += "reco_nShowers > 0";
  }
  else if(PDG == "2212"){
    histName += "proton";
    recoDef += "reco_nTracks > 0";
  }
  else if(PDG == "211"){
    histName += "piplus";
    recoDef += "reco_nTracks > 0";
  }
  else if(PDG == "-211"){
    histName += "piminus";
    recoDef += "reco_nTracks > 0";
  }

  histName += "_";
  histName += variableName;

  TH1F *hist = new TH1F(histName,title,nBins,xMin,xMax);

  tree->Draw(variable + ">>" + histName,"mc_PDG == " + PDG + " && " + recoDef);

  outputFile->cd();
  outputFile->cd(directoryName.Data());

  hist->Write(histName,TObject::kOverwrite);

  TCanvas *canvas = new TCanvas("canvas","canvas");
  hist->Draw();
  canvas->Print(Form("%s/%s.png",directoryName.Data(),histName.Data()));

  delete hist, canvas;
}
