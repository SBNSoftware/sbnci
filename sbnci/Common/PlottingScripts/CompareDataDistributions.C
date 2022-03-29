//Shamelessly stolen from Vito's script in lar_ci/test and modified to suit sbn needs.

#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include "TAxis.h"
#include "TCanvas.h"
#include "TChain.h"
#include "TEfficiency.h"
#include "TFile.h"
#include "TKey.h"
#include "TH1.h"
#include "TLegend.h"
#include "TPaveText.h"
#include "TROOT.h"
#include "TString.h"
#include "TStyle.h"
#include "TSystem.h"

void LoopHistos(vector<uint> &folder_size, vector<string> &list_folder, vector<string> &list_histo, TString InputHistoFile) {

  string word;

  TFile *f1 = TFile::Open(InputHistoFile.Data());
  f1->ls();
  cout << f1->GetListOfKeys()->GetSize() << endl;

  const Int_t nBranches=f1->GetListOfKeys()->GetSize();
  cout << endl;
  TDirectory *savdir = gDirectory;
  TIter next_key(f1->GetListOfKeys());
  TKey *key_dir;
  TKey *key_histo;
  while ((key_dir = (TKey*)next_key())) {
    if ( key_dir->IsFolder() ){
      f1->cd();
      f1->cd(key_dir->GetName());
      cout << gDirectory->GetListOfKeys()->GetSize() << endl;
      const Int_t nHistos=gDirectory->GetListOfKeys()->GetSize();
      TString s_branch_names[nHistos];
      cout << key_dir->GetName() << endl;
      list_folder.push_back(key_dir->GetName());
      TIter next_key_histo(gDirectory->GetListOfKeys());
      UInt_t idx=0;
      while ((key_histo = (TKey*)next_key_histo())) {
        TClass *cl = gROOT->GetClass(key_histo->GetClassName());
        if (! (cl->InheritsFrom("TH1") || cl->InheritsFrom("TEfficiency"))  ) continue;
        cout << idx++ << " " << key_dir->GetName() << "/" << key_histo->GetName() << endl;
        list_histo.push_back(key_histo->GetName());
      }
      folder_size.push_back(idx);
    }
    cout << endl;
  }

}

double calculateChiSqDistance(TH1D* O, TH1D* E){

    double chisq = 0;
    for (int i = 1; i < O->GetNbinsX()+1; i++){
        double O_i = O->GetBinContent(i);
        double E_i = E->GetBinContent(i);
        double O_ierr = O->GetBinError(i);
        double E_ierr = E->GetBinError(i);
        if (O_i == 0 || E_i == 0){
            chisq += 0;
        }
        else{
            chisq += std::pow(O_i - E_i,2)/(std::sqrt(std::pow(O_ierr,2) + std::pow(E_ierr,2)));
//             chisq += std::pow(O_i - E_i,2)/(E_i);
        }
    }

    if ( O->Integral()+O->GetBinContent(0)+O->GetBinContent(O->GetNbinsX()+1) > 0)
      return chisq/(O->Integral()+O->GetBinContent(0)+O->GetBinContent(O->GetNbinsX()+1));
    else
      return -1;
}

void setStyle(TH1* h, int dhc/*, TString histoYAxis*/) {

    gStyle->SetOptStat(0);
    h->GetYaxis()->SetTitleOffset(1.4);
    // nom data
    if ( dhc == 0 ) {
        h->SetLineColor(kBlack);
        h->SetMarkerStyle(20);
        h->SetMarkerColor(kBlack);
        h->SetMarkerSize(0.7);
//         h->GetYaxis()->SetTitle(histoYAxis);
        h->SetLineWidth(2);
        h->GetXaxis()->SetTitleSize(0);
        h->GetXaxis()->SetLabelSize(0);
    }
    // nom mc
    else if ( dhc == 3 ) {
        h->SetLineColor(kOrange+10);
//         h->GetYaxis()->SetTitle(histoYAxis);
        h->SetFillColorAlpha(kOrange+6, 0.5);
        h->SetMarkerColor(2);
        h->SetLineWidth(2);
        h->GetXaxis()->SetTitleSize(0);
        h->GetXaxis()->SetLabelSize(0);
    }
    // alt data
    else if ( dhc == 2 ) {
        h->SetLineColor(kOrange+10);
        h->SetMarkerStyle(20);
        h->SetMarkerColor(kOrange+10);
        h->SetMarkerSize(0.7);
//         h->GetYaxis()->SetTitle(histoYAxis);
        h->SetLineWidth(2);
        h->GetXaxis()->SetTitleSize(0);
        h->GetXaxis()->SetLabelSize(0);
    }
    // alt mc
    else if ( dhc == 1 ) {
        h->SetLineColor(kBlack);
//         h->GetYaxis()->SetTitle(histoYAxis);
        h->SetFillColorAlpha(kAzure+6, 0.5);
        h->SetMarkerColor(2);
        h->SetLineWidth(2);
        h->GetXaxis()->SetTitleSize(0);
        h->GetXaxis()->SetLabelSize(0);
    }
    gStyle->SetTitleX(0.15);

}


double getMax(TH1* hData, TH1* hMC) {
    double datamax = hData->GetMaximum() * 1.15;
    double mcmax = hMC->GetMaximum() * 1.15;
    if (datamax > mcmax) return datamax;
    else if (mcmax > datamax) return mcmax;
    else return datamax;
}

void getHistMaxMinWithError(TH1* hist, double &histMax, double &histMin) {
  double maxi = -std::numeric_limits<double>::max();
  double mini = std::numeric_limits<double>::max();

  for(int i=1; i<hist->GetNbinsX()+1; ++i)
    {
      maxi = std::max(maxi, hist->GetBinContent(i) + hist->GetBinErrorUp(i));
      mini = std::min(mini, hist->GetBinContent(i) - hist->GetBinErrorLow(i));
    }
  histMax = maxi;
  histMin = mini;
}

void setStyleRatio(TH1D* h, TString file1_label, TString file2_label){
    h->SetTitle("");
    h->SetNdivisions(510, "Y");
    h->GetYaxis()->SetTitle("#frac{("+file2_label+")-("+file1_label+")}{("+file1_label+")}");
    h->GetYaxis()->SetTitleOffset(0.57);
    h->GetYaxis()->SetTitleSize(0.07);
    h->GetYaxis()->SetLabelSize(0.075);
    h->GetXaxis()->SetTitleSize(0.09);
    h->GetXaxis()->SetLabelSize(0.075);
}

void setLegend(TH1* hfile1, int file1_dmc, TString file1_label, TH1* hfile2, int file2_dmc, TString file2_label)
{

  if (file1_dmc == 0 && file2_dmc == 1){
        TLegend* legfile1 = new TLegend(0.40, 0.91, 0.70, 1.0);
        legfile1->AddEntry(hfile1, file1_label, "lep");
        legfile1->SetBorderSize(0);
        legfile1->SetTextAlign(12);
        legfile1->SetTextSize(0.05);
        legfile1->Draw("same");
        TLegend* legfile2 = new TLegend(0.70, 0.91, 1.00, 1.0);
        legfile2->AddEntry(hfile2, file2_label);
        legfile2->SetTextAlign(12);
        legfile2->SetTextSize(0.05);
        legfile2->SetBorderSize(0);
        legfile2->Draw("same");
    }
    if (file1_dmc == 3 && file2_dmc == 1){
        TLegend* legfile1 = new TLegend(0.40, 0.91, 0.70, 1.0);
        legfile1->AddEntry(hfile1, file1_label);
        legfile1->SetBorderSize(0);
        legfile1->SetTextAlign(12);
        legfile1->SetTextSize(0.05);
        legfile1->Draw("same");
        TLegend* legfile2 = new TLegend(0.700, 0.91, 1.00, 1.0);
        legfile2->AddEntry(hfile2, file2_label);
        legfile2->SetTextAlign(12);
        legfile2->SetTextSize(0.05);
        legfile2->SetBorderSize(0);
        legfile2->Draw("same");
    }
    if (file1_dmc == 0 && file2_dmc == 2){
        TLegend* legfile1 = new TLegend(0.40, 0.91, 0.70, 1.0);
        legfile1->AddEntry(hfile1, file1_label, "lep");
        legfile1->SetBorderSize(0);
        legfile1->SetTextAlign(12);
        legfile1->SetTextSize(0.05);
        legfile1->Draw("same");
        TLegend* legfile2 = new TLegend(0.70, 0.91, 1.00, 1.0);
        legfile2->AddEntry(hfile2, file2_label, "lep");
        legfile2->SetTextAlign(12);
        legfile2->SetTextSize(0.05);
        legfile2->SetBorderSize(0);
        legfile2->Draw("same");
    }
}

int getNBins(TH1D *h){
    int nBins = 0;
    for (int i = 1; i < h->GetNbinsX()+1; i++){
        if (h->GetBinContent(i) !=0)
            nBins++;
    }
    return nBins;
}

map<string, vector<float>> thresholdmap(ifstream &infile){
    //build threshold map from file
    map<string, vector<float>> mapping;
    string line;
    string out1;

    vector<string> out;

    if (infile.is_open())
    {
        while ( getline (infile ,line))
        {
            line.erase(remove(line.begin(), line.end(), ' '), line.end());
            stringstream ss(line);
            while ( getline (ss ,out1, '#'))
                out.push_back(out1);
            if(out.size()==3 && out.at(0).size()>0){
                float v = stof(out.at(1));
                vector<float> th{stof(out.at(1)),stof(out.at(2))};
                mapping[out.at(0)] = th;
            }
            out.clear();
        }
    }

    return mapping;
}

void CompareDataDistributions(TString gCurVersion="v07_06_00", TString gRefVersion="v07_00_00", const bool useNormalised = true){

//   TString ref_xrootd_path=gSystem->Getenv("ref_dunetpc_ana_hist_xrootd_path");

  //set default chi2 thresholds
  float th1 = 0.25;
  float th2 = 0.5;
  // threshold file is retrived through env var CI_COMP_TH_FILE
  string thfile_env="CI_COMP_TH_FILE";
  char * thfile = getenv( thfile_env.c_str() );
  map<string, vector<float>> thmap;
  if (!thfile){
      cout << "*** WARNING no threshold file found with env var " << thfile_env << endl;
      cout << "\t default chi2 threshold values are used\n" << endl;
  }
  else{
    ifstream infile(thfile);
    thmap = thresholdmap(infile);
    if(thmap.size()==0){
        cout << "*** WARNING no thresholds found in " << thfile << endl;
        cout << "\t default chi2 threshold values are used\n" << endl;
    }
  }

  std::ofstream ComparisonChiSquare_file("ComparisonChiSquare.txt");

  vector<string> list_folder, list_histo;
  vector<uint> folder_size;
  LoopHistos(folder_size, list_folder, list_histo, Form("ci_validation_histos_testing_%s.root", gCurVersion.Data()));

  // For debug purpose - BEGIN
  cout << folder_size.size() << endl;
  cout << list_folder.size() << endl;
  cout << list_histo.size() << endl;

  TString histo_names="", folder_names="";
  UInt_t idx_offset_tmp=0;
  for (UInt_t idx=0; idx<folder_size.size(); idx++){
      folder_names=list_folder.at(idx);
      cout << endl << folder_names << endl;
      for (UInt_t idy=0; idy<folder_size.at(idx); idy++){
          histo_names=list_histo.at(idy+idx_offset_tmp);
          cout << histo_names << endl;
      }
      idx_offset_tmp+=folder_size.at(idx);

  }
  // For debug purpose - END


  TFile *RefFile = TFile::Open(Form("ci_validation_histos_ref_%s.root", gRefVersion.Data()));
  TFile *CurFile = TFile::Open(Form("ci_validation_histos_testing_%s.root", gCurVersion.Data()));

  UInt_t nBranches=list_folder.size();
  cout << "nBranches: " << nBranches << endl;

  UInt_t nHistos=list_histo.size();
  cout << "nHistos: " << nHistos << endl;


  TString s_histo_names="";
  TString s_branch_names="";
  TString s_local_branch_names="";
  UInt_t idx_offset=0;
  for(UInt_t idb=0; idb<nBranches; idb++){
      s_branch_names=list_folder.at(idb);
      s_local_branch_names = s_branch_names+"Compare";

      gSystem->Exec(Form("mkdir -pv %s", s_local_branch_names.Data()));

      nHistos = folder_size.at(idb);
      for (UInt_t idx=0; idx<nHistos; idx++){
          s_histo_names=list_histo.at(idx+idx_offset);
          TCanvas *c1 = new TCanvas(Form("c1_%d",idx),Form("%s_%s_ShwEff", s_branch_names.Data(), s_histo_names.Data()), 1000, 800);
          TPad *topPad = new TPad(Form("topPad_%d",idx), "", 0.005, 0.3, 0.995, 0.995);
          TPad *bottomPad = new TPad(Form("bottomPad_%d",idx), "", 0.005, 0.005, 0.995, 0.3);
          topPad->SetBottomMargin(0.02);
          bottomPad->SetTopMargin(0.0);
          bottomPad->SetBottomMargin(0.18);
          bottomPad->SetGridy();
          topPad->Draw();
          bottomPad->Draw();
          topPad->cd();

          TH1D* hRef = new TH1D();
          TH1D* hCur = new TH1D();

          if ( ! ((TString)(RefFile->Get(Form("%s/%s", s_branch_names.Data(), s_histo_names.Data())))->ClassName()).Contains("TEfficiency") ){

              hRef = (TH1D*)RefFile->Get(Form("%s/%s", s_branch_names.Data(), s_histo_names.Data()));
              hCur = (TH1D*)CurFile->Get(Form("%s/%s", s_branch_names.Data(), s_histo_names.Data()));

              if (useNormalised && hRef->Integral() > 0 && hCur->Integral() > 0)
                  hCur->Scale((hRef->Integral()+hRef->GetBinContent(0)+hRef->GetBinContent(hRef->GetNbinsX()+1))/(hCur->Integral()+hCur->GetBinContent(0)+hCur->GetBinContent(hCur->GetNbinsX()+1)));
	      
	      double maxext = getMax(hRef, hCur);
	      hRef->SetMaximum(maxext);
          }
          else{

              TEfficiency* hRefEff = (TEfficiency*)RefFile->Get(Form("%s/%s", s_branch_names.Data(), s_histo_names.Data()));
              TEfficiency* hCurEff = (TEfficiency*)CurFile->Get(Form("%s/%s", s_branch_names.Data(), s_histo_names.Data()));

              topPad->cd();

              TEfficiency *hRefEff_clone = (TEfficiency*)hRefEff->Clone("hRefEff_clone");
              TEfficiency *hCurEff_clone = (TEfficiency*)hCurEff->Clone("hCurEff_clone");

              TH1D*hp1=(TH1D*)hCurEff_clone->GetPassedHistogram();
              TH1D*ht1=(TH1D*)hCurEff_clone->GetTotalHistogram();
              hCur=(TH1D*)hCurEff_clone->GetTotalHistogram();
	      hCur->Sumw2();
              hCur->Divide(hp1, ht1, 1., 1., "");

              TH1D*hp2=(TH1D*)hRefEff_clone->GetPassedHistogram();
              TH1D*ht2=(TH1D*)hRefEff_clone->GetTotalHistogram();
              hRef=(TH1D*)hRefEff_clone->GetTotalHistogram();
	      hRef->Sumw2();
              hRef->Divide(hp2, ht2, 1., 1., "");


          }

          setStyle(hRef, 3);
          setStyle(hCur, 1);

          topPad->cd();
          // draw MC histo error bars...
          hRef->Draw("e2");
          TH1F* hRefc = (TH1F*)hRef->Clone("hRefc");
          hRefc->SetDirectory(0);
          hRefc->SetFillColor(0);
          hRefc->Draw("hist same");
          // and data
          hCur->Draw("e2same");
          // clone, and draw as histogram
          TH1F* hCurc = (TH1F*)hCur->Clone("hCurc");
          hCurc->SetDirectory(0);
          hCurc->SetFillColor(0);
          hCurc->Draw("hist same");

          TString RefLegend="Ref: "+gRefVersion;
          TString CurLegend="Cur: "+gCurVersion;
          setLegend(hRef, 3, RefLegend.Data(), hCur, 1, CurLegend.Data());

          double chisqv = calculateChiSqDistance(hCur, hRef);
          Color_t CanvasColor = kWhite;
          map<string, vector<float>>::const_iterator iter = thmap.find(Form("%s/%s", s_local_branch_names.Data(), s_histo_names.Data()));
          if(iter!=thmap.end()){
              vector<float> th = vector<float>(iter->second);
              th1 = th.at(0);
              th2 = th.at(1);
          }
          if (chisqv>th2) CanvasColor = kRed;
          else if (chisqv>th1) CanvasColor = kYellow;
          topPad->SetFillColor(CanvasColor);
          topPad->SetFrameFillColor(-1);
          bottomPad->SetFillColor(CanvasColor);
          bottomPad->SetFrameFillColor(-1);
          TString chisq = Form("#chi^{2}: %g", chisqv);
          int nBins = std::max(getNBins(hRef),getNBins(hCur));
          TString NDF = Form("No. Bins: %i", nBins);
          topPad->cd();
          TPaveText *pt = new TPaveText(0.35, 0.78, 0.75, 0.88, "NDC");
          pt->AddText(chisq);
          pt->AddText(NDF);
          pt->SetFillStyle(0);
          pt->SetBorderSize(0);
          pt->SetTextAlign(31);
          pt->Draw("same");


          bottomPad->cd();
          TH1D *ratioPlotFile2 = (TH1D*)hRef->Clone("ratioPlotFile2");
          ratioPlotFile2->Add(hRef, -1);
          ratioPlotFile2->Divide(hRef);
          TH1D* ratioPlotFile2C = (TH1D*)ratioPlotFile2->Clone("ratioPlotFile2C");
          ratioPlotFile2C->SetFillColor(0);
          TH1D *ratioPlotFile1 = (TH1D*)hCur->Clone("ratioPlotFile1");
          ratioPlotFile1->Add(hRef, -1);
          ratioPlotFile1->Divide(hRef);
	  double maxvalue, minvalue;
	  getHistMaxMinWithError(ratioPlotFile1, maxvalue, minvalue);
	  double range = 1.2 * std::max(std::abs(maxvalue), std::abs(minvalue));
          ratioPlotFile2->GetYaxis()->SetRangeUser(-range, range);
          setStyleRatio(ratioPlotFile2, gRefVersion.Data(), gCurVersion.Data());
          ratioPlotFile2->Draw("hist");
          ratioPlotFile2C->Draw("histsame");
          ratioPlotFile1->Draw("e1same");

          c1->Print(Form("%s/%s.png",s_local_branch_names.Data(), s_histo_names.Data() ));
          ComparisonChiSquare_file << Form("%s/%s %f\n", s_local_branch_names.Data(), s_histo_names.Data(), chisqv);

          delete hRef;
          delete hCur;
          delete bottomPad;
          delete topPad;
          delete c1;

      }

      idx_offset+=folder_size.at(idb);

  }

  ComparisonChiSquare_file.close();


}
