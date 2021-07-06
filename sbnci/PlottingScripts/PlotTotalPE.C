#define PlotTotalPE_cxx
#include "PlotTotalPE.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void PlotTotalPE::Loop()
{
//   In a ROOT session, you can do:
//      root> .L PlotTotalPE.C
//      root> PlotTotalPE t
//      root> t.GetEntry(12); // Fill t data members with entry number 12
//      root> t.Show();       // Show values of entry 12
//      root> t.Show(16);     // Read and show values of entry 16
//      root> t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
//
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
   }
   
   TH1F* TTotalPE = new TH1F("TotalFlashPE",";TotalFlashPE",100,0,1000);
   TH1F* TTime= new TH1F("Time",";Time",100,0,100);
   TH1F* TTimeWidth= new TH1F("TimeWidth",";TimeWidth",20,0,20);
   TH1F* TYCenter= new TH1F("YCenter",";YCenter",50,0, 500);
   TH1F* TZCenter= new TH1F("Zcenter",";ZCenter",50,0,500);
   TH1F* TYWidth= new TH1F("YWidth",";YWidth",50,0,500);
   TH1F* TZWidth= new TH1F("ZWidth",";ZWidth",50,0,500);

   TH2F* PEvsYC = new TH2F("TotalPEvsYCenter", ";YCenter;TotalFlashPE", 25, 0, 500, 50, 0, 1000);
   TH2F* PEvsZC = new TH2F("TotalPEvsZCenter", ";ZCenter;TotalFlashPE", 25, 0, 500, 50, 0, 1000);
   TH2F* PEvsYW = new TH2F("TotalPEvsYWidth", ";YWidth;TotalFlashPE", 25, 0, 500, 50, 0, 1000);
   TH2F* PEvsZW = new TH2F("TotalPEvsZWidth", ";ZWidth;TotalFlashPE", 25, 0, 500, 50, 0, 1000);
   TH2F* PEvsTW = new TH2F("TotalPEvsTimeWidth", ";TimeWidth;TotalFlashPE", 20, 0, 20, 50, 0, 1000);
   TH2F* YCvsZC = new TH2F("YCentervsZCenter", ";Zcenter;YCenter", 20, 0, 1000, 20, 0, 1000);

   PEvsYC->SetMarkerStyle(8);
   PEvsZC->SetMarkerStyle(8);
   PEvsYW->SetMarkerStyle(8);
   PEvsZW->SetMarkerStyle(8);
   PEvsTW->SetMarkerStyle(8);
   YCvsZC->SetMarkerStyle(8);

   float marksz = 1.5;

   PEvsYC->SetMarkerSize(marksz);
   PEvsZC->SetMarkerSize(marksz);
   PEvsYW->SetMarkerSize(marksz);
   PEvsZW->SetMarkerSize(marksz);
   PEvsTW->SetMarkerSize(marksz);
   YCvsZC->SetMarkerSize(marksz);
 
   for(size_t i=0; i<fChain->GetEntries(); i++){ 

      fChain->GetEntry(i); 

      for(size_t iflash=0; iflash<NFlash; iflash++){

         TTotalPE->Fill(TotalFlashPE->at(iflash));
         TTime->Fill(FlashTime->at(iflash));
         TTimeWidth->Fill(FlashTimeWidth->at(iflash));
         TYCenter->Fill(YCenter->at(iflash));
         TZCenter->Fill(ZCenter->at(iflash));
         TYWidth->Fill(YWidth->at(iflash));
         TZWidth->Fill(ZWidth->at(iflash));

         PEvsYC->Fill(YCenter->at(iflash), TotalFlashPE->at(iflash));
         PEvsZC->Fill(ZCenter->at(iflash), TotalFlashPE->at(iflash));
         PEvsYW->Fill(YWidth->at(iflash), TotalFlashPE->at(iflash));
         PEvsZW->Fill(ZWidth->at(iflash), TotalFlashPE->at(iflash));
         PEvsTW->Fill(FlashTimeWidth->at(iflash), TotalFlashPE->at(iflash));
         YCvsZC->Fill(ZCenter->at(iflash), YCenter->at(iflash));

      }
   }

   TCanvas* c1 = new TCanvas();
   TTotalPE->Draw();
   c1->SaveAs("TotalPE.png");
   TCanvas* c2 = new TCanvas();
   TTime->Draw();
   c2->SaveAs("Time.png");
   TCanvas* c3 = new TCanvas();
   TTimeWidth->Draw();
   c3->SaveAs("TimeWidth.png");
   TCanvas* c4 = new TCanvas();
   TYCenter->Draw();
   c4->SaveAs("YCenter.png");
   TCanvas* c5 = new TCanvas();
   TZCenter->Draw();
   c5->SaveAs("ZCenter.png");
   TCanvas* c6 = new TCanvas();
   TYWidth->Draw();
   c6->SaveAs("YWidth.png");
   TCanvas* c7 = new TCanvas();
   TZWidth->Draw();
   c7->SaveAs("ZWidth.png");

   TCanvas* d1 = new TCanvas();
   PEvsYC->Draw();
   TCanvas* d2 = new TCanvas();
   PEvsZC->Draw();
   TCanvas* d3 = new TCanvas();
   PEvsYW->Draw();
   TCanvas* d4 = new TCanvas();
   PEvsZW->Draw();
   TCanvas* d5 = new TCanvas();
   PEvsTW->Draw();
   TCanvas* d6 = new TCanvas();
   YCvsZC->Draw();
   d1->SaveAs("TotalPEvsYCenter.png");
   d2->SaveAs("TotalPEvsZCenter.png");
   d3->SaveAs("TotalPEvsYWidth.png");
   d4->SaveAs("TotalPEvsZWidth.png");
   d5->SaveAs("TotalPEvsTimeWidth.png");
   d6->SaveAs("YCentervsZCenter.png");

}
