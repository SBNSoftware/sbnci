//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Jul  1 14:28:25 2021 by ROOT version 6.22/06
// from TTree pdsTree/Tree with PDS validation information
// found on file: pdsValidationTree.root
//////////////////////////////////////////////////////////

#ifndef PlotTotalPE_h
#define PlotTotalPE_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "vector"

class PlotTotalPE {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Int_t           Run;
   Int_t           SubRun;
   Int_t           Event;
   Int_t           NFlash;
   vector<double>  *TotalFlashPE;
   vector<double>  *FlashTime;
   vector<double>  *FlashTimeWidth;
   vector<double>  *YCenter;
   vector<double>  *ZCenter;
   vector<double>  *YWidth;
   vector<double>  *ZWidth;

   // List of branches
   TBranch        *b_Run;   //!
   TBranch        *b_SubRun;   //!
   TBranch        *b_Event;   //!
   TBranch        *b_NFlash;   //!
   TBranch        *b_TotalFlashPE;   //!
   TBranch        *b_FlashTime;   //!
   TBranch        *b_FlashTimeWidth;   //!
   TBranch        *b_YCenter;   //!
   TBranch        *b_ZCenter;   //!
   TBranch        *b_YWidth;   //!
   TBranch        *b_ZWidth;   //!

   PlotTotalPE(TTree *tree=0);
   virtual ~PlotTotalPE();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef PlotTotalPE_cxx
PlotTotalPE::PlotTotalPE(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("pdsValidationTree.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("pdsValidationTree.root");
      }
      TDirectory * dir = (TDirectory*)f->Get("pdsValidationTree.root:/pdsvalidation");
      dir->GetObject("pdsTree",tree);

   }
   Init(tree);
}

PlotTotalPE::~PlotTotalPE()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t PlotTotalPE::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t PlotTotalPE::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void PlotTotalPE::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   TotalFlashPE = 0;
   FlashTime = 0;
   FlashTimeWidth = 0;
   YCenter = 0;
   ZCenter = 0;
   YWidth = 0;
   ZWidth = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("Run", &Run, &b_Run);
   fChain->SetBranchAddress("SubRun", &SubRun, &b_SubRun);
   fChain->SetBranchAddress("Event", &Event, &b_Event);
   fChain->SetBranchAddress("NFlash", &NFlash, &b_NFlash);
   fChain->SetBranchAddress("TotalFlashPE", &TotalFlashPE, &b_TotalFlashPE);
   fChain->SetBranchAddress("FlashTime", &FlashTime, &b_FlashTime);
   fChain->SetBranchAddress("FlashTimeWidth", &FlashTimeWidth, &b_FlashTimeWidth);
   fChain->SetBranchAddress("YCenter", &YCenter, &b_YCenter);
   fChain->SetBranchAddress("ZCenter", &ZCenter, &b_ZCenter);
   fChain->SetBranchAddress("YWidth", &YWidth, &b_YWidth);
   fChain->SetBranchAddress("ZWidth", &ZWidth, &b_ZWidth);
   Notify();
}

Bool_t PlotTotalPE::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void PlotTotalPE::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t PlotTotalPE::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef PlotTotalPE_cxx
