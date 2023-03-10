#define tdata_track_cxx

#include <iostream>
#include <iomanip>
#include <sstream>
#include <stack>
#include <vector>
#include <cmath>

#include <TError.h>
#include <TH1.h>
#include <TF1.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TPoint.h>
#include <TMath.h>
#include <Math/MinimizerOptions.h>

#include "tdata_track.h"

int tdata_track::GetWireLayer(unsigned int widx)
{
  if (widx <= 11)
  {
    return 0;
  }
  else if ((widx >= 12) && (widx <= 23))
  {
    return 1;
  }
  else if ((widx >= 48) && (widx <= 59))
  {
    return 2;
  }
  else if ((widx >= 60) && (widx <= 71))
  {
    return 3;
  }

  return -1;
}


void tdata_track::CreateWireLookupTable()
{
  for(int i=0;i<12;i++)
  {
    x_filo[i]    =trasl_bot+d0x+cdx+xc/2+spsc/2+(xc+spsc)*i;   // x-view  bottom chamber layer=1
    x_filo[i+12] =trasl_bot+d0x+cdx+xc+spsc+(xc+spsc)*i;       // x-view  bottom chamber layer=2
    x_filo[i+24] =trasl_bot+d0y+cdy+xc/2+spsc/2+(xc+spsc)*i;   // y-view  bottom chamber layer=1
    x_filo[i+36] =trasl_bot+d0y+cdy+xc+spsc+(xc+spsc)*i;       // y-view  bottom chamber layer=2
    x_filo[i+48] =trasl_top+d0y+cdy+xc/2+spsc/2+(xc+spsc)*i;   // x-view  top    chamber layer=1
    x_filo[i+60] =trasl_top+d0x+cdx+xc+spsc+(xc+spsc)*i;       // x-view  top    chamber layer=2
    x_filo[i+72] =trasl_top+d0y+cdy+xc/2+spsc/2+(xc+spsc)*i;   // y-view  top    chamber layer=1
    x_filo[i+84] =trasl_top+d0y+cdy+xc+spsc+(xc+spsc)*i;       // y-view  top    chamber layer=2

    z_filo[i]    =sbdz+sb1z+zc/2;                               //bottom
    z_filo[i+12] =sbdz+sb1z+zc/2+zc;
    z_filo[i+24] =sbdz+sb1z+zc/2+2*zc;
    z_filo[i+36] =sbdz+sb1z+zc/2+3*zc;
    z_filo[i+48] =sbdz+sb1z+zc/2+zz0;                           //top
    z_filo[i+60] =sbdz+sb1z+zc/2+zc+zz0;
    z_filo[i+72] =sbdz+sb1z+zc/2+2*zc+zz0;
    z_filo[i+84] =sbdz+sb1z+zc/2+3*zc+zz0;
  }
}

Long64_t tdata_track::GetProcessedEvents()
{
  return this->nevents;
}

Long64_t tdata_track::GetDiscardedEvents()
{
  return this->ev_tt_disc;
}

Long64_t tdata_track::GetDiscardedSCEvents()
{
  return this->ev_sc_disc;
}

Long64_t tdata_track::GetDiscardedWREvents()
{
  return this->ev_wr_disc;
}

Long64_t tdata_track::GetSCEvents()
{
  return this->ev_sc_read;
}

Long64_t tdata_track::GetWREvents()
{
  return this->ev_wr_read_4 + this->ev_wr_read_3 + this->ev_wr_read_2;
}

void tdata_track::PrintStatistics()
{
  // computing how many events have been discarded

  if (this->nevents == 0)
  {
    Info("", this->MSG_INFO_NOSTATS);
  }
  else
  {
    Info("", this->MSG_INFO_STATS,
         this->nevents,
         this->GetWREvents(),
         100.0 * ((double)this->GetWREvents()/(double)this->nevents),
         this->ev_wr_read_2,
         100.0 * ((double)this->ev_wr_read_2/(double)this->nevents),
         this->ev_wr_read_3,
         100.0 * ((double)this->ev_wr_read_3/(double)this->nevents),
         this->ev_wr_read_4,
         100.0 * ((double)this->ev_wr_read_4/(double)this->nevents),
         this->ev_tt_disc,
         100.0 * ((double)this->ev_tt_disc/(double)this->nevents),
         this->ev_sc_disc,
         100.0 * ((double)this->ev_sc_disc/(double)this->nevents),
         this->ev_wr_disc,
         100.0 * ((double)this->ev_wr_disc/(double)this->nevents));
  }
}



void tdata_track::Loop(std::list<Float_t> *hist_sizes)
{
  //   In a ROOT session, you can do:
  //      root> .L tdata.C
  //      root> tdata t
  //      root> t.GetEntry(12); // Fill t data members with entry number 12
  //      root> t.Show();       // Show values of entry 12
  //      root> t.Show(16);     // Read and show values of entry 16
  //      root> t.Loop();       // Loop on all entries
  //
  //  This is the loop skeleton where:
  //    jentry is the global entry number in the chain
  //    ientry is the entry number in the current Tree
  //
  //  Note that the argument to GetEntry must be:
  //    jentry for TChain::GetEntry
  //    ientry for TTree::GetEntry and TBranch::GetEntry
  //
  //  To read only selected branches, Insert statements like:
  //
  //    METHOD1:
  //      fChain->SetBranchStatus("*",0);  // disable all branches
  //      fChain->SetBranchStatus("branchname",1);  // activate branchname
  //    METHOD2: replace line
  //      fChain->GetEntry(jentry);       //read all branches
  //      b_branchname->GetEntry(ientry); //read only this branch

  if (fChain == 0) return;

  Long64_t nentries = fChain->GetEntriesFast();

  Long64_t nbytes = 0, nb = 0;

  this->CreateWireLookupTable();

  Long64_t event_read_scint = 0;
  Long64_t event_read_wires_4 = 0;
  Long64_t event_read_wires_3 = 0;
  Long64_t event_read_wires_2 = 0;
  Long64_t total_events = 0;

  Float_t distance, delta_t;

  TTree *tt = new TTree("dists", "distributions");

  tt->Branch("TIME", &delta_t, "TIME/F");
  tt->Branch("DIST", &distance, "DIST/F");

  for (Long64_t jentry=0; jentry<nentries;jentry++) 
  {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0)
    {
      break;
    }

    nb = fChain->GetEntry(jentry);
    nbytes += nb;

    total_events++;

    // Selecting only those events in which only two scintillators
    // have been triggered, and the scintillators belong to two
    // differente layers.
    if (SCINT_COUNT != 2)
    {
      // Here more or less that two scintillators have been triggered:
      continue; //skipping the event!
    }
    else if ((SCINT_INDEX[0] < 4) && (SCINT_INDEX[1] < 4))
    {
      // Here both scintillators belong to the BOTTOM layer
      continue; //skipping the event!
    }
    else if ((SCINT_INDEX[0] >= 4) && (SCINT_INDEX[1] >= 4))
    {
      // Here both scintillators belong to the TOP layer
      continue; //skipping the event!
    }
    else if (SCINT_INDEX[0] == 6 || SCINT_INDEX[1] == 6)
    {
      //continue;
    }
    else
    {
      // incrementing the numeber of events that match the above
      // criterias
      event_read_scint++;
    }


    // Computing the Time-Of-Flight (TOF) of the particle. The 'time_mul'
    // constant is needed to get times in ns.
    delta_t = fabs(SCINT_TIME[1] - SCINT_TIME[0]) * this->time_mul;

    // NOTE: Here we only use X drift chambers

    // This variable holds the indexes of all the wires that
    //  have been triggered
    unsigned int x_layers[4] = {0}; 

    // This variable holds the number of drift chambers that
    // have been triggered in each layer.
    std::stack<unsigned int> x_wires_id;

    // For each wire found in the current event...
    for(int i=0; i < WIRES_COUNT; i++)
    {
      // ... get the wire index ...
      Int_t wire_layer = GetWireLayer(WIRE_INDEX[i]);

      // ... and if the wire is along x axis...
      if(wire_layer < 0)
      {
        continue;
      }
      else
      {
        // ... then store the wire index into the stack ...
        x_wires_id.push(WIRE_INDEX[i]);

        // ... and increment the number of wires in the
        // respective layer.
        x_layers[wire_layer]++;
      }
    } 


    if ((x_layers[0] <= 1) && (x_layers[1] <= 1) &&
        (x_layers[2] <= 1) && (x_layers[3] <= 1))
    {
      // selecting only those events in which at most
      // one drift chamber per layer have been triggered,
      switch((x_layers[0] + x_layers[1] + x_layers[2] + x_layers[3]))
      {
        case 4:
          event_read_wires_4++;
          break;

        case 3:
          event_read_wires_3++;
          break;

        case 2:
          // An extra check: there must exactly one wire in TOP
          // and exactly one wire in BOTTOM 
          if (x_layers[0] + x_layers[1] != 1)
          {
            continue;
          }
          else
          {
            event_read_wires_2++;
          }
          break;
        default:
          continue;
      }
    }
    else
    {
      continue; // skiping the event!
    }


    tt->Fill();

  }

  tt->Write();

  // saving statistics to the object private memory space
  this->nevents = total_events;
  this->ev_sc_read = event_read_scint;
  this->ev_wr_read_4 = event_read_wires_4;
  this->ev_wr_read_3 = event_read_wires_3;
  this->ev_wr_read_2 = event_read_wires_2;
  this->ev_sc_disc = nevents - event_read_scint;
  this->ev_wr_disc = event_read_scint - this->GetWREvents();
  this->ev_tt_disc = nevents - this->GetWREvents();

  // printing the statistics
  PrintStatistics();
}




/*
 * NOTE: the following section was originally in tdata.h,
 *      I moved the class constructors into the source file to edit it in an easy way
 *       
 */
tdata_track::tdata_track(TTree *tree, Float_t my_trasl_top, Float_t my_trasl_bot) :
trasl_top(my_trasl_top), trasl_bot(my_trasl_bot), fChain(0)
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
  if (tree == 0) {
    TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("data_tree_dyn.root");
    if (!f || !f->IsOpen()) {
      f = new TFile("data_tree_dyn.root");
    }
    f->GetObject("tdata",tree);
  }
  Init(tree);
}

tdata_track::~tdata_track()
{
  if (!fChain) return;
  delete fChain->GetCurrentFile();
}

Int_t tdata_track::GetEntry(Long64_t entry)
{
  // Read contents of entry.
  if (!fChain) return 0;
  return fChain->GetEntry(entry);
}
Long64_t tdata_track::LoadTree(Long64_t entry)
{
  // Set the environment to read one entry
  if (!fChain)
    return -5;

  Long64_t centry = fChain->LoadTree(entry);
  if (centry < 0)
    return centry;

  if (fChain->GetTreeNumber() != fCurrent)
  {
    fCurrent = fChain->GetTreeNumber();
    Notify();
  }
  return centry;
}

void tdata_track::Init(TTree *tree)
{
  // The Init() function is called when the selector needs to initialize
  // a new tree or chain. Typically here the branch addresses and branch
  // pointers of the tree will be set.
  // It is normally not necessary to make changes to the generated
  // code, but the routine can be extended by the user if needed.
  // Init() will be called many times when running on PROOF
  // (once per file to be processed).

  // Set branch addresses and branch pointers
  if (!tree)
    return;

  fChain = tree;
  fCurrent = -1;
  fChain->SetMakeClass(1);

  fChain->SetBranchAddress("SCINT_COUNT", &SCINT_COUNT, &b_SCINTS_COUNT);
  fChain->SetBranchAddress("SCINT_INDEX", SCINT_INDEX, &b_SCINT_INDEX);
  fChain->SetBranchAddress("SCINT_TIME", SCINT_TIME, &b_SCINT_TIME);
  fChain->SetBranchAddress("WIRES_COUNT", &WIRES_COUNT, &b_WIRES_COUNT);
  fChain->SetBranchAddress("WIRE_INDEX", WIRE_INDEX, &b_WIRE_INDEX);
  fChain->SetBranchAddress("WIRE_TIME", WIRE_TIME, &b_WIRE_TIME);

  Notify();
}

Bool_t tdata_track::Notify()
{
  // The Notify() function is called when a new file is opened. This
  // can be either for a new TTree in a TChain or when when a new TTree
  // is started when using PROOF. It is normally not necessary to make changes
  // to the generated code, but the routine can be extended by the
  // user if needed. The return value is currently not used.

  return kTRUE;
}

void tdata_track::Show(Long64_t entry)
{
  // Print contents of entry.
  // If entry is not specified, print current entry
  if (!fChain) return;
  fChain->Show(entry);
}

Int_t tdata_track::Cut(Long64_t entry)
{
  // This function may be called from Loop.
  // returns  1 if entry is accepted.
  // returns -1 otherwise.
  return 1;
}

