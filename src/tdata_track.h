//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon Jan 18 11:50:31 2016 by ROOT version 6.06/00
// from TTree tdata/data
// found on file: data_tree_dyn.root
//////////////////////////////////////////////////////////

#ifndef tdata_track_h
#define tdata_track_h

#include <list>
#include <stack>
#include <string.h>

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TH1.h>

// Header file for the classes stored in the TTree if any.

class tdata_track {
  private:
    Float_t d0x=0;
    Float_t d0y=0;
    Float_t cdx=88.75;
    Float_t cdy=87.75;
    Float_t xc=92;
    Float_t sb1z=5;
    Float_t sbdz=276; /// DA CONTROLLARE CON GIVI
    Float_t spsc=1.5;
    Float_t zc=29;
    Float_t zz0=2857; /// DA CONTROLLARE CON GIVI
    Float_t trasl_top=0;
    Float_t trasl_bot=0;
		Float_t vdrift[4] = {0.04}; // mm/ns 1 VDRIFT FOR ALL LAYERS TO BE IMPROOVED

    Int_t n_xbot=0;
    Int_t n_xtop=0;
    Int_t i_filo_bot[96]={0};
    Int_t i_filo_top[96]={0};

		Int_t view_filo[96]={0};
    Float_t x_filo[96]={0};
    Float_t z_filo[96]={0};

    //scintillators z distance
    Float_t scint_top_z = 3245; // mm
    Float_t scint_w     = 1000; // mm

    Float_t time_mul = 0.025; // ns CAMAC TDC for SCINT
    Float_t time_res = 1.04; // ns VME TDC for DRIFT CHAMBER

    // How many events have been processed
    Long64_t nevents = 0;

    // How many events have been succesfully read
    Long64_t ev_sc_read = 0;
    Long64_t ev_wr_read_2 = 0;
    Long64_t ev_wr_read_3 = 0;
    Long64_t ev_wr_read_4 = 0;

    // How many events have been discarded
    Long64_t ev_sc_disc = 0;
    Long64_t ev_wr_disc = 0;
    Long64_t ev_tt_disc = 0;

	
    TH1F *drifttime_hist = new TH1F("DriftTime", "DriftTime", 1000, 0, 1000);
    TH1F *driftlenght_hist = new TH1F("DriftLenght", "DriftLenght", 500, 0, 100);
  

    const char *MSG_INFO_NOSTATS="No statistics available yet!\n";
    const char *MSG_INFO_STATS = "\n"
      "-------------------------------------\n"
      "Total events processed: %lld\n"
      "Total events used: %lld (%06.3f%%)\n"
      "  2 wires hit: %lld (%06.3f%%)\n"
      "  3 wires hit: %lld (%06.3f%%)\n"
      "  4 wires hit: %lld (%06.3f%%)\n"
      "Total events discarded: %lld (%06.3f%%)\n"
      "  events discarded (scint): %lld (%06.3f%%)\n"
      "  events discarded (wires): %lld (%06.3f%%)\n"
      "-------------------------------------\n";

    void    CreateWireLookupTable();
    Int_t   GetWireLayer(unsigned int widx);

  public :
    TTree          *fChain;   //!pointer to the analyzed TTree or TChain
    Int_t           fCurrent; //!current Tree number in a TChain

    // Fixed size dimensions of array or collections stored in the TTree,
    // if any.

    // Declaration of leaf types
    Int_t           SCINT_COUNT;
    Int_t           SCINT_INDEX[8];   //[SCINTS_COUNT]
    Float_t         SCINT_TIME[8];   //[SCINTS_COUNT]
    Int_t           WIRES_COUNT;
    Int_t           WIRE_INDEX[95];   //[WIRES_COUNT]
    Float_t         WIRE_TIME[95];   //[WIRES_COUNT]

    // List of branches
    TBranch        *b_SCINTS_COUNT;   //!
    TBranch        *b_SCINT_INDEX;   //!
    TBranch        *b_SCINT_TIME;   //!
    TBranch        *b_WIRES_COUNT;   //!
    TBranch        *b_WIRE_INDEX;   //!
    TBranch        *b_WIRE_TIME;   //!

    tdata_track(TTree *tree=0, Float_t my_trasl_top=0, Float_t my_trasl_bot=0);
    virtual ~tdata_track();
    virtual Int_t    Cut(Long64_t entry);
    virtual Int_t    GetEntry(Long64_t entry);
    virtual Long64_t LoadTree(Long64_t entry);
    virtual void     Init(TTree *tree);
    virtual void     Loop(std::list<Float_t> *hist_sizes=NULL);
    virtual Bool_t   Notify();
    virtual void     Show(Long64_t entry = -1);

    Long64_t GetProcessedEvents();   // Returns the number of events processed
    Long64_t GetDiscardedEvents();   // Returns the number of events discarded

    Long64_t GetDiscardedSCEvents(); // Returns the number of events discarded
                                     // due to scintillator filtering

    Long64_t GetDiscardedWREvents(); // Returns the number of events discarded
                                     // due to drift chamber filtering

    Long64_t GetSCEvents(); // Returns the number of events discarded
                            // due to scintillator filtering

    Long64_t GetWREvents(); // Returns the number of events discarded
                            // due to drift chamber filtering

    void PrintStatistics(); // Prints the statistics about the processed events
};

#endif

