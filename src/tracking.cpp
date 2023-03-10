#include <iostream>
#include <list>
#include <typeinfo>

#include <TCanvas.h>
#include <TError.h>
#include <TFitResult.h>
#include <TGraphErrors.h>
#include <TPaveStats.h>

#include "tdata_track.h"

void tracking()
{
  TChain *chain_s = new TChain("tdata");
  chain_s->Add("data/Tracker_sovrapposte_02122015_1.txt.root");

  tdata_track tt_s(chain_s,0 , 0);



  tt_s.Loop();



}

int main()
{
  tracking();
  return 0;
}
