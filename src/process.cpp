#include <iostream>
#include <list>
#include <typeinfo>

#include <TCanvas.h>
#include <TError.h>
#include <TFitResult.h>
#include <TGraphErrors.h>
#include <TPaveStats.h>

#include "tdata.h"

void process()
{
  TChain *chain_s = new TChain("tdata");
  chain_s->Add("data/Tracker_sovrapposte_02122015_1.txt.root");
  chain_s->Add("data/Tracker_sovrapposte_02122015_2.txt.root");
  chain_s->Add("data/Tracker_sovrapposte_02122015_3.txt.root");
  chain_s->Add("data/Tracker_sovrapposte_02122015_4.txt.root");
  chain_s->Add("data/Tracker_sovrapposte_02122015_5.txt.root");
  chain_s->Add("data/Tracker_sovrapposte_02122015_6.txt.root");
  chain_s->Add("data/Tracker_sovrapposte_02122015_7.txt.root");
  chain_s->Add("data/Tracker_sovrapposte_02122015_8.txt.root");
  chain_s->Add("data/Tracker_sovrapposte_02122015_9.txt.root");
  chain_s->Add("data/Tracker_sovrapposte_02122015_10.txt.root");

  TChain *chain_i = new TChain("tdata");
  chain_i->Add("data/Tracker_intermedie_03122015_1.txt.root");
  chain_i->Add("data/Tracker_intermedie_03122015_2.txt.root");
  chain_i->Add("data/Tracker_intermedie_03122015_3.txt.root");
  chain_i->Add("data/Tracker_intermedie_03122015_4.txt.root");
  chain_i->Add("data/Tracker_intermedie_03122015_5.txt.root");
  chain_i->Add("data/Tracker_intermedie_03122015_6.txt.root");
  chain_i->Add("data/Tracker_intermedie_03122015_7.txt.root");
  chain_i->Add("data/Tracker_intermedie_03122015_8.txt.root");
  chain_i->Add("data/Tracker_intermedie_03122015_9.txt.root");
  chain_i->Add("data/Tracker_intermedie_03122015_10.txt.root");

  TChain *chain_d = new TChain("tdata");
  chain_d->Add("data/Tracker_distanziate_02122015_1.txt.root");
  chain_d->Add("data/Tracker_distanziate_02122015_2.txt.root");
  chain_d->Add("data/Tracker_distanziate_02122015_3.txt.root");
  chain_d->Add("data/Tracker_distanziate_02122015_4.txt.root");
  chain_d->Add("data/Tracker_distanziate_02122015_5.txt.root");
  chain_d->Add("data/Tracker_distanziate_02122015_6.txt.root");
  chain_d->Add("data/Tracker_distanziate_02122015_7.txt.root");
  chain_d->Add("data/Tracker_distanziate_02122015_8.txt.root");
  chain_d->Add("data/Tracker_distanziate_02122015_9.txt.root");
  chain_d->Add("data/Tracker_distanziate_03122015_10.txt.root");

  tdata tt_s(chain_s, 0);
  tdata tt_i(chain_i, 1010);
  tdata tt_d(chain_d, 2350);

  //std::list<Float_t> s_hist_sizes = {0.10, 0.30, 1};
  //std::list<Float_t> i_hist_sizes = {0.10, 0.30, 1};
  //std::list<Float_t> d_hist_sizes = {0.20, 0.40, 1};

  std::list<Float_t> s_hist_sizes = {1};
  std::list<Float_t> i_hist_sizes = {1};
  std::list<Float_t> d_hist_sizes = {1};

  tt_s.Loop(&s_hist_sizes);
  tt_i.Loop(&i_hist_sizes);
  tt_d.Loop(&d_hist_sizes);

  tt_s.SaveHists("chart-sov.root");
  tt_i.SaveHists("chart-int.root");
  tt_d.SaveHists("chart-dst.root");

  tt_s.MakeFits();
  tt_i.MakeFits();
  tt_d.MakeFits();

  tt_s.DrawHists();
  tt_i.DrawHists();
  tt_d.DrawHists();

  std::list<std::vector<Float_t>> all_points;

  all_points.splice(all_points.end(), tt_s.GetPoints());
  all_points.splice(all_points.end(), tt_i.GetPoints());
  all_points.splice(all_points.end(), tt_d.GetPoints());

  std::cout << std::endl;
  Info("", "Making linear regression with %lu points...", all_points.size());

  TCanvas canv("c-graph", "td-canvas", 4000, 3000);
  canv.Divide(1);

  TGraphErrors graph(all_points.size());
  graph.SetMarkerColor(4);
  graph.SetMarkerSize(1.5);
  graph.SetMarkerStyle(21);

  std::list<std::vector<Float_t>>::iterator it;
  int i;
  for(i = 0, it = all_points.begin(); it != all_points.end(); it++, i++)
  {
    graph.SetPoint(i, (*it)[1], (*it)[0]);
    graph.SetPointError(i, (*it)[3], (*it)[2]);
  }

  TFitResultPtr fitres = graph.Fit("pol1", "S");

  if (!fitres->IsValid())
  {
    Error("", "There was a problem durng the fitting!");
  }
  else
  {
    Float_t t0  = fitres->Parameter(0);
    Float_t t0e = fitres->ParError(0);
    Float_t v  = 1.0/fitres->Parameter(1);
    Float_t ve = fitres->ParError(1)*(v*v);

    std::cout << std::endl;
    std::cout << "t_delay = "<< t0 << "±" << t0e << std::endl;
    std::cout << "v_muons =" << v << "±" << ve << std::endl;
  }

  graph.Sort();

  canv.cd(0);
  std::cout << std::endl;
  graph.Draw("AP");

  gPad->Update();
  //TPaveStats *stats =(TPaveStats*)canv.GetPrimitive("stats");
  TPaveStats *stats = (TPaveStats*) graph.GetListOfFunctions()->FindObject("stats");
  if (stats != NULL)
  {
    stats->SetName("h1stats");
    stats->SetY1NDC(.95);
    stats->SetY2NDC(.75);
    stats->SetX1NDC(.15);
    stats->SetX2NDC(.35);
    stats->SetTextColor(1);
  }

  canv.SaveAs("linear-regression.png");
}

int main()
{
  process();
  return 0;
}
