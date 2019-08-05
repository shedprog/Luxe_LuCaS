#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TProfile.h"
#include "TRandom.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TNtuple.h"
#include "TProfile.h"
#include "TBenchmark.h"
#include "TStyle.h"
#include "TPaveText.h"
#include "TFrame.h"
#include "TF1.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TInterpreter.h"
// #include "gRoot.h"
#include <vector>

void Build_Occupancy(const Char_t *datapath= "") {

  TFile *f = new TFile(datapath);

  TTree *tree = (TTree*)f->Get("Tracks_true");

  double x[500],y[500];

  int n = tree->GetEntries();
  std::cout<<"Events: "<<n<<"\n";

  TH2F* hist = new TH2F("hist","Occupancy",600,-300.,300.,200,-20.,20.);

  for(int i=0; i<n; i++){
  tree->GetEntry(i);
  y[i] = tree->GetLeaf("y")->GetValue();
  x[i] = tree->GetLeaf("x")->GetValue();
  std::cout<<x[i]<<" "<<y[i]<<"\n"; 
  hist->Fill(x[i],y[i]);
  }

  TCanvas *c1 = new TCanvas("c1","The Ntuple canvas",1500,900);
  c1->Divide(1,1);
  c1->cd();
  c1->cd(1);

  hist->Draw("colz");

  c1->Update();
  c1->SaveAs("./2D_Occupancy.png");

}

void Build_2DHIts(const Char_t *datapath= "") {

   // gStyle->
   // gROOT->SetStyle("ATLAS");

   TFile *f = new TFile(datapath);

   TTree *Events = (TTree*)f->Get("Events");
   TTree *Track = (TTree*)f->Get("Tracks");
   TTree *Hits = (TTree*)f->Get("Hits");

   double Weight = 1;

   std::vector<int>* Track_PDG = 0;

   std::vector<double>* Hits_xHit = 0;
   std::vector<double>* Hits_yHit = 0; 
   std::vector<double>* Hits_zCell = 0;
   
   // Events->SetBranchAddress("Weight",&Weight);

   Track->SetBranchAddress("Tracks_PDG",&Track_PDG);
   
   Hits->SetBranchAddress("Hits_xHit",&Hits_xHit);
   Hits->SetBranchAddress("Hits_yHit",&Hits_yHit);
   Hits->SetBranchAddress("Hits_zCell",&Hits_zCell);

   std::vector<TH2F*> Hists;

   for (Int_t i = 1; i<=20; i++)
   {
    TH2F* hist = new TH2F(Form("hist%d",i),Form("Occupancy of Hits layer %d",i),200,-300.,300.,50,-100.,100.);
    Hists.push_back(hist);
   }

   Int_t nentries = (Int_t)Events->GetEntries();
   std::cout<<"Entries: "<<nentries<<" \n";
   for (Int_t i=0;i<nentries;i++) 
   {
     Events->GetEntry(i);
     Track->GetEntry(i);
     Hits->GetEntry(i);

     // if((*Track_PDG)[0]==11 or (*Track_PDG)[0]==-11){
     //      std::cout << Weight << "\n"; 
     // }

     int Number_of_Hits = Hits_xHit->size();

     // if(Number_of_Hits!=0 and (*Track_PDG)[0]==22){
     if(Number_of_Hits!=0){

      for (Int_t j=0;j<Number_of_Hits;j++){
        int layer =  ((int) (*Hits_zCell)[j]) - 1;
        // std::cout << layer << "\n";
        Hists.at(layer)->Fill((*Hits_xHit)[j],(*Hits_yHit)[j],Weight);
      }

     }
   }


   TCanvas *c1 = new TCanvas("c1","The Ntuple canvas",1500,900);
   c1->Divide(2,2);
   c1->cd();

  for (Int_t i = 0; i<4; i++)
  {
    c1->cd(i+1);
    Hists.at(i)->Draw("colz");
  }  
  c1->Update();
  c1->SaveAs("./layers1.png");

  TCanvas *c2 = new TCanvas("c2","The Ntuple canvas",1500,900);
  c2->Divide(2,2);
  c2->cd();

  for (Int_t i = 0; i<4; i++)
  {
    c2->cd(i+1);
    Hists.at(i+4)->Draw("colz");
  }  
  c2->Update();
  c2->SaveAs("./layers2.png");

  TCanvas *c3 = new TCanvas("c3","The Ntuple canvas",1500,900);
  c3->Divide(2,2);
  c3->cd();

  for (Int_t i = 0; i<4; i++)
  {
    c3->cd(i+1);
    Hists.at(i+8)->Draw("colz");
  }  
  c3->Update();
  c3->SaveAs("./layers3.png");

     TCanvas *c4 = new TCanvas("c4","The Ntuple canvas",1500,900);
   c4->Divide(2,2);
   c4->cd();

  for (Int_t i = 0; i<4; i++)
  {
    c4->cd(i+1);
    Hists.at(i+12)->Draw("colz");
  }  
  c4->Update();
  c4->SaveAs("./layers4.png");

  TCanvas *c5 = new TCanvas("c5","The Ntuple canvas",1600,900);
  c5->Divide(2,2);
  c5->cd();

  for (Int_t i = 0; i<4; i++)
  {
    c5->cd(i+1);
    Hists.at(i+16)->Draw("colz");
  }
  c5->Update();
  c5->SaveAs("./layers5.png");
}

void EnergyOfx_reco(const Char_t *datapath= "") 
{
   TFile *f = new TFile(datapath);

   TTree *tree = (TTree*)f->Get("mc");

   double x[100],E[100];

   int n = tree->GetEntries();
   tree->GetEntry(0);
   std::cout<<"Events: "<<n<<"\n";

   Int_t cal_n_clusters = tree->GetLeaf("cal_n_clusters")->GetValue();
   std::cout<<"Clasters: "<<cal_n_clusters<<"\n";

   for(int i=0; i<cal_n_clusters; i++){
   E[i] = tree->GetLeaf("cal_cluster_energy")->GetValue(i)/1000; 
   x[i] = tree->GetLeaf("cal_cluster_x")->GetValue(i); 
   }
   TGraph* gr = new TGraph(cal_n_clusters,x,E);
   gr->SetTitle("Cluster Energy (x) - no calibration");
   gr->GetXaxis()->SetTitle("x [mm]");
   gr->GetYaxis()->SetTitle("E in silicon [GeV]");
   gr->SetMinimum(0);

   TCanvas *c1 = new TCanvas("c1","The Ntuple canvas",1500,900);
   gr->Draw("A*");

  c1->Update();
  c1->SaveAs("./EnergyOfx_reco.png");
}

void EnergyOfx_detector(const Char_t *datapath= "") 
{
   TFile *f = new TFile(datapath);

   TTree *tree = (TTree*)f->Get("Tracks_true");

   double x[200],E[200];

   int n = tree->GetEntries();
   std::cout<<"Events: "<<n<<"\n";


   for(int i=0; i<n; i++){
   tree->GetEntry(i);
   E[i] = tree->GetLeaf("E")->GetValue()/1000;
   // std::cout<<i<<" "<<E[i];
   x[i] = tree->GetLeaf("x")->GetValue();
   // std::cout<<" "<<x[i]<<"\n";
   }

   TGraph* gr = new TGraph(n,x,E);
   gr->SetTitle("Energy (x) - initial particle");
   gr->GetXaxis()->SetTitle("x [mm]");
   gr->GetYaxis()->SetTitle("E [GeV]");
   gr->SetMinimum(0);

   TCanvas *c1 = new TCanvas("c1","The Ntuple canvas",1500,900);
   gr->Draw("A*");

  c1->Update();
  c1->SaveAs("./EnergyOfx_detector.png");
}

void PlotLCOut(const Char_t *datapath= "") 
{
  // Build_2DHIts(datapath);
  // EnergyOfx_reco(datapath);
  // EnergyOfx_detector(datapath);
  Build_Occupancy(datapath);
}