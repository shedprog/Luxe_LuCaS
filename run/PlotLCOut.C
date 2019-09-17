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
#include "TMath.h"
#include "TMarker.h"
// #include "gRoot.h"
#include <vector>
#include <typeinfo>

std::string basename(const std::string &filename)
{
    if (filename.empty()) {
        return {};
    }

    auto len = filename.length();
    auto index = filename.find_last_of("/\\");

    if (index == std::string::npos) {
        return filename;
    }

    if (index + 1 >= len) {

        len--;
        index = filename.substr(0, len).find_last_of("/\\");

        if (len == 0) {
            return filename;
        }

        if (index == 0) {
            return filename.substr(1, len - 1);
        }

        if (index == std::string::npos) {
            return filename.substr(0, len);
        }

        return filename.substr(index + 1, len - index - 1);
    }

    return filename.substr(index + 1, len - index);
}

void build_tracks_flow(const Char_t *datapath= "") {

   // gStyle->
   // gROOT->SetStyle("ATLAS");

   TFile *f = new TFile(datapath);

   TTree *Events = (TTree*)f->Get("Events");
   TTree *Track = (TTree*)f->Get("Tracks");
   TTree *Hits = (TTree*)f->Get("Hits");
   TTree *Track_true = (TTree*)f->Get("Tracks_true");

   double Weight = 1;

   // std::vector<double>* Hits_yHit = 0;
   std::vector<double>* Hits_zCell = 0;
   std::vector<double>* Hits_xCell = 0;
   std::vector<double>* Hits_yCell = 0;
   std::vector<double>* Hits_eHit = 0;
   std::vector<double>* Hits_Sensor = 0;
   Hits->SetBranchAddress("Hits_xCell",&Hits_xCell);
   Hits->SetBranchAddress("Hits_yCell",&Hits_yCell);
   Hits->SetBranchAddress("Hits_zCell",&Hits_zCell);
   Hits->SetBranchAddress("Hits_eHit",&Hits_eHit);
   Hits->SetBranchAddress("Hits_Sensor",&Hits_Sensor);
   TH1F* hist_E = new TH1F("hist","Flow of initial electrons before calorimeter weighted by deposit energy in 1-st layer",200,-300.,300.);

   double xpos=0;
   Track_true->SetBranchAddress("x",&xpos);
   TH1F* hist_x = new TH1F("hist_x","Flow of initial electrons before calorimeter",200,-300.,300.);


   //Filling of Hits
   Int_t nentries = (Int_t)Hits->GetEntries();
   // std::cout<<"Entries: "<<nentries<<" \n";
   // for (Int_t i=0;i<nentries;i++)
   // {
   //   Hits->GetEntry(i);

   //   int Number_of_Hits = Hits_yCell->size();

   //   if(Number_of_Hits!=0){
   //    for (Int_t j=0;j<Number_of_Hits;j++){
   //      int layer =  ((int) (*Hits_zCell)[j]);
   //      if(layer==1){

   //        int pad = (*Hits_yCell)[j];
   //        int sector = (*Hits_xCell)[j];
   //        int sensor = (*Hits_Sensor)[j];


   //        double shift_x = 195.2 - 80.0 + 1.61;
   //        double rho = 80. + 0.9 + 1.8 * pad;
   //        double phi = - TMath::Pi()/12 + TMath::Pi()/48 + (sector - 1) * TMath::Pi()/24;

   //        if(sensor <= 2) phi = phi + TMath::Pi();

   //        double y = rho * TMath::Sin(phi);
   //        double x_ = rho * TMath::Cos(phi);

   //        if(sensor == 4) x_ = x_ + shift_x;
   //        else if(sensor == 2) x_ = x_ - shift_x;

   //      // Cordinates start not in the center of curvature
   //      // but in the begining of cordinates in Geant4
   //        double cord_shift = 195.2 - 30.0 - 120.0;

   //        if(sensor <= 2) x_ = x_ + cord_shift;
   //        else if(sensor >= 3) x_ = x_ - cord_shift;

   //        Weight = (*Hits_eHit)[j];

   //        //~~~~~Have to be changed for one BX~~~~
   //        //hist_E->Fill(x_,Weight);
   //        hist_E->Fill(x_,Weight*0.001);

   //      }
   //    }
   //    }
   // }


   //Filling of track position
   nentries = (Int_t)Track_true->GetEntries();
   // std::cout<<"Entries: "<<nentries<<" \n";
   for (Int_t i=0;i<nentries;i++)
   {
     Track_true->GetEntry(i);

     Hits->GetEntry(i);
     int Number_of_Hits = Hits_eHit->size();

     double E_weight = 0;

     if(Number_of_Hits!=0){
      for (Int_t j=0;j<Number_of_Hits;j++){
        int layer =  ((int) (*Hits_zCell)[j]);
        if(layer==1) E_weight = E_weight + (*Hits_eHit)[j];
      }
      }

     //~~~~~Have to be changed for one BX~~~~
     // hist_x->Fill(xpos);
     hist_E->Fill(xpos,E_weight);
     hist_x->Fill(xpos);
   }

  TCanvas *c1 = new TCanvas("c1","The Ntuple canvas",1500,1000);
  c1->Divide(1,2);
  c1->cd();

  c1->cd(1);
  hist_x->Draw("HIST");

  c1->cd(2);
  hist_E->Draw("HIST");

  c1->Update();
  c1->SaveAs("./Flow.png");

}
void Build_Occupancy(const Char_t *datapath= "") {


  TFile *f = new TFile(datapath);

  TTree *tree = (TTree*)f->Get("Tracks_testplane");

  double x[500],y[500];

  int n = tree->GetEntries();
  std::cout<<"Events: "<<n<<"\n";

  TH2F* hist = new TH2F("hist",basename(datapath),400,-800.,800.,400,-8.,8.);

  for(int i=0; i<n; i++){
  std::cout << i << " ";
  tree->GetEntry(i);
  // y= tree->GetLeaf("y")->GetValue();
  // x[i] = tree->GetLeaf("x")->GetValue();
  // std::cout<<x[i]<<" "<<y[i]<<"\n";
  double weight = tree->GetLeaf("Weight")->GetValue();
  hist->Fill(tree->GetLeaf("x")->GetValue(),tree->GetLeaf("y")->GetValue(),weight/1000);
  }

  TCanvas *c1 = new TCanvas("c1","The Ntuple canvas",1500,900);
  c1->Divide(1,1);
  c1->cd();
  c1->cd(1);

  hist->Draw("colz");
  hist->GetXaxis()->SetTitle("x [mm]");
  hist->GetYaxis()->SetTitle("y [mm]");

  c1->Update();
  std::cout << "Press Enter to Continue";
  std::cin.ignore();

  c1->SaveAs(Form("./Occupancy%s.png",basename(datapath)));
gROOT->SetStyle("ATLAS");
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
    TH2F* hist = new TH2F(Form("hist%d",i),Form("Occupancy of Hits layer %d",i),200,-1500.,1500.,200,-5.,5.);
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

void EnergyProfile(const Char_t *datapath= "")
{

   // gStyle->
   gROOT->SetStyle("ATLAS");

   TFile *f = new TFile(datapath);

   TTree *Events = (TTree*)f->Get("Events");
   TTree *Track = (TTree*)f->Get("Tracks");
   TTree *Track_true = (TTree*)f->Get("Tracks_true");
   TTree *Hits = (TTree*)f->Get("Hits");

   double Weight = 1;

   std::vector<int>* Track_PDG = 0;

   std::vector<double>* Hits_xCell = 0;
   std::vector<double>* Hits_yCell = 0;
   std::vector<double>* Hits_zCell = 0;
   std::vector<double>* Hits_eHit = 0;
   std::vector<double>* Hits_Sensor = 0;

   // Events->SetBranchAddress("Weight",&Weight);

   Track->SetBranchAddress("Tracks_PDG",&Track_PDG);

   Hits->SetBranchAddress("Hits_xCell",&Hits_xCell);
   Hits->SetBranchAddress("Hits_yCell",&Hits_yCell);
   Hits->SetBranchAddress("Hits_zCell",&Hits_zCell);
   Hits->SetBranchAddress("Hits_Sensor",&Hits_Sensor);
   Hits->SetBranchAddress("Hits_eHit",&Hits_eHit);

   std::vector< std::vector<TH2F*> > Hists;

   for (Int_t i = 1; i<=20; i++)
   {
    std::vector<TH2F*> one_layer;
      // for (Int_t j = 1; j<=4; j++){

      TH2F* hist = new TH2F(Form("hist%d_%d",i,1),Form("layer %d sensor 1",i),64,0.5,64.5,4,0.5,4.5);
      one_layer.push_back(hist);

      hist = new TH2F(Form("hist%d_%d",i,2),Form("layer %d sensor 2",i),64,0.5,64.5,4,0.5,4.5);
      one_layer.push_back(hist);

      hist = new TH2F(Form("hist%d_%d",i,3),Form("layer %d sensor 3",i),64,-64.5,-0.5,4,-4.5,-0.5);
      one_layer.push_back(hist);

      hist = new TH2F(Form("hist%d_%d",i,4),Form("layer %d sensor 4",i),64,-64.5,-0.5,4,-4.5,-0.5);
      one_layer.push_back(hist);


      // }
    Hists.push_back(one_layer);
   }

   Int_t nentries = (Int_t)Events->GetEntries();
   std::cout<<"Entries: "<<nentries<<" \n";

   for (Int_t i=0;i<nentries;i++)
   {
     Hits->GetEntry(i);


     int Number_of_Hits = Hits_xCell->size();

     if(Number_of_Hits!=0){
      for (Int_t j=0;j<Number_of_Hits;j++){
        int layer =  ((int) (*Hits_zCell)[j]) - 1;
        int sensor = ((int) (*Hits_Sensor)[j]) - 1;

        if(sensor==0 or sensor==1)
                Hists.at(layer).at(sensor)->Fill((*Hits_yCell)[j],
                                                 (*Hits_xCell)[j],
                                                 (*Hits_eHit)[j]*Weight);
        else if(sensor==2 or sensor==3)
                Hists.at(layer).at(sensor)->Fill(-(*Hits_yCell)[j],
                                                 -(*Hits_xCell)[j],
                                                 (*Hits_eHit)[j]*Weight);
        }


      }
     }



  for (Int_t i = 0; i<20; i++)
  {
    TCanvas *c1 = new TCanvas("c1","The Ntuple canvas",1500,900);
    c1->Divide(4,1);
    c1->cd();

    c1->cd(1);
    Hists.at(i).at(3)->Draw("colz");

    c1->cd(2);
    Hists.at(i).at(2)->Draw("colz");

    c1->cd(3);
    Hists.at(i).at(0)->Draw("colz");

    c1->cd(4);
    Hists.at(i).at(1)->Draw("colz");

    c1->Update();
    c1->SaveAs(Form("./output/layer_%d.png",i));
  }
}

void EnergyTower4(const Char_t *datapath= "")
{
   gROOT->SetStyle("ATLAS");

   TFile *f = new TFile(datapath);

   TTree *Hits = (TTree*)f->Get("Hits");

   std::vector<double>* Hits_xCell = 0;
   std::vector<double>* Hits_yCell = 0;
   // std::vector<double>* Hits_zCell = 0;
   std::vector<double>* Hits_eHit = 0;
   std::vector<double>* Hits_Sensor = 0;

   Hits->SetBranchAddress("Hits_xCell",&Hits_xCell);
   Hits->SetBranchAddress("Hits_yCell",&Hits_yCell);
   // Hits->SetBranchAddress("Hits_zCell",&Hits_zCell);
   Hits->SetBranchAddress("Hits_Sensor",&Hits_Sensor);
   Hits->SetBranchAddress("Hits_eHit",&Hits_eHit);

   TH1F* hist = new TH1F("hist","Energies (joined by 4 sectors)",276,-60,-40);

   Int_t nentries = (Int_t)Hits->GetEntries();
   std::cout<<"Entries: "<<nentries<<" \n";

   for (Int_t i=0;i<nentries;i++)
   {
     Hits->GetEntry(i);
     int Number_of_Hits = Hits_xCell->size();
     if(Number_of_Hits!=0){
      for (Int_t j=0;j<Number_of_Hits;j++){
        // int layer =  ((int) (*Hits_zCell)[j]) - 1;
        int sensor = ((int) (*Hits_Sensor)[j]);

        if(sensor==1)
        {
          hist->Fill(1+(*Hits_yCell)[j],(*Hits_eHit)[j]);
        }

        else if(sensor==2)
        {
          hist->Fill(1+(*Hits_yCell)[j]+64,(*Hits_eHit)[j]);
        }

        else if(sensor==3)
        {
          hist->Fill(-1-(*Hits_yCell)[j],(*Hits_eHit)[j]);
        }
        else if(sensor==4)
        {
          hist->Fill(-1-(*Hits_yCell)[j]-64,(*Hits_eHit)[j]);
        }

       }
      }
     }

    TFile *out = new TFile("TowerEnergy4.root","RECREATE");
    out->cd();
    TCanvas *c1 = new TCanvas("c1","The Ntuple canvas",1500,900);
    c1->Divide(1,1);
    c1->cd();
    hist->Draw("HISTO");
    hist->GetXaxis()->SetTitle("Number of pad");
    hist->GetYaxis()->SetTitle("E [MeV]");
    hist->SetLineColor(2);
    hist->SetLineWidth(3);
    c1->Update();
    int l;
    std::cin>>l;
    // c1->SaveAs("./TowerEnergy4.pdf");
    c1->Write();
    out->Close();
    f->Close();
}

void EnergyTower_andInputParticle(const Char_t *datapath= "")
{
   gROOT->SetStyle("ATLAS");

   TFile *f = new TFile(datapath);

   TTree *Tracks_true_test = (TTree*)f->Get("Tracks_testplane");

   TTree *Hits = (TTree*)f->Get("Hits");
   std::vector<double>* Hits_xCell = 0;
   std::vector<double>* Hits_yCell = 0;
   std::vector<double>* Hits_xHit = 0;
   std::vector<double>* Hits_zHit = 0;
   std::vector<double>* Hits_zCell = 0;
   std::vector<double>* Hits_eHit = 0;
   std::vector<double>* Hits_Sensor = 0;
   Hits->SetBranchAddress("Hits_xCell",&Hits_xCell);
   Hits->SetBranchAddress("Hits_yCell",&Hits_yCell);
   Hits->SetBranchAddress("Hits_xHit",&Hits_xHit);
   Hits->SetBranchAddress("Hits_zHit",&Hits_zHit);
   Hits->SetBranchAddress("Hits_zCell",&Hits_zCell);
   Hits->SetBranchAddress("Hits_Sensor",&Hits_Sensor);
   Hits->SetBranchAddress("Hits_eHit",&Hits_eHit);
   TH1F* hist_energy = new TH1F("hist_energy","Energies (joined by 4 sectors)",2600,-130,130);

   TTree *Tracks_true = (TTree*)f->Get("Tracks_true");
   std::vector<double> x_first_layer;
   std::vector<double> y_first_layer;
   std::vector<double> x;
   std::vector<double> y;
   std::vector<double> z;
   std::vector<double> Mx;
   std::vector<double> My;
   std::vector<double> Mz;
   std::vector<double> E;
   // Tracks_true->SetBranchAddress("x",&x);
   // Tracks_true->SetBranchAddress("y",&y);
   // Tracks_true->SetBranchAddress("z",&z);
   // Tracks_true->SetBranchAddress("Mx",&Mx);
   // Tracks_true->SetBranchAddress("My",&My);
   // Tracks_true->SetBranchAddress("Mz",&Mz);
   // Tracks_true->SetBranchAddress("E",&E);
   TH1F* hist_initial_particle_E = new TH1F("hist_initial_particle_E","Energy for fixed pad",2600,-130,130);

   Int_t nentries = (Int_t)Hits->GetEntries();
   std::cout<<"Entries: "<<nentries<<" \n";
   for (Int_t i=0;i<nentries;i++)
   {
     Hits->GetEntry(i);
     int Number_of_Hits = Hits_xCell->size();
     if(Number_of_Hits!=0){
      for (Int_t j=0;j<Number_of_Hits;j++){
        // int layer =  ((int) (*Hits_zCell)[j]);
        int sensor = ((int) (*Hits_Sensor)[j]);
        if(sensor==1) hist_energy->Fill(1+(*Hits_yCell)[j],(*Hits_eHit)[j]*87);
        else if(sensor==2) hist_energy->Fill(1+(*Hits_yCell)[j]+64,(*Hits_eHit)[j]*87);
        else if(sensor==3) hist_energy->Fill(-1-(*Hits_yCell)[j],(*Hits_eHit)[j]*87);
        else if(sensor==4) hist_energy->Fill(-1-(*Hits_yCell)[j]-64,(*Hits_eHit)[j]*87);
       }
      }
     }

   // Select Track coressponding to it's hit in the first layer of the callorimeter
   // Because in the first layer there will be still some back-scattered particles
   Int_t nentries_t = (Int_t)Tracks_true->GetEntries();
   std::cout<<"Entries: "<<nentries_t<<" \n";
   for (Int_t t=0;t<nentries_t;t++)
   {
     Tracks_true->GetEntry(t);
     double E = Tracks_true->GetLeaf("E")->GetValue();
     double sensor = Tracks_true->GetLeaf("Sensor")->GetValue();
     double pad = Tracks_true->GetLeaf("xCell")->GetValue();

      if(sensor==1) hist_initial_particle_E->Fill(1+pad,E);
      else if(sensor==2) hist_initial_particle_E->Fill(1+pad+64,E);
      else if(sensor==3) hist_initial_particle_E->Fill(-1-pad,E);
      else if(sensor==4) hist_initial_particle_E->Fill(-1-pad-64,E);


    }




    TFile *out = new TFile(Form("TowerEnergy_Tracks%s",basename(datapath)),"RECREATE");
    out->cd();
    TCanvas *c1 = new TCanvas("c1","The Ntuple canvas",1500,900);
    c1->Divide(1,1);
    c1->cd();

    hist_initial_particle_E->SetTitle(Form("Energy for fixed pad, tracks: %d",nentries_t));
    hist_initial_particle_E->GetXaxis()->SetTitle("Number of pad");
    hist_initial_particle_E->GetYaxis()->SetTitle("E [MeV]");
    hist_initial_particle_E->SetLineColor(3);
    hist_initial_particle_E->SetLineWidth(1);
    hist_initial_particle_E->Draw("HISTO");

    hist_energy->GetXaxis()->SetTitle("Number of pad");
    hist_energy->GetYaxis()->SetTitle("E [MeV]");
    hist_energy->SetLineColor(2);
    hist_energy->SetLineWidth(3);
    hist_energy->Draw("HISTOSame");

    auto legend = new TLegend(0.17,0.83,0.37,0.93);
    // legend->SetHeader("The Legend Title","C"); // option "C" allows to center the header
    legend->AddEntry(hist_initial_particle_E,"Energy of initial particle","l");
    legend->AddEntry(hist_energy,"Energy absorbed in tower","l");
    legend->AddEntry("",Form("number of e^{+}e^{-}: %d",(Int_t)Tracks_true_test->GetEntries()),"");
    legend->AddEntry("",Form("number of e^{+}e^{-} in calorimeter: %d",(Int_t)Tracks_true->GetEntries()),"");
    // legend->SetFillColor(0);
    // legend->SetFillStyle(0);
    legend->SetTextFont(42);
    legend->SetTextSize(0.015);
    // legend->SetBorderSize(0);
    // legend->AddEntry("gr","Graph with error bars","lep");
    legend->Draw("");

    c1->Update();
    std::cout << "Press Enter to Continue";
    std::cin.ignore();
    c1->SaveAs(Form("TowerEnergy_Tracks%s.png",basename(datapath)));
    c1->Write();
    out->Close();
    f->Close();
}

void LongitudinalProfile(const Char_t *datapath= "")
{

   // gROOT->SetStyle("ATLAS");

   TFile *f = new TFile(datapath);

   TTree *Hits = (TTree*)f->Get("Hits");
   std::vector<double>* Hits_eHit = 0;
   int EventID;
   Hits->SetBranchAddress("Hits_eHit",&Hits_eHit);
   Hits->SetBranchAddress("eventID",&EventID);

   TTree *True = (TTree*)f->Get("Tracks_true");

   TH2F* LongitudinalProfile = new TH2F("LongitudinalProfile","Energies absorbed vs angle",200,400.0,1700,200,0,0.2);

   double Energies[10000] = {0};
   double Angles[10000] = {0};

   Int_t nentries = (Int_t)Hits->GetEntries();
   std::cout<<"Entries: "<<nentries<<" \n";
   for (Int_t i=0;i<nentries;i++)
   {
     Hits->GetEntry(i);
     int Number_of_Hits = Hits_eHit->size();
     if(Number_of_Hits!=0){
      for (Int_t j=0;j<Number_of_Hits;j++){
        // std::cout<<EventID<<"\n";
        Energies[EventID] = Energies[EventID] + (*Hits_eHit)[j];
        // std::cout<<Energies[EventID]<<"\n";
       }
      }
     }

   nentries = (Int_t)True->GetEntries();
   std::cout<<"Entries: "<<nentries<<" \n";
   for (Int_t i=0;i<nentries;i++)
   {
     True->GetEntry(i);
     int ID = True->GetLeaf("eventID")->GetValue();
     double Mz = True->GetLeaf("Mz")->GetValue();
     double Mx = True->GetLeaf("Mx")->GetValue();
     Angles[ID] = TMath::Abs(TMath::ATan(Mx/Mz));
     // std::cout<<Angles[ID]<<"\n";
   }

   for(int i=0; i<10000; ++i){
    if(Angles[i]!=0.0){
      // std::cout<<"Fill: "<<Energies[i]<<Angles[i]<<"\n";
      LongitudinalProfile->Fill(Energies[i],Angles[i]);
    }
   }

    TFile *out = new TFile("LongProf.root","RECREATE");
    out->cd();
    TCanvas *c1 = new TCanvas("c1","The Ntuple canvas",1500,900);
    c1->Divide(1,1);
    c1->cd();

    LongitudinalProfile->GetXaxis()->SetTitle("E absorbed silicon [MeV]");
    LongitudinalProfile->GetYaxis()->SetTitle("#phi [rad]");
    LongitudinalProfile->Draw("colz");

    c1->Update();
    std::cout << "Press Enter to Continue";
    std::cin.ignore();
    c1->SaveAs("Phi_AbsorbedE.png");
    c1->Write();
    out->Close();
    f->Close();


}
void EnergyTower2D_Tracks(const Char_t *datapath= "")
{
   // gROOT->SetStyle("ATLAS");

   TFile *f = new TFile(datapath);

   TTree *Hits = (TTree*)f->Get("Hits");
   TTree *Tracks = (TTree*)f->Get("Tracks_true");

   std::vector<double>* Hits_xCellpos = 0;
   std::vector<double>* Hits_yCellpos = 0;
   std::vector<double>* Hits_zCellpos = 0;
   std::vector<double>* Hits_eHit = 0;
   std::vector<double>* Hits_Sensor = 0;

   Hits->SetBranchAddress("Hits_xCellpos",&Hits_xCellpos);
   Hits->SetBranchAddress("Hits_yCellpos",&Hits_yCellpos);
   Hits->SetBranchAddress("Hits_zCellpos",&Hits_zCellpos);
   Hits->SetBranchAddress("Hits_Sensor",&Hits_Sensor);
   Hits->SetBranchAddress("Hits_eHit",&Hits_eHit);

   double Track_x;
   double Track_y;

   Tracks->SetBranchAddress("x",&Track_x);
   Tracks->SetBranchAddress("y",&Track_y);

   TH2F* hist_s1 = new TH2F("hist_s1","Energies (joined by 4 sectors)",110,90,640,11,-27.5,27.5);
   // TH1F* hist_s2 = new TH1F("hist","Energies (joined by 4 sectors)",220,90,640,22,-2.75,2.75);

   Int_t nentries = (Int_t)Hits->GetEntries();
   std::cout<<"Entries: "<<nentries<<" \n";
   for (Int_t i=0;i<nentries;i++)
   {
     Hits->GetEntry(i);
     int Number_of_Hits = Hits_xCellpos->size();
     if(Number_of_Hits!=0){
      for (Int_t j=0;j<Number_of_Hits;j++){
        // int layer =  ((int) (*Hits_zCell)[j]) - 1;
        int sensor = ((int) (*Hits_Sensor)[j]);
        // std::cout<<sensor<<"\n";
        if(sensor==2)
        {
          hist_s1->Fill((*Hits_xCellpos)[j],(*Hits_yCellpos)[j],(*Hits_eHit)[j]);
          std::cout<<(*Hits_xCellpos)[j]<<"\t"<<(*Hits_yCellpos)[j]<<"\t"<<(*Hits_eHit)[j]<<"\n";
        }


       }
      }
     }



    TFile *out = new TFile("TowerEnergy2D_Tracks.root","RECREATE");
    out->cd();
    TCanvas *c1 = new TCanvas("c1","The Ntuple canvas",1500,900);
    c1->Divide(1,1);
    c1->cd();
    c1->SetLogz();
    hist_s1->SetTitle("Energy Per Tower");
    hist_s1->GetXaxis()->SetTitle("x [mm]");
    hist_s1->GetYaxis()->SetTitle("y [mm]");
    hist_s1->Draw("colz");

    int tr = 0;
    //Tracks
    std::vector<TMarker*> Marks;
    nentries = (Int_t)Tracks->GetEntries();
    std::cout<<"Entries: "<<nentries<<" \n";
    for (Int_t i=0;i<nentries;i++)
    {
      Tracks->GetEntry(i);
      std::cout<<Track_x<<"\t"<<Track_y<<"\n";
      TMarker *m = new TMarker(Track_x, Track_y, 20);
      Marks.push_back(m);
      Marks[i]->SetMarkerColor(2);
      Marks[i]->SetMarkerSize(1);
      Marks[i]->Draw("same");
      // std::cout << "Press Enter to Continue";
      // std::cin.ignore();
      if(Track_x>0.0) tr++;
    }

    auto legend = new TLegend(0.7,0.5,0.8,0.6);
    legend->AddEntry("",Form("Tracks: %d",tr),"");
    // legend->AddEntry("",Form("number of e^{+}e^{-} in calorimeter: %d",(Int_t)Tracks_true->GetEntries()),"");
    legend->SetBorderSize(0);
    legend->SetTextFont(42);
    legend->SetTextSize(0.04);
    legend->Draw("");

    c1->Update();
    std::cout << "Press Enter to Continue";
    std::cin.ignore();
    c1->SaveAs("./TowerEnergy2D_Tracks.pdf");
    c1->Write();
    out->Close();
    f->Close();
}

void PlotLCOut(const Char_t *datapath= "")
{
  // const Char_t *datapath = argv[0];
  // std::cout << argv[0];
  // Build_2DHIts(datapath);
  // EnergyOfx_reco(datapath);
  // EnergyOfx_detector(datapath);
  // Build_Occupancy(datapath);
  // build_tracks_flow(datapath);
  // EnergyProfile(datapath);
  // EnergyTower4(datapath);
  // EnergyTower_andInputParticle(datapath);
  //LongitudinalProfile(datapath);
  EnergyTower2D_Tracks(datapath);
}
