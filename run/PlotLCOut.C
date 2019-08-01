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
#include <vector>

void PlotLCOut(const Char_t *datapath= "") {

   TFile *f = new TFile(datapath);

   TTree *Events = (TTree*)f->Get("Events");
   TTree *Track = (TTree*)f->Get("Tracks");
   TTree *Hits = (TTree*)f->Get("Hits");

   double Weight = 1;

   std::vector<int>* Track_PDG = 0;

   std::vector<double>* Hits_xHit = 0;
   std::vector<double>* Hits_yHit = 0; 
   std::vector<double>* Hits_zCell = 0;
   
   Events->SetBranchAddress("Weight",&Weight);

   Track->SetBranchAddress("Tracks_PDG",&Track_PDG);
   
   Hits->SetBranchAddress("Hits_xHit",&Hits_xHit);
   Hits->SetBranchAddress("Hits_yHit",&Hits_yHit);
   Hits->SetBranchAddress("Hits_zCell",&Hits_zCell);


   TH2F *hist_1 = new TH2F("hist_1","Occupancy of Hits layer 1",700,-300.,300.,600,-100.,100.);
   TH2F *hist_2 = new TH2F("hist_2","Occupancy of Hits layer 2",700,-300.,300.,600,-100.,100.);
   TH2F *hist_3 = new TH2F("hist_3","Occupancy of Hits layer 3",700,-300.,300.,600,-100.,100.);
   TH2F *hist_4 = new TH2F("hist_4","Occupancy of Hits layer 4",700,-300.,300.,600,-100.,100.);
   TH2F *hist_5 = new TH2F("hist_5","Occupancy of Hits layer 5",700,-300.,300.,600,-100.,100.);


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

     if(Number_of_Hits!=0 and (*Track_PDG)[0]==22){
     // if(Number_of_Hits!=0){

      for (Int_t j=0;j<Number_of_Hits;j++){
         switch ((int) (*Hits_zCell)[j])
         {
            case 1:
            hist_1->Fill((*Hits_xHit)[j],(*Hits_yHit)[j],Weight);
            break;

            case 2:
            hist_2->Fill((*Hits_xHit)[j],(*Hits_yHit)[j],Weight);
            break;

            case 3:
            hist_3->Fill((*Hits_xHit)[j],(*Hits_yHit)[j],Weight);
            break;

            case 4:
            hist_4->Fill((*Hits_xHit)[j],(*Hits_yHit)[j],Weight);
            break;

            case 5:
            hist_5->Fill((*Hits_xHit)[j],(*Hits_yHit)[j],Weight);
            break;
         }
      }
     }
   }


   TCanvas *c1 = new TCanvas("c1","The Ntuple canvas",2000,2000);
   c1->Divide(2,5);
   c1->cd();

   c1->cd(1);
   hist_1->Draw("colz");

   c1->cd(2);
   hist_2->Draw("colz");

   c1->cd(3);
   hist_3->Draw("colz");

   c1->cd(4);
   hist_4->Draw("colz");

   c1->cd(5);
   hist_5->Draw("colz");
   
   c1->Update();
   c1->SaveAs("./test_photons.png");
}