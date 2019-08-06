/*
 * LCRootOut.cc
 *
 *  Created on: Jul 27, 2019
 *      Author: shchedroloisiev
 */

#include "Setup.hh"
#include "DetectorConstruction.hh"
#include "tls.hh"
#include "G4String.hh"
#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4ios.hh"

#include "LCSensitiveDetector.hh"
#include "LCHit.hh"
#include "G4AnalysisMessenger.hh"
#include "LCRootOut.hh"

#include "G4SystemOfUnits.hh"
#include "GlobalVars.hh"
#include "PrimaryGeneratorAction.hh"


LCRootOut::LCRootOut()
{ 
  RootOutFile = "Default.root";
}

LCRootOut::LCRootOut(const G4String name )
{
  RootOutFile = name;
}

LCRootOut::~LCRootOut()
{ 
  delete G4AnalysisManager::Instance();
}

void LCRootOut::Init()
{
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  //analysisManager->SetFileName(RootOutFile);
  std::cout<<"Set File Name\n";
  analysisManager->SetFileName(root_file_name+".root");
  analysisManager->SetVerboseLevel(1);
  analysisManager->SetActivation(true);    // enable inactivation of histograms

  analysisManager->CreateNtuple("Events", "Events for LumCaL in LUXE");
  analysisManager->CreateNtupleDColumn("vX"); //0
  analysisManager->CreateNtupleDColumn("vY"); //1
  analysisManager->CreateNtupleDColumn("vZ"); //2
  analysisManager->CreateNtupleIColumn("numPrim"); //3
  analysisManager->CreateNtupleIColumn("numHit"); //4
  analysisManager->CreateNtupleDColumn("Etot_1"); //5
  analysisManager->CreateNtupleDColumn("Emax"); //6
  analysisManager->CreateNtupleDColumn("Weight"); //6
  analysisManager->FinishNtuple(0);
  
  analysisManager->CreateNtuple("Tracks", "Tracks for LumCaL in LUXE");
  analysisManager->CreateNtupleDColumn(1, "Tracks_pX",Tracks_pX); //double pX; 0
  analysisManager->CreateNtupleDColumn(1, "Tracks_PY",Tracks_pY); //double pY; 1
  analysisManager->CreateNtupleDColumn(1, "Tracks_PZ",Tracks_pZ); //double pZ; 2
  analysisManager->CreateNtupleIColumn(1, "Tracks_ID",Tracks_ID); //int    ID; 3
  analysisManager->CreateNtupleIColumn(1, "Tracks_PDG",Tracks_PDG); //int   PDG; 4
  analysisManager->FinishNtuple(1);

  analysisManager->CreateNtuple("Hits", "Hits for LumCaL in LUXE");
  analysisManager->CreateNtupleIColumn(2,"Hits_cellID",Hits_cellID);	// int    cellID;
  analysisManager->CreateNtupleDColumn(2,"Hits_eHit",Hits_eHit);	//    double eHit; 0 
  analysisManager->CreateNtupleDColumn(2,"Hits_xCell",Hits_xCell);	//    double xCell; 1
  analysisManager->CreateNtupleDColumn(2,"Hits_yCell",Hits_yCell);	//    double yCell; 2
  analysisManager->CreateNtupleDColumn(2,"Hits_zCell",Hits_zCell);	//    double zCell; 3
  analysisManager->CreateNtupleDColumn(2,"Hits_xHit",Hits_xHit);	//    double xHit; 4
  analysisManager->CreateNtupleDColumn(2,"Hits_yHit",Hits_yHit);	//    double yHit; 5
  analysisManager->CreateNtupleDColumn(2,"Hits_zHit",Hits_zHit);	//    double zHit; 6
  analysisManager->CreateNtupleDColumn(2,"Hits_TOF",Hits_TOF);	//    double TOF; 7
  analysisManager->CreateNtupleIColumn(2,"Hits_Sensor",Hits_Sensor);  //    double TOF; 7
  analysisManager->FinishNtuple(2);

  analysisManager->CreateNtuple("Tracks_true", "Tracks before first absorber");
  analysisManager->CreateNtupleDColumn(3,"x");
  analysisManager->CreateNtupleDColumn(3,"y");
  analysisManager->CreateNtupleDColumn(3,"z");
  analysisManager->CreateNtupleDColumn(3,"E");
  analysisManager->CreateNtupleIColumn(3,"PDG");
  analysisManager->FinishNtuple(3);

  
  analysisManager->OpenFile();
}

void LCRootOut::ProcessEvent(const G4Event* event, LCHitsCollection *collection)
{
    G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();

    double numHits = 0; 
    double numPrim = 0;
    double Etot = 0.0 ;
    double Emax  = 0. ;

    std::vector<G4double> Tracks_pX_;
    std::vector<G4double> Tracks_pY_;
    std::vector<G4double> Tracks_pZ_;
    std::vector<G4int> Tracks_ID_;
    std::vector<G4int> Tracks_PDG_;

    std::vector<G4int> Hits_cellID_;
    std::vector<G4double> Hits_eHit_;
    std::vector<G4double> Hits_xCell_;
    std::vector<G4double> Hits_yCell_;
    std::vector<G4double> Hits_zCell_;
    std::vector<G4double> Hits_xHit_;
    std::vector<G4double> Hits_yHit_;
    std::vector<G4double> Hits_zHit_;
    std::vector<G4double> Hits_TOF_;
    std::vector<G4int> Hits_Sensor_;

    //
    // get all primary MC particles
    G4int nv = event->GetNumberOfPrimaryVertex();
    G4int k=0;

    double vX,vY,vZ;

    for (int v = 0; v < nv ; v++) {
      // Take trajectory before absorber
      // auto track_vec = event->GetTrajectoryContainer()->GetVector();
      // std::cout<<event->GetTrajectoryContainer()->GetVector()<<" "<<event->GetTrajectoryContainer()->GetVector()->size()<<"\n";
    
      G4PrimaryVertex *pv = event->GetPrimaryVertex(v);
    	vX = pv->GetX0();
    	vY = pv->GetY0();
    	vZ = pv->GetZ0();

  	  G4PrimaryParticle *pp = pv->GetPrimary();


      while (pp) {
  	 
    	  Tracks_pX_.push_back(pp->GetMomentum().getX());
    	  Tracks_pY_.push_back(pp->GetMomentum().getY());
    	  Tracks_pZ_.push_back(pp->GetMomentum().getZ());
    	  Tracks_ID_.push_back(pp->GetTrackID());
    	  Tracks_PDG_.push_back(pp->GetPDGcode());

    	  k++;
    	  pp = pp->GetNext();
            
      }
      }

    numPrim = k;

    if( collection ){

    numHits = collection->entries();
    //std::cout<<numHits<<"\n";
    G4int i = 0;

    while ( i < numHits){
      // Hit_t hit;
      int    cellID;
      double eHit;
      double xCell;
      double yCell;
      double zCell;
      double xHit;
      double yHit;
      double zHit;
      double TOF;
      double Sensor;


      cellID =(*collection)[i]->GetCellCode();
      // G4int    side = (cellID >> 24) & 0xff ;
      eHit  = (*collection)[i]->GetEnergy();
      
      //hit.xCell = (*collection)[i]->GetXcell();
      //hit.yCell = (*collection)[i]->GetYcell();
      //hit.zCell = (*collection)[i]->GetZcell();
      xCell = (cellID >> 8) & 0xFF;
      yCell = (cellID     ) & 0xFF;
      zCell = (cellID >> 16) & 0xFF;
      Sensor =  (cellID >> 24) & 0xFF;
      xHit  = (*collection)[i]->GetXhit();
      yHit  = (*collection)[i]->GetYhit();
      zHit  = (*collection)[i]->GetZhit();
      TOF   = (*collection)[i]->GetTOF();


		  Etot += eHit ;
		  if ( Emax < eHit  ) Emax = eHit;
		
      Hits_cellID_.push_back(cellID);
      Hits_eHit_.push_back(eHit);
      Hits_xCell_.push_back(xCell);
      Hits_yCell_.push_back(yCell);
      Hits_zCell_.push_back(zCell);
      Hits_xHit_.push_back(xHit);
      Hits_yHit_.push_back(yHit);
      Hits_zHit_.push_back(zHit);
      Hits_TOF_.push_back(TOF);
      Hits_Sensor_.push_back(Sensor);

      
      i++;
    }

    analysisManager->FillNtupleDColumn(0, vX);
    analysisManager->FillNtupleDColumn(1, vY);
    analysisManager->FillNtupleDColumn(2, vZ);
    analysisManager->FillNtupleIColumn(3, numPrim);
    analysisManager->FillNtupleIColumn(4, numHits);
    analysisManager->FillNtupleDColumn(5, Etot);
    analysisManager->FillNtupleDColumn(6, Emax);
    analysisManager->FillNtupleDColumn(7, weight_fromMC);
    analysisManager->AddNtupleRow(0);

    Fill_Tracks_pX(Tracks_pX_);
    Fill_Tracks_pY(Tracks_pY_);
    Fill_Tracks_pZ(Tracks_pZ_);
    Fill_Tracks_ID(Tracks_ID_);
    Fill_Tracks_PDG(Tracks_PDG_);
    analysisManager->AddNtupleRow(1);

    Fill_Hits_cellID(Hits_cellID_);
    Fill_Hits_eHit(Hits_eHit_);
    Fill_Hits_xCell(Hits_xCell_);
    Fill_Hits_yCell(Hits_yCell_);
    Fill_Hits_zCell(Hits_zCell_);
    Fill_Hits_xHit(Hits_xHit_);
    Fill_Hits_yHit(Hits_yHit_);
    Fill_Hits_zHit(Hits_zHit_);
    Fill_Hits_TOF(Hits_TOF_);
    Fill_Hits_Sensor(Hits_Sensor_);
	  analysisManager->AddNtupleRow(2);

    }

}

void LCRootOut::TestPlaneFill(G4double x ,G4double y ,G4double z,G4double E,G4int PDG){

  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  analysisManager->FillNtupleDColumn(3,0, x);
  analysisManager->FillNtupleDColumn(3,1, y);
  analysisManager->FillNtupleDColumn(3,2, z);
  analysisManager->FillNtupleDColumn(3,3, E);
  analysisManager->FillNtupleIColumn(3,4, PDG);
  analysisManager->AddNtupleRow(3);

}


void LCRootOut::End()
{
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();

	// if ( analysisManager->IsActive() ) {
	analysisManager->Write();
	analysisManager->CloseFile();
}
