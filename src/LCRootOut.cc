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
#include "G4EventManager.hh"

#include "G4UIcmdWithAString.hh"

LCRootOut::LCRootOut() : fMessenger(0), RootOutFile("Default.root")
{
  fMessenger = new LCRootMesseger(this);
}

// LCRootOut::LCRootOut(const G4String name )
// {
//   RootOutFile = name;
// }

LCRootOut::~LCRootOut()
{
  delete G4AnalysisManager::Instance();
  delete RootOutFile;
  delete fMessenger;
}

void LCRootOut::Init()
{
  // std::cout << "@@@@@@@@@@@@Create Ntuple\n";
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  analysisManager->SetFileName(RootOutFile);
  // std::cout<<"Set File Name\n";
  // analysisManager->SetFileName(root_file_name+".root");
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
  analysisManager->CreateNtupleIColumn("eventID");
  analysisManager->FinishNtuple(0);

  analysisManager->CreateNtuple("Tracks", "Tracks for LumCaL in LUXE");
  analysisManager->CreateNtupleDColumn(1, "Tracks_pX",Tracks_pX); //double pX; 0
  analysisManager->CreateNtupleDColumn(1, "Tracks_PY",Tracks_pY); //double pY; 1
  analysisManager->CreateNtupleDColumn(1, "Tracks_PZ",Tracks_pZ); //double pZ; 2
  analysisManager->CreateNtupleIColumn(1, "Tracks_ID",Tracks_ID); //int    ID; 3
  analysisManager->CreateNtupleIColumn(1, "Tracks_PDG",Tracks_PDG); //int   PDG; 4
  analysisManager->CreateNtupleIColumn(1, "eventID");
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
  analysisManager->CreateNtupleIColumn(2,"Hits_Sensor",Hits_Sensor);  //    double TOF; 8
  analysisManager->CreateNtupleIColumn(2,"eventID");
  analysisManager->CreateNtupleDColumn(2,"Hits_xCellpos",Hits_xCellpos);	//    double xCell; 1
  analysisManager->CreateNtupleDColumn(2,"Hits_yCellpos",Hits_yCellpos);	//    double yCell; 2
  analysisManager->CreateNtupleDColumn(2,"Hits_zCellpos",Hits_zCellpos);	//    double zCell; 3

  analysisManager->FinishNtuple(2);

  analysisManager->CreateNtuple("Tracks_true", "Tracks in the first sensor");
  analysisManager->CreateNtupleDColumn(3,"x");
  analysisManager->CreateNtupleDColumn(3,"y");
  analysisManager->CreateNtupleDColumn(3,"z");
  analysisManager->CreateNtupleDColumn(3,"Mx");
  analysisManager->CreateNtupleDColumn(3,"My");
  analysisManager->CreateNtupleDColumn(3,"Mz");
  analysisManager->CreateNtupleDColumn(3,"E");
  analysisManager->CreateNtupleIColumn(3,"PDG");
  analysisManager->CreateNtupleIColumn(3,"xCell");
  analysisManager->CreateNtupleIColumn(3,"yCell");
  analysisManager->CreateNtupleIColumn(3,"zCell");
  analysisManager->CreateNtupleIColumn(3,"Sensor");
  analysisManager->CreateNtupleDColumn(3,"Weight"); //6
  analysisManager->CreateNtupleIColumn(3,"eventID");
  analysisManager->FinishNtuple(3);

  analysisManager->CreateNtuple("Tracks_testplane", "Tracks in the test plane");
  analysisManager->CreateNtupleDColumn(4,"x");
  analysisManager->CreateNtupleDColumn(4,"y");
  analysisManager->CreateNtupleDColumn(4,"z");
  analysisManager->CreateNtupleDColumn(4,"Mx");
  analysisManager->CreateNtupleDColumn(4,"My");
  analysisManager->CreateNtupleDColumn(4,"Mz");
  analysisManager->CreateNtupleDColumn(4,"E");
  analysisManager->CreateNtupleIColumn(4,"PDG");
  analysisManager->CreateNtupleDColumn(4,"Weight"); //8
  analysisManager->CreateNtupleIColumn(4,"eventID");
  analysisManager->FinishNtuple(4);

  analysisManager->OpenFile();
}

void LCRootOut::ProcessEvent(const G4Event* event, LCHitsCollection *collection)
{
    G4int eventID = event->GetEventID();

    // std::cout<<"Processed ProcessEvent with eventID: "<< eventID << "\n";
    // std::cout << "Press Enter to Continue";
    // std::cin.ignore();

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

    std::vector<G4double> Hits_xCellpos_;
    std::vector<G4double> Hits_yCellpos_;
    std::vector<G4double> Hits_zCellpos_;

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
    // std::cout<<numHits<<"\n";
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

      double xCellpos;
      double yCellpos;
      double zCellpos;


      cellID =(*collection)[i]->GetCellCode();
      // G4int    side = (cellID >> 24) & 0xff ;
      eHit  = (*collection)[i]->GetEnergy();

      xCellpos = (*collection)[i]->GetXcell();
      yCellpos = (*collection)[i]->GetYcell();
      zCellpos = (*collection)[i]->GetZcell();
      // std::cout << "readout position of hit: " << xCellpos << " "
      //                                          << yCellpos << " "
      //                                          << zCellpos << "\n";

      yCell = (cellID >> 8) & 0xFF;
      xCell = (cellID     ) & 0xFF;
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

      Hits_xCellpos_.push_back(xCellpos);
      Hits_yCellpos_.push_back(yCellpos);
      Hits_zCellpos_.push_back(zCellpos);

      // std::cout<<eHit<<"\n";

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
    analysisManager->FillNtupleIColumn(8, eventID);
    analysisManager->AddNtupleRow(0);

    Fill_Tracks_pX(Tracks_pX_);
    Fill_Tracks_pY(Tracks_pY_);
    Fill_Tracks_pZ(Tracks_pZ_);
    Fill_Tracks_ID(Tracks_ID_);
    Fill_Tracks_PDG(Tracks_PDG_);
    analysisManager->FillNtupleIColumn(1,5, eventID);
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
    // analysisManager->FillNtupleDColumn(2,7, weight_fromMC);
    analysisManager->FillNtupleIColumn(2,10, eventID);

    Fill_Hits_xCellpos(Hits_xCellpos_);
    Fill_Hits_yCellpos(Hits_yCellpos_);
    Fill_Hits_zCellpos(Hits_zCellpos_);

	  analysisManager->AddNtupleRow(2);

    }

}

void LCRootOut::TestPlaneFill(G4double x ,G4double y ,G4double z,G4double Mx ,G4double My ,G4double Mz,G4double E,G4int PDG){

  G4int evtNb = G4EventManager::GetEventManager()->GetConstCurrentEvent()->GetEventID();
  // std::cout<<"Processed TestPlaneFill with eventID: "<< evtNb << "\n";
  // std::cout << "Press Enter to Continue";
  // std::cin.ignore();

  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  analysisManager->FillNtupleDColumn(4,0, x);
  analysisManager->FillNtupleDColumn(4,1, y);
  analysisManager->FillNtupleDColumn(4,2, z);
  analysisManager->FillNtupleDColumn(4,3, Mx);
  analysisManager->FillNtupleDColumn(4,4, My);
  analysisManager->FillNtupleDColumn(4,5, Mz);
  analysisManager->FillNtupleDColumn(4,6, E);
  analysisManager->FillNtupleIColumn(4,7, PDG);
  analysisManager->FillNtupleDColumn(4,8, weight_fromMC);
  analysisManager->FillNtupleIColumn(4,9, evtNb);
  analysisManager->AddNtupleRow(4);

}


void LCRootOut::End()
{
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();

	// if ( analysisManager->IsActive() ) {
	analysisManager->Write();
	analysisManager->CloseFile();
}

LCRootMesseger::LCRootMesseger(LCRootOut* histoM) :G4UImessenger(), fHistManager(histoM), fRootName()
{
  fRootName = new G4UIcmdWithAString("/analysis/filename",this);
  fRootName->SetGuidance("Name of the root file for output");
  fRootName->SetParameterName("rootname",false);
  fRootName->AvailableForStates(G4State_PreInit, G4State_Idle);
}


LCRootMesseger::~LCRootMesseger()
{
  delete fRootName;
}


void LCRootMesseger::SetNewValue(G4UIcommand* command, G4String newValue)
{
  if (command == fRootName)
    { fHistManager->SetName(newValue);}
}
