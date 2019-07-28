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

LCRootOut::LCRootOut()
{ 
  RootOutFile = "Default.root";
  std::cout<<"@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ Root init 1\n";
  // fHMessanger = new HistoMessenger(this);
   // G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
}

LCRootOut::LCRootOut(const G4String name )
{
  RootOutFile = name;
  std::cout<<"@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ Root init 2\n";
  // fHMessanger = new HistoMessenger(this);
//    G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
}

LCRootOut::~LCRootOut()
{ 
  G4cout << " LCRootOut deleted ! " << G4endl;
  delete G4AnalysisManager::Instance();
  // delete analysisManager;
}

void LCRootOut::CreateNewTree()
{
  std::cout<<"@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ Create New Tree\n"; 


  // Create or get analysis manager
  // The choice of analysis technology is done via selection of a namespace
  // in HistoManager.hh
 
  // analysisManager->SetFileName(RootOutFile);
  // analysisManager->SetVerboseLevel(1);
  // analysisManager->SetActivation(true);    // enable inactivation of histograms

  // analysisManager->CreateNtuple("Events", "Events for LumCaL in LUXE");
  // analysisManager->CreateNtupleDColumn(0,"vX/D"); //0
  // analysisManager->CreateNtupleDColumn(0,"vY/D"); //1
  // analysisManager->CreateNtupleDColumn(0,"vZ/D"); //2
  // analysisManager->CreateNtupleIColumn(0,"numPrim/I"); //3
  // analysisManager->CreateNtupleIColumn(0,"numHits/I"); //4
  // analysisManager->CreateNtupleDColumn(0,"Etot[1]/D"); //5
  // analysisManager->CreateNtupleDColumn(0,"Emax/D"); //6
  // analysisManager->FinishNtuple(0);
  
  // analysisManager->CreateNtuple("Tracks", "Tracks for LumCaL in LUXE");
  // analysisManager->CreateNtupleDColumn(1, "Tracks_pX",Tracks_pX); //double pX; 0
  // analysisManager->CreateNtupleDColumn(1, "Tracks_PY",Tracks_pY); //double pY; 1
  // analysisManager->CreateNtupleDColumn(1, "Tracks_PZ",Tracks_pZ); //double pZ; 2
  // analysisManager->CreateNtupleIColumn(1, "Tracks_ID",Tracks_ID); //int    ID; 3
  // analysisManager->CreateNtupleIColumn(1, "Tracks_PDG",Tracks_PDG); //int   PDG; 4
  // analysisManager->FinishNtuple(1);

  // analysisManager->CreateNtuple("Hits", "Hits for LumCaL in LUXE");
  // analysisManager->CreateNtupleIColumn(2,"Hits_cellID",Hits_cellID);	// int    cellID;
  // analysisManager->CreateNtupleDColumn(2,"Hits_eHit",Hits_eHit);	//    double eHit; 0 
  // analysisManager->CreateNtupleDColumn(2,"Hits_xCell",Hits_xCell);	//    double xCell; 1
  // analysisManager->CreateNtupleDColumn(2,"Hits_yCell",Hits_yCell);	//    double yCell; 2
  // analysisManager->CreateNtupleDColumn(2,"Hits_zCell",Hits_zCell);	//    double zCell; 3
  // analysisManager->CreateNtupleDColumn(2,"Hits_xHit",Hits_xHit);	//    double xHit; 4
  // analysisManager->CreateNtupleDColumn(2,"Hits_yHit",Hits_yHit);	//    double yHit; 5
  // analysisManager->CreateNtupleDColumn(2,"Hits_zHit",Hits_zHit);	//    double zHit; 6
  // analysisManager->CreateNtupleDColumn(2,"Hits_TOF",Hits_TOF);	//    double TOF; 7
  // analysisManager->FinishNtuple(2);

}

void LCRootOut::SetAddresses()
{
	std::cout<<"@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ SetAddresses\n";
}

void LCRootOut::Init()
{
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  analysisManager->SetFileName(RootOutFile);
  analysisManager->SetVerboseLevel(1);
  analysisManager->SetActivation(true);    // enable inactivation of histograms

  analysisManager->CreateNtuple("Events", "Events for LumCaL in LUXE");
  analysisManager->CreateNtupleDColumn(0,"vX/D"); //0
  analysisManager->CreateNtupleDColumn(0,"vY/D"); //1
  analysisManager->CreateNtupleDColumn(0,"vZ/D"); //2
  analysisManager->CreateNtupleIColumn(0,"numPrim/I"); //3
  analysisManager->CreateNtupleIColumn(0,"numHits/I"); //4
  analysisManager->CreateNtupleDColumn(0,"Etot[1]/D"); //5
  analysisManager->CreateNtupleDColumn(0,"Emax/D"); //6
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
  analysisManager->FinishNtuple(2);
  
  analysisManager->OpenFile();
}

void LCRootOut::ProcessEvent(const G4Event* event, LCHitsCollection *collection)
{
	// G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();

    std::cout<<"@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ ProcessEvent\n";
    numHits = 0; 
    numPrim = 0;
    Etot[0] = 0. ;
    Etot[1] = 0. ;
    Emax  = 0. ;

    //
    // get all primary MC particles
    G4int nv = event->GetNumberOfPrimaryVertex();
    G4int k=0;
    for (int v = 0; v < nv ; v++) {
    
    G4PrimaryVertex *pv = event->GetPrimaryVertex(v);
	vX = pv->GetX0();
	vY = pv->GetY0();
	vZ = pv->GetZ0();

	G4PrimaryParticle *pp = pv->GetPrimary();


    while (pp) {
	 
	  double pX = (pp->GetMomentum()).getX();
	  double pY = (pp->GetMomentum()).getY();
	  double pZ = (pp->GetMomentum()).getZ();
	  int ID =  pp->GetTrackID();
	  int PDG = pp->GetPDGcode();

	  //
	  // pTracks->push_back( t );
	  FillTracks_pX(pX);
	  FillTracks_pY(pY);
	  FillTracks_pZ(pZ);
	  FillTracks_ID(ID);
	  FillTracks_PDG(PDG);

	  k++;
	  pp = pp->GetNext();
        
    }
    }

    numPrim = k;
    //    G4cout<<" Number of primary particles: "<<numPrim << G4endl;

    if( collection ){
    numHits = collection->entries();
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


      cellID =(*collection)[i]->GetCellCode();
      G4int    side = (cellID >> 24) & 0xff ;
      eHit  = (*collection)[i]->GetEnergy();
      
      //hit.xCell = (*collection)[i]->GetXcell();
      //hit.yCell = (*collection)[i]->GetYcell();
      //hit.zCell = (*collection)[i]->GetZcell();
      xCell = (cellID>>8) & 0xFF;
      yCell = (cellID) & 0xFF;
      zCell = (cellID>>16) & 0xFF;
      xHit  = (*collection)[i]->GetXhit();
      yHit  = (*collection)[i]->GetYhit();
      zHit  = (*collection)[i]->GetZhit();
      TOF   = (*collection)[i]->GetTOF();

      bool cellInTestBeam = false; // to mask 2014 TB
      
      if ( ((xCell == 12)&&(yCell>49))  ||  ((xCell == 13)&&(yCell>45)) )
	  cellInTestBeam=true;
      
      if((Setup::Lcal_n_layers > 12) ||( (Setup::Lcal_n_layers <= 12)&&(cellInTestBeam) ) ){
		Etot[side-1] += eHit ;
		if ( Emax < eHit  ) Emax = eHit;
		
		// // pHits->push_back( hit );
		// Hits_cellID->push_back(cellID);
		// Hits_eHit->push_back(eHit);
		// Hits_xCell->push_back(xCell);
		// Hits_yCell->push_back(yCell);
		// Hits_zCell->push_back(zCell);
		// Hits_xHit->push_back(xHit);
		// Hits_yHit->push_back(yHit);
		// Hits_zHit->push_back(zHit);
		// Hits_TOF->push_back(TOF);

		// histManager->SetTracks(ftrackId);
		// histManager->SetTracksVtx(ftrackVtxPos);
		// histManager->SetTracksMomentum(ftrackMomentum);
		// histManager->SetTracksE(ftrackE);
		// histManager->SetTracksPDG(ftrackPDG);
		// histManager->SetTracksPhysProc(ftrackPProc);
		// histManager->SetTracksPTId(ftrackPTD);
		// analysisManager->AddNtupleRow(3);

		FillcellID(cellID);
		FilleHit(eHit);
		FillxCell(xCell);
		FillyCell(yCell);
		FillzCell(zCell);
		FillxHit(xHit);
		FillyHit(yHit);
		FillzHit(zHit);
		FillTOF(TOF);

      }
      //
      i++;
    }

    analysisManager->FillNtupleDColumn(0, 0, vX);
    analysisManager->FillNtupleDColumn(0, 1, vY);
    analysisManager->FillNtupleDColumn(0, 2, vZ);
    analysisManager->FillNtupleIColumn(0, 3, numPrim);
    analysisManager->FillNtupleIColumn(0, 4, numHits);
    analysisManager->FillNtupleDColumn(0, 5, Etot[0]);
    analysisManager->FillNtupleDColumn(0, 6, Emax);

    analysisManager->AddNtupleRow(0);
    analysisManager->AddNtupleRow(1);
	  analysisManager->AddNtupleRow(2);

	  // ClearData();
    // _LcalData->Fill();
    //
    // clear vectors
    // 
    // pTracks->clear();
    // pHits->clear();

	// Hits_cellID->clear();
	// Hits_eHit->clear();
	// Hits_xCell->clear();
	// Hits_yCell->clear();
	// Hits_zCell->clear();
	// Hits_xHit->clear();
	// Hits_yHit->clear();
	// Hits_zHit->clear();
	// Hits_TOF->clear();

	// Tracks_pX->clear();
	// Tracks_pY->clear();
	// Tracks_pZ->clear();
	// Tracks_ID->clear();
	// Tracks_PDG->clear();

    }


    /*    G4cout << "LCRootOut::ProcessEvent completed." << G4endl;
    G4cout << " Number of Hits "<< numHits << G4endl;
    G4cout << " Energy " << Etot1  / GeV << " [GeV]"<< G4endl; */
}

// void LCRootOut::ProcEventAccumulate( LCHitsCollection *collection )
// {
//   std::cout<<"@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ ProcEventAccumulate\n";
//   // - B.P. 29.05.2009
//   // this method is devoted to merge events from entire run
//   // into one tree. It is to be used solely for beam-background data.
//   // Primaries are not stored in the tree ( typically there is 10^5
//   // primaries per event ). The list of tracks contributing to
//   // a hit is neglected as well.  

//     if( collection ){

//     G4int       i = 0;
//     G4int   nHits = collection->entries();
 
//     while ( i < nHits){

//         G4int address = (*collection)[i]->GetCellCode();
// 	G4int    side = (address >> 24) & 0xff ;
// 	G4double eH  = (*collection)[i]->GetEnergy();
// 	unsigned int it = 0;
// 	while ( it < theUsedCells.size() ){
// 	  if( theUsedCells[it] == address ) break;
// 	  it++;
// 	}

// 	if ( it == theUsedCells.size() ) { // cellID not in the list, add

// 	  theUsedCells.push_back( address );
// 	  Hit_t hit;
// 	  hit.eHit  = eH;
// 	  hit.xCell = (*collection)[i]->GetXcell();
// 	  hit.yCell = (*collection)[i]->GetYcell();
// 	  hit.zCell = (*collection)[i]->GetZcell();
// 	  hit.xHit  = (*collection)[i]->GetXhit();
// 	  hit.yHit  = (*collection)[i]->GetYhit();
// 	  hit.zHit  = (*collection)[i]->GetZhit();
// 	  hit.TOF   = (*collection)[i]->GetTOF();
// 	  if ( Emax < eH ) Emax = eH;
// 	  pHits->push_back( hit ); 
// 	  numHits++;
// 	}else{        // cell in list, increase energy deposit only
// 	  (*pHits)[it].eHit += eH;
// 	  if ( (*pHits)[it].eHit > Emax ) Emax = (*pHits)[it].eHit;
// 	}
//         Etot[side-1] += eH ;
//         i++;
//     }

//     }

// }

void  LCRootOut::ClearData()
{
  std::cout<<"@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ Clear Data\n"; 
	
  Hits_cellID.clear();
	Hits_eHit.clear();
	Hits_xCell.clear();
	Hits_yCell.clear();
	Hits_zCell.clear();
	Hits_xHit.clear();
	Hits_yHit.clear();
	Hits_zHit.clear();
	Hits_TOF.clear();

	Tracks_pX.clear();
	Tracks_pY.clear();
	Tracks_pZ.clear();
	Tracks_ID.clear();
	Tracks_PDG.clear();
  	
  	// fNHits = 0;
}

void LCRootOut::End()
{
  std::cout<<"@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ End\n";
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  // end of run action
  // if ( Setup::AccumulateEvents ){          //fill tree, do not close the file
  //   _LcalData->Fill();                     // it will be closed by main()
  //   theUsedCells.clear();
  
  // }else{
	// _file->Write();
	// _file->Close();
	// G4cout << "LCRootOut::Closed file: " << Setup::RootFileName << G4endl;
	// delete _file;
	// _file = NULL;
	// _LcalData = NULL;
  
  // }
	// G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();  
	if ( analysisManager->IsActive() ) {    
	analysisManager->Write();
	analysisManager->CloseFile();
  }
  std::cout<<"@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ End of End\n";
}

// void FillTracks_pX(G4double pX_) {Tracks_pX.push_back(pX_)}
// void FillTracks_pY(G4double pY_) {Tracks_pY.push_back(pY_)}
// void FillTracks_pZ(G4double pZ_) {Tracks_pZ.push_back(pZ_)}
// void FillTracks_ID(G4int ID_) {Tracks_ID.push_back(ID_)}
// void FillTracks_PDG(G4int PDG_) {Tracks_PDG.push_back(PDG_)}

// void FillcellID(G4int cellID_) {Hits_cellID.push_back(cellID_)}
// void FilleHit(G4double eHit_) {Hits_eHit.push_back(eHit_)}
// void FillxCell(G4double xCell_) {Hits_xCell.push_back(xCell_)}
// void FillyCell(G4double yCell_) {Hits_yCell.push_back(yCell_)}
// void FillzCell(G4double zCell_) {Hits_zCell.push_back(zCell_)}
// void FillxHit(G4double xHit_) {Hits_xHit.push_back(xHit_)}
// void FillyHit(G4double yHit_) {Hits_yHit.push_back(yHit_)}
// void FillzHit(G4double zHit_) {Hits_zHit.push_back(zHit_)}
// void FillTOF(G4double TOF_) {Hits_TOF.push_back(TOF_)}