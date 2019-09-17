/*
 * LCSensitiveDetector.cc
 * 2ModLumiCal
 *
 *  Created on: Mar 26, 2009
 *      Author: aguilar
 */

// time optimization
#include <time.h>

#include "LCSensitiveDetector.hh"
#include "LCHit.hh"
// #include "LCEventAction.hh"
#include "Setup.hh"

#include "G4SystemOfUnits.hh"

#include "G4AnalysisMessenger.hh"
#include "g4root.hh"
#include "G4EventManager.hh"
#include "G4Event.hh"

#include "GlobalVars.hh"

LCSensitiveDetector::LCSensitiveDetector(G4String sdname,
                                         G4double CalX_half,
                                         G4double CalY_half,
                                         G4double dCellX,
                                         G4double dCellY,
                                         G4bool   cellvirtual)
    : G4VSensitiveDetector(sdname),
      HCID(-1),
      SDName(sdname), // this sets string SDName to sdname
      hitsColl(0),
      VirtualCell( cellvirtual )
{

    collName = SDName+"_HC"; // not dynamic - name becomes LumiCalSD_HC
    // collectionName.insert(collName); // a G4VSensitiveDetector element
    collectionName.insert(collName);
    origin = G4ThreeVector();
    hitMap = new LCHitMap;

    // set primary args - this is how you tell LCSensDet
    // the geometric parameters of LumiCal
		// SetRhoMin(CalRhoMin);
		// SetPhiOffset(CalPhiOffset);
		// SetRhoCellDim(cellDimRho);
		// SetPhiCellDim( cellDimPhi);
		// SetNCellRho(nCellRho);
		// SetNCellPhi(nCellPhi);

    f_CalX_half = CalX_half;
    f_CalY_half = CalY_half;

    f_dCellX = dCellX;
    f_dCellY = dCellY;

    f_nCellX = (G4int)floor(2 * CalX_half / dCellX);
    f_nCellY = (G4int)floor(2 * CalY_half / dCellY);

    G4cout << "SD created <" << sdname << ">" << G4endl;
    std::cout  << "Parametrs -->>" << "\n"
               << "f_CalX_half:    " << f_CalX_half   << "\n"
               << "f_CalY_half: "    << f_CalY_half   << "\n"
               << "f_dCellX:     " << f_dCellX    << "\n"
               << "f_dCellY:     " << f_dCellY    << "\n"
               << "f_nCellX:     " << f_nCellX    << "\n"
               << "f_nCellY:     " << f_nCellY    << "\n"
               << "cellvirtual:  " << cellvirtual << "\n";

    // std::cout << "Press Enter to Continue";
    // std::cin.ignore();
}

LCSensitiveDetector::~LCSensitiveDetector()
{
    hitMap->clear();
    delete hitMap;
}

void LCSensitiveDetector::Initialize(G4HCofThisEvent *HCE)
{
    // std::cout<<"@@@@@@@@@@@@@@@@@@@@@@@ SD was initialized!\n";
    //Create a (G4) hit collection for this event
    // hitsColl = new LCHitsCollection(SDName, collectionName[0]);
    hitsColl = new LCHitsCollection(SDName, collName);

    if (HCID < 0) {
      HCID = G4SDManager::GetSDMpointer()->GetCollectionID(collName);
    }

    // Would be sufficient to add the collection at the end of the event, but
    // doing it here suppresses a warning during compilation
    HCE->AddHitsCollection(HCID, hitsColl);
    primaryID = 0;
    //Reset the hit map
    hitMap->clear();
}

//Pass the hit collection to the other G4 hits collection
//These collections will be transformed into the LCIO Collection you
//find in the output file
void LCSensitiveDetector::EndOfEvent(G4HCofThisEvent* HCE)
{
    // std::cout<<"@@@@@@@@@@@@@@@@@@@@@@@ SD was End!\n";
    LCHitMapIterator iter;

    // the hit map stores all the hits - they are inserted during ProcessHits()
    // This loop cycles through the hit map and adds them all to hitsColl
    for (iter=hitMap->begin(); iter != hitMap->end(); ++iter) {
        hitsColl->insert(iter->second);
    }

    hitMap->clear(); // hit map has already been copied to hitsColl
                     // but what about the cell ID? is that just lost?

    // Add hit collection to HCE
    HCE->AddHitsCollection(HCID, hitsColl);

}


G4bool LCSensitiveDetector::ProcessHits(G4Step *aStep, G4TouchableHistory *)
{
    // -------------- LOCAL VARIABLES AND CHECKS --------------

    // store energy deposition and momentum change from primary particles
    G4double edep = aStep->GetTotalEnergyDeposit();

    // make sure you aren't using geantinos
    if (edep <= 0 && aStep->GetTrack()->GetDefinition()->GetParticleType() != "geantino") return true;

    // PARTICLE INFORMATION ---------------------------------

    // Get parent ID
    G4int parentID = aStep->GetTrack()->GetParentID();
    // time of flight
    G4double ToF   = aStep->GetTrack()->GetGlobalTime();

     // Set primary particle track ID/PDG
    G4int PDG=aStep->GetTrack()->GetDefinition()->GetPDGEncoding();

    if ( parentID <= 0 ){
    G4int PID = aStep->GetTrack()->GetTrackID();
    if( PID != primaryID )
       {
           primaryID = PID;
           primaryPDG= PDG;
       }
     }

     // Get location (volumes & coordinates) with PreStepPoint->TouchableHandle
    const G4StepPoint *preStepPoint = aStep->GetPreStepPoint();
    const G4StepPoint *postStepPoint = aStep->GetPostStepPoint();
    G4ThreeVector GlobalHitPos = ( (preStepPoint->GetPosition())+(postStepPoint->GetPosition())) / 2.;
    //
    // std::cout << "GlobalHitPosition " << GlobalHitPos/mm << " mm " <<G4endl;
    // std::cout << "Press Enter to Continue";
    // std::cin.ignore();

    G4int LC_num = 0;
    G4int layer_num = 0;
    G4int sector_num = 0;
    G4int cell_num = 0;
    // G4int sensor_number = 0; // 1 - closer to the center
    //                       // 2 - farther from center


    // Use the touchable handle to get the volume hierarchy for the hit
    G4TouchableHandle theTouchable = preStepPoint->GetTouchableHandle();
    // LC_num     = theTouchable->GetHistory()->GetVolume(1)->GetCopyNo();
    // layer_num  = theTouchable->GetHistory()->GetVolume(2)->GetCopyNo();
    // sector_num = theTouchable->GetReplicaNumber(1);
    // cell_num   = theTouchable->GetReplicaNumber(0);                     G4double CalXMin,

    // std::cout << "touchable: " << LC_num << " "
    //                            << layer_num << " "
    //                            << sector_num << " "/data/Projects_physics/LUXE/IPstrong/Lists/bppp_17.5GeV_5mfoiltoIP_Enelas_1.0J.out
    //                            << cell_num << "\n";

    // std::cout << "Press Enter to Continue";
    // std::cin.ignore();


    // sensor_number = theTouchable->GetHistory()->GetVolume(4)->GetCopyNo();

    // Find the unique volume path for a cell which has a hit
    // Which copy of LumiCal?
    if (Setup::LcalTBeam)
    {
      LC_num = 1;
      layer_num  = theTouchable->GetHistory()->GetVolume(1)->GetCopyNo();
    }
    else{
      LC_num     = theTouchable->GetHistory()->GetVolume(1)->GetCopyNo();
      // Which layer inside LumiCal?
      layer_num  = theTouchable->GetHistory()->GetVolume(2)->GetCopyNo();
    }


    if ( !VirtualCell ){

        // Which sector inside the tile? (sectors are replicated *after* the cells)
        sector_num = theTouchable->GetReplicaNumber(1);
        // Which cell inside sector?
        cell_num   = theTouchable->GetReplicaNumber(0);



    }else{

        G4ThreeVector LocalHitPos = theTouchable->GetHistory()->GetTopTransform().TransformPoint(GlobalHitPos) ;
        //
        //G4cout << "LocalHitPosition " << LocalHitPos/mm <<" mm" << G4endl;
        //
        //~~~~~~~~~~~~~~~~~~~~~~LocalHitPos~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ rectamgle calorimeter START

        // Pay Attention that cordinates are reversed
        G4double Y = LocalHitPos.x() + f_CalY_half;
        G4double X = LocalHitPos.y() + f_CalX_half;

        cell_num = (G4int)floor( X / f_dCellX );
        sector_num = (G4int)floor( Y / f_dCellY );

        //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ rectamgle calorimeter END
        // std::cout << "LocalHitPosition " << LocalHitPos/mm <<" mm" << G4endl;
        // std::cout << "GlobalHitPosition " << GlobalHitPos/mm <<" mm" << G4endl;
        // std::cout << "My x y:" << X <<" " << Y << G4endl;

        // std::cout << "LCAL arm : " << LC_num
        //           << " Layer: "     << layer_num
        //           << " Sector: "    << sector_num
        //           << " Cell: "      << cell_num
        //           << G4endl;
        //
        // std::cout << "Press Enter to Continue";
        // std::cin.ignore();
    }

    // ENCODE THE HIT ---------------------------------

    // Use bitwise operations to store the location of a cell in a single
    // 32-bit word
    // 4 bytes, 4 levels in the volume tree - 1 byte per volumeGlobalCellPos
    // | LC number | Layer number | Sector number | Cell number |
    // note - this only works because no volume will have a number > 255 (2^8)
    cell_code cellID;
    cellID.id0 = 0; cellID.id1 = 0;   // 32 0's in a row per member
        // LumiCal only uses id0; id1 is a legacy from the Mokka version

    cellID.id0 |= ((cell_num+1)   <<  0); // store the cell # in the lowest 8 digits
    cellID.id0 |= ((sector_num+1) <<  8); // shift the sector # to the next byte up
    cellID.id0 |= (layer_num  << 16); // shift the layer # to the next byte up
    cellID.id0 |= (LC_num     << 24); // shift the LumiCal # to the highest byte

    // get local and global position using touchable
    // check if the cell has been hit before, and check the indices
    if (!FindHit(cellID.id0, edep, primaryID, PDG) &&
        !( (cell_num > f_nCellX) || (sector_num > f_nCellY) ) )  {

        // Global hit coordinates

        // FG: the local coordinate system of the cell replica is the same as for the disk
        //     containing the the pad - rotated by phi of the phi replica, i.e.
        //     the origin is at the center of the disk
        //     hence the position of the cell center in its coordinate system is
        //     given by:
  //       G4ThreeVector localCellPos(CalRhoMin+((G4double)cell_num)*cellDimRho, 0. , 0. );
  //
	// if ( VirtualCell ) localCellPos.setPhi( (0.5000 + ((G4double)sector_num-1))*cellDimPhi );Global Cell Position (212,0.5,5021.92) mm


	// else localCellPos.setPhi(0.5000 * cellDimPhi);
G4ThreeVector LocalHitPos = theTouchable->GetHistory()->GetTopTransform().TransformPoint(GlobalHitPos);

G4ThreeVector localCellPos( f_dCellY/2.0 + ((G4double)sector_num)*f_dCellY - f_CalY_half,
                            f_dCellX/2.0 + ((G4double)cell_num)*f_dCellX - f_CalX_half,
                            LocalHitPos.z() );

	//
	//		  G4cout << "LCP1 " << localCellPos << G4endl;
	//

        // compute GlobalCellPos based on touchable with localCellPos
G4ThreeVector GlobalCellPos = theTouchable->GetHistory()->GetTopTransform().Inverse().TransformPoint(localCellPos);
	//
	//		  G4cout << "GCP1 " << GlobalCellPos << G4endl;
	//
    // std::cout << "OUTPUT:\n";
    // std::cout << f_dCellY/2 << " " << (G4double)sector_num << " " << f_dCellY << " " << f_CalY_half << "\n";
    // std::cout << "pad cell: " << sector_num << " " << cell_num << "\n";
    // std::cout << "Global Cell Position " << GlobalCellPos/mm << " mm " <<G4endl;
    // std::cout << "Local Cell Position " << localCellPos/mm <<" mm" << G4endl;
    // std::cout << "Press Enter to Continue";
    // std::cin.ignore();

        // assert id w/in valid range using bitwise ops
        // (0xff = 255 in hex)
        assert((cellID.id0 & 0xff)<=120);


        // Use all this information to create an LCHit instance
        LCHit *theHit = new LCHit(GlobalHitPos.x(),    // global hit position
                                  GlobalHitPos.y(),    //
                                  GlobalHitPos.z(),     //
                                  GlobalCellPos.x(),   // GlobalCellPosRho
                                  GlobalCellPos.y(),   // GlobalCellPosPhi
                                  GlobalCellPos.z(),   // GlobalCellPosZ
                                  edep,                // edep
                                  primaryID,           // track number
                                  primaryPDG,          // PDG encoding
                                  cellID,              // Cell ID code
                                  ToF);                // time since event begin


        //Insert the hit and cellID into the hitmap
        hitMap->insert(LCHitMapPair(cellID.id0, theHit));
    }

    // ~~~~~~~~~~~~~~~~~~Added "Tracker" as the first layer of silicon~~~~~~~~~~~~~~~~~~~~~~~~~~~
    auto event = G4EventManager::GetEventManager()->GetConstCurrentEvent();
    // auto evtNb = event->GetEventID();
    auto eventNb = event->GetEventID();

    G4Track* theTrack = aStep->GetTrack();
    // check if it is entering the test volume
    G4StepPoint* thePrePoint = aStep->GetPreStepPoint();
    G4VPhysicalVolume* thePrePV = thePrePoint->GetPhysicalVolume();
    G4LogicalVolume* thePreLV = thePrePV->GetLogicalVolume();

    G4LogicalVolume* thePostLV = 0;
    G4StepPoint* thePostPoint = aStep->GetPostStepPoint();
    if (thePostPoint) {
       G4VPhysicalVolume* thePostPV = thePostPoint->GetPhysicalVolume();
       if (thePostPV) thePostLV = thePostPV->GetLogicalVolume();
    }

    if (thePostLV != thePreLV and
        theTrack->GetParentID()==0 and
        PDG!=22 and
        layer_num==1) {

          G4double x = theTrack->GetPosition().x();
          G4double y = theTrack->GetPosition().y();
          G4double z = theTrack->GetPosition().z();

          G4double E = theTrack->GetKineticEnergy();

          G4ThreeVector theMomentumDirection = theTrack->GetDynamicParticle()->GetMomentumDirection();

          G4double Mx = theMomentumDirection.x();
          G4double My = theMomentumDirection.y();
          G4double Mz = theMomentumDirection.z();

          // std::cout << "PDG" <<PDG <<" "<< theTrack->GetParentID() <<"\n";
          // std::cout<<"SteppingAction: "<< x << " " << y << " " << z<<"\n";
          // std::cout << "Press Enter to Continue";
          // std::cin.ignore();


          G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
          analysisManager->FillNtupleDColumn(3,0, x);
          analysisManager->FillNtupleDColumn(3,1, y);
          analysisManager->FillNtupleDColumn(3,2, z);
          analysisManager->FillNtupleDColumn(3,3, Mx);
          analysisManager->FillNtupleDColumn(3,4, My);
          analysisManager->FillNtupleDColumn(3,5, Mz);
          analysisManager->FillNtupleDColumn(3,6, E);
          analysisManager->FillNtupleIColumn(3,7, PDG);
          analysisManager->FillNtupleIColumn(3,8, cell_num);
          analysisManager->FillNtupleIColumn(3,9, sector_num);
          analysisManager->FillNtupleIColumn(3,10,layer_num);
          analysisManager->FillNtupleIColumn(3,11,LC_num);
          analysisManager->FillNtupleDColumn(3,12,weight_fromMC);
          analysisManager->FillNtupleIColumn(3,13,eventNb);
          analysisManager->AddNtupleRow(3);

    }


    return true;
}

G4bool LCSensitiveDetector::FindHit(G4int    cellID,
                                    G4double edep,
                                    G4int    TID,
                                    G4int    PDG)
{
    /*Parameter explanations:
     *  cellID = unique cell identifier
     *  edep = amount of energy being added to the calorimeter cell
     *  TID = track ID of the particle (from Step)
     *  PDG = particle data group identifier code
     */

    //fg: use find method of map !
    // the find method of map returns:
        // successful   - a pointer to the appropriate pair
        // unsuccessful - a pointer to the end marker ( end() )
    LCHitMapIterator iter = hitMap->find(cellID);

    // If you find the hit already in hitMap (as indexed by the cellID),
    // increment its energy!
    if(iter != hitMap->end()) {
         (iter->second)->AddEdep(TID, PDG, edep) ;
         return true;
     }
    return false;
}
