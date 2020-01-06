//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
/// \file field/field04/src/F04SteppingAction.cc
/// \brief Implementation of the F04SteppingAction class
//

#include "G4Track.hh"

#include "SteppingAction.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "GlobalVars.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::SteppingAction(LCRootOut *RO)
 : G4UserSteppingAction(),
   fTestPlaneVolume(0)
{
   RootOut = RO;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::~SteppingAction()
{

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SteppingAction::UserSteppingAction(const G4Step* theStep)
{

  G4Track* theTrack = theStep->GetTrack();
  G4int PDG = theTrack->GetDynamicParticle()->GetPDGcode();
  // if(PDG!=-11 and PDG!=11){
  //   std::cout << "initial non electron detected!\n" << PDG << "\n";
  //   G4ThreeVector theMomentumDirection = theTrack->GetDynamicParticle()->GetMomentumDirection();
  //   G4double Mx = theMomentumDirection.x();
  //   G4double My = theMomentumDirection.y();
  //   G4double Mz = theMomentumDirection.z();
  //   std::cout<<"Moment: "<< Mx << " " << My << " " << Mz<<"\n";
  //
  //   G4ThreeVector vertex_pos = theTrack->GetVertexPosition();
  //   G4double x = vertex_pos.x();
  //   G4double y = vertex_pos.y();
  //   G4double z = vertex_pos.z();
  //   std::cout<<"Vertex: " << x << " " << y << " " << z<<"\n";
  //
  //   G4double xp = theTrack->GetPosition().x();
  //   G4double yp = theTrack->GetPosition().y();
  //   G4double zp = theTrack->GetPosition().z();
  //   std::cout<<"Pos: "<< xp << " " << yp << " " << zp <<"\n";
  //
  //   std::cout << "Press Enter to Continue";
  //   std::cin.ignore();
  // };

  if(theTrack->GetDynamicParticle()->GetMomentumDirection().z()<=0.0) return;
  // if(PDG==22) return;
  // if(theTrack->GetParentID() != 0) return;

  fTestPlaneVolume = G4LogicalVolumeStore::GetInstance()->GetVolume("TestPlane");
  G4ParticleDefinition* particleType = theTrack->GetDefinition();

  // check if it is entering the test volume
  G4StepPoint* thePrePoint = theStep->GetPreStepPoint();
  G4VPhysicalVolume* thePrePV = thePrePoint->GetPhysicalVolume();
  G4LogicalVolume* thePreLV = thePrePV->GetLogicalVolume();

  G4LogicalVolume* thePostLV = 0;
  G4StepPoint* thePostPoint = theStep->GetPostStepPoint();
  if (thePostPoint) {
     G4VPhysicalVolume* thePostPV = thePostPoint->GetPhysicalVolume();
     if (thePostPV) thePostLV = thePostPV->GetLogicalVolume();
  }

  if (thePostLV == fTestPlaneVolume and
      thePreLV  != fTestPlaneVolume) {

    G4double x = theTrack->GetPosition().x();
    G4double y = theTrack->GetPosition().y();
    G4double z = theTrack->GetPosition().z();

    G4double E = theTrack->GetKineticEnergy();

    G4ThreeVector theMomentumDirection = theTrack->GetDynamicParticle()->GetMomentumDirection();

    G4double Mx = theMomentumDirection.x();
    G4double My = theMomentumDirection.y();
    G4double Mz = theMomentumDirection.z();

    // if(PDG!=11 and PDG!=-11 and Mz>=0.0){
    // std::cout << "\n PDG: " <<PDG <<" "<< theTrack->GetParentID() <<"\n";
    // std::cout<<"Moment: "<< Mx << " " << My << " " << Mz<<"\n";
    // std::cout<<"SteppingAction: "<< x << " " << y << " " << z<<"\n";
    // std::cout << "Press Enter to Continue";
    // std::cin.ignore();
    // }

    RootOut->TestPlaneFill(x,y,z,Mx,My,Mz,E,PDG);
  }


  auto event = G4EventManager::GetEventManager()->GetConstCurrentEvent();
  // auto evtNb = event->GetEventID();
  auto eventNb = event->GetEventID();

  auto SensorTestPlaneVolume = G4LogicalVolumeStore::GetInstance()->GetVolume("logicBaseUnit_T");

  if (thePostLV == SensorTestPlaneVolume and
      thePreLV  != SensorTestPlaneVolume) {

        // std::cout<<"Taken logicBaseUnit_T Track with eventID: "<< eventNb << "\n";
        // std::cout << "Press Enter to Continue";
        // std::cin.ignore();

        G4double x = theTrack->GetPosition().x();
        G4double y = theTrack->GetPosition().y();
        G4double z = theTrack->GetPosition().z();

        G4double E = theTrack->GetKineticEnergy();

        // G4ThreeVector theMomentumDirection = theTrack->GetDynamicParticle()->GetMomentumDirection();
        //
        // G4double Mx = theMomentumDirection.x();
        // G4double My = theMomentumDirection.y();
        // G4double Mz = theMomentumDirection.z();

        G4ThreeVector theMomentumDirection = theTrack->GetDynamicParticle()->GetMomentumDirection();
        G4double Mx = theMomentumDirection.x();
        G4double My = theMomentumDirection.y();
        G4double Mz = theMomentumDirection.z();


        // std::cout << "\n PDG: " <<PDG <<" "<< theTrack->GetParentID() <<"\n";
        // std::cout<<"Moment: "<< Mx << " " << My << " " << Mz<<"\n";
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
        // analysisManager->FillNtupleIColumn(3,8, cell_num);
        // analysisManager->FillNtupleIColumn(3,9, sector_num);
        // analysisManager->FillNtupleIColumn(3,10,layer_num);
        // analysisManager->FillNtupleIColumn(3,11,LC_num);
        analysisManager->FillNtupleDColumn(3,12,weight_fromMC);
        analysisManager->FillNtupleIColumn(3,13,eventNb);
        analysisManager->AddNtupleRow(3);

  }

}
