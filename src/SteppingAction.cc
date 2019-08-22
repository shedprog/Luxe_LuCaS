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

  if(PDG==22) return;
  if(theTrack->GetParentID() != 0) return;

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
    // std::cout << "PDG" <<PDG <<" "<< theTrack->GetParentID() <<"\n";
    // std::cout<<"SteppingAction: "<< x << " " << y << " " << z;

    RootOut->TestPlaneFill(x,y,z,Mx,My,Mz,E,PDG);

    return;

  }

}
