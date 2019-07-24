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
/// \file electromagnetic/TestEm5/src/DetectorMessenger.cc
/// \brief Implementation of the DetectorMessenger class
//
// $Id: DetectorMessenger.cc 77083 2013-11-21 10:35:55Z gcosmo $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "DetectorMessenger.hh"

#include "DetectorConstruction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithoutParameter.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorMessenger::DetectorMessenger(DetectorConstruction * Det)
:G4UImessenger(),fDetector(Det),
 fTestemDir(0),
 fDetDir(0),
 fAbsMaterCmd(0),
 fAbsThickCmd(0),
 fAbsSizXYCmd(0),
 fAbsZposCmd(0),
 fWorldMaterCmd(0),
 fWorldZCmd(0),
 fWorldXYCmd(0),
 fAbsTypeCmd(0)
{ 
  fTestemDir = new G4UIdirectory("/lxphoton/");
  fTestemDir->SetGuidance("UI commands specific to this example.");
  
  fDetDir = new G4UIdirectory("/lxphoton/det/");
  fDetDir->SetGuidance("detector construction commands");
      
  fAbsMaterCmd = new G4UIcmdWithAString("/lxphoton/det/setAbsMat",this);
  fAbsMaterCmd->SetGuidance("Select Material of the Absorber.");
  fAbsMaterCmd->SetParameterName("choice",false);
  fAbsMaterCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  fAbsMaterCmd->SetToBeBroadcasted(false);
  
  fAbsTypeCmd = new G4UIcmdWithAString("/lxphoton/det/setAbsType",this);
  fAbsTypeCmd->SetGuidance("Set tareget shape: (foil/wire).");
  fAbsTypeCmd->SetParameterName("TargetType",false);
  fAbsTypeCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  fAbsTypeCmd->SetToBeBroadcasted(false);

  fWorldMaterCmd = new G4UIcmdWithAString("/lxphoton/det/setWorldMat",this);
  fWorldMaterCmd->SetGuidance("Select Material of the World.");
  fWorldMaterCmd->SetParameterName("wchoice",false);
  fWorldMaterCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  fWorldMaterCmd->SetToBeBroadcasted(false);
    
  fAbsThickCmd = new G4UIcmdWithADoubleAndUnit("/lxphoton/det/setAbsThick",this);
  fAbsThickCmd->SetGuidance("Set Thickness of the Absorber");
  fAbsThickCmd->SetParameterName("SizeZ",false);  
  fAbsThickCmd->SetRange("SizeZ>0.");
  fAbsThickCmd->SetUnitCategory("Length");  
  fAbsThickCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  fAbsThickCmd->SetToBeBroadcasted(false);
  
  fAbsSizXYCmd = new G4UIcmdWithADoubleAndUnit("/lxphoton/det/setAbsXY",this);
  fAbsSizXYCmd->SetGuidance("Set sizeXY of the Absorber");
  fAbsSizXYCmd->SetParameterName("SizeXY",false);
  fAbsSizXYCmd->SetRange("SizeXY>0.");
  fAbsSizXYCmd->SetUnitCategory("Length");
  fAbsSizXYCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  fAbsSizXYCmd->SetToBeBroadcasted(false);
  
  fAbsZposCmd = new G4UIcmdWithADoubleAndUnit("/lxphoton/det/setAbsZpos",this);
  fAbsZposCmd->SetGuidance("Set X pos. of the Absorber");
  fAbsZposCmd->SetParameterName("Xpos",false);
  fAbsZposCmd->SetUnitCategory("Length");
  fAbsZposCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  fAbsZposCmd->SetToBeBroadcasted(false);
  
  fWorldZCmd = new G4UIcmdWithADoubleAndUnit("/lxphoton/det/setWorldZ",this);
  fWorldZCmd->SetGuidance("Set X size of the World");
  fWorldZCmd->SetParameterName("WSizeX",false);
  fWorldZCmd->SetRange("WSizeX>0.");
  fWorldZCmd->SetUnitCategory("Length");
  fWorldZCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  fWorldZCmd->SetToBeBroadcasted(false);
  
  fWorldXYCmd = new G4UIcmdWithADoubleAndUnit("/lxphoton/det/setWorldXY",this);
  fWorldXYCmd->SetGuidance("Set sizeXY of the World");
  fWorldXYCmd->SetParameterName("WSizeXY",false);
  fWorldXYCmd->SetRange("WSizeXY>0.");
  fWorldXYCmd->SetUnitCategory("Length");
  fWorldXYCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  fWorldXYCmd->SetToBeBroadcasted(false);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorMessenger::~DetectorMessenger()
{
  delete fAbsMaterCmd; 
  delete fAbsThickCmd; 
  delete fAbsSizXYCmd;  
  delete fAbsZposCmd; 
  delete fWorldMaterCmd;
  delete fWorldZCmd;
  delete fWorldXYCmd;
  delete fDetDir;  
  delete fTestemDir;
  delete fAbsTypeCmd;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{ 
  if ( command == fAbsMaterCmd )
   {fDetector->SetAbsorberMaterial(newValue);}
   
  if ( command == fAbsTypeCmd )
   {fDetector->SetAbsorberType(newValue);}
   
  if ( command == fWorldMaterCmd )
   {fDetector->SetWorldMaterial(newValue);}
   
  if ( command == fAbsThickCmd )
  {fDetector->SetAbsorberThickness(fAbsThickCmd->GetNewDoubleValue(newValue));}
   
  if ( command == fAbsSizXYCmd )
   {fDetector->SetAbsorberSizeXY(fAbsSizXYCmd->GetNewDoubleValue(newValue));}
   
  if ( command == fAbsZposCmd )
   {fDetector->SetAbsorberZpos(fAbsZposCmd->GetNewDoubleValue(newValue));}
   
  if ( command == fWorldZCmd )
   {fDetector->SetWorldSizeZ(fWorldZCmd->GetNewDoubleValue(newValue));}
   
  if ( command == fWorldXYCmd )
   {fDetector->SetWorldSizeXY(fWorldXYCmd->GetNewDoubleValue(newValue));}
   
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
