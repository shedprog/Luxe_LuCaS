//
/// \brief Implementation of the DetectorConstruction class
//

#include "DetectorConstruction.hh"

#include "G4Material.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4UniformMagField.hh"
#include "G4FieldManager.hh"

#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"

#include "G4UnitsTable.hh"
#include "G4NistManager.hh"
#include "G4RunManager.hh"

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

#include "G4GlobalMagFieldMessenger.hh"
#include "G4AutoDelete.hh"

#include "G4UserLimits.hh"

DetectorConstruction::DetectorConstruction()

:G4VUserDetectorConstruction(),
 fWorldMaterial(0),fDefaultWorld(true),
 fSolidWorld(0),fLogicWorld(0),fPhysiWorld(0),
 fMagnetFieldValue(0.0), fSolidMagnet(0), fLogicMagnet(0), fPhysiMagnet(0),
 fMagnetSizeX(0), fMagnetSizeY(0), fMagnetSizeZ(0), fMagnetZPos(0),

 // LumiCal part
  logicWorld(0),
  physiWorld(0),
  logicWholeLC(0),
  logicSensor(0),
  logicCell(0)

{

  fWorldSizeZ = 11.0*m;
  fWorldSizeXY= 3.0*m;

  DefineMaterials();
  SetWorldMaterial("Galactic");


  // create commands for interactive definition of the calorimeter
  // LumiCal = new LCDetectorConstruction;
}

DetectorConstruction::~DetectorConstruction(){ }


G4VPhysicalVolume* DetectorConstruction::Construct()
{
  return ConstructCalorimeter();
}



void DetectorConstruction::DefineMaterials()
{
  //This function illustrates the possible ways to define materials

  G4String symbol;             //a=mass of a mole;
  G4double a, z, density;      //z=mean number of protons;

  G4int ncomponents, natoms;
  G4double fractionmass;
  G4double temperature, pressure;

  // //
  // // example of vacuum
  // //

  density     = universe_mean_density;    //from PhysicalConstants.h
  pressure    = 3.e-18*pascal;
  temperature = 2.73*kelvin;
  new G4Material("Galactic", z=1, a=1.01*g/mole,density,
                 kStateGas,temperature,pressure);
}


G4VPhysicalVolume* DetectorConstruction::ConstructCalorimeter()
{
  // Cleanup old geometry
  //
  G4GeometryManager::GetInstance()->OpenGeometry();
  G4PhysicalVolumeStore::GetInstance()->Clean();
  G4LogicalVolumeStore::GetInstance()->Clean();
  G4SolidStore::GetInstance()->Clean();

  // World
  //
  fSolidWorld = new G4Box("World",                                //its name
                   fWorldSizeXY/2,fWorldSizeXY/2,fWorldSizeZ/2);   //its size

  fLogicWorld = new G4LogicalVolume(fSolidWorld,          //its solid
                                    fWorldMaterial,        //its material
                                    "World");              //its name

  fPhysiWorld = new G4PVPlacement(0,                      //no rotation
                                  G4ThreeVector(),       //at (0,0,0)
                                  fLogicWorld,             //its logical volume
                                  "World",                 //its name
                                  0,                       //its mother  volume
                                  false,                   //no boolean operation
                                  0);                        //copy number


  // std::cout<<"@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ DetectorConstruction::ConstructCalorimeter\n";
  InitDetectorParameters();

  BuildTBeamPT16();

  ConstructMagnet();


  // TestPlane
  G4Box* solidTestPlane = new G4Box("TestPlane",fWorldSizeXY/2,fWorldSizeXY/2,1*mm);   //its size

  G4LogicalVolume* logicTestPlane = new G4LogicalVolume(solidTestPlane,
                                                        Vacuum,
                                                        "TestPlane");

  G4ThreeVector testPlaneCenter = G4ThreeVector(0.,0.,4.3*m);
  new G4PVPlacement(0,
                    testPlaneCenter,
                    "TestPlane",
                    logicTestPlane,
                    fPhysiWorld,
                    false,
                    0);


  // std::cout<<"@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ DetectorConstruction::End\n";
  return fPhysiWorld;
}


void DetectorConstruction::SetWorldMaterial(G4String materialChoice)
{
  // search the material by its name
  G4Material* pttoMaterial =
    G4NistManager::Instance()->FindOrBuildMaterial(materialChoice);

  if (pttoMaterial && fWorldMaterial != pttoMaterial) {
    fWorldMaterial = pttoMaterial;
    if(fLogicWorld) fLogicWorld->SetMaterial(fWorldMaterial);
    if(fLogicMagnet) fLogicMagnet->SetMaterial(fWorldMaterial);
    G4RunManager::GetRunManager()->PhysicsHasBeenModified();
  }
}



void DetectorConstruction::SetWorldSizeZ(G4double val)
{
  fWorldSizeZ = val;
  fDefaultWorld = false;
  // ComputeCalorParameters();
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}



void DetectorConstruction::SetWorldSizeXY(G4double val)
{
  fWorldSizeXY = val;
  fDefaultWorld = false;
  // ComputeCalorParameters();
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}


void DetectorConstruction::ConstructMagnet()
{
  //Pay Attention:
  //Here only active aria of the magnet simulated
  //Do not put 1.5m magnet - because the active aria only 1.0m

   fMagnetFieldValue = 1.4*tesla;
//  fMagnetFieldValue = 0.0*gauss;
  fMagnetSizeX = 60*cm;
  fMagnetSizeY = 20*cm;
  fMagnetSizeZ = 100*cm;
  // fMagnetZPos = fZposAbs + 30*cm + fMagnetSizeZ/2.0 + fAbsorberThickness/2.0;
  fMagnetZPos = 100*cm + 150*cm/2;

  G4ThreeVector  fieldVector( 0.0, fMagnetFieldValue, 0.0);
  G4MagneticField *magField = new G4UniformMagField( fieldVector );
  G4FieldManager  *localFieldMgr = new G4FieldManager (magField);
  localFieldMgr->CreateChordFinder(magField);

  fSolidMagnet = new G4Box("Magnet", fMagnetSizeX/2.0, fMagnetSizeY/2.0, fMagnetSizeZ/2.0);
  fLogicMagnet = new G4LogicalVolume(fSolidMagnet, fWorldMaterial, "Magnet");
  fPhysiMagnet = new G4PVPlacement(0, G4ThreeVector(0.,0., fMagnetZPos), fLogicMagnet, "Magnet", fLogicWorld, false, 0);

  fLogicMagnet->SetFieldManager(localFieldMgr, true);
}
