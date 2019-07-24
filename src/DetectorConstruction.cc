//
/// \brief Implementation of the DetectorConstruction class
//

#include "DetectorConstruction.hh"
#include "DetectorMessenger.hh"

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

// #include "LCDetectorConstruction.hh"

// DetectorConstruction::DetectorConstruction

DetectorConstruction::DetectorConstruction()

:G4VUserDetectorConstruction(),
 fAbsorberMaterial(0),fWorldMaterial(0),fDefaultWorld(true),
 fSolidWorld(0),fLogicWorld(0),fPhysiWorld(0),
 fSolidAbsorber(0),fLogicAbsorber(0),fPhysiAbsorber(0),
 fDetectorMessenger(0),
 fMagnetFieldValue(0.0), fSolidMagnet(0), fLogicMagnet(0), fPhysiMagnet(0),
 fMagnetSizeX(0), fMagnetSizeY(0), fMagnetSizeZ(0), fMagnetZPos(0), fTargetType(tfoil),

 // LumiCal part
  logicWorld(0),      physiWorld(0),           // World
  logicWholeLC(0),                             // WholeLC
  logicSensor(0),
  logicCell(0)                                 // Cell

{
  // default parameter values of the calorimeter
  // fAbsorberThickness = 5.0*mm;
  // fAbsorberSizeXY    = 10.0*cm;
  // fZposAbs           = 0.*cm;
  ComputeCalorParameters();
  
  // materials  
  DefineMaterials();
  SetWorldMaterial   ("Galactic");
//   SetWorldMaterial   ("Air20");
  SetAbsorberMaterial("Copper");
//   SetAbsorberMaterial("Air20");
 
  // create commands for interactive definition of the calorimeter  
  fDetectorMessenger = new DetectorMessenger(this);

  // LumiCal = new LCDetectorConstruction;
}



DetectorConstruction::~DetectorConstruction()
{ 
  delete fDetectorMessenger;
  // delete LumiCal;
}



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
  
  //
  // define Elements
  //

  G4Element* H  = new G4Element("Hydrogen",symbol="H",  z= 1, a=   1.01*g/mole);
  G4Element* C  = new G4Element("Carbon",  symbol="C",  z= 6, a=  12.01*g/mole);
  G4Element* N  = new G4Element("Nitrogen",symbol="N",  z= 7, a=  14.01*g/mole);
  G4Element* O  = new G4Element("Oxygen",  symbol="O",  z= 8, a=  16.00*g/mole);
  G4Element* Na = new G4Element("Sodium",  symbol="Na", z=11, a=  22.99*g/mole);
  G4Element* Ar = new G4Element("Argon",   symbol="Ar", z=18, a=  39.95*g/mole);
  G4Element* I  = new G4Element("Iodine",  symbol="I" , z=53, a= 126.90*g/mole);
  G4Element* Xe = new G4Element("Xenon",   symbol="Xe", z=54, a= 131.29*g/mole);

  //
  // define simple materials
  //

  new G4Material("H2Liq"    , z= 1, a= 1.01*g/mole, density= 70.8*mg/cm3);
  new G4Material("Beryllium", z= 4, a= 9.01*g/mole, density= 1.848*g/cm3);
  new G4Material("Aluminium", z=13, a=26.98*g/mole, density= 2.700*g/cm3);
  new G4Material("Silicon"  , z=14, a=28.09*g/mole, density= 2.330*g/cm3);

  G4Material* lAr = 
    new G4Material("liquidArgon", density= 1.390*g/cm3, ncomponents=1);
  lAr->AddElement(Ar, natoms=1);

  new G4Material("Iron",     z=26, a= 55.85*g/mole, density= 7.870*g/cm3);
  new G4Material("Copper",   z=29, a= 63.55*g/mole, density= 8.960*g/cm3);
  new G4Material("Germanium",z=32, a= 72.61*g/mole, density= 5.323*g/cm3);
  new G4Material("Silver",   z=47, a=107.87*g/mole, density= 10.50*g/cm3);
  new G4Material("Tungsten", z=74, a=183.85*g/mole, density= 19.30*g/cm3);
  new G4Material("Gold",     z=79, a=196.97*g/mole, density= 19.32*g/cm3);
  new G4Material("Lead",     z=82, a=207.19*g/mole, density= 11.35*g/cm3);
  new G4Material("Nickel",   z=28, a=58.693*g/mole, density= 8.9*g/cm3);
 

  //
  // define a material from elements.   case 1: chemical molecule
  //

  G4Material* H2O = new G4Material("Water",density= 1.000*g/cm3,ncomponents=2);
  H2O->AddElement(H, natoms=2);
  H2O->AddElement(O, natoms=1);
  H2O->GetIonisation()->SetMeanExcitationEnergy(78*eV);

  G4Material* CH = new G4Material("Plastic",density= 1.04*g/cm3,ncomponents=2);
  CH->AddElement(C, natoms=1);
  CH->AddElement(H, natoms=1);

  G4Material* NaI = new G4Material("NaI", density= 3.67*g/cm3, ncomponents=2);
  NaI->AddElement(Na, natoms=1);
  NaI->AddElement(I , natoms=1);
  NaI->GetIonisation()->SetMeanExcitationEnergy(452*eV);

  //
  // define a material from elements.   case 2: mixture by fractional mass
  //

  G4Material* Air = new G4Material("Air", density= 1.290*mg/cm3, ncomponents=2);
  Air->AddElement(N, fractionmass=0.7);
  Air->AddElement(O, fractionmass=0.3);

  G4Material* Air20 = 
    new G4Material("Air20", density= 1.205*mg/cm3, ncomponents=2,
                   kStateGas, 293.*kelvin, 1.*atmosphere);
  Air20->AddElement(N, fractionmass=0.7);
  Air20->AddElement(O, fractionmass=0.3);

  //Graphite
  //
  G4Material* Graphite = 
    new G4Material("Graphite", density= 1.7*g/cm3, ncomponents=1);
  Graphite->AddElement(C, fractionmass=1.);

  //Havar
  //
  G4Element* Cr = new G4Element("Chrome", "Cr", z=25, a=  51.996*g/mole);
  G4Element* Fe = new G4Element("Iron"  , "Fe", z=26, a=  55.845*g/mole);
  G4Element* Co = new G4Element("Cobalt", "Co", z=27, a=  58.933*g/mole);
  G4Element* Ni = new G4Element("Nickel", "Ni", z=28, a=  58.693*g/mole);
  G4Element* W  = new G4Element("Tungsten","W", z=74, a= 183.850*g/mole);

  G4Material* Havar = 
    new G4Material("Havar", density= 8.3*g/cm3, ncomponents=5);
  Havar->AddElement(Cr, fractionmass=0.1785);
  Havar->AddElement(Fe, fractionmass=0.1822);
  Havar->AddElement(Co, fractionmass=0.4452);
  Havar->AddElement(Ni, fractionmass=0.1310);
  Havar->AddElement(W , fractionmass=0.0631);

  //
  // examples of gas
  //  
  new G4Material("ArgonGas", z=18, a=39.948*g/mole, density= 1.782*mg/cm3,
                 kStateGas, 273.15*kelvin, 1*atmosphere);
                           
  new G4Material("XenonGas", z=54, a=131.29*g/mole, density= 5.458*mg/cm3,
                 kStateGas, 293.15*kelvin, 1*atmosphere);
                           
  G4Material* CO2 =
    new G4Material("CarbonicGas", density= 1.977*mg/cm3, ncomponents=2);
  CO2->AddElement(C, natoms=1);
  CO2->AddElement(O, natoms=2);

  G4Material* ArCO2 =
    new G4Material("ArgonCO2",   density= 1.8223*mg/cm3, ncomponents=2);
  ArCO2->AddElement (Ar,  fractionmass=0.7844);
  ArCO2->AddMaterial(CO2, fractionmass=0.2156);

  //another way to define mixture of gas per volume
  G4Material* NewArCO2 =
    new G4Material("NewArgonCO2", density= 1.8223*mg/cm3, ncomponents=3);
  NewArCO2->AddElement (Ar, natoms=8);
  NewArCO2->AddElement (C,  natoms=2);
  NewArCO2->AddElement (O,  natoms=4);

  G4Material* ArCH4 = 
    new G4Material("ArgonCH4",    density= 1.709*mg/cm3,  ncomponents=3);
  ArCH4->AddElement (Ar, natoms=93);
  ArCH4->AddElement (C,  natoms=7);
  ArCH4->AddElement (H,  natoms=28);

  G4Material* XeCH = 
    new G4Material("XenonMethanePropane", density= 4.9196*mg/cm3, ncomponents=3,
                   kStateGas, 293.15*kelvin, 1*atmosphere);
  XeCH->AddElement (Xe, natoms=875);
  XeCH->AddElement (C,  natoms=225);
  XeCH->AddElement (H,  natoms=700);

  G4Material* steam = 
    new G4Material("WaterSteam", density= 1.0*mg/cm3, ncomponents=1);
  steam->AddMaterial(H2O, fractionmass=1.);
  steam->GetIonisation()->SetMeanExcitationEnergy(71.6*eV);  

  //
  // example of vacuum
  //

  density     = universe_mean_density;    //from PhysicalConstants.h
  pressure    = 3.e-18*pascal;
  temperature = 2.73*kelvin;
  new G4Material("Galactic", z=1, a=1.01*g/mole,density,
                 kStateGas,temperature,pressure);
}



void DetectorConstruction::ComputeCalorParameters()
{
  // Compute derived parameters of the calorimeter
  // fZstartAbs = fZposAbs-0.5*fAbsorberThickness; 
  // fZendAbs   = fZposAbs+0.5*fAbsorberThickness;
  // fZstartAbs = 300*cm;
  // fZendAbs = 350*cm;

  if (fDefaultWorld) {
//      fWorldSizeZ = 1.5*fAbsorberThickness; fWorldSizeXY= 1.2*fAbsorberSizeXY;
     fWorldSizeZ = 10.0*m; 
//      fWorldSizeZ = 400.0*cm; 
     fWorldSizeXY= 2.0*m;
  }         
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
                                 
  // Absorber
  // // 
  // if (fTargetType == twire) {
  //   fSolidAbsorber = new G4Tubs("Absorber",        
  //                    0, fAbsorberThickness/2.0 , fWorldSizeXY/2.0, 0, 2.0*M_PI );
                          
  //   fLogicAbsorber = new G4LogicalVolume(fSolidAbsorber,    //its solid
  //                                           fAbsorberMaterial, //its material
  //                                         "Absorber");       //its name
                                                
  //   fPhysiAbsorber = new G4PVPlacement(new G4RotationMatrix(0.0, M_PI/2.0, 0.0),         //rotation
  //                       G4ThreeVector(0.,0., fZposAbs),    //its position
  //                               fLogicAbsorber,     //its logical volume
  //                               "Absorber",         //its name
  //                               fLogicWorld,        //its mother
  //                               false,             //no boulean operat
  //                               0);                //copy number
  // } else { 
  //   fSolidAbsorber = new G4Box("Absorber",        
  //                     fAbsorberSizeXY/2.,fAbsorberSizeXY/2.,fAbsorberThickness/2.);
                          
  //   fLogicAbsorber = new G4LogicalVolume(fSolidAbsorber,    //its solid
  //                                           fAbsorberMaterial, //its material
  //                                         "Absorber");       //its name
                                                
  //   fPhysiAbsorber = new G4PVPlacement(0,                   //no rotation
  //                       G4ThreeVector(0.,0., fZposAbs),    //its position
  //                               fLogicAbsorber,     //its logical volume
  //                               "Absorber",         //its name
  //                               fLogicWorld,        //its mother
  //                               false,             //no boulean operat
  //                               0);                //copy number
  // }              


  InitDetectorParameters();                   
  // ConstructLumiCal();
  BuildTBeamPT16();

  ConstructMagnet();

//   ConstructTelescope();
//   G4double maxStep = 10.0*nm;
//   G4double maxStep = 1.0*um;
//   fStepLimit = new G4UserLimits(maxStep);
//   fStepLimit->SetMaxAllowedStep(maxStep);
//   fLogicAbsorber->SetUserLimits(fStepLimit);
  
  PrintCalorParameters();         
  
  return fPhysiWorld;
}

// G4VPhysicalVolume* DetectorConstruction::ConstructCalorimeter()
// {
//   G4GeometryManager::GetInstance()->OpenGeometry();
//   G4PhysicalVolumeStore::GetInstance()->Clean();
//   G4LogicalVolumeStore::GetInstance()->Clean();
//   G4SolidStore::GetInstance()->Clean();

//   fPhysiWorld = LumiCal->Construct();

//   ConstructMagnet();

//   return fPhysiWorld;
// }



void DetectorConstruction::PrintCalorParameters()
{
  G4cout << "\n" << fWorldMaterial    << G4endl;
  G4cout << "\n" << fAbsorberMaterial << G4endl;
    
  G4cout << "\n The  WORLD   is made of "  << G4BestUnit(fWorldSizeZ,"Length")
         << " of " << fWorldMaterial->GetName();
  G4cout << ". The transverse size (XY) of the world is " 
         << G4BestUnit(fWorldSizeXY,"Length") << G4endl;
  G4cout << " The ABSORBER is made of " 
         <<G4BestUnit(fAbsorberThickness,"Length")
         << " of " << fAbsorberMaterial->GetName() << (fTargetType == tfoil ? " foil." : " wire.");
  G4cout << ". The transverse size (XY) is " 
         << G4BestUnit(fAbsorberSizeXY,"Length") << G4endl;
  G4cout << " Z position of the middle of the absorber "
         << G4BestUnit(fZposAbs,"Length");
  G4cout << G4endl;
}



void DetectorConstruction::SetAbsorberMaterial(G4String materialChoice)
{
  // search the material by its name
  G4Material* pttoMaterial = G4NistManager::Instance()->FindOrBuildMaterial(materialChoice);

  if (pttoMaterial && fAbsorberMaterial != pttoMaterial) {
    fAbsorberMaterial = pttoMaterial;                  
    if(fLogicAbsorber) fLogicAbsorber->SetMaterial(fAbsorberMaterial);
    G4RunManager::GetRunManager()->PhysicsHasBeenModified();
  }
}


void DetectorConstruction::SetAbsorberType(G4String val)
{
  if (val == "foil") {
    fTargetType = tfoil;
  } else if (val == "wire") {
    fTargetType = twire;
  } else {
    G4cout << "DetectorConstruction::SetAbsorberType: <" << val << ">"
           << " is not supported!"
           << G4endl;
  }
  G4RunManager::GetRunManager()->ReinitializeGeometry();
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
    


void DetectorConstruction::SetAbsorberThickness(G4double val)
{
  fAbsorberThickness = val;
  ComputeCalorParameters();
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}



void DetectorConstruction::SetAbsorberSizeXY(G4double val)
{
  fAbsorberSizeXY = val;
  ComputeCalorParameters();
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}



void DetectorConstruction::SetWorldSizeZ(G4double val)
{
  fWorldSizeZ = val;
  fDefaultWorld = false;
  ComputeCalorParameters();
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}



void DetectorConstruction::SetWorldSizeXY(G4double val)
{
  fWorldSizeXY = val;
  fDefaultWorld = false;
  ComputeCalorParameters();
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}



void DetectorConstruction::SetAbsorberZpos(G4double val)
{
  fZposAbs  = val;
  ComputeCalorParameters();
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}  

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....




void DetectorConstruction::ConstructSDandField()
{
    if ( fFieldMessenger.Get() == 0 ) {
        // Create global magnetic field messenger.
        // Uniform magnetic field is then created automatically if
        // the field value is not zero.
        G4ThreeVector fieldValue = G4ThreeVector();
        G4GlobalMagFieldMessenger* msg = new G4GlobalMagFieldMessenger(fieldValue);
        //msg->SetVerboseLevel(1);
        G4AutoDelete::Register(msg);
        fFieldMessenger.Put( msg );
        
    }
}



void DetectorConstruction::ConstructMagnet()
{
   fMagnetFieldValue = 1.4*tesla;
//  fMagnetFieldValue = 0.0*gauss;
  fMagnetSizeX = 50*cm;
  fMagnetSizeY = 20*cm;
  fMagnetSizeZ = 100*cm;
  // fMagnetZPos = fZposAbs + 30*cm + fMagnetSizeZ/2.0 + fAbsorberThickness/2.0;
  fMagnetZPos = 150*cm;
  
  G4ThreeVector  fieldVector( 0.0, fMagnetFieldValue, 0.0);  
  G4MagneticField *magField = new G4UniformMagField( fieldVector );    
  G4FieldManager  *localFieldMgr = new G4FieldManager (magField); 
  localFieldMgr->CreateChordFinder(magField);
  
  fSolidMagnet = new G4Box("Magnet", fMagnetSizeX/2.0, fMagnetSizeY/2.0, fMagnetSizeZ/2.0); 
  fLogicMagnet = new G4LogicalVolume(fSolidMagnet, fWorldMaterial, "Magnet"); 
  fPhysiMagnet = new G4PVPlacement(0, G4ThreeVector(0.,0., fMagnetZPos), fLogicMagnet, "Magnet", fLogicWorld, false, 0); 
  
  fLogicMagnet->SetFieldManager(localFieldMgr, true);   
}   



// void DetectorConstruction::ConstructTelescope()
// {
//   fTelescopeNPlanes = 6;  
//   fTelscopeSensorSizeX = 2.0*cm;
//   fTelscopeSensorSizeY = 2.0*cm;
//   fTelscopeSensorSizeZ = 0.05*mm;
//   fTelscopePlanesDz = 15.0*cm;
//   fTelscopeFirstPlaneZPos = fMagnetZPos + 50*cm + fMagnetSizeZ/2.0 + fTelscopeSensorSizeZ/2.0;
  
//   G4Material* TelescopeSensorMaterial = G4NistManager::Instance()->FindOrBuildMaterial("Silicon");
//   fSolidTelescopeSensor = new G4Box("TelescopeSensor", fTelscopeSensorSizeX/2.0, fTelscopeSensorSizeY/2.0, fTelscopeSensorSizeZ/2.0); 
//   fLogicTelescopeSensor = new G4LogicalVolume(fSolidTelescopeSensor, TelescopeSensorMaterial, "TelescopeSensor"); 

  
//   for(int i = 0; i < fTelescopeNPlanes; ++i) {
//     G4double plane_z_pos = fTelscopeFirstPlaneZPos + i*fTelscopePlanesDz;
//     if (i > 2) plane_z_pos += 10.0*cm;
//     new G4PVPlacement(0, G4ThreeVector(5.0*mm, 0., plane_z_pos), fLogicTelescopeSensor, "TelescopeSensor", fLogicWorld, false, i, true); 
//   }
// }



