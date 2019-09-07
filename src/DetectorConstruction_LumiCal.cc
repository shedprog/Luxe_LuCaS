

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


#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Polycone.hh"
#include "G4Polyhedra.hh"
#include "G4EllipticalTube.hh"
#include "G4LogicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"
#include "G4PVParameterised.hh"
#include "G4SubtractionSolid.hh"
#include "G4UnionSolid.hh"

#include "G4Material.hh"
#include "G4MaterialTable.hh"
#include "G4Element.hh"
#include "G4ElementTable.hh"

#include "G4TransportationManager.hh"
#include "G4FieldManager.hh"
#include "G4SDManager.hh"
#include "G4Region.hh"
#include "G4ProductionCuts.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "G4ios.hh"

#include "G4GeometryTolerance.hh"
#include "G4GeometryManager.hh"

#include "G4UserLimits.hh"

#include "G4SystemOfUnits.hh"

//Add mag
#include "G4UniformMagField.hh"
#include "G4FieldManager.hh"
#include "G4GlobalMagFieldMessenger.hh"

#include <G4PVDivision.hh>

// Taken from LUCAS part  
#include "Setup.hh"
#include <iostream>
#include <sstream>
#include <string>
#include "globals.hh"

#include "LCSensitiveDetector.hh"
#include "DetectorConstruction.hh"
// #include "DetectorMessenger.hh"

//----------------------------------------
//// DESY 2016 Prototype test beam
//----------------------------------------

G4double SensorGap(G4double R, G4double r, G4double alpha)
{
    G4double b = -R * sqrt(4.0-2.0*r*r/(R*R)*(1.0-cos(alpha)));
    G4double a = 1.0;
    G4double c = R*R - r*r;
    G4double dh = (-b-sqrt(b*b-4.0*a*c))/(2.0*a);
    G4double res = r-R+dh;
    return res;
}

void DetectorConstruction::BuildTBeamPT16(){

// Detector parametrs:

bool overlap_check = true;

double base_airx = 63.*mm;
double base_airy = 120.*mm;
double base_airz = 100.*mm;



//---------------
// WORLD
//---------------
// Top level
G4cout<< "DetectorConstrucion::Construct(): creating World ....";
G4cout<< "The paramters : dx dy dz WorldMat = " <<Setup::world_hdx <<"- "  << Setup::world_hdy<<" -  "  << Setup::world_hdz << "- .... " << WorldMat ;

// G4Box *solidWorld = new G4Box("World", 3.*m, 3.*m, 3.*m);
G4Box *solidWorld = new G4Box("World", base_airx, base_airy, base_airz);
logicWorld = new G4LogicalVolume(solidWorld, WorldMat, "World", 0, 0, 0);
physiWorld = new G4PVPlacement(0,               // no rotation
                           G4ThreeVector(), // origin
                           logicWorld,      // its logical volume
                           "World",         // name
                           0,               // no mother log; null ptr
                           true,           // no boolean ops
                           overlap_check);              // copy #

logicWorld->SetVisAttributes(G4VisAttributes::Invisible);

G4cout<< ".......... done! " << G4endl;



//==============================================================================================================
//base units for 2016 TB  at desy : 
//		1. absorber and sensor - "AS"  / "AS_PL", 1 mm gup, 1 slot, total =4.5mm
//		2. absorber only base unit - "A" / "A_PL"  1 mm gup, 1 slot, total =4.5mm with Wabsorber_MGS (93% W )or Wabsorber_PL (95 % W) 
//		3. For Tracker base unit with only sensor muodul and air gup  - "JSM:X",  total wide= 1 + X mm
//              4. To add the plastic scintilator inftont use SC
//
//
//==============================================================================================================


 G4cout<< " Building Test Beam 2016..." << G4endl;

//---------------------------
//   x, y  direction (beam) sizes 
//---------------------------
  // G4double airhx = 300.0 *mm;
  // G4double airhy = 300.0 *mm;
//---------------------------
//   z direction (beam) sizes 
//---------------------------
//need to add z distences of elementts : 
// 1 mm /2 for SensorMoudule and air gups. 
//4.5 mm /2 for base unite with absorber.  
//
  G4double airhz =  hAbsorberPitchDZ +  hTungstenDZ ;  // should be 2.25 mm half a slut 
  G4double airhz1mm = hAbsorberPitchDZ; // should be 0.5 mm 
 std::cout << "$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ pitchhz =  "<< 2* airhz *mm  <<" mm &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&"<< std::endl;
 std::cout << "$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ pitchhz1mm =  "<< 2* airhz1mm *mm <<" mm &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&"<< std::endl;

 G4double airhz_PSC =  10.0*mm ;
 G4double airhz_TR_PSC =  3.75*mm ;

//---------------------------
//   building air base units (6 0f them ) 
//---------------------------
 std::cout << "-------------------------------------------------" << std::endl;
 std::cout << "-----------    building air base units           : "    << std::endl;
 std::cout << "-------------------------------------------------" << std::endl;

  G4Box *solidBaseUnit_S = new G4Box ( "solidBaseUnit_P", base_airx, base_airy, airhz ); 
  G4LogicalVolume *logicBaseUnit_S = new G4LogicalVolume (solidBaseUnit_S, Air, "logicBaseUnit", 0, 0, 0);//1
  logicBaseUnit_S->SetVisAttributes( G4VisAttributes::Invisible );

  G4Box *solidBaseUnit = new G4Box ( "solidBaseUnit", base_airx, base_airy, airhz );
  G4Box *solidBaseUnit_PL = new G4Box ( "solidBaseUnit_PL", base_airx, base_airy, airhz );  
  G4LogicalVolume *logicBaseUnit = new G4LogicalVolume (solidBaseUnit, Air, "logicBaseUnit", 0, 0, 0);//1
  G4LogicalVolume *logicBaseUnit_PL = new G4LogicalVolume (solidBaseUnit_PL, Air, "logicBaseUnit_PL", 0, 0, 0);//2
  logicBaseUnit->SetVisAttributes( G4VisAttributes::Invisible );
  logicBaseUnit_PL->SetVisAttributes( G4VisAttributes::Invisible );
 
  G4Box *solidBaseUnitAB = new G4Box ( "solidBaseUnitAB", base_airx, base_airy, airhz );
  G4Box *solidBaseUnitAB_PL = new G4Box ( "solidBaseUnitAB_PL", base_airx, base_airy, airhz );  
  G4LogicalVolume *logicBaseUnitONLYAB = new G4LogicalVolume (solidBaseUnitAB, Air, "logicBaseUnitONLYAB", 0, 0, 0);//3
   G4LogicalVolume *logicBaseUnitONLYAB_PL = new G4LogicalVolume (solidBaseUnitAB_PL, Air, "logicBaseUnitONLYAB_PL", 0, 0, 0);//4
  logicBaseUnitONLYAB->SetVisAttributes( G4VisAttributes::Invisible );
  logicBaseUnitONLYAB_PL->SetVisAttributes( G4VisAttributes::Invisible );

  G4Box *solidBaseUnitFor1mm = new G4Box ( "solidBaseUnitFor1mm", base_airx, base_airy, airhz1mm );
  G4Box *solidBaseUnitFor1mm_PL = new G4Box ( "solidBaseUnitFor1mm_NoS", base_airx, base_airy,airhz1mm );  
  G4LogicalVolume *logicBaseUnitFor1mm = new G4LogicalVolume (solidBaseUnitFor1mm, Air, "logicBaseUnitFor1mm", 0, 0, 0);//5
  G4LogicalVolume *logicBaseUnitFor1mm_NoS = new G4LogicalVolume (solidBaseUnitFor1mm_PL, Air, "logicBaseUnitFor1mm_NoS", 0, 0, 0);//6
  logicBaseUnitFor1mm->SetVisAttributes( G4VisAttributes::Invisible );
  logicBaseUnitFor1mm_NoS->SetVisAttributes( G4VisAttributes::Invisible );
  
  G4Colour *SC_Color  = new	G4Colour(1., 0., 1., 1.);

  ///for Plastic SC
  G4Box *solidBase_PSC = new G4Box ( "solidBase_PSC", base_airx, base_airy, airhz_PSC );
  G4LogicalVolume *logicBaseUnit_PSC = new G4LogicalVolume (solidBase_PSC, PLASTIC_SC , "logicBaseUnit_PSC", 0, 0, 0);//5
  logicBaseUnit_PSC->SetVisAttributes( SC_Color );

  //for trigger SC
  G4Box *solidBase_TR_PSC = new G4Box ( "solidBase_TR_PSC", base_airx, base_airy, airhz_TR_PSC );
  G4LogicalVolume *logicBaseUnit_TR_PSC = new G4LogicalVolume (solidBase_TR_PSC, PLASTIC_SC , "logicBaseUnit_TR_PSC", 0, 0, 0);//5
  logicBaseUnit_TR_PSC->SetVisAttributes(SC_Color);
  
  //--------------------------------------------
  //add absorber (single X0) in front sensor 
  //---------------------------------------------
  std::cout << "-------------------------------------------------" << std::endl;
  std::cout << "-----------    building  absorber               : "<< std::endl;
  std::cout << "-------------------------------------------------" << std::endl;

  G4double zposAbs0 = -airhz + hTungstenDZ;
  G4double zposAbs0_MSG = -airhz + hTungstenDZ*1.02;

  G4Box *solidAbs0 = new G4Box("solidAbs0", base_airx, base_airy, hTungstenDZ );
  G4Box *solidAbs0_MGS = new G4Box("solidAbs0_MGS", base_airx, base_airy, hTungstenDZ*1.02 ); // nominal thiknes of MSG plate is ~3.57
  G4Box *solidAbs0_PL = new G4Box("solidAbs0_PL", base_airx, base_airy, hTungstenDZ );

  G4LogicalVolume *logicAbs0 = new G4LogicalVolume(solidAbs0, Tungsten, "logicAbs0", 0, 0, 0);
  G4LogicalVolume *logicAbs0_MGS = new G4LogicalVolume(solidAbs0_MGS,Setup::Wabsorber_MGS , "logicAbs0_MGS", 0, 0, 0); // 2 kind of absorber plate we use in 2014 -2015 TB
  G4LogicalVolume *logicAbs0_PL = new G4LogicalVolume(solidAbs0_PL,Setup::Wabsorber_PL , "logicAbs0_PL", 0, 0, 0); 
 
 //G4double zpos1mmAbs0= -airhz1slot + hTungstenDZ;


 std::cout << "-------------------------------------------------" << std::endl;
 std::cout << "-----------  Placement  absorber                : "    << std::endl;
 std::cout << "-------------------------------------------------" << std::endl;
  //G4Colour *Abs0Color  = new	G4Colour(1., 0., 1., 1.);
  //logicAbs0->SetVisAttributes( G4VisAttributes::Abs0Color );
 new G4PVPlacement ( 0, G4ThreeVector( 0., 0.*mm, zposAbs0_MSG ), logicAbs0_MGS, "Absorber0_MGS"  ,logicBaseUnit , false, 0,1); // absorber in base unit 1 
 new G4PVPlacement ( 0, G4ThreeVector( 0., 0.*mm, zposAbs0 )    , logicAbs0_PL,  "Absorber0_PL"   ,logicBaseUnit_PL , false, 0,1); // absorber in base unit 2 
 new G4PVPlacement ( 0, G4ThreeVector( 0., 0.*mm, zposAbs0_MSG ), logicAbs0_MGS, "Absorber0_MGS"  ,logicBaseUnitONLYAB , false, 0, 1); // absorber in base unit 3 
 new G4PVPlacement ( 0, G4ThreeVector( 0., 0.*mm, zposAbs0 )    , logicAbs0_PL,  "Absorber0_PL"   ,logicBaseUnitONLYAB_PL , false, 0, 1); // absorber in base unit 4 

  G4VisAttributes *Abs0Att = new G4VisAttributes(G4Colour(0.0,0.0,1.0));
	//Abs0Att->SetForceWireframe ( true );
  logicAbs0->SetVisAttributes(Abs0Att);
  logicAbs0_MGS->SetVisAttributes(Abs0Att);
  logicAbs0_PL->SetVisAttributes(Abs0Att);  

//----------------------------------------------
//	carbon fiber suppurt 
//---------------------------------------------
 
  // G4double CFhx = 70.0 *mm;
  // G4double CFhy = 132.0 *mm;
  G4double CFhz = 0.395 *mm;

 std::cout << "-------------------------------------------------" << std::endl;
 std::cout << "-----------   building carbon fiber suppurt      : "    << std::endl;
 std::cout << "-------------------------------------------------" << std::endl;

  // G4Box *solidCF = new G4Box ( "solidCF", CFhx, CFhy, CFhz );
  G4Box *solidCF = new G4Box ( "solidCF", base_airx, base_airy, CFhz );
  G4LogicalVolume *logicCF = new G4LogicalVolume (solidCF, Setup::C_fiber, "logicCF", 0, 0, 0);
  
  G4VisAttributes *CF0Att = new G4VisAttributes(G4Colour(1.0,1.0,1.0));
  logicCF->SetVisAttributes(CF0Att);

//---------------------------------------------
//		sensor stracture from front to back  :
//			1. 0.150 kapton front.
//			1.2. 0.010 epoxy glue
//			2.1. 0.020 AL on sensor
//			2.2. 0.320 Si sensor
//			2.3. 0.020 AL on sensor
//			3.1. 0.040 epoxy (condactive - not take into acount)
//			3.2. 0.025 cupper (on kapton)
//			3.3. 0.065 kapton back
// 			3.4. 0.020 epoxy glue
//			total 0.670  
//---------------------------------------------
 

  G4double sPhi = 75. *deg;
  G4double dPhi = 30. *deg;
  // G4double sPhi = 60.*deg;
  // G4double dPhi = 30.*deg;

 std::cout << "-------------------------------------------------" << std::endl;
 std::cout << "-----------   building sensor partrs             : "    << std::endl;
 std::cout << "-------------------------------------------------" << std::endl;

  G4Tubs *solidSensV = new G4Tubs("solidSensorV", SensRadMin + deadPhi ,SensRadMax - deadPhi, hSiliconDZ, sPhi, dPhi);
  G4Tubs *solidMetal = new G4Tubs("solidMetal", SensRadMin + deadPhi ,SensRadMax - deadPhi, hMetalDZ, sPhi, dPhi);
  G4Tubs *solidFanoutFrnt  = new G4Tubs("solidFanoutFrnt",SensRadMin, SensRadMax, hFanoutFrontDZ, sPhi,dPhi);
  G4Tubs *solidFanoutBack  = new G4Tubs("solidFanoutBack",SensRadMin, SensRadMax, hFanoutBackDZ, sPhi,dPhi);

  logicSensorV = new G4LogicalVolume(solidSensV, Silicon, "logicSensorV", 0,0,0);
  logicMetalV = new G4LogicalVolume(solidMetal, Aluminium, "logicMetalV", 0,0,0);
  logicFanoutFrnt = new G4LogicalVolume(solidFanoutFrnt, FanoutMatF,"logicFanoutFront", 0, 0, 0); 
  logicFanoutBack = new G4LogicalVolume(solidFanoutBack, FanoutMatB,"logicFanoutFront", 0, 0, 0); 
  
  // logicSensorV->SetVisAttributes( G4VisAttributes::Invisible );
  // logicMetalV->SetVisAttributes( G4VisAttributes::Invisible );
  G4VisAttributes *FanoutFrntAtt = new G4VisAttributes(G4Colour(0.5,0.1,0.0));
  logicFanoutFrnt->SetVisAttributes(FanoutFrntAtt);

  G4VisAttributes *VisFanoutBack = new G4VisAttributes(G4Colour(0.5,0.1,0.3));
  logicFanoutBack->SetVisAttributes(VisFanoutBack);

  G4VisAttributes *VisSensorV = new G4VisAttributes(G4Colour(0.0,1.0,0.0));
  logicSensorV->SetVisAttributes(VisSensorV);

  G4VisAttributes *VisMetalV = new G4VisAttributes(G4Colour(0.0,1.0,0.7));
  logicMetalV->SetVisAttributes(VisMetalV );


//-------------------------------------------------------------
//		placment of sensor part in CF
//-------------------------------------------------------------
//------placment of first sensor
G4double ypos = - SensRadMax;
// G4double ypos = 0.0;
std::cout << "-------------------------------------------------" << std::endl;
std::cout << "-----------  placment  sensor partrs  in CF suppurt : "    << std::endl;
std::cout << "-------------------------------------------------" << std::endl;
std::cout << "  some importent sizes " << std::endl;
std::cout << " ypos            : " << ypos <<  std::endl;
std::cout << " CFhz            : " << CFhz <<  std::endl;
std::cout << " hFanoutFrontDZ  : " << hFanoutFrontDZ <<  std::endl;
std::cout << " hMetalDZ        : " <<  hMetalDZ<<  std::endl;
std::cout << " hSiliconDZ      : " <<  hSiliconDZ<<  std::endl;
std::cout << " hFanoutBackDZ   : " <<  hFanoutBackDZ<<  std::endl;


G4double SensorAtCF = -CFhz +  hFanoutFrontDZ ; 
std::cout << " 1.SensorAtCF : " << SensorAtCF <<  std::endl;
new G4PVPlacement ( 0, G4ThreeVector(0., ypos, ( SensorAtCF)),logicFanoutFrnt , "FunOut0", logicCF, false,1, overlap_check);

SensorAtCF = SensorAtCF + hFanoutFrontDZ +  hMetalDZ ; 
std::cout << " 2.SensorAtCF : " << SensorAtCF <<  std::endl;
new G4PVPlacement ( 0, G4ThreeVector(0., ypos, ( SensorAtCF)),logicMetalV , "PadMetal0", logicCF, false,1, overlap_check);

SensorAtCF = SensorAtCF + hMetalDZ + hSiliconDZ; 
std::cout << " 3.SensorAtCF : " << SensorAtCF <<  std::endl;
new G4PVPlacement ( 0, G4ThreeVector(0.,ypos , ( SensorAtCF)),logicSensorV , "SensorV0", logicCF, false,1, overlap_check);

SensorAtCF = SensorAtCF + hSiliconDZ + hMetalDZ;
std::cout << " 4.SensorAtCF : " << SensorAtCF <<  std::endl;
new G4PVPlacement ( 0, G4ThreeVector(0., ypos, ( SensorAtCF)),logicMetalV , "PadMetal1", logicCF, false,1, overlap_check);

SensorAtCF =  SensorAtCF + hMetalDZ + hFanoutBackDZ ; 
std::cout << " 5.SensorAtCF : " << SensorAtCF <<  std::endl;
new G4PVPlacement ( 0, G4ThreeVector(0.,ypos , ( SensorAtCF)),logicFanoutBack , "FunOut1", logicCF, false,1, overlap_check);

// ------placement of second sensor
double delta_NoOverlap = SensorGap(SensRadMax,SensRadMin,dPhi);
ypos = - SensRadMin + delta_NoOverlap;

std::cout << "-------------------------------------------------" << std::endl;
std::cout << "-----------  placment  sensor partrs  in CF suppurt : "    << std::endl;
std::cout << "-------------------------------------------------" << std::endl;
std::cout << "  some importent sizes " << std::endl;
std::cout << " ypos            : " << ypos <<  std::endl;
std::cout << " CFhz            : " << CFhz <<  std::endl;
std::cout << " hFanoutFrontDZ  : " << hFanoutFrontDZ <<  std::endl;
std::cout << " hMetalDZ        : " <<  hMetalDZ<<  std::endl;
std::cout << " hSiliconDZ      : " <<  hSiliconDZ<<  std::endl;
std::cout << " hFanoutBackDZ   : " <<  hFanoutBackDZ<<  std::endl;


SensorAtCF = -CFhz +  hFanoutFrontDZ ; 
std::cout << " 1.SensorAtCF : " << SensorAtCF <<  std::endl;
new G4PVPlacement ( 0, G4ThreeVector(0., ypos, ( SensorAtCF)),logicFanoutFrnt , "FunOut0", logicCF, false,2, overlap_check);

SensorAtCF = SensorAtCF + hFanoutFrontDZ +  hMetalDZ ; 
std::cout << " 2.SensorAtCF : " << SensorAtCF <<  std::endl;
new G4PVPlacement ( 0, G4ThreeVector(0., ypos, ( SensorAtCF)),logicMetalV , "PadMetal0", logicCF, false,2, overlap_check);

SensorAtCF = SensorAtCF + hMetalDZ + hSiliconDZ; 
std::cout << " 3.SensorAtCF : " << SensorAtCF <<  std::endl;
new G4PVPlacement ( 0, G4ThreeVector(0.,ypos , ( SensorAtCF)),logicSensorV , "SensorV0", logicCF, false,2, overlap_check);

SensorAtCF = SensorAtCF + hSiliconDZ + hMetalDZ;
std::cout << " 4.SensorAtCF : " << SensorAtCF <<  std::endl;
new G4PVPlacement ( 0, G4ThreeVector(0., ypos, ( SensorAtCF)),logicMetalV , "PadMetal1", logicCF, false,2, overlap_check);

SensorAtCF =  SensorAtCF + hMetalDZ + hFanoutBackDZ ; 
std::cout << " 5.SensorAtCF : " << SensorAtCF <<  std::endl;
new G4PVPlacement ( 0, G4ThreeVector(0.,ypos , ( SensorAtCF)),logicFanoutBack , "FunOut1", logicCF, false,2, overlap_check);


//-------------------------------------------------------------
//		placment of sensor & CF in base units 1,2,5
//-------------------------------------------------------------
 std::cout << "-------------------------------------------------" << std::endl;
 std::cout << "-----------  placment of CF suppurt in air base units : "    << std::endl;
 std::cout << "-------------------------------------------------" << std::endl;
  new G4PVPlacement ( 0, G4ThreeVector(0., 0.,airhz - CFhz ),logicCF , "CF0",logicBaseUnit , false, 0,overlap_check); //1
  new G4PVPlacement ( 0, G4ThreeVector(0., 0.,airhz - CFhz ),logicCF , "CF0",logicBaseUnit_PL , false, 0,overlap_check); //2
  new G4PVPlacement ( 0, G4ThreeVector(0., 0.,airhz1mm - CFhz ),logicCF , "CF0",logicBaseUnitFor1mm , false, 0,overlap_check); //5

  new G4PVPlacement ( 0, G4ThreeVector(0., 0.,airhz - CFhz ),logicCF , "CF0",logicBaseUnit_S , false, 0,overlap_check); //1

     // put LCAL to world
     //
  //G4double zposLC = 1000. *mm;
  
  G4double zposLC =  0.0 *mm; // for 2016 from 4th Telescope plane to box 
  // G4double zyposLC = -(22./64.)*(SensRadMax-SensRadMin) *mm;
  G4double zyposLC = 0.0 *mm;
  // G4double DUTairhx = 100.0 *mm;
  // G4double DUTairhy = 100.0 *mm;
  // G4double DUTairhz = 54 *mm; // need to defint that  as Lcal_tungsten_hdz + C (spacer) 
  G4double DUTextrahz = 0.002 *mm; 
  // G4double ShiftForSigleMIP = 10 *mm; 

	G4double zpos_PSC = 1130. *mm; // for 2016 from 4th Telescope plane to Plastic scintilator  
	G4double zpos_TR_PSC = 20. *mm; // for 2016 from 4th Telescope plane to TRiger Plastic scintilator  
	G4int n_layers = Setup::Lcal_n_layers;
//-----------------
// for general use 
//-----------------

 std::cout << "-------------------------------------------------" << std::endl;
 std::cout << "-----------plaising layers : "       << n_layers  << std::endl;
 std::cout << "-------------------------------------------------" << std::endl;

int Is_SC = 0; 
int Is_TSC = 0; 
int Is_stag = 0; 
G4double ypos_stag[30] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}; 
std::istringstream iss((std::string)Setup::TBeam_senrio);
G4int iplacelayer = 0; // for sensors 
G4int iplaceElements = 0; // for all 
std::string delimiter = ":";
    while ((iss)&& (iplacelayer < n_layers))
    {
    	std::cout<<"!!!!!!!!!!!!!!!Put layer!!!!!!!!!!!!!!\n";
		std::stringstream placement_name("");
		std::string Osub;
		iss >> Osub;
		std::cout << "Substring:" << Osub << std::endl;
		std::string delimiter = ":";
		std::string tmpsub = Osub.substr(0, Osub.find(delimiter));
		std::string sub;
		int i_air = 0 ;

	if ( tmpsub.length() < Osub.length()) // check if the ":" for air gup size for tracker is there 
	{
		Osub.erase(0, 1 + tmpsub.length());
		i_air = std::atoi(Osub.c_str()) *mm;
		sub = tmpsub;
	}

	else {
		sub = tmpsub;
	}

  if (sub == "S") {
      placement_name << "DUTAS" << iplacelayer;
      new  G4PVPlacement ( 0, G4ThreeVector( 0., zyposLC+ypos_stag[iplacelayer], zposLC + airhz),logicBaseUnit_S, placement_name.str().c_str(), logicWorld, 0, iplacelayer+1, overlap_check);
      G4cout<< " placed "<<iplaceElements << " as sensor number  : " << iplacelayer<< " at z = " << zposLC<<" layer with name "<< placement_name.str().c_str() << G4endl;
      zposLC= zposLC+ 2*airhz + DUTextrahz;
      iplacelayer++;
      }
	else if (sub == "AS") {
			placement_name << "DUTAS" << iplacelayer;
			new  G4PVPlacement ( 0, G4ThreeVector( 0., zyposLC+ypos_stag[iplacelayer], zposLC + airhz),logicBaseUnit, placement_name.str().c_str(), logicWorld, 0, iplacelayer+1, overlap_check);
			G4cout<< " placed "<<iplaceElements << " as sensor number  : " << iplacelayer<< " at z = " << zposLC<<" layer with name "<< placement_name.str().c_str() << G4endl;
			zposLC= zposLC+ 2*airhz + DUTextrahz;
			iplacelayer++;
			}

	else if(sub == "AS_PL"){
			placement_name << "DUTASP" << iplacelayer;
			new  G4PVPlacement ( 0, G4ThreeVector( 0., zyposLC+ypos_stag[iplacelayer], zposLC + airhz  ),logicBaseUnit_PL, placement_name.str().c_str(), logicWorld, 0, iplacelayer+1,overlap_check);
			G4cout<< " placed "<<iplaceElements << " as sensor number  : " << iplacelayer<<" at z = " << zposLC<<" layer with name "<< placement_name.str().c_str() << G4endl;
			zposLC= zposLC+2*airhz + DUTextrahz;			
			iplacelayer++;
			}		

	else if (sub == "A_PL") {
			placement_name << "ABP" << iplaceElements;
			new  G4PVPlacement ( 0, G4ThreeVector( 0., zyposLC, zposLC +airhz ),logicBaseUnitONLYAB_PL, placement_name.str().c_str(), logicWorld, 0,iplaceElements +1, overlap_check);
			G4cout<< " placed "<<iplaceElements <<" at z = " << zposLC<<" layer with name "<< placement_name.str().c_str() << G4endl;
			zposLC= zposLC+ 2* airhz + DUTextrahz;
			//iplaceA++;
			}

	else if(sub == "A"){
			placement_name << "AB" << iplaceElements;
			new  G4PVPlacement ( 0, G4ThreeVector( 0., zyposLC, zposLC + airhz ),logicBaseUnitONLYAB, placement_name.str().c_str(), logicWorld, 0, iplaceElements+1,overlap_check);
			G4cout<< " placed "<<iplaceElements <<" at z = " << zposLC<<" layer with name "<< placement_name.str().c_str() << G4endl;
			zposLC= zposLC+ 2*airhz + DUTextrahz;			
//iplaceA++;
			}		

	else if(sub == "JSM"){
			placement_name << "DUTJSM" << iplacelayer;	     		
			new  G4PVPlacement ( 0, G4ThreeVector( 0., zyposLC+ypos_stag[iplacelayer], zposLC+ airhz1mm ),logicBaseUnitFor1mm, placement_name.str().c_str(), logicWorld, 0, iplacelayer+1, overlap_check);
			G4cout<< " placed "<<iplaceElements << " as sensor number  : " << iplacelayer<<" at z = " << zposLC<<" layer with name "<< placement_name.str().c_str() << G4endl;
			zposLC= zposLC+ 2*airhz1mm + DUTextrahz + i_air; //add the air gup after 
			iplacelayer++;
			}
	else if(sub == "SC"){
			if(Is_SC == 0 ){
				placement_name << "SC" ;	     		
				new  G4PVPlacement ( 0, G4ThreeVector( 0., zyposLC, zpos_PSC ),logicBaseUnit_PSC, placement_name.str().c_str(), logicWorld, 0, 1, overlap_check);
				G4cout<< " placed SC" <<" at z = " << zpos_PSC <<" layer with name "<< placement_name.str().c_str() << G4endl;
				Is_SC =1; // can by placed only 1 time  
			}
			}
	else if(sub == "TSC"){
			if(Is_TSC == 0 ){
				placement_name << "TR_SC" ;	     		
				new  G4PVPlacement ( 0, G4ThreeVector( 0., zyposLC, zpos_TR_PSC ),logicBaseUnit_TR_PSC, placement_name.str().c_str(), logicWorld, 0, 1, overlap_check);
				G4cout<< " placed SC" <<" at z = " << zpos_TR_PSC <<" layer with name "<< placement_name.str().c_str() << G4endl;
				Is_TSC =1; // can by placed only 1 time  
			}
			}
	else if(sub == "stag"){
			if(Is_stag == 0 ){
				G4cout<< "Stag should applay before all sensor layers " << G4endl;
                                double TB_stag[8] = {0.0, 0.0, 0.2*mm, -0.7 *mm, 1.5 *mm, -1.0 *mm, 0.0 ,0.0};// asked by marina
				//double TB_stag[8] = {-0.11 * mm, -1.26* mm, 0.46*mm, -0.275644 *mm, 1.87705 *mm, -1.2183 *mm, 0.53323 *mm,0 *mm};// proper stagging 
                               // double TB_stag[8] = {0.11 * mm, 1.26* mm, -0.46*mm, 0.275644 *mm, -1.87705 *mm, 1.2183 *mm, -0.53323 *mm,0 *mm};
				for(int istag = 0 ; istag < 8 ; istag++){
					ypos_stag[istag] = TB_stag[istag] ;
                                        G4cout << "Stugging: layer " << istag <<" movment :  " << ypos_stag[istag] << G4endl;

				}	     		

				Is_stag =1; // can by placed only 1 time  
			}
			}			
 	else {
		G4cout << "Substring: " << sub <<" in layer " << iplaceElements <<" is not recognized we will brack"<< G4endl;
		break;
		}
	iplaceElements++;
	
    } 

   // ---------------
   //  SENSITIVE DETECTOR
   //  ---------------
    G4SDManager* SDman = G4SDManager::GetSDMpointer();

    // Initialize the Sensitive Detector
    LCSensitiveDetector *SensDet = new LCSensitiveDetector("LumiCalSD",  // name
                                                            Cell0_Rad,    // inner LC radius
                                                            startPhi,     // start angle
                                                            CellPitch,    // radial cell size
                                                            sectorPhi,    // angular cell width
                                                            nCells,       // # cells in the rad dir
                                                            nSectors,     // # cells in the phi dir
                      				                              VirtualCell); // cell type real/virtual =  false/true
        
    SDman->AddNewDetector(SensDet);
    // the Cells are the sensitive detectors
    logicSensorV->SetSensitiveDetector(SensDet);

    // if ( VirtualCell )  logicSensorV->SetSensitiveDetector( SensDet );
    //   //logicCell->SetSensitiveDetector(SensDet);
    // else
    //   //logicSensorV->SetSensitiveDetector( SensDet );
    //   G4cout << "  there is no VirtualCell.... " << G4endl;


	// G4cout <<  " Test Beam setup done !  "  << G4endl;

    G4double fromedge_to_center = 14.0*cm;
    G4double shift_y = base_airy + fromedge_to_center;

    G4Transform3D tr1( G4RotationMatrix().rotateZ( 90.0*deg ),
                       G4ThreeVector( -shift_y, 0.0, 5.0*m));
    new G4PVPlacement(tr1,
    logicWorld,
    "LumiCal", // an updated string
    fLogicWorld,
    false,
    1,
    overlap_check); // copy number

    G4Transform3D tr2( G4RotationMatrix().rotateZ( -90.0*deg ),
                       G4ThreeVector( shift_y, 0.0, 5.0*m));
    new G4PVPlacement(tr2,
    logicWorld,
    "LumiCal", // an updated string
    fLogicWorld,
    false,
    2,
    overlap_check); // copy number
    
  // SDman = G4SDManager::GetSDMpointer();
  // std::cout<<"@@@@@@@@@@@@@@@@@@@@@@@@@ Pointer from detector Construction: "<<SDman->GetCollectionID("LumiCalSD")<<"\n";

   std::cout << "LumiCal Done! Size from the center to absorber: " << base_airy;

}

void DetectorConstruction::InitDetectorParameters()
{

  // Initialize LCAL local parameters

  Vacuum      = Setup::Vacuum;                // beam pipe filler
  Air         = Setup::Air;                   // world filler
  PLASTIC_SC  = Setup::PLASTIC_SC;                   // world filler
  Silicon     = Setup::Silicon;               // Sensor Cell
  Iron        = Setup::Iron;                  // material for lateral beam tubes
  Aluminium   = Setup::Alu;                   // material for LCAL mechanical parts
  Tungsten    = Setup::Wabsorber;              // Absorber plate
  Graphite    = Setup::Graphite;
  FanoutMatF  = Setup::FanoutMatF;            // Front Fanout material
  FanoutMatB  = Setup::FanoutMatB;            // Back Fanout material
  BeamPipeMat = Setup::Beryllium;
  LHcalMat    = Setup::Tungsten;
  BCalMat     = Setup::Tungsten;
  Mask_Mat    = Setup::Tungsten;
  if ( Setup::Build_Beampipe == "Yes" ) WorldMat  = Air;
  else WorldMat = Vacuum;
  if(Setup::LcalTBeam >0 ) WorldMat  = Air;
  rotAng      = Setup::Beam_Crossing_Angle / 2.;
  rotAng1     = 180.*deg - rotAng;
  rotAng2     = rotAng;
  // geometric parameters
  // LCAL
  VirtualCell = Setup::Lcal_virtual_cells;                // cell type 
  nLayers    = Setup::Lcal_n_layers;
  nSectors   = Setup::Lcal_n_sectors;
  nTiles     = Setup::Lcal_n_tiles;;
  nCells     = Setup::Lcal_n_rings;

  Lcal_zend  = Setup::Lcal_z_end;
  innerRad   = Setup::Lcal_inner_radius;
  Cell0_Rad  = Setup::Lcal_Cell0_radius;
  SensRadMin = Setup::Lcal_SensRadMin;        // silicon sensor r-min including dead edges
  SensRadMax = Setup::Lcal_SensRadMax;        // silicon sensor r-max

  deadPhi    = Setup::Lcal_sector_dead_gap;
  phi_offset = Setup::Lcal_layers_phi_offset;
 // derived
  hLumiCalDZ     = Setup::Lcal_hdz;          // half length of LCAL
  Lcal_zbegin    = Lcal_zend - 2.*hLumiCalDZ ;
  SensDZ         = Setup::Lcal_sensor_dz;
  CellPitch      = Setup::Lcal_CellPitch;
  layer_gap      = Setup::Lcal_layer_gap;         // air gap between layers
  hSiliconDZ     = Setup::Lcal_silicon_hdz;       // half thickness of the silicon
  hMetalDZ      = Setup::Lcal_pad_metal_thickness/2.; // half thickness of the pad metallization
  hFanoutFrontDZ = Setup::Lcal_fanoutF_hdz;       // half thickness fanout front
  hFanoutBackDZ  = Setup::Lcal_fanoutB_hdz;       // half thickness fanout back 
  hTungstenDZ    = Setup::Lcal_tungsten_hdz;      // half thickness absorber
  hLayerDZ       = Setup::Lcal_layer_hdz;         // half thickness of the layer
  hSensorDZ      = hSiliconDZ + hMetalDZ;        // sensor half thickness including metallization
  hAbsorberPitchDZ = Setup::Lcal_absorber_pitch / 2.; // half pitch between absorber layer (defulte 1 mm) 
  startPhi  = Setup::Lcal_start_phi;
  endPhi    = Setup::Lcal_end_phi;
  sectorPhi = Setup::Lcal_sector_dphi;
  tilePhi   = endPhi / G4double( nTiles );
  assert ( SensRadMin >= (innerRad / cos ( tilePhi /2. )));
  // FE - space
  FECave_hDZ =   Setup::Lcal_ChipCaveDepth/2.;
  FECave_rmin =  Setup::Lcal_SensRadMax ;
  FECave_rmax =  Setup::Lcal_FEChip_rmax;
  PCB_hDZ     =  Setup::Lcal_PCB_thickness / 2.;
  FEChip_hDZ  =  Setup::Lcal_FEChip_space / 2. ; 
  G4double  FE_size = FECave_rmax - FECave_rmin;
  Lcal_extra_size = (  FE_size > Setup::Lcal_space_for_ears ) ? FE_size : Setup::Lcal_space_for_ears;
  outerRad   = SensRadMax + Lcal_extra_size;
  assert ( FECave_hDZ >= ( PCB_hDZ + FEChip_hDZ )); 
  // LHcal
  LHcalToLCalDist = 45. *mm;
  LHcal_zbegin = Lcal_zend + LHcalToLCalDist ;
  LHcal_rmin  =   93.0 *mm;
  LHcal_rmax  =  330.0 *mm;
  LHcal_hDZ   = (525.0 /2.) *mm;
  LHcal_zend = LHcal_zbegin + 2.*LHcal_hDZ ;
  // BCal
  BCalToLHcalDist = 390. *mm;
  BCal_rmin       =  25. *mm;
  BCal_rmax       = 150. *mm;
  BCal_zbegin     = LHcal_zend + BCalToLHcalDist;
  BCal_dPairMoni  =   1. *mm;   
  BCal_dgraphite  = 100. *mm;   
  BCal_zlength    = 117. *mm;
  BCal_sphi       = 200. *deg;    
  BCal_dphi       = 320. *deg;    
  
  // beam
  pipe_th       = Setup::Beam_pipe_thickness;
  Lcal_inner_pipe_zend     = BCal_zbegin - 5.0*mm; 
  LcalToBeamTol = Setup::Lcal_to_BeamPipe_clearance; 
  BCal_inner_outube_rmax = BCal_rmin - 2.*mm;
  BCal_inner_intube_rmax = 15. *mm;
  // mask
  Mask_zstart = LHcal_zend + 5.*mm;
  Mask_thickness = 50.*mm;
  Mask_hDX  = 290.*mm;
  Mask_hDY  = 290.*mm;
  Mask_rmin_at_zstart = Setup::Lcal_inner_radius + Setup::Lcal_to_BeamPipe_clearance;
}


// =========== CELL PARAMETERIZATION ============
LCCellParam::LCCellParam(G4int    NoCells,
                         G4double startR,
                         G4double endR,
			 G4double SensHalfZ,
                         G4double SihalfZ,
                         G4double AlhalfZ,
                         G4double clipSize,
                         G4double startPhi,
                         G4double deltaPhi)
{
    lNoCells  = NoCells;
    lstartR   = startR + clipSize;
    lendR     = endR - clipSize;
    lSenshalfZ= SensHalfZ;
    lSihalfZ  = SihalfZ;
    lAlhalfZ  = AlhalfZ;
    lstartPhi = startPhi;
    ldeltaPhi = deltaPhi;
    lclipSize = ( clipSize > 0. ) ?  clipSize + 0.05: 0.;
    ldeltaR   = (lendR - lstartR)/(G4double)NoCells;
}

LCCellParam::~LCCellParam() {}

void LCCellParam::ComputeTransformation(const G4int, G4VPhysicalVolume *physVol ) const
{
  physVol->SetTranslation(  G4ThreeVector(0., 0., 0));
  physVol->SetRotation(0);
}

void LCCellParam::ComputeDimensions(G4Tubs &Cell, const G4int copyNo,
                                    const G4VPhysicalVolume* physVol ) const
{
    G4double innerRad = lstartR + copyNo * ldeltaR;
    G4double midRad   = innerRad + ldeltaR/2;
    G4double outerRad = innerRad + ldeltaR;
    G4double cutPhi   = atan(lclipSize / midRad) *rad ;
    G4double delPhi   = ldeltaPhi - cutPhi;
    G4double startPhi = lstartPhi;
    G4double metPhi   = atan(0.05/midRad)*rad;
    G4double halfZ    = lSihalfZ;

    G4String MotherLogName = physVol->GetMotherLogical()->GetName();

    if ( MotherLogName.contains("MetSector") ){
	   innerRad +=  0.05;
	   outerRad -=  0.05;
           delPhi   -= metPhi;
           startPhi += metPhi;
           halfZ     = lAlhalfZ; 
	 }

    Cell.SetInnerRadius(innerRad);
    Cell.SetOuterRadius(outerRad);
    Cell.SetZHalfLength(halfZ); 

    if ( MotherLogName.contains("Sector1") )
      {      Cell.SetStartPhiAngle( lstartPhi + cutPhi);
	     Cell.SetDeltaPhiAngle( delPhi    ); }
    else if ( MotherLogName.contains("Sector4") )
      {      Cell.SetStartPhiAngle( lstartPhi );
	     Cell.SetDeltaPhiAngle( delPhi    ); }
    else 
      {      Cell.SetStartPhiAngle( lstartPhi );
	     Cell.SetDeltaPhiAngle( ldeltaPhi ); }

}
