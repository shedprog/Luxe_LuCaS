//
/// \brief Definition of the DetectorConstruction class
//

#ifndef DetectorConstruction_h
#define DetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "G4LogicalVolume.hh"
#include "globals.hh"
#include "G4Cache.hh"

#include "G4PVDivision.hh"
#include "G4VPVParameterisation.hh"
#include "G4NistManager.hh"


class G4Box;
class G4VSolid;
class G4VPhysicalVolume;
class G4Material;
class G4MaterialCutsCouple;
class G4UniformMagField;
class G4UserLimits;


class DetectorConstruction : public G4VUserDetectorConstruction
{
  public:

    DetectorConstruction();
   ~DetectorConstruction();

  public:

     void SetWorldMaterial(G4String);
     void SetWorldSizeZ   (G4double);
     void SetWorldSizeXY  (G4double);

     void SetMagField(G4double);
     void DefineMaterials();

     virtual G4VPhysicalVolume* Construct();
     // virtual void ConstructSDandField();



  public:

     // void PrintCalorParameters();
     G4VPhysicalVolume* ConstructCalorimeter();
     G4Material* GetWorldMaterial()     {return fWorldMaterial;};
     G4double    GetWorldSizeZ()    const {return fWorldSizeZ;};
     G4double    GetWorldSizeXY()   const {return fWorldSizeXY;};

     const G4VPhysicalVolume* GetphysiWorld() {return fPhysiWorld;};

     G4double GetMagnetZend() { return fMagnetZPos + fMagnetSizeZ/2.0; };

  protected:
     enum tTargetType : G4int {tfoil, twire};


  private:


     G4Material*        fWorldMaterial;
     G4double           fWorldSizeZ;
     G4double           fWorldSizeXY;

     G4bool             fDefaultWorld;

     G4Box*             fSolidWorld;
     G4LogicalVolume*   fLogicWorld;
     G4VPhysicalVolume* fPhysiWorld;


     G4double           fMagnetFieldValue; //2000 - 13000 Gauss
     G4Box*             fSolidMagnet;
     G4LogicalVolume*   fLogicMagnet;
     G4VPhysicalVolume* fPhysiMagnet;
     G4double           fMagnetSizeX, fMagnetSizeY, fMagnetSizeZ;
     G4double           fMagnetZPos;

  private:

     void ConstructMagnet();

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
//~~~~~~~~~~~LumiCal Part~~~~~~~~~~~~~~~!
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!

  private:

    void ConstructLumiCal();
    void InitDetectorParameters();
    void BuildTBeamPT16();

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ LUMUCAL WORLD~~~~~~~~~~~~~~~~~! START
    // world
    G4LogicalVolume *logicWorld;
    G4VPhysicalVolume *physiWorld;

    // Beam Pipe
    G4LogicalVolume *logicLCalInnerTube;
    G4LogicalVolume *logEndVac;

    // LumiCal Shell
    G4LogicalVolume *logicWholeLC;

    // layer = SENSORS + TUNGSTEN ABSORBER + SUPPORT
    // This volume is sensor, tungsten plate and connection between them.
    // It is itterated inside LumiCal volume (logicWholeLC)
    G4LogicalVolume *logicLayer;

    // TUNGSTEN ABSORBER PLATE
    G4LogicalVolume *logicAbsorber;

    // Front and back FANOUT LAYERS
    G4LogicalVolume *logicFanoutFrnt;
    G4LogicalVolume *logicFanoutBack;

    // SENSOR
    G4LogicalVolume *logicSensor;
    G4LogicalVolume *logicSensorV;
    G4LogicalVolume *logicSensorEnv;

    // Metallization mother volume
    G4LogicalVolume *logicMetSec;
    G4LogicalVolume *logicMetalV;
    G4LogicalVolume *MetSector1;
    G4LogicalVolume *MetSector2;
    G4LogicalVolume *MetSector4;

    // SECTOR
    // 48 sectors per sensor; contains cells
    G4LogicalVolume *logicFESector;
    G4LogicalVolume *logicChip;
    G4LogicalVolume *logicPCB;
    G4LogicalVolume *logicSector1;
    G4LogicalVolume *logicSector2;
    G4LogicalVolume *logicSector4;

    // Cell and pad metalization volumes
    // 64 cells per sector
    G4LogicalVolume *logicCell;
    G4LogicalVolume *logicCellMet;

    // FE chips electronics mother volume
    G4LogicalVolume *logicFEmother;

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ LUMUCAL WORLD~~~~~~~~~~~~~~~~~! END

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ LUMUCAL Matirial~~~~~~~~~~~~~~! start
    G4Material  *Vacuum,                // World
    *Air,                   // Filler
    *PLASTIC_SC,
    *Silicon,               // Cell
    *Aluminium,             // Cell metallization
    *BeamPipeMat,           // central Beam pipe mat
    *Iron,                  // lateral Beam pipe mat
    *Tungsten,              // Absorber plate
    *Copper,                // Fanout component
    *Graphite,              // BCal front shield
    *Kapton,                // Fanout component
    *Epoxy,                 // Fanout component
    *FanoutMatF,            // Front Fanout material
    *FanoutMatB,            // Back Fanout material
    *LHcalMat,
    *BCalMat,
    *Mask_Mat,
    *WorldMat;
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ LUMUCAL Matirial~~~~~~~~~~~~~~! End
    //New calorimeter setups
    G4double pix_x_size;
    G4double pix_y_size;
    // G4 regions created to set production cuts
    //
    G4Region *regionLCal, *regionBCal, *regionLHcal, *regionMask;
    G4bool VirtualCell;
    // rotation angle
    G4double rotAng, rotAng1, rotAng2;
    // geometric parameters
    G4double Lcal_zbegin;
    G4double Lcal_zend;
    G4double innerRad;
    G4double outerRad;
    G4double deadPhi;
    G4int    nLayers;
    G4int    nTiles;
    G4int    nSectors;
    G4int    nCells;
    G4double startPhi;
    G4double endPhi;
    G4double phi_offset;
    G4double Zcave;
    G4double CaveDepth;
    G4double hLumiCalDZ;          // half length of LCAL
    G4double Cell0_Rad;
    G4double SensRadMin;        // silicon sensor r-min
    G4double SensRadMax;        // silicon sensor r-max
    G4double SensZ0;
    G4double SensDZ ;           // z distance between sensors
    G4double CellPitch;
    G4double hLayerDZ;         // half thickness of the layer
    G4double layer_gap;         // air gap between layers
    G4double hSiliconDZ;       // half thickness of the siliconZ
    G4double hMetalDZ;         // half thickness of pad metallization
    G4double hSensorDZ;
    G4double hFanoutFrontDZ;       // half thickness fanout front
    G4double hFanoutBackDZ;       // half thickness fanout back
    G4double hTungstenDZ;      // half thickness absorber
    G4double hAbsorberPitchDZ;// pitch between absorber layer defulte 1 mm
    G4double sectorPhi;
    G4double tilePhi;
    G4double FECave_hDZ;
    G4double FECave_rmin ;
    G4double FECave_rmax ;
    G4double FEChip_hDZ;
    G4double PCB_hDZ;
    G4double Lcal_extra_size;
    // Beam pipe
    G4double Lcal_inner_pipe_zend;
    G4double pipe_th;
    G4double LcalToBeamTol;
    G4double BCal_inner_outube_rmax;
    G4double BCal_inner_intube_rmax;
    // LHcal
    G4double LHcalToLCalDist;
    G4double LHcal_zbegin;
    G4double LHcal_zend;
    G4double LHcal_hDZ;
    G4double LHcal_rmin;
    G4double LHcal_rmax;
    // BCal
    G4double BCalToLHcalDist;
    G4double BCal_rmin;
    G4double BCal_rmax;
    G4double BCal_zbegin;
    G4double BCal_dPairMoni;
    G4double BCal_dgraphite;
    G4double BCal_zlength;
    G4double BCal_sphi;
    G4double BCal_dphi;
    // mask
    G4double Mask_zstart;
    G4double Mask_thickness;
    G4double Mask_hDX;
    G4double Mask_hDY;
    G4double Mask_rmin_at_zstart;
    //

};

class LCCellParam : public G4VPVParameterisation
{
public:

    LCCellParam(G4int    NoCells,
                G4double startR,
                G4double endR,
        G4double SensHalfZ,
                G4double SihalfZ,
                G4double AlhalfZ,
                G4double clipSize,
                G4double startPhi,
                G4double deltaPhi);
    virtual ~LCCellParam();

public:

virtual  void ComputeTransformation(const G4int repNo, G4VPhysicalVolume *physVol) const;
virtual  void ComputeDimensions(G4Tubs &Cell, const G4int copyNo, const G4VPhysicalVolume* physvol) const;

private: // Dummy declarations to get rid of warnings
    void ComputeDimensions (G4Trd&,const G4int,const G4VPhysicalVolume*) const {}
    void ComputeDimensions (G4Trap&,const G4int,const G4VPhysicalVolume*) const {}
    void ComputeDimensions (G4Cons&,const G4int,const G4VPhysicalVolume*) const {}
    void ComputeDimensions (G4Sphere&,const G4int,const G4VPhysicalVolume*) const {}
    void ComputeDimensions (G4Orb&,const G4int,const G4VPhysicalVolume*) const {}
    void ComputeDimensions (G4Torus&,const G4int,const G4VPhysicalVolume*) const {}
    void ComputeDimensions (G4Para&,const G4int,const G4VPhysicalVolume*) const {}
    void ComputeDimensions (G4Hype&,const G4int,const G4VPhysicalVolume*) const {}
    void ComputeDimensions (G4Box&,const G4int,const G4VPhysicalVolume*) const {}
    void ComputeDimensions (G4Polycone&,const G4int,const G4VPhysicalVolume*) const {}
    void ComputeDimensions (G4Polyhedra&,const G4int,const G4VPhysicalVolume*) const {}

private:

    G4int    lNoCells;
    G4int    lNoLayers;
    G4int    lNoSectors;
    G4double lstartR;
    G4double lendR;
    G4double lSenshalfZ;
    G4double lSihalfZ;
    G4double lAlhalfZ;
    G4double lstartPhi;
    G4double ldeltaPhi;
    G4double lclipSize;
    G4double ldeltaR;
};


#endif
