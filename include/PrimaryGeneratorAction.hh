//
/// \brief Definition of the PrimaryGeneratorAction class
//


#ifndef PrimaryGeneratorAction_h
#define PrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4ParticleGun.hh"
#include "globals.hh"

class G4Event;
class DetectorConstruction;
class PrimaryGeneratorMessenger;
class PrimarySpectra;
class LuxeTestGenerator;

class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
  public:
    PrimaryGeneratorAction(DetectorConstruction*);    
   ~PrimaryGeneratorAction();

  public:
    void SetDefaultKinematic();
    void SetBeamType(G4String val);
    void SetSpectraFile(G4String val); 
    void SetMCParticleFile(G4String val, const G4bool list = false);
    virtual void GeneratePrimaries(G4Event*);
    G4ParticleGun* GetParticleGun() const {return fParticleGun;}

    void GenerateGaussian(G4Event* anEvent);
    void GenerateMono(G4Event* anEvent);
    void GeneratefromMC(G4Event* anEvent);

    G4double GetX0() const {return fx0;}
    G4double GetY0() const {return fy0;}
    G4double GetZ0() const {return fz0;}
    G4double GetSigmaX() const {return fsigmax;}
    G4double GetBetaX()  const {return fbetax;}
    G4double GetEmittanceX() const {return femittancex;}
    
    G4double GetSigmaY() const {return fsigmay;}
    G4double GetBetaY()  const {return fbetay;}
    G4double GetEmittanceY() const {return femittancey;}
    
    void SetSigmaX(const G4double sigma);
    void SetSigmaY(const G4double sigma);

    G4double GetSigmaZ() const {return fsigmaz;}
    
  protected:
    enum tBeamType : G4int {beamGauss, beamMono, beamMC};
    
    G4double TestHitTarget(const std::vector <double> &pp, const double *vtx);
    
  private:
    G4ParticleGun*         fParticleGun;
    DetectorConstruction*  fDetector;
    G4int                  fBeamType;
    G4double               fx0, fy0, fz0;
    G4double               fsigmax, fsigmay, fsigmaz;
    G4double               femittancex, femittancey;
    G4double               fbetax, fbetay;

    PrimarySpectra        *fSpectra;
    LuxeTestGenerator     *lxgen;
    G4String               fMCfile;
    G4bool                 fMCList;
    G4int                  fnfixparticles;   
    
    PrimaryGeneratorMessenger* fGunMessenger;     
};



#endif


