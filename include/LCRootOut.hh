#ifndef LCROOTOUT_HH_
#define LCROOTOUT_HH_ 1

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <sstream>

#include "tls.hh"
#include "G4String.hh"

#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4ios.hh"

#include "LCSensitiveDetector.hh"
#include "LCHit.hh"
#include "g4root.hh"

#include "G4UImessenger.hh"

#include "G4UIcmdWithAString.hh"
// namespaces
using namespace std;
// using namespace TMath;

class LCRootMesseger;

class LCRootOut
{
public:
    LCRootOut();
    // LCRootOut(const G4String fName );
   ~LCRootOut();

public:
 
  void Init();                                                         // opens a file, creates a Tree 
  void ProcessEvent(const G4Event* event, LCHitsCollection *HitsColl);
  void SetName(G4String newValue) {RootOutFile = newValue;}

  void ClearData();
  void End();                                                         // writes to file and closes it
  void TestPlaneFill(G4double ,G4double ,G4double,G4double ,G4double ,G4double,G4double,G4int);



private:
  LCRootMesseger* fMessenger;
  G4String RootOutFile;
 
private:

  std::vector<G4double> Tracks_pX;
  std::vector<G4double> Tracks_pY;
  std::vector<G4double> Tracks_pZ;
  std::vector<G4int> Tracks_ID;
  std::vector<G4int> Tracks_PDG;

  std::vector<G4int> Hits_cellID;
  std::vector<G4double> Hits_eHit;
  std::vector<G4double> Hits_xCell;
  std::vector<G4double> Hits_yCell;
  std::vector<G4double> Hits_zCell;
  std::vector<G4double> Hits_xHit;
  std::vector<G4double> Hits_yHit;
  std::vector<G4double> Hits_zHit;
  std::vector<G4double> Hits_TOF;
  std::vector<G4int> Hits_Sensor;


  void Fill_Tracks_pX(const std::vector<G4double> &v){ Tracks_pX = v;};
  void Fill_Tracks_pY(const std::vector<G4double> &v){ Tracks_pY = v;};
  void Fill_Tracks_pZ(const std::vector<G4double> &v){ Tracks_pZ = v;};
  void Fill_Tracks_ID(const std::vector<G4int> &v){ Tracks_ID = v;};
  void Fill_Tracks_PDG(const std::vector<G4int> &v){ Tracks_PDG = v;};

  void Fill_Hits_cellID(const std::vector<G4int> &v){ Hits_cellID = v;};
  void Fill_Hits_eHit(const std::vector<G4double> &v){ Hits_eHit = v;};
  void Fill_Hits_xCell(const std::vector<G4double> &v){ Hits_xCell = v;};
  void Fill_Hits_yCell(const std::vector<G4double> &v){ Hits_yCell = v;};
  void Fill_Hits_zCell(const std::vector<G4double> &v){ Hits_zCell = v;};
  void Fill_Hits_xHit(const std::vector<G4double> &v){ Hits_xHit = v;};
  void Fill_Hits_yHit(const std::vector<G4double> &v){ Hits_yHit = v;};
  void Fill_Hits_zHit(const std::vector<G4double> &v){ Hits_zHit = v;};
  void Fill_Hits_TOF(const std::vector<G4double> &v){ Hits_TOF = v;};
  void Fill_Hits_Sensor(const std::vector<G4int> &v){ Hits_Sensor = v;};

};

class LCRootMesseger : public G4UImessenger
{

  public:
    LCRootMesseger(LCRootOut* root);
   ~LCRootMesseger();
    
    virtual void SetNewValue(G4UIcommand*, G4String);
    
  private:
    LCRootOut* fHistManager;
    G4UIcmdWithAString* fRootName;
    
};

#endif;