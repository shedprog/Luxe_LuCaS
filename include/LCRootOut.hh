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
// namespaces
using namespace std;
// using namespace TMath;

class LCRootOut
{
public:
    LCRootOut();
    LCRootOut(const G4String fName );
   ~LCRootOut();

public:
 
  void Init();                                                         // opens a file, creates a Tree 
  void ProcessEvent(const G4Event* event, LCHitsCollection *HitsColl);
  // void ProcEventAccumulate( LCHitsCollection *HitsColl);
  void ClearData();
  void End();                                                         // writes to file and closes it
  void SetAddresses();                                                // sets branch addresses in "UPDATE" mode
  void CreateNewTree();                                               // creates new Tree


private:
  // root output file name 
  G4String RootOutFile;
  G4AnalysisManager* analysisManager;
 

private:
  G4double vX, vY, vZ;
  G4int numPrim;       // number of primary particles
  G4int numHits;     // total number of hits
  //  caloHit
  G4double Etot[2];       // total energy deposit in arm per arm
  G4double Emax;          // max  energy deposit in cell

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


void FillTracks_pX(G4double pX_) {Tracks_pX.push_back(pX_);}
void FillTracks_pY(G4double pY_) {Tracks_pY.push_back(pY_);}
void FillTracks_pZ(G4double pZ_) {Tracks_pZ.push_back(pZ_);}
void FillTracks_ID(G4int ID_) {Tracks_ID.push_back(ID_);}
void FillTracks_PDG(G4int PDG_) {Tracks_PDG.push_back(PDG_);}

void FillcellID(G4int cellID_) {Hits_cellID.push_back(cellID_);}
void FilleHit(G4double eHit_) {Hits_eHit.push_back(eHit_);}
void FillxCell(G4double xCell_) {Hits_xCell.push_back(xCell_);}
void FillyCell(G4double yCell_) {Hits_yCell.push_back(yCell_);}
void FillzCell(G4double zCell_) {Hits_zCell.push_back(zCell_);}
void FillxHit(G4double xHit_) {Hits_xHit.push_back(xHit_);}
void FillyHit(G4double yHit_) {Hits_yHit.push_back(yHit_);}
void FillzHit(G4double zHit_) {Hits_zHit.push_back(zHit_);}
void FillTOF(G4double TOF_) {Hits_TOF.push_back(TOF_);}

//------------------------------------------------------------------------
};

#endif;