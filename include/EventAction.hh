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
/// \file electromagnetic/TestEm5/include/EventAction.hh
/// \brief Definition of the EventAction class
//
// $Id: EventAction.hh 76464 2013-11-11 10:22:56Z gcosmo $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef EventAction_h
#define EventAction_h 1

#include "G4UserEventAction.hh"
#include "globals.hh"

#include <vector>

#include "LCRootOut.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class EventAction : public G4UserEventAction
{
  public:
    EventAction();
   ~EventAction();

  public:
    virtual void BeginOfEventAction(const G4Event*);
    virtual void   EndOfEventAction(const G4Event*);
    
    void AddEnergy      (G4double edep)   {fEnergyDeposit  += edep;};
    void AddTrakLenCharg(G4double length) {fTrakLenCharged += length;};
    void AddTrakLenNeutr(G4double length) {fTrakLenNeutral += length;};
    
    void CountStepsCharg ()               {fNbStepsCharged++ ;};
    void CountStepsNeutr ()               {fNbStepsNeutral++ ;};
    
    void SetTransmitFlag (G4int flag) 
                           {if (flag > fTransmitFlag) fTransmitFlag = flag;};
    void SetReflectFlag  (G4int flag) 
                           {if (flag > fReflectFlag)   fReflectFlag = flag;};
                           
    void AddKinEnegyCharged (const G4double ek) { fKinetikEnergyCharged.push_back(ek); };                 
    void AddKinEnegyGamma (const G4double ek) { fKinetikEnergyGamma.push_back(ek); };                 
                                             
        
  private:
    G4double fEnergyDeposit;
    G4double fTrakLenCharged, fTrakLenNeutral;
    G4int    fNbStepsCharged, fNbStepsNeutral;
    G4int    fTransmitFlag,   fReflectFlag; 
    
    std::vector<G4double> fKinetikEnergyGamma;
    std::vector<G4double> fKinetikEnergyCharged;

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//~~~~~~~~~~Lumi Cal modifications~~~~~~~
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
      EventAction(LCRootOut *RO);
      // LCEventAction();
      // ~LCEventAction();

  public:
         // void BeginOfEventAction(const G4Event* );
         // void EndOfEventAction(const G4Event* );
         void AddLeakEnergy ( G4double elost );
         G4double GetLeakEnergy();

  private:
      G4int collID;       // ID number for a collection of calorimeter hits
      G4String collName;  // hits collection name
      LCRootOut *RootOut; // Handles writing the ROOT tree
      G4SDManager *SDman;
      time_t  _start, _end;
      G4int noTrackKilled;
      G4double LeakEnergy;  // Accumulates energy leak

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

    
