//
/// \brief Implementation of the EventAction class
//

#include "EventAction.hh"
#include "Run.hh"
// #include "HistoManager.hh"

#include "G4RunManager.hh"
#include "G4Event.hh"

#include <sys/times.h>
#include <fstream>
#include "Setup.hh"

#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4HCofThisEvent.hh"
#include "G4VHitsCollection.hh"
#include "G4TrajectoryContainer.hh"
#include "G4VVisManager.hh"
#include "G4SDManager.hh"
#include "G4UImanager.hh"
#include "G4RunManager.hh"
#include "G4ios.hh"
#include "G4PrimaryVertex.hh"
#include "G4PrimaryParticle.hh"

#include "G4SystemOfUnits.hh"




EventAction::EventAction()
:G4UserEventAction(),
 fEnergyDeposit(0.),
 fTrakLenCharged(0.), fTrakLenNeutral(0.),
 fNbStepsCharged(0), fNbStepsNeutral(0),
 fTransmitFlag(0), fReflectFlag(0)
{ }

// EventAction::EventAction(LCRootOut * RO)
// :G4UserEventAction(),
//  fEnergyDeposit(0.),
//  fTrakLenCharged(0.), fTrakLenNeutral(0.),
//  fNbStepsCharged(0), fNbStepsNeutral(0),
//  fTransmitFlag(0), fReflectFlag(0)
// {    
// 	  collID = -1;
//     RootOut = RO;
// }

EventAction::~EventAction()
{ }


void EventAction::BeginOfEventAction(const G4Event* )
{
 // initialisation per event
 fEnergyDeposit  = 0.;
 fTrakLenCharged = fTrakLenNeutral = 0.; 
 fNbStepsCharged = fNbStepsNeutral = 0;
 fTransmitFlag   = fReflectFlag    = 0;    
 
 fKinetikEnergyGamma.clear();
 fKinetikEnergyCharged.clear();


}



void EventAction::EndOfEventAction(const G4Event* event)
{
  Run* run = static_cast<Run*>(G4RunManager::GetRunManager()->GetNonConstCurrentRun());
              
 run->AddEnergy(fEnergyDeposit);
 run->AddTrakLenCharg(fTrakLenCharged);
 run->AddTrakLenNeutr(fTrakLenNeutral);

 run->CountStepsCharg(fNbStepsCharged);
 run->CountStepsNeutr(fNbStepsNeutral);

 run->CountTransmit (fTransmitFlag);
 run->CountReflect  (fReflectFlag);


}



