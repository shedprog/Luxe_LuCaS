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
{
	collID = -1;
	RootOut = 0;
}

EventAction::EventAction(LCRootOut * RO)
:G4UserEventAction(),
 fEnergyDeposit(0.),
 fTrakLenCharged(0.), fTrakLenNeutral(0.),
 fNbStepsCharged(0), fNbStepsNeutral(0),
 fTransmitFlag(0), fReflectFlag(0)
{    
	collID = -1;
	RootOut = RO;
}

EventAction::~EventAction()
{ }

void EventAction::AddLeakEnergy ( G4double eleak )
{
  noTrackKilled++;
  LeakEnergy += eleak;
}
G4double EventAction::GetLeakEnergy (){ return LeakEnergy ;}

void EventAction::BeginOfEventAction(const G4Event* )
{
 // initialisation per event
 fEnergyDeposit  = 0.;
 fTrakLenCharged = fTrakLenNeutral = 0.; 
 fNbStepsCharged = fNbStepsNeutral = 0;
 fTransmitFlag   = fReflectFlag    = 0;    
 
 fKinetikEnergyGamma.clear();
 fKinetikEnergyCharged.clear();


// reset energy counter
LeakEnergy =0.;
noTrackKilled = 0;

// Use a Sensitive Detector manager to assign an ID #
// to the hits collections

SDman = G4SDManager::GetSDMpointer();
time(&_start);
if (collID < 0) {
// there is only 1 hits collection, so the name is constant
collID = SDman->GetCollectionID("LumiCalSD_HC");
}


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


// Check that a collection ID has been assigned
if (collID < 0) { return; }
//
G4int evtnum = (event->GetEventID())+1;

//
// report on track killed
//
if ( noTrackKilled > 0 ){
		G4cout << " Event: "<<evtnum + Setup::EventStartNumber 
			   <<" Back Energy Leak : "<<LeakEnergy / GeV <<" GeV"<<G4endl;
}


//
G4HCofThisEvent *HCE = event->GetHCofThisEvent();

LCHitsCollection *HitsColl = 0;

if (HCE) { 
HitsColl = (LCHitsCollection*)(HCE->GetHC(collID)); 
// std::cout << HitsColl->entries() << "\n";
// std::cout << "Press Enter to Continue";
// std::cin.ignore();
}

if (HitsColl) { // fill the ROOT Tree

if( RootOut ) {
	// if ( Setup::AccumulateEvents ){
	// 	RootOut->ProcEventAccumulate( HitsColl );
	// }
	// else {
		RootOut->ProcessEvent( event, HitsColl );
	// }
}
}


}



