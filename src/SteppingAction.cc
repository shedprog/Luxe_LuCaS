// 
/// \brief Implementation of the SteppingAction class
//



#include "SteppingAction.hh"

#include "DetectorConstruction.hh"
#include "RunAction.hh"
#include "EventAction.hh"
// #include "HistoManager.hh"

#include "G4Step.hh"



SteppingAction::SteppingAction(DetectorConstruction* DET, EventAction* EA)
:G4UserSteppingAction(),fDetector(DET), fEventAction(EA)
{ }



SteppingAction::~SteppingAction()
{ }



void SteppingAction::UserSteppingAction(const G4Step* aStep)
{
 if (aStep->GetPreStepPoint()->GetTouchableHandle()->GetVolume() 
     != fDetector->GetAbsorber()) return;

 fEventAction->AddEnergy (aStep->GetTotalEnergyDeposit());
   
 G4double charge = aStep->GetTrack()->GetDefinition()->GetPDGCharge();
 if (charge != 0.) { 
   fEventAction->AddTrakLenCharg(aStep->GetStepLength());
   fEventAction->CountStepsCharg();
 } else {
   fEventAction->AddTrakLenNeutr(aStep->GetStepLength());
   fEventAction->CountStepsNeutr();
 }
}



