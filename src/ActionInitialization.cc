//
/// \brief Implementation of the ActionInitialization class

#include "ActionInitialization.hh"
#include "DetectorConstruction.hh"
#include "PrimaryGeneratorAction.hh"
#include "RunAction.hh"
#include "EventAction.hh"
#include "TrackingAction.hh"
#include "SteppingAction.hh"
#include "StackingAction.hh"
#include "SteppingVerbose.hh"



ActionInitialization::ActionInitialization(DetectorConstruction* det)
 : G4VUserActionInitialization(),fDetector(det)
{ 
    RootOut = 0;
}

ActionInitialization::ActionInitialization(DetectorConstruction* det, LCRootOut *RO)
 : G4VUserActionInitialization(),fDetector(det)
{ 
    RootOut = RO;
}

ActionInitialization::~ActionInitialization()
{ }



void ActionInitialization::BuildForMaster() const
{
 SetUserAction(new RunAction(fDetector, RootOut));
}




void ActionInitialization::Build() const
{
  PrimaryGeneratorAction* primary = new PrimaryGeneratorAction(fDetector);
  
  SetUserAction(primary);
 
  RunAction* runaction = new RunAction(fDetector,RootOut,primary);
  SetUserAction(runaction); 
  
  EventAction* eventaction = new EventAction(RootOut);
  SetUserAction(eventaction);

  SetUserAction(new TrackingAction(fDetector,eventaction));

  SetUserAction(new SteppingAction(fDetector,eventaction));

  SetUserAction(new StackingAction(eventaction));
}  



G4VSteppingVerbose* ActionInitialization::InitializeSteppingVerbose() const
{
  return new SteppingVerbose();
}  


