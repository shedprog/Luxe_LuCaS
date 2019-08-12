
/// \brief Main program of the electromagnetic/TestEm5 example
//
// #define G4MULTITHREADED

// #ifdef G4MULTITHREADED
// #include "G4MTRunManager.hh"
// #else
#include "G4RunManager.hh"
// #endif

#include "G4UImanager.hh"
#include "Randomize.hh"

#include "DetectorConstruction.hh"


#include "G4StepLimiterPhysics.hh"
#include "QGSP_BERT.hh"
#include "G4EmStandardPhysics.hh"

#include "Setup.hh"


#ifdef G4VIS_USE
#include "G4VisExecutive.hh"
#endif

#ifdef G4UI_USE
#include "G4UIExecutive.hh"
#endif

#include "PrimaryGeneratorAction.hh"
#include "RunAction.hh"
#include "EventAction.hh"

#include "SteppingAction.hh"

#include "LCRootOut.hh"
bool isRoot = 1;

int main(int argc,char** argv) {

  //choose the Random engine
  G4Random::setTheEngine(new CLHEP::RanecuEngine);
  long seeds[2] = {12345L, 5L};
  G4Random::setTheSeeds(seeds, 3);

  Setup *theSetup = Setup::GetSetup();
  theSetup->SetupInit(argc, argv);    
  
  // Construct the default run manager
// #ifdef G4MULTITHREADED
//   G4MTRunManager* runManager = new G4MTRunManager;
//   G4int nThreads = G4Threading::G4GetNumberOfCores();
//   if (argc==3) nThreads = G4UIcommand::ConvertToInt(argv[2]);
//   runManager->SetNumberOfThreads(nThreads);
//   G4cout << "===== lxbeamsim is started with " 
//          <<  runManager->GetNumberOfThreads() << " threads =====" << G4endl;
// #else
// G4VSteppingVerbose::SetInstance(new SteppingVerbose);
  G4RunManager* runManager = new G4RunManager;
// #endif



  // set mandatory initialization classes
  DetectorConstruction* detector = new DetectorConstruction;
  // LCDetectorConstruction* detector = new LCDetectorConstruction;
  runManager->SetUserInitialization(detector);

  //==================== default physics
  // PhysicsList *plist = new PhysicsList();
  //==================== G4 physics
  G4VUserPhysicsList *plist = new QGSP_BERT;
  plist->SetDefaultCutValue(Setup::rangeCut);
//   G4VModularPhysicsList *plist = new QGSP_BERT;
//   plist->RegisterPhysics(new G4StepLimiterPhysics());
  runManager->SetUserInitialization(plist);
  

  if(!isRoot){
    
    PrimaryGeneratorAction* primary = new PrimaryGeneratorAction(detector);
    runManager->SetUserAction(primary);

    RunAction* runaction = new RunAction(detector,primary);
    runManager->SetUserAction(runaction); 

    EventAction* eventaction =  new EventAction();
    runManager->SetUserAction(eventaction);

  }else{

    LCRootOut *theRootOut = new LCRootOut();

    PrimaryGeneratorAction* primary = new PrimaryGeneratorAction(detector);
    runManager->SetUserAction(primary);
 
    RunAction* runaction = new RunAction(detector,theRootOut,primary);
    runManager->SetUserAction(runaction); 

    EventAction* eventaction =  new EventAction(theRootOut);
    runManager->SetUserAction(eventaction);
    
    SteppingAction* steps = new SteppingAction(theRootOut);
    runManager->SetUserAction(steps);

  }

  // get the pointer to the User Interface manager 
  G4UImanager* UI = G4UImanager::GetUIpointer();  
  if (argc!=1)   // batch mode  
    {
      G4String command = "/control/execute ";
      G4String fileName = argv[1];
      UI->ApplyCommand(command+fileName);
    }
    
  else           //define visualization and UI terminal for interactive mode
    { 
#ifdef G4VIS_USE
      G4VisManager* visManager = new G4VisExecutive;
      visManager->Initialize();
#endif    
     
#ifdef G4UI_USE
      G4UIExecutive * ui = new G4UIExecutive(argc,argv);      
      ui->SessionStart();
      delete ui;
#endif
     
#ifdef G4VIS_USE
      delete visManager;
#endif     
    }
    
  // job termination
  //  
  delete runManager;

  return 0;
}



