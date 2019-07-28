//
/// \brief Implementation of the RunAction class
//

#include "RunAction.hh"
#include "DetectorConstruction.hh"
#include "PrimaryGeneratorAction.hh"
// #include "HistoManager.hh"
#include "G4UImessenger.hh"

#include "G4AnalysisMessenger.hh"

#include "Run.hh"

#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4UnitsTable.hh"
#include "G4EmCalculator.hh"

#include "Randomize.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include <iomanip>

#include "GlobalVars.hh"


// //static const double     pi  = 3.14159265358979323846;

// RunAction::RunAction(DetectorConstruction* det, LCRootOut *RO, PrimaryGeneratorAction* kin)
// :G4UserRunAction(),fDetector(det), fPrimary(kin), fRun(0), fFileName("lxphoton_out_vg1")
// { 
//   // Book predefined histograms
//   //fHistoManager = new HistoManager();
  
//   // G4AnalysisManager* man = G4AnalysisManager::Instance();
//   // // Open an output file
//   // man->SetFileName(fFileName);
//   // man->SetVerboseLevel(1);
//   // man->SetActivation(true);    // enabl

//   RootOut = RO;
// }

RunAction::RunAction(DetectorConstruction* det, PrimaryGeneratorAction* kin)
:G4UserRunAction(),fDetector(det), fPrimary(kin), fRun(0)
{ 
  // Book predefined histograms
  //fHistoManager = new HistoManager();
  
  // G4AnalysisManager* man = G4AnalysisManager::Instance();
  // // Open an output file
  // man->SetFileName(fFileName);
  // man->SetVerboseLevel(1);
  // man->SetActivation(true);    // enabl

  // RootOut = 0;
}


RunAction::~RunAction()
{ 
  //delete fHistoManager;
}



G4Run* RunAction::GenerateRun()
{ 
  fRun = new Run(fDetector);
  //fRun->SetHistoManager(fHistoManager);
  return fRun;
}



void RunAction::BeginOfRunAction(const G4Run* Run)
{
  // save Rndm status
  ////  G4RunManager::GetRunManager()->SetRandomNumberStore(true);
  if (isMaster) G4Random::showEngineStatus();
     
  // keep run condition
  if ( fPrimary ) { 
    G4ParticleDefinition* particle = fPrimary->GetParticleGun()->GetParticleDefinition();
    G4double energy = fPrimary->GetParticleGun()->GetParticleEnergy();
    fRun->SetPrimary(particle, energy);
    std::cout<<"~~~~~~~~~~~~~~~~~~~~~~~~~~~ E: "<<energy<<" "<<particle<<"\n";
  }

}



void RunAction::EndOfRunAction(const G4Run* Run)
{  
  // print Run summary
  //
  if (isMaster) fRun->EndOfRun();    
      

  if (isMaster) G4Random::showEngineStatus();


}