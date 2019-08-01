//
/// \brief Implementation of the RunAction class
//

#include "RunAction.hh"
#include "DetectorConstruction.hh"
#include "PrimaryGeneratorAction.hh"
// #include "HistoManager.hh"
#include "G4UImessenger.hh"

#include "G4AnalysisMessenger.hh"

// #include "Run.hh"

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

RunAction::RunAction(DetectorConstruction* det, LCRootOut *RO, PrimaryGeneratorAction* kin)
:G4UserRunAction(), fDetector(det), fPrimary(kin)
{
  // fDetector = det;
  // fPrimary = kin;
  std::cout<<"@@@@@@@@@@@@@@@@@@@@ RunAction Init\n";
  RootOut = RO;
}

RunAction::RunAction(DetectorConstruction* det, PrimaryGeneratorAction* kin)
:G4UserRunAction(),fDetector(det), fPrimary(kin)

{
  RootOut = 0;
}


RunAction::~RunAction()
{ 
  //delete fHistoManager;
}



// G4Run* RunAction::GenerateRun()
// { 
//   // fRun = new Run(fDetector);
//   //fRun->SetHistoManager(fHistoManager);
//   // return fRun;
// }



void RunAction::BeginOfRunAction(const G4Run* Run)
{
    // std::cout<<"@@@@@@@@@@@@@@@@@@@@ RunAction BeginOfRunAction\n";
  // save Rndm status
  ////  G4RunManager::GetRunManager()->SetRandomNumberStore(true);
  if (isMaster) G4Random::showEngineStatus();
  
  // keep run condition
  // if ( fPrimary ) { 
  //   // std::cout<<"LOL !";
  //   G4ParticleDefinition* particle = fPrimary->GetParticleGun()->GetParticleDefinition();
  //   G4double energy = fPrimary->GetParticleGun()->GetParticleEnergy();
  //   fRun->SetPrimary(particle, energy);
  //   std::cout<<"~~~~~~~~~~~~~~~~~~~~~~~~~~~ E: "<<energy<<" "<<particle<<"\n";
  // }


  tms fTimeNow;

  G4cout << "Run " << Run->GetRunID() << " start." << G4endl;
  G4cout << " You may stop the run gently, any time you wish." << G4endl;
  G4cout << "  To do it, create in the current directory semaphore file " << G4endl;
  G4cout << "  named <aStopRun> , file content is ignored "<< G4endl; 
  G4cout << "  ( eg. shell> touch aStopRun  )"<< G4endl; 

  Setup::StartTime = times( &fTimeNow );
  Setup::NoOfEventsToProcess = Run->GetNumberOfEventToBeProcessed() - Setup::EventStartNumber;
  if( Setup::Particle_Generator.contains("pairs.") ) Setup::AccumulateEvents = true ;
 //
  // set random seed manually
  //
  G4long seed = Setup::StartTime;
  G4RandGauss::setTheSeed(seed);
  G4cout << " Starting with random seed " << seed << G4endl;
  Setup::StartRandomSeed = seed;

  if( RootOut ) RootOut->Init();

  Print("RUN_BEGIN", Run);

}



void RunAction::EndOfRunAction(const G4Run* Run)
{
  // std::cout<<"@@@@@@@@@@@@@@@@@@@@ RunAction EndOfRunAction\n";
    // std::cout<<"@@@@@@@@@@@@@@@@@@ END RUN ACTION\n";
  // print Run summary
  //
  // if (isMaster) fRun->EndOfRun();         

  if (isMaster) G4Random::showEngineStatus();

  Print("END_OF_RUN", Run);
  if( RootOut )  RootOut->End();

}

void RunAction::Print(G4String now, const G4Run* Run)
{
  time_t tnow = time(NULL);
  G4cout << "======================================================================"<<G4endl;
  if ( now == "END_OF_RUN")
    {
      tms fTimeNow;
      clock_t EndRunT = times( &fTimeNow );
      G4double diff = 10.*( EndRunT - Setup::StartTime ) *ms ;
      G4int EventsProcessed = Run->GetNumberOfEvent(); 
  G4cout << "|                End of Run  :  "<< Run->GetRunID()<< G4endl;
  G4cout << "|                  time now  :  "<< ctime(&tnow) ;
  G4cout << "|            Events Processed:  "<< EventsProcessed<< G4endl;
  if ( Setup::batchMode )
    {
  G4cout << "|             written to file:  "<< Setup::RootFileName << G4endl;
    }
  G4cout << "|                Time elapsed:  "<< diff / s << " seconds."<< G4endl;
  G4cout << "|    Time to process an event:  "<< diff/G4double(EventsProcessed) / s << " seconds."<< G4endl;
    }
  else 
    {
  G4cout << "|                     Begin of Run  : "<< Run->GetRunID()<< G4endl;
  G4cout << "|                         time now  : "<< ctime(&tnow) ;
  G4cout << "|                      Random Seed  : "<< Setup::StartRandomSeed << G4endl;
  G4cout << "| Global Parameters for the Run are : "           <<G4endl;
  G4cout << "---------------------------------------------------------------------"<<G4endl;
  G4cout << "|                   batchMode:  "<< Setup::batchMode << G4endl;
  G4cout << "|                   macroName:  "<< Setup::macroName << G4endl;
  G4cout << "|                  PrintLevel:  "<< Setup::PrintLevel << G4endl;
  G4cout << "|           Logging frequency:  "<< Setup::LogFreq << G4endl;
  G4cout << "|             PhysicsListName:  "<< Setup::PhysicsListName << G4endl;
  if ( Setup::batchMode )
    {
  G4cout << "|       ROOT output file name:  "<< Setup::RootFileName << G4endl;
  G4cout << "|       ROOT output open mode:  "<< Setup::RootFileMode << G4endl;
  G4cout << "|           accumulate events:  "<< Setup::AccumulateEvents << G4endl;
    }
    if (Setup::LcalTBeam){
  G4cout << "|             The senrio is  :  "<< Setup::TBeam_senrio << G4endl;
    }
  G4cout << "|                   SetupFile:  "<< Setup::SetupFile << G4endl;
  G4cout << "|         Beam_Crossing_Angle:  "<< Setup::Beam_Crossing_Angle / mrad << " [mrad]" << G4endl;
  G4cout << "|         Nominal field value:  "<< Setup::Nominal_Field_value / tesla << " [T]"<< G4endl;
  G4cout << "|          Particle generator:  "<< Setup::Particle_Generator<< G4endl;
  G4cout << "| Number of events to process:  "<< Setup::NoOfEventsToProcess << G4endl; 
  G4cout << "|   Detector components build:  "<< G4endl;
  G4cout << "|                   Beam Tube:  "<< Setup:: Build_Beampipe << G4endl;
  G4cout << "|                        LCAL:  "<< Setup:: Build_LCal << G4endl;
  G4cout << "|                       LHCAL:  "<< Setup:: Build_LHcal << G4endl;
  G4cout << "|                        BCAL:  "<< Setup:: Build_BCal << G4endl;
  G4cout << "|                        MASK:  "<< Setup:: Build_Mask << G4endl;
  G4cout << "|                    rangeCut:  "<< Setup::rangeCut/ mm << " [mm]"<<G4endl;
  G4cout << "|   Region Production Cuts:     "<< G4endl;
  G4cout << "|                        LCAL:  "<<  Setup::LCal_Region_Cut / mm <<" [mm]"<< G4endl;
  G4cout << "|                       LHCAL:  "<<  Setup::LHcal_Region_Cut / mm <<" [mm]"<< G4endl;
  G4cout << "|                        BCAL:  "<<  Setup::BCal_Region_Cut / mm <<" [mm]"<< G4endl;
  G4cout << "|                        MASK:  "<<  Setup::Mask_Region_Cut / mm <<" [mm]"<< G4endl;
   }
  G4cout << "========================================================================"<<G4endl;
}
