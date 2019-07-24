//
/// \brief Implementation of the RunAction class
//

#include "RunAction.hh"
#include "DetectorConstruction.hh"
#include "PrimaryGeneratorAction.hh"
#include "HistoManager.hh"
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

//static const double     pi  = 3.14159265358979323846;

RunAction::RunAction(DetectorConstruction* det, PrimaryGeneratorAction* kin)
:G4UserRunAction(),fDetector(det), fPrimary(kin), fRun(0), fFileName("lxphoton_out_vg1")
{ 
  // Book predefined histograms
  //fHistoManager = new HistoManager();
  
  G4AnalysisManager* man = G4AnalysisManager::Instance();
  // Open an output file
  man->SetFileName(fFileName);
  man->SetVerboseLevel(1);
  man->SetActivation(true);    // enabl
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



void RunAction::BeginOfRunAction(const G4Run*)
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
  // My histograms 14.06.2019
  
  // Get analysis manager
  G4AnalysisManager* man = G4AnalysisManager::Instance();
  // Open an output file
  // man->SetFileName(fFileName);
  man->OpenFile();
  // Create histogram(s)
  //extern double weight_fromMC;
  //G4cout<<"###############Weights of histograms######################: -> "<<weight_fromMC<<"\n";
  // // Initial particles
  // man->CreateH1( "1" , "Initial Gamma Energy [GeV]", 100, 0, 200*MeV);
  // man->CreateH1( "2" , "Initial e- Energy [GeV]", 100, 0, 200*MeV);
  // man->CreateH1( "3" , "Initial e+ Energy [GeV]", 100, 0, 200*MeV);
  // man->CreateH1( "4" ,"Initial Phi angle", 120,-pi,pi);
  // man->CreateH1( "5" ,"Initial Phi angle", 120,-pi,pi);
  // man->CreateH1( "6" ,"Initial Phi angle", 120,-pi,pi);
  // man->CreateH1( "7" ,"Initial Thetha angle", 60,0,0.0005);
  // man->CreateH1( "8" ,"Initial Thetha angle", 60,0,0.0005);
  // man->CreateH1( "9" ,"Initial Thetha angle", 60,0,0.0005);

  // man->CreateH1( "10", "Initial (vt) Gamma x pos", 1000, -100.0*cm,100.0*cm);
  // man->CreateH1( "11", "Initial (vt) Gamma y pos", 1000, -100.0*cm,100.0*cm);
  // man->CreateH1( "12", "Initial (vt) Gamma z pos", 1000, -100.0*cm,100.0*cm);

  // man->CreateH1( "13", "Initial (vt) e- x pos", 1000, -100.0*cm,100.0*cm);
  // man->CreateH1( "14", "Initial (vt) e- y pos", 1000, -100.0*cm,100.0*cm);
  // man->CreateH1( "15", "Initial (vt) e- z pos", 1000, -100.0*cm,100.0*cm);

  // man->CreateH1( "16", "Initial (vt) e+ x pos", 1000, -100.0*cm,100.0*cm);
  // man->CreateH1( "17", "Initial (vt) e+ y pos", 1000, -100.0*cm,100.0*cm);
  // man->CreateH1( "18", "Initial (vt) e+ z pos", 1000, -100.0*cm,100.0*cm);

  // // Final parametrs
  // man->CreateH1( "19" , "Fin Gamma Energy [GeV]", 100, 0, 200*MeV);
  // man->CreateH1( "20" , "Fin e- Energy [GeV]", 100, 0, 200*MeV);
  // man->CreateH1( "21" , "Fin e+ Energy [GeV]", 100, 0, 200*MeV);
  
  // man->CreateH1( "22" ,"Fin  Gamma Phi angle", 120,-pi,pi);
  // man->CreateH1( "23" ,"Fin e- Phi angle", 120,-pi,pi);
  // man->CreateH1( "24" ,"Fin e+ Phi angle", 120,-pi,pi);
  
  // man->CreateH1( "25" ,"Fin Gamma Thetha angle", 60,0,0.0005);
  // man->CreateH1( "26" ,"Fin e- Thetha angle", 60,0,0.0005);
  // man->CreateH1( "27" ,"Fin e+ Thetha angle", 60,0,0.0005);

  // man->CreateH1( "28", "Final (vt) Gamma x pos", 1000, -100.0*cm,100.0*cm);
  // man->CreateH1( "29", "Final (vt) Gamma y pos", 1000, -100.0*cm,100.0*cm);
  // man->CreateH1( "30", "Final (vt) Gamma z pos", 1000, -100.0*cm,100.0*cm);

  // man->CreateH1( "31", "Final (vt) e- x pos", 1000, -100.0*cm,100.0*cm);
  // man->CreateH1( "32", "Final (vt) e- y pos", 1000, -100.0*cm,100.0*cm);
  // man->CreateH1( "33", "Final (vt) e- z pos", 1000, -100.0*cm,100.0*cm);

  // man->CreateH1( "34", "Final (vt) e+ x pos", 1000, -100.0*cm,100.0*cm);
  // man->CreateH1( "35", "Final (vt) e+ y pos", 1000, -100.0*cm,100.0*cm);
  // man->CreateH1( "36", "Final (vt) e+ z pos", 1000, -100.0*cm,100.0*cm);

  // man->CreateH1( "37", "PDG number", 100, -50, 50);

    // Initial particles
  man->CreateH1( "0" , "Initial Gamma Energy [GeV]", 100, 0, 20*GeV);
  man->CreateH1( "1" , "Initial e- Energy [GeV]", 100, 0, 20*GeV);
  man->CreateH1( "2" , "Initial e+ Energy [GeV]", 100, 0, 20*GeV);

  man->CreateH1( "3" ,"Initial gamma Phi angle", 120,-pi,pi);
  man->CreateH1( "4" ,"Initial e- Phi angle", 120,-pi,pi);
  man->CreateH1( "5" ,"Initial e+ Phi angle", 120,-pi,pi);

  man->CreateH1( "6" ,"Initial gamma Thetha angle", 60,0,0.0005);
  man->CreateH1( "7" ,"Initial e- Thetha angle", 60,0,0.0005);
  man->CreateH1( "8" ,"Initial e+ Thetha angle", 60,0,0.0005);

  man->CreateH1( "9", "Initial (vt) Gamma x pos", 200, -0.1*cm,0.1*cm);
  man->CreateH1( "10", "Initial (vt) Gamma y pos", 200, -0.1*cm,0.1*cm);
  man->CreateH1( "11", "Initial (vt) Gamma z pos", 200, -0.1*cm,0.1*cm);

  man->CreateH1( "12", "Initial (vt) e- x pos", 200, -0.1*cm,0.1*cm);
  man->CreateH1( "13", "Initial (vt) e- y pos", 200, -0.1*cm,0.1*cm);
  man->CreateH1( "14", "Initial (vt) e- z pos", 200, -0.1*cm,0.1*cm);

  man->CreateH1( "15", "Initial (vt) e+ x pos", 200, -0.1*cm,0.1*cm);
  man->CreateH1( "16", "Initial (vt) e+ y pos", 200, -0.1*cm,0.1*cm);
  man->CreateH1( "17", "Initial (vt) e+ z pos", 200, -0.1*cm,0.1*cm);

  // Final parametrs
  man->CreateH1( "18" , "Fin Gamma Energy [GeV]", 100, 0, 200*MeV);
  man->CreateH1( "19" , "Fin e- Energy [GeV]", 100, 0, 200*MeV);
  man->CreateH1( "20" , "Fin e+ Energy [GeV]", 100, 0, 200*MeV);
  
  man->CreateH1( "21" ,"Fin  Gamma Phi angle", 120,-pi,pi);
  man->CreateH1( "22" ,"Fin e- Phi angle", 120,-pi,pi);
  man->CreateH1( "23" ,"Fin e+ Phi angle", 120,-pi,pi);
  
  man->CreateH1( "24" ,"Fin Gamma Thetha angle", 60,0,0.0005);
  man->CreateH1( "25" ,"Fin e- Thetha angle", 60,0,0.0005);
  man->CreateH1( "26" ,"Fin e+ Thetha angle", 60,0,0.0005);

  G4double pixel_size = 0.0050; //#cm
  G4double x_size = 30.0; //#cm
  G4double y_size = 1.0; //#cm

  G4int bin_x = (G4int) (x_size/pixel_size);
  G4int bin_y = (G4int) (y_size/pixel_size);

  man->CreateH1( "27", "Final (vt) Gamma x pos", 200, -20.0*cm,  20.0*cm);
  man->CreateH1( "28", "Final (vt) Gamma y pos", 200, -20.0*cm,  20.0*cm);
  man->CreateH1( "29", "Final (vt) Gamma z pos", 200, 300.0*cm, 310.0*cm);

  man->CreateH1( "30", "Final (vt) e- x pos", 200,   0.0*cm,   20.0*cm);
  man->CreateH1( "31", "Final (vt) e- y pos",   200,  -0.5*cm,   0.5*cm);
  man->CreateH1( "32", "Final (vt) e- z pos",   200, 300.0*cm, 310.0*cm);

  man->CreateH1( "33", "Final (vt) e+ x pos", 200,  -30.0*cm,  0.0*cm);
  man->CreateH1( "34", "Final (vt) e+ y pos", 200,   -0.5*cm,  0.5*cm);
  man->CreateH1( "35", "Final (vt) e+ z pos", 200,  300.0*cm,310.0*cm);

  man->CreateH1( "36", "PDG number", 100, -50, 50);

  man->CreateH2( "h2_0","x vs. y at tracker",bin_x, -30.0*cm, 0.0*cm, bin_y, -0.5*cm, 0.5*cm);


  /* 
  //histograms
  //        
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  if ( analysisManager->IsActive() ) {
    PrimaryGeneratorAction *gen = new PrimaryGeneratorAction(fDetector);
    if ( gen ) { 
      G4double zmax = fabs(0.5 * fDetector->GetWorldSizeZ() - gen->GetZ0()); 
   
      G4double beta = gen->GetBetaX();
      G4double range = 5.0*gen->GetSigmaX() * sqrt(1.0 + pow(zmax/beta, 2.0) ); 
      G4double drange = range*zmax/(zmax*zmax + beta*beta) + 5.0*sqrt((gen->GetEmittanceX()*beta)/(zmax*zmax + beta*beta));
      analysisManager->SetH2Activation(10, true);
      analysisManager->SetH2(10, 500, -range, range, 500, -drange, drange);  
      analysisManager->SetH1Activation(53, true);
      analysisManager->SetH1(53, 1000, -range, range);  
     
      G4cout << "x range : " << range << "   " << drange << G4endl;

      beta = gen->GetBetaY();
      range = 5.0*gen->GetSigmaY() * sqrt(1.0 + pow(zmax/beta, 2.0) ); 
      drange = range*zmax/(zmax*zmax + beta*beta) + 5.0*sqrt((gen->GetEmittanceY()*beta)/(zmax*zmax + beta*beta));
      analysisManager->SetH2Activation(11, true);
      analysisManager->SetH2(11, 500, -range, range, 500, -drange, drange);  
      analysisManager->SetH1Activation(54, true);
      analysisManager->SetH1(54, 1000, -range, range);  

      G4cout << "y range : " << range << "   " << drange << G4endl;
      
      range = zmax/c_light;
      drange = 10.0*gen->GetSigmaZ()/c_light;
      analysisManager->SetH1Activation(61, true);
      analysisManager->SetH1(61, 1000, range - drange, range + drange);  
      analysisManager->SetH1Activation(62, true);
      analysisManager->SetH1(62, 1000, range - drange, range+drange);  
      G4cout << "t range : " << range - drange << "   " << range + drange << G4endl;
    }
    delete gen;
    std::cout<<"before OpenFile\n";
    analysisManager->OpenFile();
    std::cout<<"file Opened\n";
  }
  */ 
}



void RunAction::EndOfRunAction(const G4Run*)
{  
  // print Run summary
  //
  if (isMaster) fRun->EndOfRun();    
      
  G4AnalysisManager* man = G4AnalysisManager::Instance();
  man->Write();
  man->CloseFile();

  // save histograms
  // G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();  
  // if ( analysisManager->IsActive() ) {    
  //   analysisManager->Write();
  //   analysisManager->CloseFile();
  // }  

  // show Rndm status
  if (isMaster) G4Random::showEngineStatus();
}


