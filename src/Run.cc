/// \brief Implementation of the Run class
//

#include "Run.hh"
#include "DetectorConstruction.hh"
#include "PrimaryGeneratorAction.hh"
// #include "HistoManager.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"

#include "G4EmCalculator.hh"
#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"

#include <iomanip>


Run::Run(DetectorConstruction* det)
: G4Run(),
  fDetector(det), 
  fParticle(0), fEkin(0.), 
  fTreeCutX(1.0*mm), fTreeCutY(1.0*mm),
  fTreeParticle(0)
{
  fEnergyDeposit  = fEnergyDeposit2  = 0.;
  fTrakLenCharged = fTrakLenCharged2 = 0.;
  fTrakLenNeutral = fTrakLenNeutral2 = 0.;
  fNbStepsCharged = fNbStepsCharged2 = 0.;
  fNbStepsNeutral = fNbStepsNeutral2 = 0.;
  fMscProjecTheta = fMscProjecTheta2 = 0.;
  fMscThetaCentral = 0.;

  fNbGamma = fNbElect = fNbPosit = 0;

  fTransmit[0] = fTransmit[1] = fReflect[0] = fReflect[1] = 0;
  
  fMscEntryCentral = 0;
  
  fEnergyLeak[0] = fEnergyLeak[1] = fEnergyLeak2[0] = fEnergyLeak2[1] = 0.;
 }



Run::~Run()
{ }



void Run::SetPrimary(G4ParticleDefinition* particle, G4double energy)
{ 
  fParticle = particle;
  fEkin = energy;
  
  //compute theta0
  // fMscThetaCentral = 3*ComputeMscHighland();
}
 


void Run::Merge(const G4Run* run)
{ 
  const Run* localRun = static_cast<const Run*>(run);

  // pass information about primary particle
  fParticle = localRun->fParticle;
  fEkin     = localRun->fEkin;
  
  fMscThetaCentral = localRun->fMscThetaCentral;

  // accumulate sums
  //
  fEnergyDeposit   += localRun->fEnergyDeposit;  
  fEnergyDeposit2  += localRun->fEnergyDeposit2;  
  fTrakLenCharged  += localRun->fTrakLenCharged;
  fTrakLenCharged2 += localRun->fTrakLenCharged2;   
  fTrakLenNeutral  += localRun->fTrakLenNeutral;  
  fTrakLenNeutral2 += localRun->fTrakLenNeutral2;
  fNbStepsCharged  += localRun->fNbStepsCharged;
  fNbStepsCharged2 += localRun->fNbStepsCharged2;
  fNbStepsNeutral  += localRun->fNbStepsNeutral;
  fNbStepsNeutral2 += localRun->fNbStepsNeutral2;
  fMscProjecTheta  += localRun->fMscProjecTheta;
  fMscProjecTheta2 += localRun->fMscProjecTheta2;

    
  fNbGamma += localRun->fNbGamma;
  fNbElect += localRun->fNbElect;      
  fNbPosit += localRun->fNbPosit;
  
  fTransmit[0] += localRun->fTransmit[0];  
  fTransmit[1] += localRun->fTransmit[1];
  fReflect[0]  += localRun->fReflect[0];
  fReflect[1]  += localRun->fReflect[1];
  
  fMscEntryCentral += localRun->fMscEntryCentral;
  
  fEnergyLeak[0]  += localRun->fEnergyLeak[0];
  fEnergyLeak[1]  += localRun->fEnergyLeak[1];
  fEnergyLeak2[0] += localRun->fEnergyLeak2[0];
  fEnergyLeak2[1] += localRun->fEnergyLeak2[1];
                  
  G4Run::Merge(run); 
} 



void Run::EndOfRun()
{
  // std::cout<<"@@@@@@@@@@@@@@@@@@ AND OF RUN\n";
  // compute mean and rms
  //
 
  /*
  G4int TotNbofEvents = numberOfEvent;
  if (TotNbofEvents == 0) return;
  
  G4double EnergyBalance = fEnergyDeposit + fEnergyLeak[0] + fEnergyLeak[1];
  EnergyBalance /= TotNbofEvents;

  fEnergyDeposit /= TotNbofEvents; fEnergyDeposit2 /= TotNbofEvents;
  G4double rmsEdep = fEnergyDeposit2 - fEnergyDeposit*fEnergyDeposit;
  if (rmsEdep>0.) rmsEdep = std::sqrt(rmsEdep/TotNbofEvents);
  else            rmsEdep = 0.;

  fTrakLenCharged /= TotNbofEvents; fTrakLenCharged2 /= TotNbofEvents;
  G4double rmsTLCh = fTrakLenCharged2 - fTrakLenCharged*fTrakLenCharged;
  if (rmsTLCh>0.) rmsTLCh = std::sqrt(rmsTLCh/TotNbofEvents);
  else            rmsTLCh = 0.;

  fTrakLenNeutral /= TotNbofEvents; fTrakLenNeutral2 /= TotNbofEvents;
  G4double rmsTLNe = fTrakLenNeutral2 - fTrakLenNeutral*fTrakLenNeutral;
  if (rmsTLNe>0.) rmsTLNe = std::sqrt(rmsTLNe/TotNbofEvents);
  else            rmsTLNe = 0.;

  fNbStepsCharged /= TotNbofEvents; fNbStepsCharged2 /= TotNbofEvents;
  G4double rmsStCh = fNbStepsCharged2 - fNbStepsCharged*fNbStepsCharged;
  if (rmsStCh>0.) rmsStCh = std::sqrt(rmsTLCh/TotNbofEvents);
  else            rmsStCh = 0.;

  fNbStepsNeutral /= TotNbofEvents; fNbStepsNeutral2 /= TotNbofEvents;
  G4double rmsStNe = fNbStepsNeutral2 - fNbStepsNeutral*fNbStepsNeutral;
  if (rmsStNe>0.) rmsStNe = std::sqrt(rmsTLCh/TotNbofEvents);
  else            rmsStNe = 0.;


  G4double Gamma = (G4double)fNbGamma/TotNbofEvents;
  G4double Elect = (G4double)fNbElect/TotNbofEvents;
  G4double Posit = (G4double)fNbPosit/TotNbofEvents;

  G4double transmit[2];
  transmit[0] = 100.*fTransmit[0]/TotNbofEvents;
  transmit[1] = 100.*fTransmit[1]/TotNbofEvents;

  G4double reflect[2];
  reflect[0] = 100.*fReflect[0]/TotNbofEvents;
  reflect[1] = 100.*fReflect[1]/TotNbofEvents;


  G4double rmsMsc = 0., tailMsc = 0.;
  if (fMscEntryCentral > 0) {
    fMscProjecTheta /= fMscEntryCentral; fMscProjecTheta2 /= fMscEntryCentral;
    rmsMsc = fMscProjecTheta2 - fMscProjecTheta*fMscProjecTheta;
    if (rmsMsc > 0.) { rmsMsc = std::sqrt(rmsMsc); }
    if(fTransmit[1] > 0.0) {
      tailMsc = 100.- (100.*fMscEntryCentral)/(2*fTransmit[1]);
    }    
  }
  

  fEnergyLeak[0] /= TotNbofEvents; fEnergyLeak2[0] /= TotNbofEvents;
  G4double rmsEl0 = fEnergyLeak2[0] - fEnergyLeak[0]*fEnergyLeak[0];
  if (rmsEl0>0.) rmsEl0 = std::sqrt(rmsEl0/TotNbofEvents);
  else           rmsEl0 = 0.;
  
  fEnergyLeak[1] /= TotNbofEvents; fEnergyLeak2[1] /= TotNbofEvents;
  G4double rmsEl1 = fEnergyLeak2[1] - fEnergyLeak[1]*fEnergyLeak[1];
  if (rmsEl1>0.) rmsEl1 = std::sqrt(rmsEl1/TotNbofEvents);
  else           rmsEl1 = 0.;    
  
  //Stopping Power from input Table.
  //
  G4Material* material = fDetector->GetAbsorberMaterial();
  G4double length      = fDetector->GetAbsorberThickness();
  G4double density     = material->GetDensity();   
  G4String partName    = fParticle->GetParticleName();

 
  G4EmCalculator emCalculator;
  G4double dEdxTable = 0., dEdxFull = 0.;
  if (fParticle->GetPDGCharge()!= 0.) { 
    dEdxTable = emCalculator.GetDEDX(fEkin,fParticle,material);
    dEdxFull  = emCalculator.ComputeTotalDEDX(fEkin,fParticle,material);    
  }
  G4double stopTable = dEdxTable/density;
  G4double stopFull  = dEdxFull /density; 
   
  //Stopping Power from simulation.
  //    
  G4double meandEdx  = fEnergyDeposit/length;
  G4double stopPower = meandEdx/density;  
 

  G4cout << "\n ======================== run summary ======================\n";

  G4int prec = G4cout.precision(3);
  
  G4cout << "\n The run was " << TotNbofEvents << " " << partName << " of "
         << G4BestUnit(fEkin,"Energy") << " through " 
         << G4BestUnit(length,"Length") << " of "
         << material->GetName() << " (density: " 
         << G4BestUnit(density,"Volumic Mass") << ")" << G4endl;
  
  G4cout.precision(4);
  
  G4cout << "\n Total energy deposit in absorber per event = "
         << G4BestUnit(fEnergyDeposit,"Energy") << " +- "
         << G4BestUnit(rmsEdep,      "Energy") 
         << G4endl;
         
  G4cout << "\n -----> Mean dE/dx = " << meandEdx/(MeV/cm) << " MeV/cm"
         << "\t(" << stopPower/(MeV*cm2/g) << " MeV*cm2/g)"
         << G4endl;
         
  G4cout << "\n From formulas :" << G4endl; 
  G4cout << "   restricted dEdx = " << dEdxTable/(MeV/cm) << " MeV/cm"
         << "\t(" << stopTable/(MeV*cm2/g) << " MeV*cm2/g)"
         << G4endl;
         
  G4cout << "   full dEdx       = " << dEdxFull/(MeV/cm) << " MeV/cm"
         << "\t(" << stopFull/(MeV*cm2/g) << " MeV*cm2/g)"
         << G4endl;
         
  G4cout << "\n Leakage :  primary = "
         << G4BestUnit(fEnergyLeak[0],"Energy") << " +- "
         << G4BestUnit(rmsEl0,       "Energy")
         << "   secondaries = "
         << G4BestUnit(fEnergyLeak[1],"Energy") << " +- "
         << G4BestUnit(rmsEl1,       "Energy")          
         << G4endl;
         
  G4cout << " Energy balance :  edep + eleak = "
         << G4BestUnit(EnergyBalance,"Energy")
         << G4endl;         
                           
  G4cout << "\n Total track length (charged) in absorber per event = "
         << G4BestUnit(fTrakLenCharged,"Length") << " +- "
         << G4BestUnit(rmsTLCh,       "Length") << G4endl;

  G4cout << " Total track length (neutral) in absorber per event = "
         << G4BestUnit(fTrakLenNeutral,"Length") << " +- "
         << G4BestUnit(rmsTLNe,       "Length") << G4endl;

  G4cout << "\n Number of steps (charged) in absorber per event = "
         << fNbStepsCharged << " +- " << rmsStCh << G4endl;

  G4cout << " Number of steps (neutral) in absorber per event = "
         << fNbStepsNeutral << " +- " << rmsStNe << G4endl;

  G4cout << "\n Number of secondaries per event : Gammas = " << Gamma
         << ";   electrons = " << Elect
           << ";   positrons = " << Posit << G4endl;

  G4cout << "\n Number of events with the primary particle transmitted = "
         << transmit[1] << " %" << G4endl;

  G4cout << " Number of events with at least  1 particle transmitted "
         << "(same charge as primary) = " << transmit[0] << " %" << G4endl;

  G4cout << "\n Number of events with the primary particle reflected = "
         << reflect[1] << " %" << G4endl;

  G4cout << " Number of events with at least  1 particle reflected "
         << "(same charge as primary) = " << reflect[0] << " %" << G4endl;

  // compute width of the Gaussian central part of the MultipleScattering
  //
  G4cout << "\n MultipleScattering:" 
         << "\n  rms proj angle of transmit primary particle = "
         << rmsMsc/mrad << " mrad (central part only)" << G4endl;

  G4cout << "  computed theta0 (Highland formula)          = "
         << ComputeMscHighland()/mrad << " mrad" << G4endl;
           
  G4cout << "  central part defined as +- "
         << fMscThetaCentral/mrad << " mrad; " 
         << "  Tail ratio = " << tailMsc << " %" << G4endl;
         
  /*
  // normalize histograms
  //
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  
  G4int ih = 1;
  G4double binWidth = analysisManager->GetH1Width(ih);
  G4double unit     = analysisManager->GetH1Unit(ih);  
  G4double fac = unit/(TotNbofEvents*binWidth);
  analysisManager->ScaleH1(ih,fac);

  ih = 10;
  binWidth = analysisManager->GetH1Width(ih);
  unit     = analysisManager->GetH1Unit(ih);  
  fac = unit/(TotNbofEvents*binWidth);
  analysisManager->ScaleH1(ih,fac);

//< ob
  ih = 20;
  binWidth = analysisManager->GetH1Width(ih);
  unit     = analysisManager->GetH1Unit(ih);  
  fac = unit/(TotNbofEvents*binWidth);
  analysisManager->ScaleH1(ih,fac);
//> ob
  
  ih = 12;
  analysisManager->ScaleH1(ih,1./TotNbofEvents);
                    
//< ob
  ih = 22;
  analysisManager->ScaleH1(ih,1./TotNbofEvents);
//> ob

  // reset default precision
  G4cout.precision(prec);
  */
  
}   



G4double Run::ComputeMscHighland()
{
  // //compute the width of the Gaussian central part of the MultipleScattering
  // //projected angular distribution.
  // //Eur. Phys. Jour. C15 (2000) page 166, formule 23.9

  // G4double t = (fDetector->GetAbsorberThickness())
  //             /(fDetector->GetAbsorberMaterial()->GetRadlen());
  // if (t < DBL_MIN) return 0.;

  // G4double T = fEkin;
  // G4double M = fParticle->GetPDGMass();
  // G4double z = std::abs(fParticle->GetPDGCharge()/eplus);

  // G4double bpc = T*(T+2*M)/(T+M);
  // G4double teta0 = 13.6*MeV*z*std::sqrt(t)*(1.+0.038*std::log(t))/bpc;
  // return teta0;
}

/*
void Run::SetHistoManager (HistoManager *hmng) 
{ 
  fHistoManager = hmng; 
  if (fHistoManager) {
    fTreeCutX = fHistoManager->GetTreeCutX();
    fTreeCutY = fHistoManager->GetTreeCutY();

    G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
    fTreeParticle = particleTable->FindParticle(fHistoManager->GetTreeParticle());
    
    std::cout << "TTree X cut: " << G4BestUnit(fTreeCutX, "Length") << std::endl; 
    std::cout << "TTree Y cut: " << G4BestUnit(fTreeCutY, "Length") << std::endl;
    std::cout << "TTree particles: " << fHistoManager->GetTreeParticle() << std::endl;
    if (fTreeParticle) std::cout << "   it was found in G4ParticleTable with PDG: " << fTreeParticle->GetPDGEncoding() << std::endl;
    else std::cout << "   it was not found in G4ParticleTable, so all particles will be saved to the TTree\n";
 }
}
*/

