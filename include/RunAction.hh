//
/// \brief Definition of the RunAction class
//

#ifndef RunAction_h
#define RunAction_h 1

#include "G4UserRunAction.hh"
#include "globals.hh"



class Run;
class DetectorConstruction;
class PrimaryGeneratorAction;
class HistoManager;



class RunAction : public G4UserRunAction
{

public:

    RunAction(DetectorConstruction* det, PrimaryGeneratorAction* prim=0);
   ~RunAction();
   
    virtual G4Run* GenerateRun();    
    virtual void BeginOfRunAction(const G4Run*);
    virtual void   EndOfRunAction(const G4Run*);

  private:
    DetectorConstruction*   fDetector;
    PrimaryGeneratorAction* fPrimary;
    Run*                    fRun;        
    HistoManager*           fHistoManager;

    G4String  fFileName;
};



#endif

