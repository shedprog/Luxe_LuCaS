//
/// \brief Implementation of the PrimaryGeneratorAction class
//


#include "PrimaryGeneratorAction.hh"

#include "DetectorConstruction.hh"
#include "PrimaryGeneratorMessenger.hh"

#include "G4Event.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include "G4UnitsTable.hh"
#include "PrimarySpectra.hh"
#include "LuxeTestGenerator.hh"
#include "G4RunManager.hh"
#include "G4String.hh"
#include "GlobalVars.hh"
#include <string>
#include <cmath>

#include <iostream>
//#include <filesystem>
//namespace fs = std::filesystem;
#include <libgen.h>
#include "Setup.hh"
// std::string root_file_name;
double weight_fromMC;


PrimaryGeneratorAction::PrimaryGeneratorAction(DetectorConstruction* DC)
:G4VUserPrimaryGeneratorAction(),
 fParticleGun(0),fDetector(DC),fBeamType(beamGauss),
 fx0(0.0), fy0(0.0), fz0(0.0), fsigmax(5.0*um), fsigmay(5.0*um), fsigmaz(24.0*um),
 femittancex(0.0), femittancey(0.0), fbetax(0.0), fbetay(0.0), fSpectra(0), lxgen(0),
 fMCfile(""), fMCList(false), fnfixparticles(0), fGunMessenger(0)

{
  G4int n_particle = 1;
  fParticleGun  = new G4ParticleGun(n_particle);
  SetDefaultKinematic();

  //create a messenger for this class
  fGunMessenger = new PrimaryGeneratorMessenger(this);
}



PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
  delete fParticleGun;
  delete fGunMessenger;
  if (fSpectra) delete fSpectra;
  if (lxgen) delete lxgen;
}



void PrimaryGeneratorAction::SetDefaultKinematic()
{
  // default particle kinematic
  // defined by the beam parameters: emittancex, emittancey, fsigmax, fsigmay, fsigmaz and drift distance to IP fz0.
  //
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4ParticleDefinition* particle = particleTable->FindParticle("e-");
  fParticleGun->SetParticleDefinition(particle);
  fParticleGun->SetParticleEnergy(17.5*GeV);

  fz0 = -20.0*cm;

  fsigmax = 5.0*um;
  fsigmay = 5.0*um;
  fsigmaz = 24.0*um;

  G4double lf = fParticleGun->GetParticleEnergy() / particle->GetPDGMass();

  femittancex = 1.4e-3 * mm / lf;
  femittancey = 1.4e-3 * mm / lf;

  fbetax = fsigmax*fsigmax/femittancex;
  fbetay = fsigmay*fsigmay/femittancey;

  fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0.,0.,1.));
  fParticleGun->SetParticlePosition(G4ThreeVector(fx0, fy0, fz0));
}


void PrimaryGeneratorAction::SetBeamType(G4String val)
{
  if (val == "gaussian") {
    fBeamType = beamGauss;
  } else if (val == "mono") {
    fBeamType = beamMono;
  } else if (val == "mc") {
    fBeamType = beamMC;
  } else if (val == "monolim") {
    fBeamType = beamMonoLimits;
  } else {
    G4cout << "PrimaryGeneratorAction::SetBeamType: <" << val << ">"
           << " is not defined"
           << G4endl;
  }
}



void PrimaryGeneratorAction::SetSigmaX(const G4double sigma)
{
  fsigmax = sigma;
  fbetax = fsigmax*fsigmax/femittancex;
  std::cout << "PrimaryGeneratorAction: Set beam sigmaX at IP to : " << G4BestUnit(fsigmax, "Length") << std::endl;
//G4cout does not print anything...   G4cout << "Set beam sigmaX at IP to : " << G4BestUnit(fsigmax, "Length") << G4endl;
}



void PrimaryGeneratorAction::SetSigmaY(const G4double sigma)
{
  fsigmay = sigma;
  fbetay = fsigmay*fsigmay/femittancey;
  std::cout << "PrimaryGeneratorAction: Set beam sigmaY at IP to : " << G4BestUnit(fsigmay, "Length") << std::endl;
//   G4cout << "Set beam sigmaY at IP to : " << G4BestUnit(fsigmay, "Length") << G4endl;
}

void PrimaryGeneratorAction::SetSpectraFile(G4String val)
{
  if (fSpectra) delete fSpectra;
  fSpectra = new PrimarySpectra();
//  fSpectra->SetVerbose(1);
  int ndata = fSpectra->LoadData(val);
  fSpectra->SetScale(GeV);
  std::cout << "PrimaryGeneratorAction: Primary spectra with " << ndata << " data points was loaded from the file " << val << std::endl;
}


void PrimaryGeneratorAction::SetMCParticleFile(G4String val, const G4bool list)
{
  if (lxgen) { delete lxgen;  lxgen = 0; }
  fMCfile = val;
  fMCList = list;
  if (fMCList) G4cout << "File with the list of files with primary MC particles" << fMCfile << G4endl;
  else         G4cout << "File with primary MC particles" << fMCfile << G4endl;

  // char *cstr = new char[val.length() + 1];
  // strcpy(cstr, val.c_str());
  // root_file_name = basename(cstr);
  // delete [] cstr;
  // std::cout<<"root file name: "<<root_file_name<<"\n";

}


void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  // std::cout<<"@@@@@@@@@@@@@@@@@@@@ GeneratePrimaries\n";
  // this function is called at the begining of event
  if (fSpectra) fParticleGun->SetParticleEnergy(fSpectra->GetRandom());
  // fParticleGun->SetParticleEnergy(fSpectra->GetRandom());

  if (fBeamType == beamGauss) {
    GenerateGaussian(anEvent);
  } else if (fBeamType == beamMono) {
    GenerateMono(anEvent);
  }
    else if (fBeamType == beamMC) {
    GeneratefromMC(anEvent);
  }
    else if (fBeamType == beamMonoLimits) {
    GenerateMonoLimits(anEvent);
    }

  // std::cout<<"@@@@@@@@@@@@@@@@@@@@@@ GenerateMono \n";
}



void PrimaryGeneratorAction::GenerateGaussian(G4Event* anEvent)
{
  /*
  // it generate particle according to defined emittance drifting from fz0 towards an IP
  // where the beam size is fsigmax, fsigmay, fsigmaz.
  G4double zshift = 0.0*cm;
  // if (zshift > fDetector->GetzstartAbs() - 100.0*fsigmaz) {
  //   zshift = 0.5 * (fDetector->GetzstartAbs() - 0.5 * fDetector->GetWorldSizeZ());
  //   if (zshift > fDetector->GetzstartAbs() - 100.0*fsigmaz) {
  //     G4String msgstr("Error setting initial position!\n");
  //     G4Exception("PrimaryGeneratorAction::", "GenerateGaussian(Event)", FatalException, msgstr.c_str());
  //   }
  // }

  fz0 = zshift;  // in SetDefaultKinematic fDetector is before update from the config file.

  G4double z0 = G4RandGauss::shoot(fz0, fsigmaz);
  z0 -= 0.5 * fDetector->GetWorldSizeZ();  // This is needed to have correct drift distance for x, y distribution.

  G4double sigmax = fsigmax * sqrt(1.0 + pow(z0/fbetax, 2.0));
  G4double x0 = G4RandGauss::shoot(fx0, sigmax);
  G4double meandx = x0*z0 / (z0*z0 + fbetax*fbetax);
  G4double sigmadx = sqrt( femittancex*fbetax / (z0*z0 + fbetax*fbetax) );
  G4double dx0 = G4RandGauss::shoot(meandx, sigmadx);

  G4double sigmay = fsigmay * sqrt(1.0 + pow(z0/fbetay, 2.0));
  G4double y0 = G4RandGauss::shoot(fy0, sigmay);
  G4double meandy = y0*z0 / (z0*z0 + fbetay*fbetay);
  G4double sigmady = sqrt( femittancey*fbetay / (z0*z0 + fbetay*fbetay) );
  G4double dy0 = G4RandGauss::shoot(meandy, sigmady);

  G4double mass = fParticleGun->GetParticleDefinition()->GetPDGMass();
  G4double E = fParticleGun->GetParticleEnergy();

  G4double pz = sqrt( (E*E - mass*mass)/ (dx0*dx0 + dy0*dy0 + 1.0) );

  fParticleGun->SetParticlePosition(G4ThreeVector(x0, y0, z0 + 0.5 * fDetector->GetWorldSizeZ() ));
  fParticleGun->SetParticleMomentumDirection(G4ThreeVector(dx0*pz, dy0*pz, pz));
  fParticleGun->GeneratePrimaryVertex(anEvent);*/

  G4double x0=0.0;
  G4double y0=0.0;
  G4double z0=0.0;

  G4double mass = fParticleGun->GetParticleDefinition()->GetPDGMass();
  // G4double E = fParticleGun->GetParticleEnergy();
  G4double E = G4RandFlat::shoot(4000.,13000);

  // G4double pz = sqrt( (E*E - mass*mass) );

  fParticleGun->SetParticlePosition(G4ThreeVector(x0, y0, z0));
  fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0, 0, 1));
  fParticleGun-> SetParticleEnergy(E);
  fParticleGun->GeneratePrimaryVertex(anEvent);

}



void PrimaryGeneratorAction::GenerateMono(G4Event* anEvent)
{

  // fz0 = 0.0;
  // fy0 = 0.0*cm + G4RandFlat::shoot(-Setup::pix_y_size/2.0, Setup::pix_y_size/2.0);
  // fx0 = 370.0*mm + G4RandFlat::shoot(-Setup::pix_x_size/2.0, Setup::pix_x_size/2.0);

  // fy0 = 0.0*mm;
  // fx0 = 370.0*mm;
  // fz0 = 0.0*mm;

  G4double mass = fParticleGun->GetParticleDefinition()->GetPDGMass();
  // G4double E = fParticleGun->GetParticleEnergy();
  // G4ThreeVector position = fParticleGun->GetParticlePosition();
  G4ThreeVector position = G4ThreeVector(480.*mm,0.0,0.0);
  // G4double E = G4RandFlat::shoot(1.*GeV,13.*GeV);
  // G4double pz = sqrt(E*E - mass*mass);

  G4double E = 0.714011463377/std::sinh(5.22758393475e-07*position.x());
  // G4double E = (269.0/4.0-position.x()*23.0/200.0)*GeV;
  G4double pz = sqrt(E*E - mass*mass);

  std::cout << "pos: "<< position << " E: " << E <<" "<<"pz: "<<pz<<"\n";

  fParticleGun->SetParticleEnergy(E);
  fParticleGun->SetParticlePosition(position);
  fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0.0, 0.0, pz));
  fParticleGun->GeneratePrimaryVertex(anEvent);

}

void PrimaryGeneratorAction::GenerateMonoLimits(G4Event* anEvent)
{
  G4double mass = fParticleGun->GetParticleDefinition()->GetPDGMass();
  G4double E = fParticleGun->GetParticleEnergy();
  G4double pz = sqrt(E*E - mass*mass);

  G4double x_ = G4RandFlat::shoot(monolim_xmin, monolim_xmax)*cm;
  G4double y_ = G4RandFlat::shoot(monolim_ymin, monolim_ymax)*cm;
  G4double z_ = 0.0;

  std::cout<<"pos: "<<x_<<" "<<y_<<" E: "<<E<<" "<<"pz: "<<pz<<"\n";

  fParticleGun->SetParticleEnergy(E);
  fParticleGun->SetParticlePosition(G4ThreeVector(x_,y_,z_));
  fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0.0, 0.0, pz));
  fParticleGun->GeneratePrimaryVertex(anEvent);
}


void PrimaryGeneratorAction::GeneratefromMC(G4Event* anEvent)
{


  if (!lxgen) {
    lxgen = new LuxeTestGenerator();
    if (fMCList)  lxgen->SetFileList(fMCfile);
    else lxgen->AddEventFile(fMCfile);
    std::cout << "Processing file: " << fMCfile << std::endl;

    lxgen->SetFileType("out", 9, 9);
    std::cout << "File type is set "  << std::endl;

  }

  if (fnfixparticles > 0) {
    std::cout << "Generating wrong\n";
    --fnfixparticles;
    fParticleGun->GeneratePrimaryVertex(anEvent);
    return;
  }

  std::vector < std::vector <double> > ptcls;
  int nscat = lxgen->GetEventFromFile(ptcls);
  // std::cout << "GetEventFromFile "  << nscat << std::endl;
  if (nscat > 1) {
      G4String msgstr("Error reading particle from a file! More than one were read, it is not supported!\n");
      G4Exception("PrimaryGeneratorAction::", "GeneratefromMC(Event)", FatalException, msgstr.c_str());
  }
  if (nscat <= 0) {
    delete lxgen;
    lxgen = 0;
    if (nscat == -1) {
      G4cout << "All particle from the file " << fMCfile << " are processed." << G4endl;
      G4RunManager::GetRunManager()->AbortRun();
    } else {
      G4String msgstr("Error reading MC particle from the file!\n");
      G4Exception("PrimaryGeneratorAction::", "GeneratefromMC(Event)", FatalException, msgstr.c_str());
    }
  }
// consider weight
// vertex units

  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4ParticleDefinition* particle = 0;
  // std::cout << nscat << "\n";
  for (int ip = 0; ip < nscat; ++ip) {
    // std::cout<<ip<<"\n";
    // if (ip%1000==0 and ip!=0) std::cout << "processed: " << ip << "\n";
    std::vector <double> pdata = ptcls[ip];
    std::vector <double> pp(4);
    int pid = static_cast<int>(pdata[7]);

    double wght = 1.0;
    double vtxtomm = 1.0;
    double  vtx[3];

//    int hindx = 1;1000000

    if (pid == 11)  { particle = particleTable->FindParticle("e-"); }
    else if (pid == -11) { particle = particleTable->FindParticle("e+"); }
    else if (pid == 22)  { particle = particleTable->FindParticle("gamma"); }
    // else if (pid == 22)  { return; }
    else {
      G4String msgstr("Error setting initial particle from a file!\n");
      G4Exception("PrimaryGeneratorAction::", "GeneratefromMC(Event)", FatalException, msgstr.c_str());
    }

//     nfpart[hindx] += wght*pdata[8];
    if (pdata[8] >= 1)  wght = wght*pdata[8];
    //extern double weight_fromMC;
    weight_fromMC = pdata[8];
    //weight_fromMC = 1.0;

    std::copy(pdata.begin()+4, pdata.begin()+7, pp.begin());
    std::copy(pdata.begin()+1, pdata.begin()+4, vtx);
    std::for_each(vtx, vtx+3, [=](double &x){x *= um;});
    //vtx[2] -= 7.0*m;
    std::for_each(vtx,vtx+3, [=](double &x){x*=vtxtomm;});
    pp[3] = pdata[0];

    fParticleGun->SetParticleDefinition(particle);
    fParticleGun->SetParticleEnergy(pp[3]*GeV);
    fParticleGun->SetParticleMomentumDirection(G4ThreeVector(pp[0], pp[1], pp[2]));
    fParticleGun->SetParticlePosition(G4ThreeVector(vtx[0], vtx[1], vtx[2]));

    //if (pid == 22) {
     // if (TestHitTarget(pp, vtx) < 1.2) {
     //   fParticleGun->SetNumberOfParticles(1);
     //   fnfixparticles = static_cast<int>(wght/100.0)-1;
     //   std::cout << "Generating " << static_cast<int>(wght) << " particle with following mC data:\n";
     // } else {
     //   std::cout << "Particle does not hit the target, " << static_cast<int>(wght) << " particles with following mC data:\n";
     //   fParticleGun->SetNumberOfParticles(0);
     //   fnfixparticles = 0;
     // }
     // std::for_each(pdata.begin(), pdata.end(), [](const double x){std::cout << x << "  ";});
     // std::cout << std::endl;
     // std::cout << "TestHitTarget: " << TestHitTarget(pp, vtx) << std::endl;
    //} else {
     // fParticleGun->SetNumberOfParticles(0);
      //fnfixparticles = 0;
    //}
//     fParticleGun->SetNumberOfParticles(static_cast<int>(wght));
  //}


  // if (pid==11 or pid==-11){
  //   //fnfixparticles = static_cast<int>(wght/100.0)-1;
  //   fParticleGun->SetNumberOfParticles(10);
  //   fnfixparticles = 10;
  // }
  // else
  // {
    fParticleGun->SetNumberOfParticles(1);
    fnfixparticles = 0;
  // }

  std::cout << "Generating " <<pid<<"\n";
  fParticleGun->GeneratePrimaryVertex(anEvent);
  std::cout << "Generated particle: " << pid << "\n";
  //std::for_each(pdata.begin(), pdata.end(), [](const double x){std::cout << x << "  ";});
  //std::cout << std::endl;
}
}
