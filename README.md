# Luxe_LuCaS
Simple Geant4 simulation of PCAL and FCAL for LUXE experiment
Simulation is based on LUCAS software for LumiCal calorim.

# Current software requirements:
Geant4 4.10.05
Root 6.14

# build
mkdir -p build_lucas
cd build_lucas
cmake <path-to-dir>/Luxe_LuCas
make

# Simple way to run visualization
./lxbeamsim
and than in the session terminal write:
/control/execute run_vis.mac

# to run
./lxbeamsim run_luxe.mac
