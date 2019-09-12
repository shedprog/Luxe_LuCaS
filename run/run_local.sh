Energies=(2.0 3.0 4.0 5.0 6.0 7.0 8.0 9.0 10.0 11.0 12.0)
for E in $Energies
do
  echo "Energy:" $E
  sed "s|/lxphoton/gun/MCParticlesFile|#/lxphoton/gun/MCParticlesFile|g;\
       s|/lxphoton/gun/beamType.*|/lxphoton/gun/beamType mono|g;\
       s|/analysis/filename.*|/analysis/filename mono_${E}_GeV.root|g;\
       s|/gun/energy.*|/gun/energy $E GeV|g;\
       s|/run/beamOn.*|/run/beamOn 5000|g\
       " run_luxe.mac > run_luxe_update.mac
  ./lxbeamsim run_luxe_update.mac

done
