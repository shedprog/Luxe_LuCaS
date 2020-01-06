
# LIST_FILE=/data/Projects_physics/LUXE/IPstrong/Lists/bppp_17.5GeV_5mfoiltoIP_Enelas_0.2J.out
# LIST_FILE=/data/Projects_physics/LUXE/IPstrong/Lists/bppp_17.5GeV_5mfoiltoIP_Enelas_0.35J.out
# LIST_FILE=/data/Projects_physics/LUXE/IPstrong/Lists/bppp_17.5GeV_5mfoiltoIP_Enelas_0.5J.out
# LIST_FILE=/data/Projects_physics/LUXE/IPstrong/Lists/bppp_17.5GeV_5mfoiltoIP_Enelas_0.7J.out
#LIST_FILE=/data/Projects_physics/LUXE/IPstrong/Lists/bppp_17.5GeV_5mfoiltoIP_Enelas_0.85J.out
# LIST_FILE=/data/Projects_physics/LUXE/IPstrong/Lists/bppp_17.5GeV_5mfoiltoIP_Enelas_1.0J.out
#LIST=(`cat $LIST_FILE`)

E_LIST=(1.0 2.0 3.0 4.0 5.0 6.0 7.0 8.0 9.0 10.0 11.0 12.0 13.0 14.0 15.0)

number=15
# for n in {1..$number}
# do
#   file=${LIST[$(($n))]}
#   name=$(basename $file)
#   echo "File: " $file
#   echo "File Name: " $name
#   sed "s|/lxphoton/gun/MCParticlesFile.*|/lxphoton/gun/MCParticlesFile $file|g;\
#        s|/lxphoton/gun/beamType.*|/lxphoton/gun/beamType mc|g;\
#        s|/analysis/filename.*|/analysis/filename ${name}.root|g;\
#        s|/gun/energy.*|#/gun/energy $E GeV|g;\
#        s|/run/beamOn.*|/run/beamOn 1000000000|g\
#        " run_luxe.mac > run_luxe_update.mac
#   ./lxbeamsim run_luxe_update.mac
#
# done
for n in {1..$number}
do
  # file=${LIST[$(($n))]}
  # name=$(basename $file)
  echo "File: " $file
  echo "File Name: " $name
  E=${E_LIST[$(($n))]}
  sed "s|/lxphoton/gun/MCParticlesFile.*|#/lxphoton/gun/MCParticlesFile|g;\
       s|/lxphoton/gun/beamType.*|/lxphoton/gun/beamType mono|g;\
       s|/analysis/filename.*|/analysis/filename ${E}GeV.root|g;\
       s|/gun/energy.*|/gun/energy $E GeV|g;\
       s|/run/beamOn.*|/run/beamOn 5000|g\
       " run_luxe.mac > run_luxe_update.mac
  ./lxbeamsim run_luxe_update.mac

done
