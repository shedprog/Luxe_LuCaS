
# LIST_FILE=/data/Projects_physics/LUXE/IPstrong/Lists/bppp_17.5GeV_5mfoiltoIP_Enelas_0.2J.out
# LIST_FILE=/data/Projects_physics/LUXE/IPstrong/Lists/bppp_17.5GeV_5mfoiltoIP_Enelas_0.35J.out
# LIST_FILE=/data/Projects_physics/LUXE/IPstrong/Lists/bppp_17.5GeV_5mfoiltoIP_Enelas_0.5J.out
# LIST_FILE=/data/Projects_physics/LUXE/IPstrong/Lists/bppp_17.5GeV_5mfoiltoIP_Enelas_0.7J.out
LIST_FILE=/data/Projects_physics/LUXE/IPstrong/Lists/bppp_17.5GeV_5mfoiltoIP_Enelas_0.85J.out
# LIST_FILE=/data/Projects_physics/LUXE/IPstrong/Lists/bppp_17.5GeV_5mfoiltoIP_Enelas_1.0J.out
LIST=(`cat $LIST_FILE`)

number=$1
for n in {1..$number}
do
  file=${LIST[$(($n))]}
  name=$(basename $file)
  echo "File: " $file
  echo "File Name: " $name
  sed "s|/lxphoton/gun/MCParticlesFile.*|/lxphoton/gun/MCParticlesFile $file|g;\
       s|/lxphoton/gun/beamType.*|/lxphoton/gun/beamType mc|g;\
       s|/analysis/filename.*|/analysis/filename ${name}.root|g;\
       s|/gun/energy.*|#/gun/energy $E GeV|g;\
       s|/run/beamOn.*|/run/beamOn 1000000000|g\
       " run_luxe.mac > run_luxe_update.mac
  ./lxbeamsim run_luxe_update.mac

done
