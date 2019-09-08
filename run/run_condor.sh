#!/bin/zsh

#Setting env
source /cvmfs/sft.cern.ch/lcg/views/setupViews.sh LCG_91 x86_64-slc6-gcc7-opt
WORKDIR=

option=mono
#option=true

if [ $option == 'true' ]
then

  LIST_FILE=/nfs/dust/zeus/group/mykytaua/LUXE/IPstrong/Lists/bppp_17.5GeV_5mfoiltoIP_Enelas_1.0J.out
  LIST=(`cat $LIST_FILE`)

  echo $TMP
  cd $TMP
  pwd

  echo "File for DATA: " ${LIST[$(($1))]}

  file=${LIST[$(($1))]}

  echo $file
  sed "s|/lxphoton/gun/MCParticlesFile.*|/lxphoton/gun/MCParticlesFile $file|g\
       " ${WORKDIR}/run_luxe.mac > ${TMP}/run_luxe.mac

  # mkdir -p ${WORKDIR}/rootOUT
  root_n=$((($1-1)/100))

  mkdir -p ${WORKDIR}/output/out_${root_n}

  cd ${WORKDIR}/output/out_${root_n}

  ${WORKDIR}/lxbeamsim ${TMP}/run_luxe.mac
  echo "Run true MC done!"

elif [ $option == 'mono' ]
then
  #statements

  N_NOW=$1
  NUMBER_OF_STEPS=10
  E_START=5.0 #GeV
  E_END=10.0 #GeV

  MC_NUMBER=1000

  echo "State: " $E_START $N_NOW $E_END $E_START $NUMBER_OF_STEPS 
  #E_current=$(( $E_START + $N * ($E_END-$E_START) / ($NUMBER_OF_STEPS-2) ))
  E_current=$( echo $E_START $N_NOW $E_END $E_START $NUMBER_OF_STEPS | awk '{print $1+$2*($3-$4)/($5-1)}')

  echo $TMP
  cd $TMP
  pwd

  sed "s|/lxphoton/gun/MCParticlesFile|#/lxphoton/gun/MCParticlesFile|g;\
       s|/lxphoton/gun/beamType.*|/lxphoton/gun/beamType mono|g;\
       s|/analysis/filename.*|/analysis/filename mono_${E_current}_GeV.root|g;\
       s|/gun/energy.*|/gun/energy $E_current|g;\
       s|/run/beamOn.*|/run/beamOn $MC_NUMBER|g\
       " ${WORKDIR}/run_luxe.mac > ${TMP}/run_luxe.mac

  cd ${TMP}
  ${WORKDIR}/lxbeamsim ${TMP}/run_luxe.mac

  cp ${TMP}/mono_${E_current}_GeV.root ${WORKDIR}/output/ 
  
  #cd ${WORKDIR}/output
  #${WORKDIR}/lxbeamsim ${TMP}/run_luxe.mac

  echo $E_current "Run mono done!"
else
	echo "No such option!"
fi
