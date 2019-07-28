#!/bin/zsh

#Setting env
source /cvmfs/sft.cern.ch/lcg/views/setupViews.sh LCG_91 x86_64-slc6-gcc7-opt
WORKDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
LIST_DIR=/nfs/dust/zeus/group/mykytaua/LUXE/IPstrong/Lists
LIST=($LIST_DIR/bppp*.out)

echo $TMP
cd $TMP
pwd

echo "File for DATA: " ${LIST[$(($1))]}

file=${LIST[$(($1))]}

sed "s|/lxphoton/gun/MCParticlesFileList.*|/lxphoton/gun/MCParticlesFileList $file|g;\
	s|/analysis/setFileName.*|/analysis/setFileName  $(basename $file)\.root|g \
     " ${WORKDIR}/luxe_gamma_new.mac > ${TMP}/luxe_gamma_new.mac

cp -rp ${WORKDIR}/hist_settings.mac $TMP

cd ${WORKDIR}
pwd

./lxbeamsim ${TMP}/luxe_gamma_new.mac 1

echo "Condor done!"
