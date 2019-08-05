#!/bin/zsh

#Setting env
source /cvmfs/sft.cern.ch/lcg/views/setupViews.sh LCG_91 x86_64-slc6-gcc7-opt
WORKDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

#LIST_DIR=/nfs/dust/zeus/group/mykytaua/LUXE/IPstrong/Lists
#LIST=($LIST_DIR/bppp*.out)

LIST_FILE=/nfs/dust/zeus/group/mykytaua/LUXE/IPstrong/Lists/bppp_17.5GeV_5mfoiltoIP_Enelas_1.0J.out
LIST=(`cat $LIST_FILE`)

echo $TMP
cd $TMP
pwd

echo "File for DATA: " ${LIST[$(($1))]}

file=${LIST[$(($1))]}

# sed "s|/lxphoton/gun/MCParticlesFileList.*|/lxphoton/gun/MCParticlesFileList $file|g;\
# 	s|/analysis/setFileName.*|/analysis/setFileName  $(basename $file)\.root|g \
#      " ${WORKDIR}/luxe_gamma_new.mac > ${TMP}/luxe_gamma_new.mac

echo $file
sed "s|/lxphoton/gun/MCParticlesFile.*|/lxphoton/gun/MCParticlesFile $file|g\
     " ${WORKDIR}/luxe_gamma_new.mac > ${TMP}/luxe_gamma_new.mac

# cp -rp ${WORKDIR}/hist_settings.mac $TMP

# cd ${WORKDIR}
pwd

${WORKDIR}/lxbeamsim luxe_gamma_new.mac
echo "Run done!"

mkdir -p ${WORKDIR}/rootOUT
 
cp -rf ./test.root ${WORKDIR}/rootOUT/$(basename $file).root

echo "Copy done!"
