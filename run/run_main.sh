#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
LIST=/nfs/dust/zeus/group/mykytaua/LUXE/IPstrong/Lists/bppp_14GeV_5mfoiltoIP_Enelas_0.2J.out
# LIST=/nfs/dust/zeus/group/mykytaua/LUXE/IPstrong/Lists/bppp_14GeV_5mfoiltoIP_Enelas_0.35J.out
# LIST=/nfs/dust/zeus/group/mykytaua/LUXE/IPstrong/Lists/bppp_14GeV_5mfoiltoIP_Enelas_0.5J.out
# LIST=/nfs/dust/zeus/group/mykytaua/LUXE/IPstrong/Lists/bppp_14GeV_5mfoiltoIP_Enelas_0.7J.out
# LIST=/nfs/dust/zeus/group/mykytaua/LUXE/IPstrong/Lists/bppp_14GeV_5mfoiltoIP_Enelas_0.85J.out
# LIST=/nfs/dust/zeus/group/mykytaua/LUXE/IPstrong/Lists/bppp_14GeV_5mfoiltoIP_Enelas_1.0J.out
# LIST=/nfs/dust/zeus/group/mykytaua/LUXE/IPstrong/Lists/bppp_17.5GeV_5mfoiltoIP_Enelas_0.2J.out
# LIST=/nfs/dust/zeus/group/mykytaua/LUXE/IPstrong/Lists/bppp_17.5GeV_5mfoiltoIP_Enelas_0.35J.out
# LIST=/nfs/dust/zeus/group/mykytaua/LUXE/IPstrong/Lists/bppp_17.5GeV_5mfoiltoIP_Enelas_0.5J.out
# LIST=/nfs/dust/zeus/group/mykytaua/LUXE/IPstrong/Lists/bppp_17.5GeV_5mfoiltoIP_Enelas_0.7J.out
# LIST=/nfs/dust/zeus/group/mykytaua/LUXE/IPstrong/Lists/bppp_17.5GeV_5mfoiltoIP_Enelas_0.85J.out
# LIST=/nfs/dust/zeus/group/mykytaua/LUXE/IPstrong/Lists/bppp_17.5GeV_5mfoiltoIP_Enelas_1.0J.out
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


WORKDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
OUTDIR=/nfs/dust/zeus/group/mykytaua/LUXE/2019.09.07

#Some stearing variables
E_start=-999
E_end=-999
Number=-999

while [ ! $# -eq 0 ]
do
    case "$1" in

        --run| -r)

            filename=$(basename $LIST)

            OUTDIR=${OUTDIR}/${filename}
            mkdir -p $OUTDIR

            mkdir -p $OUTDIR/output

            mkdir -p $OUTDIR/condor
            mkdir -p $OUTDIR/condor/log
            mkdir -p $OUTDIR/condor/err
            mkdir -p $OUTDIR/condor/out

            number=`cat $LIST | wc -l`
            echo 'number of file: '$number

            sed "s|JOBS_NUMBER|$number|g;\
            	   s|OUTDIR|$OUTDIR|g\
            	  " $WORKDIR/QA.job > $OUTDIR/QA.job

                sed "s|option=.*|option=true|g;\
                 s|LIST_FILE=.*|LIST_FILE=$LIST|g;\
                 s|WORKDIR=.*|WORKDIR=$OUTDIR|g\
            	  " $WORKDIR/run_condor.sh > $OUTDIR/run_condor.sh

	        ln -s $WORKDIR/lxbeamsim $OUTDIR/
	        cp $WORKDIR/luxe_gamma_new.mac $OUTDIR/

            cd $OUTDIR
           	condor_submit $OUTDIR/QA.job
            cd -
            shift
            ;;

        --E_start | -E1)
            shift
            E_start=$1
            shift
            ;;

        --E_end | -E2)
            shift
            E_end=$1
            shift
            ;;

        --Number | -n)
            shift
            Number=$1
            shift
            ;;

        --run_mono | -m)
            #This part of the code runs mono-energetical spectrum from E_start to E_end

            mkdir -p $OUTDIR/output
            mkdir -p $OUTDIR/condor
            mkdir -p $OUTDIR/condor/log
            mkdir -p $OUTDIR/condor/err
            mkdir -p $OUTDIR/condor/out


            sed "s|JOBS_NUMBER|$Number|g;\
            	   s|OUTDIR|$OUTDIR|g\
            	  " $WORKDIR/QA.job > $OUTDIR/QA.job

            sed "s|option=.*|option=mono|g;\
                 s|E_START=.*|E_START=$E_start|g;\
                 s|E_END=.*|E_END=$E_end|g;\  
                 s|WORKDIR=.*|WORKDIR=$OUTDIR|g\
            	  " $WORKDIR/run_condor.sh > $OUTDIR/run_condor.sh

  	        ln -s $WORKDIR/lxbeamsim $OUTDIR/
  	        cp $WORKDIR/run_luxe.mac $OUTDIR/

            cd $OUTDIR
           	condor_submit $OUTDIR/QA.job
            cd -
            shift
            ;;


        *)
            echo 'Not available flag!'
            break
            shift
            ;;
    esac
    shift
done