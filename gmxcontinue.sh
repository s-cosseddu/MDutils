#!/bin/bash
# Quickly restart gromacs jobs in a tidy way. 

function usage () {
    echo "
options:
	-i)
	    Inpref
	-t)
	    Additional time
	-o)
	    outpref
	-r)
	    run number
	-d)
	    dependency (jobid)
	-c)
	    N. Cores
	-w)
	    Walltime	
        -h)
	    usage
    "
exit 0
}

# ---------- Read CLI ------
if [ $# == 0 ]; then usage; fi

while [[ $# > 0 ]]
do
    arg="$1"

    case $arg in
	-i)
	    inpref="$2"
	    shift
	    ;;
	-t)
	    Extime="$2"
	    shift
	    ;;
	-o)
	    outpref="$2"
	    shift
	    ;;
	-r)
	    run=".$2"
	    shift
	    ;;
	-d)
	    SBATCHFLAGS="$SBATCHFLAGS --dependency=afterany:$2"
	    shift
	    ;;
	-c)
	    NCORES="$2"
	    shift # past argument
	    ;;
	-w)
	    WALLTIME="$2"
	    shift # past argument
	    ;;
	-h)
	    usage
	    exit 1
	    ;;
	*)
	    # unknown option
	    echo "wrong options"
	    usage
	    exit 1
	    ;;
    esac
    shift # past argument or value
done

if [ -z $inpref ]
then 
    echo "tpr/cpt prefix required" >&2
fi
if [ -z $Extime ];
then 
    echo "Extended time required" >&2
fi
if [ -z $outpref ];
then 
    echo "new prefix required" >&2
fi

echo "tpr/cpt prefix :: $inpref" > gmxcontinue.log
echo "Extended time  :: $Extime" >> gmxcontinue.log 
echo "new prefix     :: $outpref">> gmxcontinue.log
echo "outdir         :: $outdir" >> gmxcontinue.log
echo "gmx commands:              
gmx tpbconv -s $inpref.tpr -extend $Extime -o $outpref.tpr
gmx mdrun -s $outpref.tpr -cpi $inpref.cpt -deffnm $outpref" >> gmxcontinue.log

# ==================================================


gmxdir=gmxinput
outdir=output$run

mkdir -p gmxdir $outdir

# prepare restart
gmx tpbconv -s $inpref.tpr -extend $Extime -o ${gmxdir}/$outpref.tpr

# submit

if [ $(hostname | cut -d. -f2) == cartesius ]; then

sbatch $SBATCHFLAGS <<EOF
#!/bin/bash
#SBATCH --job-name=QDMD${run}
#SBATCH --output=QDMD.txt
#
#SBATCH -N $nnodes
#SBATCH --tasks-per-node $nc_nodes
#SBATCH -t ${walltime}:00:00

# for test, scratch not supported
srun gmx mdrun -maxh $walltime -s $gmxdir/$outname.tpr -cpi $inpref.cpt -deffnm ${outdir}/$outpref &> $outdir/${outpref}.log.mdrun

exit 0

EOF


else
    # local machine
    mpirun -np 3 gmx mdrun -s gmxinput/$outpref.tpr -cpi $inpref.cpt -deffnm ${outdir}/$outpref &> $outdir/${outpref}.log.mdrun
fi

exit 0


