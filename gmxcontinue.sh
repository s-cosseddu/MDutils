#!/bin/bash
# Quickly restart gromacs jobs in a tidy way. 

# gmx bin name (allow moving among versions)
GMXCONVERTTPR="gmx convert-tpr"
GMXMD="gmx mdrun"

function usage () {
    echo "
options:
	-i)
	    Checkpoint input file (.cpt)
	-b)
	    Bin input file (.tpr)
	-t)
	    Additional time (ps)
	-o)
	    outpref
	-r)
	    run number
	-d)
	    dependency (jobid)
        -n)
            N. nodes
	-c)
	    N. cores per node
	-w)
	    Walltime (hours)	
        -h)
	    usage

es. 
gmxcontinue -i inname.cpt -b inname.tpr -t 50000 -r 2 -o outname
    "
exit 0
}

# ---------- Read CLI ------
if [ $# == 0 ]; then usage; fi

gmxdir=gmxinput
outdir=output
MDINIT="$GMXCONVERTTPR"
MDRUN="$GMXMD"

while [[ $# > 0 ]]
do
    arg="$1"

    case $arg in
	-r)
	    run=".$2"
	    outdir=output$run
	    shift
	    ;;
	-i)
	    cpt="$2"
	    MDRUN="$MDRUN -cpi $cpt"
	    echo "cpt            :: $cpt" > gmxcontinue.log
	    shift
	    ;;
	-b)
	    tpr="$2"
	    MDINIT="$MDINIT -s $tpr"
	    echo "tpr            :: $tpr" > gmxcontinue.log
	    shift
	    ;;
	-t)
	    Extime="$2"
	    MDINIT="$MDINIT -extend $Extime"
	    echo "Extended time  :: $Extime" >> gmxcontinue.log 
	    shift
	    ;;
	-o)
	    outpref="$2"
	    MDINIT="$MDINIT -o ${gmxdir}/$outpref.tpr"
	    MDRUN="$MDRUN -s ${gmxdir}/$outpref.tpr -deffnm ${outdir}/$outpref"
	    echo "new prefix     :: $outpref">> gmxcontinue.log
	    shift
	    ;;
	-d)
	    SBATCHFLAGS="$SBATCHFLAGS --dependency=afterany:$2"
	    echo "Starting after     :: $2">> gmxcontinue.log
	    shift
	    ;;
	-n)
	    nnodes="$2"
	    shift
	    ;;
	-c)
	    nc_ncores="$2"
	    shift # past argument
	    ;;
	-w)
	    walltime="$2"
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

if [[ -z $cpt || -z $tpr ]]
then 
    echo "tpr/cpt required" >&2
    usage
    exit 1
fi
if [ -z $Extime ];
then 
    echo "Extended time required" >&2
    usage
    exit 1
fi
if [ -z $outpref ];
then 
    echo "new prefix required" >&2
    usage
    exit 1
fi
# ==================================================


# echo commands
echo "gmx commands:              
$MDINIT
$MDRUN
" >> gmxcontinue.log


# ==================================================
#  DOING the JOB
mkdir -p $gmxdir $outdir
# prepare restart
$MDINIT

# submit
if [ $(hostname | cut -d. -f2) == cartesius ]; then

    WT=${walltime:=120}
    MDRUN="$MDRUN -maxh $WT"

cat > slurm.sub <<EOF
#!/bin/bash
#SBATCH --job-name=QDMD${run}
#SBATCH --output=QDMD.txt
#
#SBATCH -N ${nnodes:=4}
#SBATCH --tasks-per-node ${nc_nodes:=24}
#SBATCH -t ${WT:=120}:00:00

# for test, scratch not supported
srun $MDRUN &> $outdir/${outpref}.log.mdrun

exit 0

EOF

sbatch $SBATCHFLAGS slurm.sub

else
    # local machine
    mpirun -np 3 $MDRUN &> $outdir/${outpref}.log.mdrun
fi

exit 0


