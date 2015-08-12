#!/bin/bash
# 

# File names
inpdir=../input
outdir=output.MD
gmxdir=gmxinput
struc_file=${inpdir}/qd_lig_s_rel.gro
top_file=${inpdir}/qd_lig_s.top
outname=CdSe_oleateMD

# queue system variable
walltime=1
nnodes=4
nc_nodes=24
if [ -z $1 ];
then
    # dependency supported only for slurm (cartesius)
    SBATCHFLAGS="--dependency=afterok:$1"
fi
    
# simulation ----------
#min=5000			# overprint runstep
runsteps=1000000
dt=0.002                        # needed to runstep (no control implemented...)

# includes
include="include = -I$FFIELD"
#define="define = -DPOSRESIONS"


# ---------
T=300
genT=300			# needed only to generate initial velocities
p=1.01325
# barostat=berendsen 		# default Parrinello-Rahman
# thermostat=berendsen            # default  v-scale
# taup=0.8
# taut=0.5


# restraints
constraints=h-bonds 	        # default none 

# AdditionalFlags="-n ${inpdir}/qd_lig_s.ndx" # flag -n necessary
# AdditionalControls="
# refcoord-scaling = all
# "

#outputs
out_positions=2500
#out_velocities=2500
out_log=1000

# ====================================================================================================
#                                        Input finished 

# -------------------------------
#  INTEGRATOR CONTROL 
if [ $min ]; then
    Runcontrol=" 
$include
$define

;Run control:
integrator	= steep		; Algorithm (steep = steepest descent minimization)
emtol		= 10.0		; Stop minimization when the maximum force < 1.0 kJ/mol
nsteps		= $min 		; Maximum number of (minimization) steps to perform
pbc          = xyz		; Periodic Boundary Conditions
"
elif [ $runsteps ] ;then
    Runcontrol="
$include
$define

;Run control: 
integrator		 = md  ; leap-frog
dt			 = $dt
nsteps  		 = $runsteps ; Maximum number of  steps to perform
pbc          = xyz		; Periodic Boundary Conditions


"
else
    echo "Put either min or runsteps"
    exit 1
fi

# ------------------------------
# init velocities
if [ $genT ]; then
    velocitycontrol="
;Velocity generation
gen_vel 		 = yes
gen_temp                 = $genT
"
else
    velocitycontrol="
;Velocity generation
gen_vel 		 = no
"
fi

# -------------------------
#  GENERAL PARAMETERS

# To be modified only rarely
GENpars="
;treatment of van der waals interactions
vdwtype			 = Cut-off
cutoff-scheme		 = Verlet
vdw-modifier             = Potential-shift-Verlet
dispcorr                 = EnerPres
coulombtype		 = PME ; PME-Switch 
nstcalclr                = 10

; CHARMM cut-off 
rlist        = 1
rvdw         = 1.2
rcoulomb     = 1.2
rvdw-switch  = 1
rcoulomb-switch  = 1

;Frequency to update the neighbor list (and the long-range forces, 
;when using twin-range cut-off's). 
nstlist 		 = 30
ns-type 		 = grid
;cut-off distance for the short-range neighbor list

"

# -------------------------
# OUTPUT FREQUENCIES 
output_freqs="
;frequency to write coordinates to output trajectory file
nstxout 		 = $out_positions   
nstvout 		 = ${out_velocities:=0}
;frequency to write energies to log file
nstlog  		 = $out_log
;frequency to write energies to energy file
nstenergy		 = $out_log
;group(s) to write to energy file 
energygrps		 = System  
"

# -------------------------
#  THERMOSTAT
if [ $T ]; then
    tempcontrol="
;Temperature coupling
tcoupl  		 = ${thermostat:=v-rescale}
tc-grps 		 = System    
tau-t			 = ${taut:=2}    
ref-t			 = $T
"
else 
    tempcontrol="
verlet-buffer-tolerance  = 5.2e-08
"
fi

# -------------------------
#  BAROSTAT

if [ $p ]; then
    presscontrol="
; ;Pressure coupling
; pcoupl  		 = berendsen
pcoupl  		 = ${barostat:=Parrinello-Rahman}
pcoupltype		 = isotropic
tau-p			 = ${taup:=5}
compressibility 	 = 10e-5 # dichloromethane
ref-p			 = $p
"
fi

if [ $constraints ]; then
    constcontrol="
; ---------------------
; constraints
constraints		 = ${constraints:=none}
; constraint_algorithm	 = shake
; shake_tol                = 0.0001
lincs-order = 6
lincs-iter = 2
lincs-warnangle = 90
"
fi

conffile=$gmxdir/$outname.mdp
ncores=$((nc_nodes*nnodes))

mkdir -p $gmxdir

cat > $conffile <<EOF
$Runcontrol
$GENpars
$output_freqs
$velocitycontrol
$tempcontrol
$presscontrol
$constcontrol
$AdditionalControls

EOF

if [ $(hostname) == bazis-h0 ]; then
# --- for VU cluster ---

module load gromacs
sbatch <<EOF
#!/bin/bash
#SBATCH --job-name=MD-run
#SBATCH --output=adf.txt
#
#SBATCH --nodes=$nnodes
#SBATCH --tasks=$nc_nodes
#SBATCH --time=${walltime}:00:00
#SBATCH --mem=60000

export PATH=$PATH:\$HOME/bin
export SCM_TMPDIR=/local/datastore0/\$SLURM_JOBID

# preparing run
mkdir -p $outdir
gmx grompp -c $struc_file -p $top_file -o $gmxdir/$outname $AdditionalFlags -f $conffile &> $outdir/$outname.log.grompp
echo "Running on $ncores" 


# preparing tmp directory
mkdir -p \$SCM_TMPDIR
cd \$SCM_TMPDIR
cp -f \$SLURM_SUBMIT_DIR/../gmxinput/$outname.tpr . 

# run
mpirun -np $ncores gmx mdrun -maxh $walltime -s $outname.tpr -deffnm $outname &> ${outname}.log.mdrun

# finalisation
rm \$SLURM_SUBMIT_DIR/$outdir/$outname.tpr
mv -f * \$SLURM_SUBMIT_DIR/$outdir


exit 0

EOF

elif [ $(hostname | cut -d. -f2) == cartesius ]; then

module load mpi
mkdir -p $outdir
gmx grompp -c $struc_file -p $top_file -o $gmxdir/$outname $AdditionalFlags -f $conffile &> $outdir/$outname.log.grompp

sbatch $SBATCHFLAGS <<EOF
#!/bin/bash
#SBATCH --job-name=QDMD
#SBATCH --output=QDMD.txt
#
#SBATCH -N $nnodes
#SBATCH --tasks-per-node $nc_nodes
#SBATCH -t ${walltime}:00:00

# for test, scratch not supported
srun gmx mdrun -maxh $walltime -s $gmxdir/$outname.tpr -deffnm $outdir/$outname &> $outdir/${outname}.log.mdrun


exit 0

EOF


else
    # local machine
    mkdir -p $outdir
    gmx grompp -c $struc_file -p $top_file -o $gmxdir/$outname $AdditionalFlags -f $conffile &> $outdir/$outname.log.grompp
    mpirun -np 3 gmx mdrun -s $gmxdir/$outname.tpr -deffnm $outdir/$outname &> $outdir/${outname}.log.mdrun
fi

exit 0
