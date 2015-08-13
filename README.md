MDutils
=======

General scripts to manage Molecular Dynamics calculations

Author
------
Salvatore M Cosseddu 2015

Description
-----------

### gmxmd ###
Quickly generate input file for MD simulations with gromacs and run simulations. Three system are supported:
- local machine
- Cartesius (National Dutch Supercomputer)
- Bazis     (VU University Cluster)
Use the supported environments as template to add your system.

#### Usage ####

The first argument passed from the CLI (if any) will be used to set
the dependency.

Required:
* General:
  - outdir = directory were outputs are saved
  - gmxdir = directory were gmx input files are saved
  - outname = prefix of the output file passed to -deffnm

* System definition : 
  - struc_file = Structure of the simulated system
  - top_file = topology file

Optional :
* Queue system 
  - walltime
  - nnodes = number of nodes
  - nc_nodes = cores per node
 
* Simulation 
  - min = minimisation steps
  or
  - runsteps = number of MD steps 
  - dt = time step size (ps)
* Includes
  - include = includes for grompp (-I, see gmx manual)
* Ensemble
  - T = Temperature
  - thermostat (default  v-scale)
  - genT = Temperature to generate initial velocities
  - taut (default 2)
  - p = Pressure
  - barostat (default Parrinello-Rahman)
  - taup (default 5)
  - compressibility (default 10e-5)
* Outputs
  - out_positions = interval to write positions (defalult 1000)
  - out_velocities = interval to write velocities (default 0)
  - out_log = interval to writing logs (defalult 1000)
* Others
  - constraints = es. h-bonds (default none) 
  - AdditionalFlags = for mdrun (es \"-n ${inpdir}/qd_lig_s.ndx\") 
  - AdditionalControls= additional lines to be included in the mdp file (es.
AdditionalControls="
   refcoord-scaling = all
"

____________________________________________________________________

### gmxcontinue ###

Provide a quick way to restart gmx run in a tidy manner (directories
gmxinput and output.$run are created for input and outputs
respectively) Two system are supported:
- local machine
- Cartesius (National Dutch Supercomputer, SLURM)
Use the supported environments as template to add your system.

#### Usage ####

Options:
* -i : Checkpoint input file (.cpt) from the previous run
* -b : Binary input file (.tpr) from the previous run
* -t : Additional time (ps
* -o : Outpref
* -r : Run number
* Queue options:
  + -d : Dependency (jobid)
  + -n : N. nodes
  + -c : N. Cores per node
  + -w : Walltime	
  + -h : Usage

es. 
gmxcontinue -i inname.cpt -b inname.tpr -t 50000 -r 2 -o outname



