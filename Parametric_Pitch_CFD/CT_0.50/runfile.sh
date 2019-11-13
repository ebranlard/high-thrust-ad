#!/bin/bash
#SBATCH --job-name=alm                # Job name
#SBATCH --time=30:00:00
#SBATCH --nodes=5                              # Number of nodes
#SBATCH --ntasks-per-node=36                    # Number of processors per node
#SBATCH --account=wakedynamics                            # Allocation
##SBATCH --account=hfm                            # Allocation
##SBATCH --qos=high     # Priority
#SBATCH --mail-user luis.martinez@nrel.gov     # E-mail adres
#SBATCH --mail-type BEGIN,END,FAIL              # Send e-mail when job begins, ends or fails

source $HOME/.bash_profile
nalu_env intel


srun -n 1 -c 1 --cpu_bind=cores /projects/hfm/shreyas/exawind/install/intel/wind-utils/bin/abl_mesh -i alm_preprocess.yaml &> log.mesh
srun -n 1 -c 1 --cpu_bind=cores /projects/hfm/shreyas/exawind/install/intel/wind-utils/bin/abl_mesh -i sampling_mesh.yaml &> log.sampling_mesh

############### NEW VERSION OF NALU!!!!!!!!!!!!!!
srun -n 180 -c 1 --cpu_bind=cores /home/lmartine/nalu-wind/build/naluX -i alm_simulation.yaml &> log
