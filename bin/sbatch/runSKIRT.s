#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=40
#SBATCH --time=10:00:00
#SBATCH --mem=80GB
#SBATCH --job-name=runSKIRT
#SBATCH --mail-type=END
#SBATCH --output=slurm_out/slurm_%x.out
#SBATCH --mail-user=ntf229@nyu.edu

module purge

cd /home/ntf229/containers

singularity exec --overlay overlay-15GB-500K.ext3:ro \
	    /scratch/work/public/singularity/cuda11.0-cudnn8-devel-ubuntu18.04.sif \
/bin/bash -c "source /ext3/env.sh; 
python /home/ntf229/CCA/880Msun_split/bin/runSKIRT.py"



