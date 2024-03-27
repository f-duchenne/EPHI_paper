#!/bin/bash

#SBATCH --job-name=jags_trait            # job name
#SBATCH --ntasks=1                        # ask for 1 task
#SBATCH --mem=16G 
#SBATCH --array=4-6
#SBATCH --cpus-per-task=1
#SBATCH --time=710:00:00                   # time limit in hrs:min:sec
#SBATCH --output=test_%j.out       # output and error log


SCRIPT="/home/duchenne/EPHI_paper/jags_model.r"
singularity exec /home/duchenne/ubuntu2.sif R --vanilla -f $SCRIPT --args $SLURM_ARRAY_TASK_ID

