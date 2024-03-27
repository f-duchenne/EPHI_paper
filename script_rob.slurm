#!/bin/bash

#SBATCH --job-name=rob            # job name
#SBATCH --ntasks=1                        # ask for 1 task
#SBATCH --mem=2G 
#SBATCH --array=1-300:1
#SBATCH --cpus-per-task=1
#SBATCH --time=10:00:00                   # time limit in hrs:min:sec
#SBATCH --output=serial_rob__%j.out       # output and error log


SCRIPT="/home/duchenne/EPHI_paper/robustness.r"
singularity exec /home/duchenne/ubuntu2.sif R --vanilla -f $SCRIPT --args $SLURM_ARRAY_TASK_ID

