#!/bin/sh
#SBATCH -N 64 -c 64
#SBATCH -p debug

#SBATCH --error="slurm.err"
#SBATCH --output="slurm.out"

#SBATCH --account=m2650
#SBATCH -t 00:30:00
#SBATCH -C haswell
#SBATCH -L project 
module load taskfarmer
#export PATH=$PATH:/usr/common/tig/taskfarmer/1.5/bin:$(pwd)
export THREADS=252
pwd
runcommands.sh taskfile.sh
