#!/bin/sh
#PBS -W group_list=ku_10024 -A ku_10024
#PBS -N dry-run-snakemake_transfer
#PBS -e dry-run-snakemake_transfer.err
#PBS -o dry-run-snakemake_transfer.out
#PBS -l nodes=1:ppn=4:thinnode
#PBS -l mem=1gb
#PBS -l walltime=5:00

# Go to the directory from where the job was submitted (initial directory is $HOME)
echo Working directory is $PBS_O_WORKDIR
cd $PBS_O_WORKDIR

exec 1>${PBS_JOBID}.out
exec 2>${PBS_JOBID}.err

module purge
module load tools
module load anaconda3/4.4.0

# I'll just run everything on one thin node instead of scheduling multiple jobs (don't know how to do that)
snakemake -R all --dry-run 
