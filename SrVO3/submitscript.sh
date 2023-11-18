#!/bin/bash
#
#PBS -N PBE_5_TZ
#PBS -m a
#PBS -l walltime=01:00:00
#PBS -l nodes=1:ppn=1
#PBS -l mem=48GB
#

STARTDIR=$PBS_O_WORKDIR
export I_MPI_COMPATIBILITY=4


#For doduo cluster
module purge
module load pyWannier90/2021-12-07-foss-2021a
module load matplotlib/3.4.2-foss-2021a


cd $STARTDIR
echo "PBS: $PBS_ID"

ls

echo "Job started at : "`date`
python main_5_TZ_MLWF.py
echo "Job ended at : "`date`
