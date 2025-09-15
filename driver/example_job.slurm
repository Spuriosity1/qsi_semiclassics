#!/bin/bash

#! Which project should be charged:
#SBATCH -A CASTELNOVO-SL2-CPU
#SBATCH -p icelake
#! How many whole nodes should be allocated?
#SBATCH --nodes=1
#! How many (MPI) tasks will there be in total? (<= nodes*76)
#! The Ice Lake (icelake) nodes have 76 CPUs (cores) each and
#! 3380 MiB of memory per CPU.
#SBATCH --ntasks=1
#! How much wallclock time will be required?
#SBATCH --time=2:00:00
#! What types of email messages do you wish to receive?
#SBATCH --mail-type=START,END
#!
#SBATCH -J static_corr
#SBATCH --output=track/static_corr.%A_%a.out
#SBATCH --error=track/static_corr.%A_%a.error
#SBATCH --cpus-per-task=1
#!
#SBATCH --array=1-80 

. /etc/profile.d/modules.sh                # Leave this line (enables the module command)
module purge                               # Removes all modules still loaded
module load rhel8/default-icl              # REQUIRED - loads the basic environment

export OPENBLAS_NUM_THREADS=$SLURM_CPUS_PER_TASK
export MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK


JOBID=$SLURM_JOB_ID
TASKID=$SLURM_ARRAY_TASK_ID

echo -e "TaskID: $TASKID\n======"
echo "Time: `date`"
echo "Running on master node: `hostname`"
echo "Current directory: `pwd`"

#if [ "$SLURM_JOB_NODELIST" ]; then
#        #! Create a machine file:
#        export NODEFILE=`generate_pbs_nodefile`
#        cat $NODEFILE | uniq > machine.file.$JOBID
#        echo -e "\nNodes allocated:\n================"
#        echo `cat machine.file.$JOBID | sed -e 's/\..*$//g'`
#fi


workdir=/home/als217/rds/hpc-work/tmp/qsi_semiclassics/driver

cd $workdir
echo -e "Changed directory to `pwd`.\n"


# Load the command from the file
COMMAND=$(sed "${SLURM_ARRAY_TASK_ID}q;d" "static_corr.plan")

echo "Running: $COMMAND"
eval "$COMMAND"

