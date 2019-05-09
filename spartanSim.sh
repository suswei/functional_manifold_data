#!/usr/bin/env bash

# example from https://stackoverflow.com/questions/50653937/how-to-run-a-job-array-in-r-using-the-rscript-command-from-the-command-line

# Partition for the job:
#SBATCH --partition=cloud


# The name of the job:
#SBATCH --job-name="manifoldFDAgeodesic"

# Maximum number of tasks/CPU cores used by the job:
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8

# Send yourself an email when the job:
# aborts abnormally (fails)
#SBATCH --mail-type=FAIL
# begins
#SBATCH --mail-type=BEGIN
# ends successfully
#SBATCH --mail-type=END

# Use this email address:
#SBATCH --mail-user=susan.wei@unimelb.edu.au

# The maximum running time of the job in days-hours:mins:sec
#SBATCH --time=0-2:00:00

#SBATCH --array=1-2400%60         ### Array index | %50: number of simultaneously tasks

# Get Array ID
i=${SLURM_ARRAY_TASK_ID}

# check that the script is launched with sbatch
if [ "x$SLURM_JOB_ID" == "x" ]; then
   echo "You need to submit your job to the queuing system with sbatch"
   exit 1
fi

# Run the job from the directory where it was launched (default)

# The modules to load:
module load MATLAB

# The job command(s):
# R --vanilla < spartanSim.R ${i} #this isn't necessary since we get the environment
R --vanilla < spartanSim.R
