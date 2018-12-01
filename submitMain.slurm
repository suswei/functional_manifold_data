#!/bin/bash
# Created by the University of Melbourne job script generator for SLURM
# Sun Dec 02 2018 08:11:14 GMT+1100 (Australian Eastern Daylight Time)

# Partition for the job:
#SBATCH --partition=cloud

# Multithreaded (SMP) job: must run on one node and the cloud partition
#SBATCH --nodes=1

# The name of the job:
#SBATCH --job-name="manifoldFDAgeodesic"

# Maximum number of tasks/CPU cores used by the job:
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1

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
#SBATCH --time=0-1:0:00

# check that the script is launched with sbatch
if [ "x$SLURM_JOB_ID" == "x" ]; then
   echo "You need to submit your job to the queuing system with sbatch"
   exit 1
fi

# Run the job from the directory where it was launched (default)

# The modules to load:
module load MATLAB

# The job command(s):
R main.R