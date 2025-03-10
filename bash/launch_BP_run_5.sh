#!/bin/bash
#
# ****************************************
# Confounder Handling Simulation Study
#
# BluePebble Launch Bash Script
# Defines and runs the R simulation procedure on the BluePebble HPC
#
# Emma Tarmey
#
# Started:          11/02/2025
# Most Recent Edit: 10/03/2025
# ****************************************
#
#SBATCH --partition=compute
#SBATCH --job-name=conf_sim_study_run_5
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --time=4-00:00:00
#SBATCH --mem-per-cpu=20G
#SBATCH --account=MATH033344
#SBATCH --mail-type=ALL
#SBATCH --mail-user=aa22294@bristol.ac.uk


# Change into working directory
cd ${SLURM_SUBMIT_DIR}
cd ..
cd R

# Record info
echo ""
echo "***** START *****"
echo "***** Confounder Handling Simulation Study - Simulation *****"
echo Start Time:        $(date)
echo Working Directory: $(pwd)
echo JOB ID:            ${SLURM_JOBID}
echo SLURM ARRAY ID:    ${SLURM_ARRAY_TASK_ID}
echo ""

# Import R
module load languages/R/4.4.1

# Execute code
# NB: we pass in the following:
# (#scenario) (#total confounders) (#measured) (#unmeasured)
# Rscript simplified_simulation_run_5.R 1  16  16  0
# Rscript simplified_simulation_run_5.R 2  16  12  4
# Rscript simplified_simulation_run_5.R 3  32  32  0
# Rscript simplified_simulation_run_5.R 4  32  24  8
Rscript simplified_simulation_run_5.R 5  64  64  0
Rscript simplified_simulation_run_5.R 6  64  48 16
Rscript simplified_simulation_run_5.R 7 128 128  0
Rscript simplified_simulation_run_5.R 8 128  96 32

echo ""
echo End Time: $(date)
echo "***** END *****"
echo ""

