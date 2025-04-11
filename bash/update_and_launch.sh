#!/bin/bash
#
# ****************************************
# Confounder Handling Simulation Study
#
# BluePebble Automation Script
# This bash scripts automates updating to the most recent
# code version from GitHub and submitting the job to BP
#
# Emma Tarmey
#
# Started:          11/02/2025
# Most Recent Edit: 11/04/2025
# ****************************************

echo ""

# delete older version
rm -f -r Confounder-Handling-Sim-Study

# clone most recent version
git clone https://github.com/RaspberryEmma/Confounder-Handling-Sim-Study

# change wd
cd Confounder-Handling-Sim-Study
cd bash

# import python
module load languages/python/3.12.3

# submit simulation to BP HPC
for i in 1 2 3 4 5 6 7 8 9;
do
	echo "Submitting job: launch_BP_run_"$i"_step_0.sh"
	sbatch "launch_BP_run_"$i"_step_0.sh"

	echo "Submitting job: launch_BP_run_"$i"_step_1.sh"
	sbatch "launch_BP_run_"$i"_step_1.sh"
done

# submit simulation to BP HPC
# null simulations
for i in 1 2 3 4 5 6 7 8 9;
do
	echo "Submitting job: launch_BP_run_null_"$i"_step_0.sh"
	sbatch "launch_BP_run_null_"$i"_step_0.sh"

	echo "Submitting job: launch_BP_run_null_"$i"_step_1.sh"
	sbatch "launch_BP_run_null_"$i"_step_1.sh"
done

# check jobs submitted correctly
sleep 5.0
echo ""
sacct -X
echo ""

