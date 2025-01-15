#!/bin/bash
# Define radius (um)
r_mean=2.00
# Define permeability range and interval
permmax=0.020
permmin=0.000
permstep=0.001
# Loop through permeability values and run one simulation at each value 
for perm in $(seq ${permmin} ${permstep} ${permmax})
	do
	# use fsl_sub to submit a job described by the following arguments to a computing cluster, could be replaced with more general slurm commands
	fsl_sub -q long "julia --project=. simulate_permeability.jl --radius ${r_mean} --perm ${perm}"  # use the --project argument to specify the correct environment 
done

