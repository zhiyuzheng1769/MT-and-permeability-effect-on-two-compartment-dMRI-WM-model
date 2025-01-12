#!/bin/bash
# Define radius range and interval (um)
rmax=5.50
rmin=0.30
rstep=0.10
# Define permeability range and interval
permmax=0.020
permmin=0.000
permstep=0.001
# Loop through radius and permeability values and run one simulation at each value pairs
for r_mean in $(seq ${rmin} ${rstep} ${rmax})
do
    for perm in $(seq ${permmin} ${permstep} ${permmax})
        do
        # use fsl_sub to submit a job described by the following arguments to a computing cluster, could be replaced with more general slurm commands
        fsl_sub  "julia --project=. simulate_perm.jl --radius ${r_mean} --perm ${perm}" # use the --project argument to specify the correct environment to use
    done
done
