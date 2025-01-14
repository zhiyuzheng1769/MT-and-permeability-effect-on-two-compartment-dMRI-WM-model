#!/bin/bash
# Define radius range and interval (um)
rmax=5.50
rmin=0.30
rstep=0.10
# Define MT induced effective T2 range and interval (ms)
mtmax=150
mtmin=10
mtstep=10
# Loop through radius and MT values and run one simulation at each value pairs
for r_mean in $(seq ${rmin} ${rstep} ${rmax})
do
    for mt in $(seq ${mtmin} ${mtstep} ${mtmax})
        do
        # use fsl_sub to submit a job described by the following arguments to a computing cluster, could be replaced with more general slurm commands
        fsl_sub "julia  --project=. simulate_mt.jl --radius ${r_mean} --MT ${mt}" # use the --project argument to specify the correct environment
    done
done
