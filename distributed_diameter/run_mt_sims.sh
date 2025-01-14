#!/bin/bash
# Define radius (um)
r_mean=2.00
# Define MT induced effective T2 range and interval (ms)
mtmax=300
mtmin=10
mtstep=10
# Loop through MT values and run one simulation at each value
for mt in $(seq ${mtmin} ${mtstep} ${mtmax})
	do
	# use fsl_sub to submit a job described by the following arguments to a computing cluster, could be replaced with more general slurm commands
	fsl_sub -q long "julia --project=. simulate_mt.jl --radius ${r_mean} --MT ${mt}"   # use the --project argument to specify the correct environment
done

