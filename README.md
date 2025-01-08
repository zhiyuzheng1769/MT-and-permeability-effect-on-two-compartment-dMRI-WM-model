# MT-and-permeability-effect-on-two-compartment-dMRI-WM-model
Repo for code that reproduces the results and analysis in the paper Investigating the sensitivity of diffusion MRI measurements to magnetization transfer and permeability via Monte-Carlo simulations &lt;add link>

## Pre-requisite
1. Install Julia from the [official website](https://julialang.org/downloads/).
2. Install MCMRSimulator version 0.9.0 following its [documentation](https://open.win.ox.ac.uk/pages/ndcn0236/mcmrsimulator.jl/v0.9/installation/).
3. Install Python 3.13 and above with jupyter notebook, numpy 2.2.1.

## Example operation
The notebook [Example_simulation](https://github.com/zhiyuzheng1769/MT-and-permeability-effect-on-two-compartment-dMRI-WM-model/blob/main/Example_simulation.ipynb) demonstrates how to set up a single simulation using parallel cylinder substrates and diffusion-weighted spin echo sequences. The simulation demonstrated there represents the most elementary unit of the mass simulations we ran for the whole project. 

## Reproduce results in paper
The [fixed_diameter](https://github.com/zhiyuzheng1769/MT-and-permeability-effect-on-two-compartment-dMRI-WM-model/tree/main/fixed_diameter) and [distributed_diameter](https://github.com/zhiyuzheng1769/MT-and-permeability-effect-on-two-compartment-dMRI-WM-model/tree/main/distributed_diameter) directories contain the files to reproduce the results from the fixed diameter case and distributed diameter case in the paper. Follow the README.md in them to reproduce the plots.
