# MT-and-permeability-effect-on-two-compartment-dMRI-WM-model
Repo for code that reproduces the results and analysis in the paperInvestigating the sensitivity of the diffusion MRI signal to magnetization transfer and permeability via Monte-Carlo simulations &lt;add link>. It is still work-in-progress and the content may be subject to substantial changes.

## Pre-requisite
1. Install Julia from the [official website](https://julialang.org/downloads/).
2. Download/clone this repo and keep the file hierachy.
3. Install pixi following the [official website](https://pixi.sh/latest/installation/) if you don't have it.
4. Open a terminal, change directory to where the repo is. Then run `pixi shell`.
5. In the same directory, run
   ```
   Julia --project=.
   Pkg.instantiate()
   ```
   
## Example operation
The notebook [Example_simulation](https://github.com/zhiyuzheng1769/MT-and-permeability-effect-on-two-compartment-dMRI-WM-model/blob/main/Example_simulation.ipynb) demonstrates how to set up a single simulation using parallel cylinder substrates and diffusion-weighted spin echo sequences. The simulation demonstrated there represents the most elementary unit of the mass simulations we ran for the whole project. Note it uses v0.11.0 of MCMRSimulator, which is the most up-to-date as of 2025-1. this version is however considerably different from the v0.9.0 which we used for producing the results in paper and in the next section. Run `pixi run julia-notebook` in terminal to start jupyter notebook with correct kernel and open the notebook.

## Reproduce results in paper
The [fixed_diameter](https://github.com/zhiyuzheng1769/MT-and-permeability-effect-on-two-compartment-dMRI-WM-model/tree/main/fixed_diameter) and [distributed_diameter](https://github.com/zhiyuzheng1769/MT-and-permeability-effect-on-two-compartment-dMRI-WM-model/tree/main/distributed_diameter) directories contain the files to reproduce the results from the fixed diameter case and distributed diameter case in the paper. Follow the README.md in them to reproduce the plots.
