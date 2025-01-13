# MT-and-permeability-effect-on-two-compartment-dMRI-WM-model
Repo for code that reproduces the results and analysis in the paper Investigating the sensitivity of diffusion MRI measurements to magnetization transfer and permeability via Monte-Carlo simulations &lt;add link>

## Pre-requisite
1. Install Julia from the [official website](https://julialang.org/downloads/).
2. Download/clone this repo and keep the file hierachy.
3. Install anaconda following the [official website](https://docs.anaconda.com/anaconda/install/) if you don't have it.
4. Open anaconda prompt or terminal, create a new environment if you want by `conda create --name dmipy` (dmipy can be replaced by any name) and activate the environment by `conda activate dmipy`. Then enter the following command to install necessary packages:
   ```
   conda install anaconda::jupyter
   conda install anaconda::pandas
   conda install conda-forge::matplotlib
   conda install numpy=1.26.4
   ```
5. Clone/download the dmipy package from [github](https://github.com/AthenaEPI/dmipy?tab=readme-ov-file) and run `python setup.py install` to install it. The pip install unfortunately no longer works reliably.

## Example operation
The notebook [Example_simulation](https://github.com/zhiyuzheng1769/MT-and-permeability-effect-on-two-compartment-dMRI-WM-model/blob/main/Example_simulation.ipynb) demonstrates how to set up a single simulation using parallel cylinder substrates and diffusion-weighted spin echo sequences. The simulation demonstrated there represents the most elementary unit of the mass simulations we ran for the whole project. Note it uses v0.11.0 of MCMRSimulator, which is the most up-to-date as of 2025-1. this version is however considerably different from the v0.9.0 which we used for producing the results in paper and in the next section.

## Reproduce results in paper
The [fixed_diameter](https://github.com/zhiyuzheng1769/MT-and-permeability-effect-on-two-compartment-dMRI-WM-model/tree/main/fixed_diameter) and [distributed_diameter](https://github.com/zhiyuzheng1769/MT-and-permeability-effect-on-two-compartment-dMRI-WM-model/tree/main/distributed_diameter) directories contain the files to reproduce the results from the fixed diameter case and distributed diameter case in the paper. Follow the README.md in them to reproduce the plots.
