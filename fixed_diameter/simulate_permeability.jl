""" 
This script runs the sequences used in the AxCaliber paper (https://doi.org/10.1002/mrm.21577) on a single cylinder substrate we created before, it should be called with two arguments --radius and --perm (permeability) provided.
"""
# Import necessary libraries for argument parsing, simulation, output formatting, file handling
using ArgParse
using MCMRSimulator
using Printf
using DelimitedFiles
using FileIO

# Function to parse command-line arguments
function parse_commandline()
    # Initialize argument parser settings
    s = ArgParseSettings()

    # Define expected command-line arguments
    @add_arg_table! s begin
        "--radius"      # The radius argument for the cylinders
        help = "radius" # Key for the radius argument
        "--perm"        # The permeability argument
        help = "perm"   # Key for the permeability argument
    end

    # Parse the arguments and return them as a dictionary
    return parse_args(s)
end

# Main function to run the simulation
function main()
    # Parse the command-line arguments
    args = parse_commandline()

    # Extract and round the radius and permeability arguments
    radius = round(parse(Float64, args["radius"]), digits=3)  # Radius rounded to 3 decimal places (um)
    perm = round(parse(Float64, args["perm"]), digits=4)      # Permeability (rounded to 4 decimal places)

    # Set MRI sequence parameters
    TR = 3000          # Repetition time (ms)
    TE = 166           # Echo time (ms)
    delta = 2.5        # Gradient duration (ms)
    Gmax = 1200        # Maximum gradient strength (mT/m)
    gamma = 0.00004257638476  # Gyromagnetic ratio (kHz/mT ̇1e-6), needed as dwi() takes gradient strength with kHz/μm unit
    
    # Initialize an empty array to store the diffusion-weighted MRI sequences
    seqs = Sequence[]

    # Recreate the sequences used in the original AxCaliber paper: https://doi.org/10.1002/mrm.21577
    # Loop through the relevant diffusion time values
    for Delta in [10, 15, 20, 30, 40, 50, 60, 80] 
        # Loop over gradient strengths in steps of 80, up to Gmax
        for g in 0:80:Gmax
            # Create a DWSE sequence with the current parameters
            seq = dwi(gradient_duration=delta, gradient_strength=g*gamma, diffusion_time=Delta, TE=TE, TR=TR)
            # Add the generated sequence to the sequences vector
            push!(seqs, seq)
        end
    end

    # Calculate the size of the bounding box based on the radius
    sz = radius * 200  # Scaling factor for the bounding box

    # Read the corresponding cylinder substrates from a JSON file
    # The file name is dynamically created based on the perm and radius arguments
    geom = MCMRSimulator.read_geometry("./perm_cyl/cylinders_MT_0_sus_0_perm_" * @sprintf("%.3f", perm) * "_rmean_" * @sprintf("%.2f", radius) * "_density_0.65.json")
    
    # Set up the simulation with the created MRI sequences and substrate geometry
    sim = Simulation(seqs, diffusivity=2.3, geometry=geom)

    # Readout the signal from intra-axonal, extra-axonal compartments and the whole volume
    # the number of isochromats is set to 250000, feel free to decrease it to speed up simulation or increase it to reduce noise floor. The bounding box defines the region of interest, and subsets define different regions
    sig = readout(250000, sim, bounding_box=BoundingBox([-sz, -sz, -sz], [sz, sz, sz]), subset=[Subset(inside=true), Subset(inside=false), Subset()])

    # Write the resulting signal data to a CSV file
    # The file name is dynamically created based on the perm and radius arguments
    writedlm("./perm/signal_MT_0_sus_0_perm_" * @sprintf("%.3f", perm) * "_rmean_" * @sprintf("%.2f", radius) * "_density_0.65.csv", transverse.(sig), ',')
end

# Run the main function
main()

