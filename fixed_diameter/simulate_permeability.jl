# Run mcmr on identical randomly distributed cylinders
using ArgParse
using MCMRSimulator
using Printf
using DelimitedFiles
using FileIO

function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table! s begin
        "--radius"
        help = "radius"
        "--perm"
        help = "perm"

    end

    return parse_args(s)
end


function main()
    args = parse_commandline()

    # Now you can access the arguments with args["argname"]
    # For example:
    radius = round(parse(Float64, args["radius"]), digits=3)
    perm = round(parse(Float64, args["perm"]), digits=4)

    # Perform operations...
    TR=3000
    TE=166
    delta=2.5
    Gmax=1200
    gamma=0.00004257638476
    seqs = Sequence[]
    for Delta in [10, 15, 20, 30, 40, 50, 60, 80]

        for g in 0:80:Gmax
            seq = dwi(gradient_duration=delta, gradient_strength=g*gamma, diffusion_time=Delta, TE=TE, TR=TR)
	    push!(seqs, seq)
        end
    end
    sz = radius*200
    # bbox = MCMRSimulator.BoundingBox([-sz/2, -sz/2, -500], [sz/2, sz/2, 500])
    geom = MCMRSimulator.read_geometry("cylinders_rep100r/permeability/cylinders_MT_0_sus_0_perm_"*@sprintf("%.3f",perm)*"_rmean_"*@sprintf("%.2f",radius)*"_density_0.65.json")
    sim = Simulation(seqs, diffusivity=3., geometry=geom)
    sig = readout(250000, sim, bounding_box=BoundingBox([-sz, -sz, -sz],[sz,sz,sz]), subset=[Subset(inside=true), Subset(inside=false), Subset()])
    writedlm("./perm/signal_MT_0_sus_0_perm_"*@sprintf("%.3f",perm)*"_rmean_"*@sprintf("%.2f",radius)*"_density_0.65.csv", transverse.(sig), ',')
end

main()
