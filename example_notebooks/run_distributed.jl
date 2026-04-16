using Distributed 
using NeutralAtoms  

num_proc = 8
addprocs(num_proc)   
@everywhere using NeutralAtoms
@everywhere include("../conf/default.jl")

@everywhere function compute_ρ(i)
    _, cfg_CZ = get_default_configs();
    c_xy = get_ideal_cxy() #get_default_cxy()

    cfg_CZ.blue_laser_params["type"] = "flattop_HG" 
    cfg_CZ.blue_laser_params["coeffs_xy"] = c_xy; 
    cfg_CZ.spontaneous_decay_intermediate = false; #true
    cfg_CZ.spontaneous_decay_rydberg = false;
    cfg_CZ.laser_noise = false;

    T = 20.0 #cfg_CZ.atom_motion = false;
    cfg_CZ.atom_params = [86.9092, T] 
    cfg_CZ.n_samples = 1 #00 #50 # 50
    
    ρ_end = NeutralAtoms.simulation_czlp(cfg_CZ)[1][end]
    println(NeutralAtoms.get_fidelity_max(ρ_end, cfg_CZ.ψ0))
    return ρ_end 
end 

function main()
    #n = nworkers()
    ρ_computed = pmap(compute_ρ, 1:num_proc)

    ρ = sum(ρ_computed) ./ num_proc

    _, cfg_CZ = get_default_configs();
    println(NeutralAtoms.get_fidelity_max(ρ, cfg_CZ.ψ0))
    #println(ρ_computed)
end

main()
# Optionally, remove the worker processes after use
rmprocs(workers())