using Distributed 
using NeutralAtoms  
using JLD2

num_proc = 7
addprocs(num_proc)   
@everywhere using NeutralAtoms
#@everywhere include("../conf/configs_gen.jl") # include("../conf/config6P_2026_04.jl") 
@everywhere include("../conf/default.jl")

@everywhere function compute_ρ(i)
    _, cfg_CZ = get_default_configs() #get_6P_config();
    
    cfg_CZ.n_samples = 1 #40
    cfg_CZ.error_options = Dict("laser_noise" => false,"spontaneous_decay_intermediate" => false,
    "spontaneous_decay_rydberg" => false,"atom_motion" => false,"free_motion" => false,
    "xy_motion" => false,"z_motion" => false,"Doppler" => true,"blockade" => false)

    ρ_end = NeutralAtoms.simulation_czlp(cfg_CZ)[1][end]
    return ρ_end 
end 

function main()
    ρ_computed = pmap(compute_ρ, 1:num_proc)

    ρ = sum(ρ_computed) ./ num_proc
    
    NeutralAtoms.save_with_JLD2("data\\rho_distributed.jld2", ρ)

    include("../conf/default.jl")
    _, cfg_CZ = get_default_configs() #get_6P_config();
    Fids = [NeutralAtoms.get_fidelity_with_rz_phi(ρ, cfg_CZ.ψ0, ϕ) for ϕ in 0:0.001:2π] 
    println(maximum(Fids))
end

main()
# Optionally, remove the worker processes after use
rmprocs(workers())