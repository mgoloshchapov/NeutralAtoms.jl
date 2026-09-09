using Distributed 
using NeutralAtoms  

include("../conf/default.jl")  
#include("../conf/configs_gen.jl") 
cfg_rydb, cfg_CZ = get_default_configs() # get_6P_config(); #

cfg = deepcopy(cfg_CZ) # deepcopy(cfg_CZ)

cfg.n_samples = 20
configs = get_rydberg_fidelity_configs(cfg, cfg.n_samples; int_prob=true) 
considered_keys = ["Total","blockade","z_motion","xy_motion","Doppler","Intermediate probability","Intermediate state decay","Rydberg state decay"]
configs_list = [[key, configs[key]] for key in considered_keys ]; # configs.keys]; # 
 
num_proc = length(configs_list)
addprocs(num_proc)   
@everywhere using NeutralAtoms
@everywhere include("../conf/default.jl") 
#@everywhere include("../conf/configs_gen.jl") 

@everywhere function compute_rydb(cnfg)    #for rydberg oscillations
    cfg = deepcopy(cnfg[2])
    cfg.ψ0 = ket_1
    T0 = T_twophoton(cfg.first_laser_params["Ω"], cfg.second_laser_params["Ω"], cfg.detuning_params[1])
    cfg.tspan = [0.0:T0/50:T0;]; 

    ρ_real = simulation(cfg)[1][end]
    fid = real(dagger(ket_1) * ρ_real * ket_1) #- calibration_error
    return [cnfg[1], 1.0 - fid]
end 

@everywhere function compute_CZ(cnfg) #for CZ
    cfg = deepcopy(cnfg[2])
    
    cfg.ψ0 = (ket_0 + ket_1) ⊗ (ket_0 + ket_1) / 2;
    ρ = NeutralAtoms.simulation_czlp(cfg)[1][end]
    Fids = [NeutralAtoms.get_fidelity_with_rz_phi(ρ, cfg.ψ0, ϕ) for ϕ in 0:0.001:2π] #    fid = NeutralAtoms.get_fidelity_with_rz_phi(ρ, cfg.ψ0, cfg.ϕ_RZ)
    fid = maximum(Fids)
    return [cnfg[1], 1.0 - fid]
end 

function main()
    err_computed = pmap(compute_CZ, configs_list) #  pmap(compute_CZ, configs_list)
    
    budget = Dict()
    for err in err_computed
        budget[err[1]] = err[2]
    end

    for key in keys(budget)
        println(key, " ", budget[key])
    end

    NeutralAtoms.save_with_JLD2("data\\err_bdgt.jld2", budget)
end

main()
# Optionally, remove the worker processes after use
rmprocs(workers())