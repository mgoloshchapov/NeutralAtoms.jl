
"""
    get_fidelity_with_rz_phi(ρ, state, ϕ_rz)

Return the fidelity with a target state after applying a global two-qubit
`RZ(ϕ_rz) ⊗ RZ(ϕ_rz)` correction.

This helper is used when calibrating the phase of the blockade-mediated CZ
sequence.
"""
function get_fidelity_with_rz_phi(ρ, state, ϕ_rz)
    ones = ket_1 ⊗ ket_1
    CZ = Id ⊗ Id - 2*(ones ⊗ dagger(ones));
    state_tr = CZ * state
    global_RZ = ϕ -> RZ(ϕ) ⊗ RZ(ϕ); 
    
    return real(dagger(state_tr) * global_RZ(ϕ_rz) * ρ * dagger( global_RZ(ϕ_rz)) * state_tr)
end

"""
    CZ_calibration_by_fidelity_oscillation(cfg::CZLPConfig; ode_kwargs...)

    Scan a global `RZ` phase and return the Bell-state fidelity curve used to
    calibrate the CZ protocol.

    # Returns
    - `(ϕ_list, F_list, ϕ_opt)`, where `ϕ_opt` maximizes the fidelity proxy.
"""
function PhiPlus_fidelity_osc(cfg::CZLPConfig; ode_kwargs...)
    cfg_parity = deepcopy(cfg)

    ket_pos = (ket_0 + ket_1) / sqrt(2) #ket_ipos = (ket_0 + 1.0im * ket_1) / sqrt(2)
    cfg_parity.ψ0 = ket_pos ⊗ ket_pos #ket_ipos ⊗ ket_pos
    ρ = simulation_czlp(cfg_parity; ode_kwargs...)[1][end]

    Had = Id ⊗ Hadamard # ρ1,2 .=  Had * ρ1,2 * dagger(Had)
    Phi_p = (ket_0 ⊗ ket_0 + ket_1 ⊗ ket_1)/sqrt(2) #Phi_ip = (ket_0 ⊗ ket_0 + 1.0im * ket_1 ⊗ ket_1)/sqrt(2)

    ϕ_list = [0.0:0.0001:2π;];
    #global_RZ = ϕ -> RZ(ϕ) ⊗ RZ(ϕ);
    #F_list = [real(dagger(Phi_p)  * Had * global_RZ(ϕ) * ρ1 * dagger(Had * global_RZ(ϕ)) * Phi_p) for ϕ in ϕ_list];
    F_list = [ get_fidelity_with_rz_phi(ρ, cfg_parity.ψ0, ϕ) for ϕ in ϕ_list]; #Had * Phi_p

    Pr_p = ket_p ⊗ dagger(ket_p)
    Proj = Pr_p ⊗ Id + Id ⊗ Pr_p - Pr_p ⊗ Pr_p 

    #use plot(ϕ_list, F_list) in notebook
    return ϕ_list, F_list, ϕ_list[argmax(F_list)],  real(expect(Proj, ρ))
end

"""
    get_parity_osc(ρ, ϕ_cal)

    Compute the parity oscillation expected after applying the CZ sequence and the
    phase correction `ϕ_cal`.
"""
function get_parity_osc(ρ, ϕ_cal)
    S_ZZ = Z ⊗ Z;

    ϕ_list =  [0.0:0.001:2π;]; #-ϕ1 .+ π/2 .+
    global_RZ = ϕ -> RZ(ϕ) ⊗ RZ(ϕ); #Had * global_RZ(ϕ) * ρ * dagger(Had * global_RZ(ϕ)) #(Phi_p ⊗ dagger(Phi_p))
    global_RX = x -> RX(x) ⊗ RX(x);

    θ = ϕ_cal; # - cfg_parity.ϕ_RZ + ϕ_cal - π
    U = a -> global_RX(π/2) * global_RZ(a) * global_RX(5*π/4)

    Par_list = [real(expect(S_ZZ, U(ϕ) * global_RZ(θ) * ρ * dagger(U(ϕ) * global_RZ(θ)))) for ϕ in ϕ_list];
    #use plot(ϕ_list, Par_list) in notebook 
    return ϕ_list, Par_list 
end 

basis_fidelity_states = [
    ket_0, 
    ket_1,
    (ket_0 + ket_1)/sqrt(2),
    (ket_0 - ket_1)/sqrt(2),
    (ket_0 + 1.0im * ket_1)/sqrt(2),
    (ket_0 - 1.0im * ket_1)/sqrt(2)
    ]

"""
    get_rydberg_fidelity_configs(cfg, n_samples=20)

    Construct a set of derived configurations used for single-atom error-budget
    analysis.

    Each returned configuration isolates one decoherence mechanism or the combined
    error budget.
"""
function get_rydberg_fidelity_configs(cfg, n_samples=20; int_prob=false)
    configs = OrderedDict()
     
    # Config to measure error from temperature
    cfg_t = deepcopy(cfg)
    cfg_t.error_options["spontaneous_decay_intermediate"] = false
    cfg_t.error_options["spontaneous_decay_rydberg"] = false
    cfg_t.error_options["laser_noise"] = false
    cfg_t.error_options["free_motion"] = true
    cfg_t.error_options["atom_motion"] = true
    cfg_t.error_options["xy_motion"] = true 
    cfg_t.error_options["z_motion"] = true
    cfg_t.error_options["Doppler"] = true
    cfg_t.n_samples = n_samples
    configs["Atom motion"] = cfg_t

    # Config to measure error from laser_noise
    cfg_t = deepcopy(cfg)
    cfg_t.error_options["atom_motion"] = false #    cfg_t.atom_params[2] = 0.1
    cfg_t.error_options["xy_motion"] = false
    cfg_t.error_options["z_motion"] = false
    cfg_t.error_options["Doppler"] = false
    cfg_t.error_options["spontaneous_decay_intermediate"] = false
    cfg_t.error_options["spontaneous_decay_rydberg"]      = false
    cfg_t.error_options["laser_noise"] = true
    cfg_t.n_samples = n_samples
    configs["Laser noise"] = cfg_t
    
    # Config to measure error from xyz_motion
    cfg_t = deepcopy(cfg)
    cfg_t.error_options["spontaneous_decay_intermediate"] = false
    cfg_t.error_options["spontaneous_decay_rydberg"] = false
    cfg_t.error_options["laser_noise"] = false
    cfg_t.error_options["free_motion"] = true
    cfg_t.error_options["atom_motion"] = true
    cfg_t.error_options["xy_motion"] = true
    cfg_t.error_options["z_motion"] = true
    cfg_t.error_options["Doppler"] = false
    cfg_t.n_samples = n_samples
    configs["xyz_motion"] = cfg_t
    
    # Config to measure error from Doppler effet
    cfg_t = deepcopy(cfg)
    cfg_t.error_options["spontaneous_decay_intermediate"] = false
    cfg_t.error_options["spontaneous_decay_rydberg"] = false
    cfg_t.error_options["laser_noise"] = false
    cfg_t.error_options["free_motion"] = true
    cfg_t.error_options["atom_motion"] = true
    cfg_t.error_options["xy_motion"] = false
    cfg_t.error_options["z_motion"] = false
    cfg_t.error_options["Doppler"] = true
    cfg_t.n_samples = n_samples
    configs["Doppler"] = cfg_t

    # Config to measure error from intermediate state decay
    cfg_t = deepcopy(cfg)
    cfg_t.error_options["blockade"] = false 
    cfg_t.error_options["atom_motion"] = false #    cfg_t.atom_params[2] = 0.1
    cfg_t.error_options["xy_motion"] = false
    cfg_t.error_options["z_motion"] = false
    cfg_t.error_options["Doppler"] = false
    cfg_t.error_options["spontaneous_decay_intermediate"] = true
    cfg_t.error_options["spontaneous_decay_rydberg"]      = false
    cfg_t.error_options["laser_noise"] = false
    cfg_t.n_samples = 1
    configs["Intermediate state decay"] = cfg_t

    # Config to measure error from rydberg state decay
    cfg_t = deepcopy(cfg)
    cfg_t.error_options["blockade"] = false 
    cfg_t.error_options["atom_motion"] = false #    cfg_t.atom_params[2] = 0.1
    cfg_t.error_options["xy_motion"] = false
    cfg_t.error_options["z_motion"] = false
    cfg_t.error_options["Doppler"] = false
    cfg_t.error_options["spontaneous_decay_intermediate"] = false
    cfg_t.error_options["spontaneous_decay_rydberg"]      = true
    cfg_t.error_options["laser_noise"] = false
    cfg_t.n_samples = 1
    configs["Rydberg state decay"] = cfg_t

    # Config to measure error from xy_motion
    cfg_t = deepcopy(cfg)
    cfg_t.error_options["blockade"] = false 
    cfg_t.error_options["spontaneous_decay_intermediate"] = false
    cfg_t.error_options["spontaneous_decay_rydberg"] = false
    cfg_t.error_options["laser_noise"] = false
    cfg_t.error_options["free_motion"] = true
    cfg_t.error_options["atom_motion"] = true
    cfg_t.error_options["xy_motion"] = true
    cfg_t.error_options["z_motion"] = false
    cfg_t.error_options["Doppler"] = false
    cfg_t.n_samples = n_samples
    configs["xy_motion"] = cfg_t
    
    # Config to measure error from z_motion
    cfg_t = deepcopy(cfg)
    cfg_t.error_options["blockade"] = false 
    cfg_t.error_options["spontaneous_decay_intermediate"] = false
    cfg_t.error_options["spontaneous_decay_rydberg"] = false
    cfg_t.error_options["laser_noise"] = false 
    cfg_t.error_options["free_motion"] = true
    cfg_t.error_options["atom_motion"] = true
    cfg_t.error_options["xy_motion"] = false
    cfg_t.error_options["z_motion"] = true
    cfg_t.error_options["Doppler"] = false
    cfg_t.n_samples = n_samples
    configs["z_motion"] = cfg_t

    # Config to measure blockade error
    cfg_t = deepcopy(cfg)
    cfg_t.error_options["blockade"] = true 
    cfg_t.error_options["spontaneous_decay_intermediate"] = false
    cfg_t.error_options["spontaneous_decay_rydberg"] = false
    cfg_t.error_options["laser_noise"] = false
    cfg_t.error_options["free_motion"] = true
    cfg_t.error_options["atom_motion"] = true
    cfg_t.error_options["xy_motion"] = false
    cfg_t.error_options["z_motion"] = false
    cfg_t.error_options["Doppler"] = false
    cfg_t.n_samples = n_samples
    configs["blockade"] = cfg_t

    # Config to measure total error
    cfg_t = deepcopy(cfg)
    cfg_t.error_options["blockade"] = true 
    cfg_t.error_options["spontaneous_decay_intermediate"] = true
    cfg_t.error_options["spontaneous_decay_rydberg"] = true
    cfg_t.error_options["laser_noise"] = true
    cfg_t.error_options["free_motion"] = true
    cfg_t.error_options["atom_motion"] = true
    cfg_t.error_options["xy_motion"] = true
    cfg_t.error_options["z_motion"] = true
    cfg_t.error_options["Doppler"] = true
    cfg_t.n_samples = n_samples
    configs["Total"] = cfg_t

    # Config to measure intermediate state probability
    if int_prob
        cfg_t = deepcopy(cfg)
        cfg_t.error_options["blockade"] = false 
        cfg_t.error_options["spontaneous_decay_intermediate"] = false
        cfg_t.error_options["spontaneous_decay_rydberg"] = false
        cfg_t.error_options["laser_noise"] = false
        cfg_t.error_options["free_motion"] = false
        cfg_t.error_options["atom_motion"] = false
        cfg_t.error_options["xy_motion"] = false
        cfg_t.error_options["z_motion"] = false
        cfg_t.error_options["Doppler"] = false
        cfg_t.n_samples = 1
        configs["Intermediate probability"] = cfg_t
    end

    return configs
end 

"""
    get_rydberg_infidelity(cfg::RydbergConfig; U=dense(identityoperator(basis)),
        states=basis_fidelity_states, n_samples=100, ode_kwargs...)

    Estimate a single-atom gate infidelity budget by averaging over input states and
    error channels.

    The returned dictionary separates motion, laser noise, spontaneous decay, and
    total error contributions for the effective two-photon model.
"""
function get_rydberg_error_budget(
    cfg::RydbergConfig;
    U=dense(identityoperator(basis)), 
    states=basis_fidelity_states, 
    n_samples=100,
    ode_kwargs...)

    configs = get_rydberg_fidelity_configs(cfg, n_samples)
    names = collect(keys(configs))
    infidelities = Dict()
    
    println("Measuring intermediate state probability...")

    cfg_t = deepcopy(cfg)
    cfg_t.error_options["atom_motion"] = false
    cfg_t.error_options["xy_motion"] = false
    cfg_t.error_options["z_motion"] = false
    cfg_t.error_options["Doppler"] = false
    cfg_t.error_options["spontaneous_decay_intermediate"]    = false
    cfg_t.error_options["spontaneous_decay_rydberg"]         = false
    cfg_t.error_options["laser_noise"]                       = false
    cfg_t.n_samples                         = 1
    T0 = T_twophoton(cfg_t.first_laser_params["Ω"],cfg_t.second_laser_params["Ω"], cfg_t.detuning_params[1])
    cfg_t.tspan = [0.0:T0/50:T0;];
    ψ_ideal = ket_1
    cfg_t.ψ0 = ψ_ideal
    
    ρ_real = simulation(cfg_t)[1][end];
    calibration_error = 1.0 - real(dagger(ψ_ideal) * ρ_real * ψ_ideal) 
    println("Intermediate state probability: $(round(100*calibration_error; digits=4)) %")
    infidelities["Intermediate state probability"] = calibration_error
    
    for name in ProgressBar(names)
        cfg_t = deepcopy(configs[name])
        T0 = T_twophoton(cfg_t.first_laser_params["Ω"],cfg_t.second_laser_params["Ω"], cfg_t.detuning_params[1])
        cfg_t.tspan = [0.0:T0/50:T0;]; #RZ(π)
        println("Measuring error from $(name)...")
        infidelity_avg = 0.0
        for state in [ket_1] #states
            ψ_ideal = state #RZ(π) * state #U * state;
            cfg_t.ψ0 = state
            ρ_real = simulation(cfg_t)[1][end]
            infidelity_avg += 1.0 - real(dagger(ψ_ideal) * ρ_real * ψ_ideal) - calibration_error
        end
        infidelities[name] = infidelity_avg / length(states)

        println()
        println("Infidelity from $(name): $(round(100.0*infidelities[name]; digits=4)) %")
    end

    return infidelities
end

"""
    get_cz_infidelity(cfg::CZLPConfig; n_samples=1, ode_kwargs...)

    Estimate a CZ-gate infidelity budget using the parity-calibration workflow.

    The result separates calibration error from the remaining decoherence channels.
"""
function get_cz_error_budget(
    cfg::CZLPConfig;
    n_samples=1,
    ode_kwargs...)

    configs = get_rydberg_fidelity_configs(cfg, n_samples)
    names = collect(keys(configs))
    infidelities = Dict()

    ket_pos = (ket_0 + ket_1) / sqrt(2)
    Φp = (ket_0 ⊗ ket_0 + ket_1 ⊗ ket_1)/sqrt(2);
    Had = Id ⊗ Hadamard
    
    cfg_t = deepcopy(cfg)
    cfg_t.error_options["atom_motion"] = false
    cfg_t.error_options["xy_motion"] = false
    cfg_t.error_options["z_motion"] = false
    cfg_t.error_options["Doppler"] = false
    cfg_t.error_options["spontaneous_decay_intermediate"]    = false
    cfg_t.error_options["spontaneous_decay_rydberg"]         = false
    cfg_t.error_options["laser_noise"]                       = false
    cfg_t.n_samples                         = 1

    println("Measuring intermediate state probability...")
    phi_list, F_list, ϕ_RZ = CZ_calibration_by_fidelity_oscillation(cfg_t);
    global_RZ = RZ(ϕ_RZ) ⊗ RZ(ϕ_RZ);

    calibration_error = 1.0 - maximum(F_list) #real(dagger(Φp) * ρ_real * Φp)

    println("Intermediate state probability: $(round(100.0*calibration_error; digits=4)) %")
    
    for name in ProgressBar(names)
        cfg_t = deepcopy(configs[name])
        println("Measuring error from $(name)...")

        cfg_t.ψ0 = ket_pos ⊗ ket_pos
        ρ_real = simulation_czlp(cfg_t)[1][end]
        ρ_real .= Had * global_RZ * ρ_real * dagger(Had * global_RZ)
        #infidelities[name] = maximum([(1.0 - real(dagger(Φp) * ρ_real * Φp) - calibration_error), 0.0])
        infidelities[name] = 1.0 - real(dagger(Φp) * ρ_real * Φp) - calibration_error
        #infidelities[name] = 1.0 - maximum(CZ_calibration_by_fidelity_oscillation(cfg_t)[2]);

        println("Infidelity from $(name): $(round(100.0*infidelities[name]; digits=4)) %")
    end

    return infidelities, calibration_error
end
