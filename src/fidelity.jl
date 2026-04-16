function get_fidelity_with_rz_phi(ρ, state, ϕ_rz)
    ones = ket_1 ⊗ ket_1
    CZ = Id ⊗ Id - 2*(ones ⊗ dagger(ones));
    state_tr = CZ * state
    global_RZ = ϕ -> RZ(ϕ) ⊗ RZ(ϕ); 
    
    return real(dagger(state_tr) * global_RZ(ϕ_rz) * ρ * dagger( global_RZ(ϕ_rz)) * state_tr)
end

function CZ_calibration_by_fidelity_oscillation(cfg::CZLPConfig; ode_kwargs...)
    cfg_parity = deepcopy(cfg)

    ket_pos = (ket_0 + ket_1) / sqrt(2) #ket_ipos = (ket_0 + 1.0im * ket_1) / sqrt(2)
    cfg_parity.ψ0 = ket_pos ⊗ ket_pos #ket_ipos ⊗ ket_pos
    ρ = simulation_czlp(cfg_parity; ode_kwargs...)[1][end]

    Had = Id ⊗ Hadamard # ρ1,2 .=  Had * ρ1,2 * dagger(Had)
    Phi_p = (ket_0 ⊗ ket_0 + ket_1 ⊗ ket_1)/sqrt(2) #Phi_ip = (ket_0 ⊗ ket_0 + 1.0im * ket_1 ⊗ ket_1)/sqrt(2)

    ϕ_list = [0.0:0.0001:2π;];
    #global_RZ = ϕ -> RZ(ϕ) ⊗ RZ(ϕ);
    #F_list = [real(dagger(Phi_p)  * Had * global_RZ(ϕ) * ρ1 * dagger(Had * global_RZ(ϕ)) * Phi_p) for ϕ in ϕ_list];
    F_list = [ get_fidelity_with_rz_phi(ρ, Had * Phi_p, ϕ) for ϕ in ϕ_list];

    #use plot(ϕ_list, F_list) in notebook
    return ϕ_list, F_list, ϕ_list[argmax(F_list)]
end

function get_parity_osc(ρ, ϕ_cal)
    S_ZZ = Z ⊗ Z;

    ϕ_list =  [0.0:0.001:2π;]; #-ϕ1 .+ π/2 .+
    global_RZ = ϕ -> RZ(ϕ) ⊗ RZ(ϕ); #Had * global_RZ(ϕ) * ρ * dagger(Had * global_RZ(ϕ)) #(Phi_p ⊗ dagger(Phi_p))
    global_RX = x -> RX(x) ⊗ RX(x);

    θ = ϕ_cal; # - cfg_parity.ϕ_RZ + ϕ_cal - π
    U = a -> global_RX(π/2) * global_RZ(a) * global_RX(5*π/4)

    Par_list = [real(expect(S_ZZ , U(ϕ) * (θ) * ρ * dagger(U(ϕ) * global_RZ(θ)) ) ) for ϕ in ϕ_list];
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

function get_rydberg_fidelity_configs(cfg, n_samples=20)
    configs = OrderedDict()

    # Config to measure error from intermediate state decay
    cfg_t = deepcopy(cfg)
    cfg_t.atom_params[2] = 0.1
    cfg_t.spontaneous_decay_intermediate = true
    cfg_t.spontaneous_decay_rydberg      = false
    cfg_t.laser_noise = false
    cfg_t.free_motion = true
    cfg_t.n_samples = 1
    configs["Intermdeiate state decay"] = cfg_t

    # Config to measure error from rydberg state decay
    cfg_t = deepcopy(cfg)
    cfg_t.atom_params[2] = 0.1
    cfg_t.spontaneous_decay_intermediate = false
    cfg_t.spontaneous_decay_rydberg      = true
    cfg_t.laser_noise = false
    cfg_t.free_motion = true
    cfg_t.n_samples = 1
    configs["Rydberg state decay"] = cfg_t

    # Config to measure error from laser_noise
    cfg_t = deepcopy(cfg)
    cfg_t.atom_params[2] = 0.1
    cfg_t.spontaneous_decay_intermediate = false
    cfg_t.spontaneous_decay_rydberg = false
    cfg_t.laser_noise = true
    cfg_t.free_motion = false
    cfg_t.n_samples = n_samples
    configs["Laser noise"] = cfg_t

    # Config to measure error from temperature
    cfg_t = deepcopy(cfg)
    cfg_t.spontaneous_decay_intermediate = false
    cfg_t.spontaneous_decay_rydberg = false
    cfg_t.laser_noise = false
    cfg_t.free_motion = true
    cfg_t.n_samples = n_samples
    configs["Atom motion"] = cfg_t

    # Config to measure total error
    cfg_t = deepcopy(cfg)
    cfg_t.spontaneous_decay_intermediate = true
    cfg_t.spontaneous_decay_rydberg = true
    cfg_t.laser_noise = true
    cfg_t.free_motion = true
    cfg_t.n_samples = n_samples
    configs["Total"] = cfg_t

    return configs
end 

function get_rydberg_infidelity(
    cfg::RydbergConfig;
    U=dense(identityoperator(basis)), 
    states=basis_fidelity_states, 
    n_samples=100,
    ode_kwargs...)

    configs = get_rydberg_fidelity_configs(cfg, n_samples)
    names = collect(keys(configs))
    infidelities = Dict()

    for name in ProgressBar(names)
        cfg_t = deepcopy(configs[name])
        println("Measuring error from $(name)...")
        infidelity_avg = 0.0
        for state in states
            ψ_ideal = U * state;
            cfg_t.ψ0 = state
            ρ_real = simulation(cfg_t)[1][end]
            infidelity_avg += 1.0 - real(dagger(ψ_ideal) * ρ_real * ψ_ideal)
        end
        infidelities[name] = infidelity_avg / length(states)

        println()
        println("Infidelity from $(name): $(round(100.0*infidelities[name]; digits=4)) %")
    end

    return infidelities
end

function get_cz_infidelity(
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
    cfg_t.spontaneous_decay_intermediate    = false
    cfg_t.spontaneous_decay_rydberg         = false
    cfg_t.laser_noise                       = false
    cfg_t.atom_params[2]                    = 0.1
    cfg_t.n_samples                         = 1

    println("Measuring error from calibration...")
    ϕ_RZ = get_parity_fidelity(cfg_t)[3];
    global_RZ = RZ(ϕ_RZ) ⊗ RZ(ϕ_RZ);

    cfg_t.ψ0 = ket_pos ⊗ ket_pos
    ρ_real = simulation_czlp(cfg_t)[1][end]
    ρ_real .= Had * global_RZ * ρ_real * dagger(Had * global_RZ)
    calibration_error = 1.0 - real(dagger(Φp) * ρ_real * Φp)

    println()
    println("Infidelity from calibration error: $(round(100.0*calibration_error; digits=4)) %")

    
    for name in ProgressBar(names)
        cfg_t = deepcopy(configs[name])
        println("Measuring error from $(name)...")

        cfg_t.ψ0 = ket_pos ⊗ ket_pos
        ρ_real = simulation_czlp(cfg_t)[1][end]
        ρ_real .= Had * global_RZ * ρ_real * dagger(Had * global_RZ)
        infidelities[name] = maximum([(1.0 - real(dagger(Φp) * ρ_real * Φp) - calibration_error), 0.0])

        println()
        println("Infidelity from $(name): $(round(100.0*infidelities[name]; digits=4)) %")
    end

    return infidelities, calibration_error
end
