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
    cfg_t.atom_params[2] = 1.0
    cfg_t.spontaneous_decay_intermediate = true
    cfg_t.spontaneous_decay_rydberg      = false
    cfg_t.laser_noise = false
    cfg_t.free_motion = true
    cfg_t.n_samples = 1
    configs["Intermdeiate state decay"] = cfg_t

    # Config to measure error from rydberg state decay
    cfg_t = deepcopy(cfg)
    cfg_t.atom_params[2] = 1.0
    cfg_t.spontaneous_decay_intermediate = false
    cfg_t.spontaneous_decay_rydberg      = true
    cfg_t.laser_noise = false
    cfg_t.free_motion = true
    cfg_t.n_samples = 1
    configs["Rydberg state decay"] = cfg_t

    # Config to measure error from laser_noise
    cfg_t = deepcopy(cfg)
    cfg_t.atom_params[2] = 1.0
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


function plot_rydberg_infidelity(
    infidelities; 
    dir_name="/Users/goloshch/ColdAtoms_test/experiments/23_07_2025/results/", 
    file_name="plot.png",
    title="Error budget for 2π pulse",
    blue=true)

    red_color  = RGBA(207.0/255, 71.0/255, 80.0/255, 1.0)
    blue_color = RGBA(135.0/255,203.0/255,230.0/255,1.0)
    if blue
        color = blue_color 
    else 
        color = red_color
    end

    keys_iF   = collect(keys(infidelities))
    keys_ordered = [
        "Total", 
        "Atom motion", 
        "Intermdeiate state decay", 
        "Rydberg state decay",
        "Laser noise"
        ]

    keys_final = [key for key in keys_ordered if key in keys_iF]
    values_final = [100*infidelities[key] for key in keys_ordered if key in keys_iF]

    p = bar(keys_final, values_final;
    xrotation=45,
    margin=10Plots.mm,
    ylabel="Infidelity, %", 
    title=title,
    label=nothing,
    dpi=300,
    size=(600, 600),
    xguidefontsize = 14,
    yguidefontsize = 14,
    xtickfontsize = 14,
    ytickfontsize = 14,
    color=color
    )

    savefig("$(dir_name)$(file_name)")

    display(p)

    return keys_final, values_final
end

function get_parity_fidelity(cfg::CZLPConfig; ode_kwargs...)
    cfg_parity = deepcopy(cfg)

    ket_pos = (ket_0 + ket_1) / sqrt(2)
    ket_ipos = (ket_0 + 1.0im * ket_1) / sqrt(2)


    cfg_parity.ψ0 = ket_pos ⊗ ket_pos
    ρ1 = simulation_czlp(cfg_parity; ode_kwargs...)[1][end]
    # cfg_parity.ψ0 = ket_ipos ⊗ ket_pos
    # ρ2 = simulation_czlp(cfg_parity; ode_kwargs...)[1][end]

    Had = Id ⊗ Hadamard
    # ρ1 .=  Had * ρ1 * dagger(Had)
    # ρ2 .=  Had * ρ2 * dagger(Had)

    Phi_p = (ket_0 ⊗ ket_0 + ket_1 ⊗ ket_1)/sqrt(2)
    Phi_ip = (ket_0 ⊗ ket_0 + 1.0im * ket_1 ⊗ ket_1)/sqrt(2)

    ϕ_list = [0.0:0.0001:2π;];
    global_RZ = ϕ -> RZ(ϕ) ⊗ RZ(ϕ);

    F_list_1 = [real(dagger(Phi_p)  * Had * global_RZ(ϕ) * ρ1 * dagger(Had * global_RZ(ϕ)) * Phi_p) for ϕ in ϕ_list];
    # F_list_2 = [real(dagger(Phi_ip) * Had * global_RZ(ϕ) * ρ2 * dagger(Had * global_RZ(ϕ)) * Phi_ip) for ϕ in ϕ_list];

    # plot(ϕ_list, [F_list_1, F_list_2])
    plot(ϕ_list, F_list_1)
    return ϕ_list, F_list_1, ϕ_list[argmax(F_list_1)]
    # return ϕ_list, F_list_1, F_list_2, ϕ_list[argmax(F_list_1)]
end

function get_parity(cfg::CZLPConfig, ϕ_cal; ode_kwargs...)
    cfg_parity = deepcopy(cfg)
    ket_ineg = (ket_0 - 1.0im * ket_1) / sqrt(2) #ket_pos = (ket_0 + ket_1) / sqrt(2)
    S_ZZ = Z ⊗ Z;

    cfg_parity.ψ0 = ket_ineg ⊗ ket_ineg #ket_pos ⊗ ket_pos
    ρ1 = simulation_czlp(cfg_parity; ode_kwargs...)[1][end] #    Had = Id ⊗ Hadamard

    ϕ_list =  [0.0:0.001:2π;]; #-ϕ1 .+ π/2 .+
    global_RZ = ϕ -> RZ(ϕ) ⊗ RZ(ϕ); #Had * global_RZ(ϕ) * ρ1 * dagger(Had * global_RZ(ϕ)) #(Phi_p ⊗ dagger(Phi_p))
    global_RX = x -> RX(x) ⊗ RX(x);
    θ = - cfg_parity.ϕ_RZ + ϕ_cal - π
    U = a -> global_RX(π/2) * global_RZ(a) * global_RX(5*π/4)

    Par_list = [real(expect(S_ZZ , U(ϕ) * global_RZ(θ) * ρ1 * dagger(U(ϕ) * global_RZ(θ)) ) ) for ϕ in ϕ_list];
    plot(ϕ_list, Par_list)
    return ϕ_list, Par_list #, ϕ_list[argmax(Par_list)]
end
#unused
function get_parity_fidelity_temp(ρ, ϕ_RZ)
    Had = Id ⊗ Hadamard

    global_RZ = RZ(ϕ_RZ) ⊗ RZ(ϕ_RZ);
    ρt = global_RZ * ρ * dagger(global_RZ);
    ρt = Had * ρt * dagger(Had)

    Phi_p = (ket_0 ⊗ ket_0 + ket_1 ⊗ ket_1)/sqrt(2)

    F = real(dagger(Phi_p) * ρt * Phi_p)
    return F, ρt
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
    cfg_t.atom_params[2]                    = 1.0
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


function plot_cz_infidelity(infidelities;
    dir_name="/Users/goloshch/ColdAtoms_test/experiments/23_07_2025/results/", 
    file_name="plot_cz.png",
    title="Error budget for 2π pulse")

    keys_iF   = collect(keys(infidelities))
    keys_ordered = [
        "Total", 
        "Atom motion", 
        "Intermdeiate state decay", 
        "Rydberg state decay",
        "Laser noise"
        ]

    keys_final = [key for key in keys_ordered if key in keys_iF]
    values_final = [100*infidelities[key] for key in keys_ordered if key in keys_iF]

    p = bar(keys_final, values_final;
    xrotation=45,
    margin=10Plots.mm,
    ylabel="Infidelity, %", 
    title=title,
    label=nothing,
    dpi=300,
    size=(600, 600),
    xguidefontsize = 14,
    yguidefontsize = 14,
    xtickfontsize = 14,
    ytickfontsize = 14,
    color=RGBA(135.0/255,203.0/255,230.0/255,1.0)
    )

    savefig("$(dir_name)$(file_name)")

    display(p)

    return keys_final, values_final
end