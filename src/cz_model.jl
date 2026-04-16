#Jump operators for master equation 
@inline function JumpOperatorsTwo(decay_params)
    Γ0, Γ1, Γl, Γr = decay_params;
    operators = [
        sqrt(Γ0)*σ0p ⊗ Id,  sqrt(Γ1)*σ1p ⊗ Id, sqrt(Γl)*σlp ⊗ Id, sqrt(Γr)*σlr ⊗ Id,
        sqrt(Γ0)*Id  ⊗ σ0p, sqrt(Γ1)*Id  ⊗ σ1p, sqrt(Γl)*Id ⊗ σlp, sqrt(Γr)*Id ⊗ σlr
        ]
    return operators
end;


@inline function get_V(sample1, sample2,  ωr, ωz, atom_motion, free_motion, c6, eps=1e-18)
    X1, Y1, Z1 = get_atom_trajectories(sample1, ωr, ωz, atom_motion, free_motion)[1:3]
    X2, Y2, Z2 = get_atom_trajectories(sample2, ωr, ωz, atom_motion, free_motion)[1:3]
    V = t -> (c6 / (eps + ((X1(t) - X2(t))^2 + (Y1(t) - Y2(t))^2 + (Z1(t) - Z2(t))^2)^3))
    return V
end


@inline function GenerateHamiltonianTwo(
    sample1, sample2,
    ωr, ωz,
    free_motion, atom_motion,
    tspan_noise, f, nodes,
    first_laser_phase_amplitudes, second_laser_phase_amplitudes,
    first_laser_params, second_laser_params,
    ϕ_first, ϕ_sec,
    Δ0, δ0,
    c6)

    operators_two = [operators .⊗ [Id]; [Id] .⊗  operators; nr ⊗ nr];
    coefficients_two = Vector{Function}();
    # coefficients_two = [];
    samples = [sample1, sample2]

    # Trajectories
    for i in 1:2
        X, Y, Z, Vx, Vy, Vz = get_atom_trajectories(samples[i], ωr, ωz, atom_motion, free_motion);

        # Generate phase noise traces for red and blue lasers
        ϕ_red_res  = ϕ(tspan_noise, f, first_laser_phase_amplitudes);
        ϕ_blue_res = ϕ(tspan_noise, f, second_laser_phase_amplitudes);

        # Interpolate phase noise traces to pass to hamiltonian
        ϕ_1  = interpolate(nodes, ϕ_red_res, Gridded(Linear()));
        ϕ_2 = interpolate(nodes, ϕ_blue_res, Gridded(Linear()));


        # Hamiltonian params trajectories
        Ω1 = t -> exp(1.0im * (ϕ_1(t) + ϕ_first(t))) * Ω(X(t), Y(t), Z(t), first_laser_params );
        Ω2 = t -> exp(1.0im * (ϕ_2(t) + ϕ_sec(t))) * Ω(X(t), Y(t), Z(t), second_laser_params);

        coefficients_two = [coefficients_two; 
            [
                t -> Δ(Vx(t), Vz(t), first_laser_params) - Δ0,
                t -> δ(Vx(t), Vz(t), first_laser_params, second_laser_params) - δ0,
                t -> Ω1(t)       / 2.0,
                t -> conj(Ω1(t)) / 2.0,
                t -> Ω2(t)       / 2.0,
                t -> conj(Ω2(t)) / 2.0,
            ]
        ];
    end;
    V = get_V(sample1, sample2, ωr, ωz, atom_motion, free_motion, c6)
    push!(coefficients_two, V)
    H = TimeDependentSum(coefficients_two, operators_two);

    return H
end;


"""
To compensate for non-ideal blockade we have to correct ΔtoΩ as ΔtoΩ - Ω/2V

Here we average V over atom temperature and get 
"""
function get_blockade_stark_shift_factor(
    trap_params,
    atom_params,
    atom_centers,
    Ω,
    c6,
    n_samples=10000
    )
    shift1, shift2 = [[atom_centers[1];zeros(3)]], [[atom_centers[2];zeros(3)]]
    samples1 = samples_generate(
        trap_params,
        atom_params,
        n_samples;
        harmonic=true
        )[1]
    samples2 = samples_generate(
        trap_params,
        atom_params,
        n_samples;
        harmonic=true
        )[1]
    samples1 .+= shift1
    samples2 .+= shift2

    Rm6 = mean(map((s1, s2) -> 1.0 / (1e-18 + sum((s1[1:3] - s2[1:3]).^2)^3), samples1, samples2))

    return - Ω / (2.0 * c6 * Rm6)
end

function simulation_czlp(
    cfg::CZLPConfig;
    ode_kwargs...)
    # Generate samples of atoms and shift their centers
    shift1, shift2 = [[cfg.atom_centers[1];zeros(3)]], [[cfg.atom_centers[2];zeros(3)]]
    samples1 = samples_generate(
        cfg.trap_params,
        cfg.atom_params,
        cfg.n_samples;
        harmonic=true
        )[1]
    samples2 = samples_generate(
        cfg.trap_params,
        cfg.atom_params,
        cfg.n_samples;
        harmonic=true
        )[1]
    samples1 .+= shift1
    samples2 .+= shift2
    
    # Unpack all parameters
    ωr, ωz = trap_frequencies(cfg.atom_params, cfg.trap_params);
    Δ0, δ0 = cfg.detuning_params;
    τ = cfg.tspan[end] / 2.0;

    tspan_noise = [0.0:cfg.tspan[end]/1000:cfg.tspan[end];];
    nodes = (tspan_noise, );
    first_laser_phase_amplitudes  = cfg.laser_noise ? cfg.first_laser_phase_amplitudes  : zero(cfg.first_laser_phase_amplitudes);
    second_laser_phase_amplitudes = cfg.laser_noise ? cfg.second_laser_phase_amplitudes : zero(cfg.second_laser_phase_amplitudes);
    ϕ_sec = t -> 0.0;
    ϕ_first = t -> t < τ ? 0.0 : cfg.ξ;

    Γ0, Γ1, Γl   = cfg.spontaneous_decay_intermediate ? cfg.decay_params[1:3] : zeros(3)
    Γr           = cfg.spontaneous_decay_rydberg      ? cfg.decay_params[4]   :  0.0
    decay_params = [Γ0, Γ1, Γl, Γr]
    J = JumpOperatorsTwo(decay_params)

    ρ0  = cfg.ψ0 ⊗ dagger(cfg.ψ0);
    #Density matrix averaged over realizations of laser noise and atom dynamics.
    ρ   = [zero(ρ0) for _ in 1:length(cfg.tspan)];
    ρt  = [zero(ρ0) for _ in 1:length(cfg.tspan)];
    #Second moment for error estimation of level populations. 
    ρ2  = [zero(ρ0) for _ in 1:length(cfg.tspan)];

    for i in ProgressBars.ProgressBar(1:cfg.n_samples)
       H = GenerateHamiltonianTwo(
                samples1[i], samples2[i],
                ωr, ωz,
                cfg.free_motion, cfg.atom_motion,
                tspan_noise, cfg.f, nodes,
                first_laser_phase_amplitudes, second_laser_phase_amplitudes,
                cfg.first_laser_params, cfg.second_laser_params,
                ϕ_first, ϕ_sec,
                Δ0, δ0,
                cfg.c6)

        ρt = timeevolution.master_dynamic(cfg.tspan, ρ0, H, J; ode_kwargs...)[2];

        ρ  .+= ρt
        ρ2 .+= ρt .^ 2
    end;

    return ρ ./ cfg.n_samples, ρ2 ./ cfg.n_samples;
end