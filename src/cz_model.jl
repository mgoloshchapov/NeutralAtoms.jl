"""
    JumpOperatorsTwo(decay_params)

Construct the Lindblad jump operators for the two-atom CZ model.
"""
@inline function JumpOperatorsTwo(decay_params)
    Γ0, Γ1, Γl, Γr = decay_params;
    operators = [
        sqrt(Γ0)*σ0p ⊗ Id,  sqrt(Γ1)*σ1p ⊗ Id, sqrt(Γl)*σlp ⊗ Id, sqrt(Γr)*σlr ⊗ Id,
        sqrt(Γ0)*Id  ⊗ σ0p, sqrt(Γ1)*Id  ⊗ σ1p, sqrt(Γl)*Id ⊗ σlp, sqrt(Γr)*Id ⊗ σlr
        ]
    return operators
end;


"""
    get_V(sample1, sample2, ωr, ωz, atom_motion, free_motion, c6, eps=1e-18)

Return the blockade interaction `V(t)` between two sampled atoms.

The interaction follows the van der Waals scaling `c6 / R(t)^6`.
"""
@inline function get_V(sample1, sample2,  ωr, ωz, atom_motion, free_motion, c6, eps=1e-18)
    X1, Y1, Z1 = get_atom_trajectories(sample1, ωr, ωz, atom_motion, free_motion)[1:3]
    X2, Y2, Z2 = get_atom_trajectories(sample2, ωr, ωz, atom_motion, free_motion)[1:3]
    V = t -> (c6 / (eps + ((X1(t) - X2(t))^2 + (Y1(t) - Y2(t))^2 + (Z1(t) - Z2(t))^2)^3))
    return V
end


"""
    GenerateHamiltonianTwo(sample1, sample2, ωr, ωz, free_motion, atom_motion,
        tspan_noise, f, nodes, red_laser_phase_amplitudes,
        blue_laser_phase_amplitudes, red_laser_params, blue_laser_params, ϕr,
        ϕb, Δ0, δ0, c6)

Assemble the time-dependent two-atom Hamiltonian for the blockade-mediated CZ
simulation.
"""
@inline function GenerateHamiltonianTwo(
    sample1, sample2,
    ωr, ωz,
    free_motion, atom_motion,
    tspan_noise, f, nodes,
    red_laser_phase_amplitudes, blue_laser_phase_amplitudes,
    red_laser_params, blue_laser_params,
    ϕr, ϕb,
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
        ϕ_red_res  = ϕ(tspan_noise, f, red_laser_phase_amplitudes);
        ϕ_blue_res = ϕ(tspan_noise, f, blue_laser_phase_amplitudes);

        # Interpolate phase noise traces to pass to hamiltonian
        ϕ_red  = interpolate(nodes, ϕ_red_res, Gridded(Linear()));
        ϕ_blue = interpolate(nodes, ϕ_blue_res, Gridded(Linear()));


        # Hamiltonian params trajectories
        Ωr = t -> exp(1.0im * (ϕ_red(t)  + ϕr(t))) * Ω(X(t), Y(t), Z(t), red_laser_params );
        Ωb = t -> exp(1.0im * (ϕ_blue(t) + ϕb(t))) * Ω(X(t), Y(t), Z(t), blue_laser_params);

        coefficients_two = [coefficients_two; 
            [
                t -> Δ(Vx(t), Vz(t), red_laser_params) - Δ0,
                t -> δ(Vx(t), Vz(t), red_laser_params, blue_laser_params) - δ0,
                t -> Ωr(t)       / 2.0,
                t -> conj(Ωr(t)) / 2.0,
                t -> Ωb(t)       / 2.0,
                t -> conj(Ωb(t)) / 2.0,
            ]
        ];
    end;
    V = get_V(sample1, sample2, ωr, ωz, atom_motion, free_motion, c6)
    push!(coefficients_two, V)
    H = TimeDependentSum(coefficients_two, operators_two);

    return H
end;


"""
    get_blockade_stark_shift_factor(trap_params, atom_params, atom_centers, Ω,
        c6, n_samples=10000)

Estimate the finite-temperature blockade correction used when calibrating the
CZ pulse.

This helper averages the inverse sixth power of the atom separation over thermal
sampling and returns the resulting Stark-shift factor.
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
   

"""
    simulation_czlp(cfg::CZLPConfig; ode_kwargs...)

Simulate the two-atom global-pulse controlled-phase protocol.

The model follows the blockade-based CZ logic highlighted in
[arXiv:1908.06101](https://arxiv.org/abs/1908.06101): two global Rydberg pulses,
an inter-pulse phase step `ξ`, finite-temperature motion, and optional laser
noise and spontaneous decay.

# Arguments
- `cfg::CZLPConfig`: two-atom phase-gate configuration.

# Keywords
- `ode_kwargs...`: keyword arguments forwarded to
  `timeevolution.master_dynamic`.

# Returns
- `(ρ, ρ2)`, the first and second moments of the two-atom density-matrix
  trajectory.
"""
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
    red_laser_phase_amplitudes  = cfg.laser_noise ? cfg.red_laser_phase_amplitudes  : zero(cfg.red_laser_phase_amplitudes);
    blue_laser_phase_amplitudes = cfg.laser_noise ? cfg.blue_laser_phase_amplitudes : zero(cfg.blue_laser_phase_amplitudes);
    ϕb = t -> 0.0;
    ϕr = t -> t < τ ? 0.0 : cfg.ξ;

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
                red_laser_phase_amplitudes, blue_laser_phase_amplitudes,
                cfg.red_laser_params, cfg.blue_laser_params,
                ϕr, ϕb,
                Δ0, δ0,
                cfg.c6)

        ρt = timeevolution.master_dynamic(cfg.tspan, ρ0, H, J; ode_kwargs...)[2];

        ρ  .+= ρt
        ρ2 .+= ρt .^ 2
    end;


    return ρ ./ cfg.n_samples, ρ2 ./ cfg.n_samples;
end
