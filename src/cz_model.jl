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
    get_V(sample1, sample2, center1, center2, ωr, ωz, error_options, c6, eps=1e-18)

Return the blockade interaction `V(t)` between two sampled atoms.

The interaction follows the van der Waals scaling `c6 / R(t)^6`.
"""
@inline function get_V(sample1, sample2, center1, center2, ωr, ωz, error_options, c6, eps=1e-18)
    errs = deepcopy(error_options)
    if errs["blockade"]
        errs["xy_motion"] = true
        errs["z_motion"] = true 
        X1, Y1, Z1 = get_atom_trajectories(sample1, center1, ωr, ωz, errs)[1:3] 
        X2, Y2, Z2 = get_atom_trajectories(sample2, center2, ωr, ωz, errs)[1:3]
        V = t -> (c6 / (eps + ((X1(t) - X2(t))^2 + (Y1(t) - Y2(t))^2 + (Z1(t) - Z2(t))^2)^3))
        return V
    else
        cx1, cy1, cz1 = center1
        cx2, cy2, cz2 = center2
        V = t -> (c6 / (eps + ((cx1-cx2)^2 + (cy1-cy2)^2 + (cz1-cz2)^2)^3))
        return V    
    end
end

"""
    GenerateHamiltonianTwo(sample1, sample2, center1, center2, ωr, ωz,
        tspan_noise, f, nodes, red_laser_phase_amplitudes,
        blue_laser_phase_amplitudes, red_laser_params, blue_laser_params, ϕr,
        ϕb, Δ0, δ0, c6)

Assemble the time-dependent two-atom Hamiltonian for the blockade-mediated CZ
simulation.
"""
@inline function GenerateHamiltonianTwo(
    sample1, sample2,
    center1, center2,
    ωr, ωz,
    error_options,
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
    centers = [center1, center2]

    # Trajectories
    for i in 1:2
        X, Y, Z, Vx, Vy, Vz = get_atom_trajectories(samples[i], centers[i], ωr, ωz, error_options);

        coefficients_two = [coefficients_two; [
            t -> Δ(Vx(t), Vz(t), first_laser_params) - Δ0,
            t -> δ(Vx(t), Vz(t), first_laser_params, second_laser_params) - δ0,
        ]]

        # Hamiltonian params trajectories
        Ω1 = t -> exp(1.0im * (ϕ_1(t) + ϕ_first(t))) * Ω(X(t), Y(t), Z(t), first_laser_params);
        Ω2 = t -> exp(1.0im * (ϕ_2(t) + ϕ_sec(t))) * Ω(X(t), Y(t), Z(t), second_laser_params);

        # Generate phase noise traces for red and blue lasers
        ϕ_red_res  = ϕ(tspan_noise, f, first_laser_phase_amplitudes);
        ϕ_blue_res = ϕ(tspan_noise, f, second_laser_phase_amplitudes);

        # Interpolate phase noise traces to pass to hamiltonian
        ϕ_1  = interpolate(nodes, ϕ_red_res, Gridded(Linear()));
        ϕ_2 = interpolate(nodes, ϕ_blue_res, Gridded(Linear()));
        
        coefficients_two = [coefficients_two; 
            [
                t -> Ω1(t)       / 2.0,
                t -> conj(Ω1(t)) / 2.0,
                t -> Ω2(t)       / 2.0,
                t -> conj(Ω2(t)) / 2.0,
            ]
        ];
    end;
    V = get_V(sample1, sample2, center1, center2, ωr, ωz, error_options, c6)
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

function pure_simulation_czlp(cfg::CZLPConfig;ode_kwargs...)
    # Generate thermal offsets around each atom center.
    samples12 = samples_generate(cfg.trap_params,cfg.atom_params,2;harmonic=true)[1]
    
    # Unpack all parameters
    ωr, ωz = trap_frequencies(cfg.atom_params, cfg.trap_params);
    Δ0, δ0 = cfg.detuning_params;
    τ = cfg.tspan[end] / 2.0;

    tspan_noise = [0.0:cfg.tspan[end]/1000:cfg.tspan[end];];
    nodes = (tspan_noise, );
    first_laser_phase_amplitudes  = cfg.error_options["laser_noise"] ? cfg.first_laser_phase_amplitudes  : zero(cfg.first_laser_phase_amplitudes);
    second_laser_phase_amplitudes = cfg.error_options["laser_noise"] ? cfg.second_laser_phase_amplitudes : zero(cfg.second_laser_phase_amplitudes);
    ϕ_sec = t -> 0.0;
    ϕ_first = t -> t < τ ? 0.0 : cfg.ξ;

    ψ0  = cfg.ψ0 
    #Density matrix averaged over realizations of laser noise and atom dynamics.
    ψt  = [zero(ψ0) for _ in 1:length(cfg.tspan)];

    H = GenerateHamiltonianTwo(samples12[1], samples12[2],
            cfg.atom_centers[1], cfg.atom_centers[2],
            ωr, ωz, cfg.error_options, tspan_noise, cfg.f, nodes,
            first_laser_phase_amplitudes, second_laser_phase_amplitudes,
            cfg.first_laser_params, cfg.second_laser_params,
            ϕ_first, ϕ_sec, Δ0, δ0, cfg.c6)

    #if  !(cfg.error_options["spontaneous_decay_rydberg"]) & !(cfg.error_options["spontaneous_decay_intermediate"])
    ψt = timeevolution.schroedinger_dynamic(cfg.tspan, cfg.ψ0, H; ode_kwargs...)[2];       

    return ψt
end 

function CZ_caliration(cfg::CZLPConfig;ode_kwargs...)
    cfg_CZ = deepcopy(cfg)

    cfg_CZ.ψ0 = (ket_0 + ket_1) ⊗ (ket_0 + ket_1) / 2
    cfg_CZ.n_samples = 1;
    cfg_CZ.atom_params[2] = 0.1; #temperature  
    cfg_CZ.error_options = Dict("laser_noise" => false,"spontaneous_decay_intermediate" => false,"spontaneous_decay_rydberg" => false,
    "atom_motion" => false,"free_motion" => false,"xy_motion" => false,"z_motion" => false,"Doppler" => false, "blockade"=>false);
    
    println("Δ = $(round(cfg_CZ.ΔtoΩ; digits=6)), ξ = $(round(cfg_CZ.ξ; digits=6))")

    ψ = pure_simulation_czlp(cfg_CZ)[end]; 

    basis = [ket_0 ⊗ ket_0, ket_0 ⊗ ket_1, ket_1 ⊗ ket_0, ket_1 ⊗ ket_1]
    state = [dagger(st) * ψ for st in basis] 

    println("Ampls difference: ", abs.(state) .- ones(4)./2)
    println("Phase on state 01: ", angle(state[2]), ", on state 10: ", angle(state[3])) 
    # " ", angle(state[3]), " ", angle(state[4]) - 2 * angle(state[2]))
    #arr = [1+0im, exp(angle(state[2])*1.0im), exp(angle(state[2])*1.0im), -exp(2.0im*angle(state[2]))*0.992] ./ 2 # norm(state - arr)
    println("avg ϕ_RZ = ", -(angle(state[2])+angle(state[3]))/2, "; phase on |11> = ", angle(state[4])-angle(state[2])-angle(state[3]), "; err_norm = ",norm(abs.(state) .- ones(4)./2))
    return -angle(state[2])
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
    # Generate thermal offsets around each atom center.
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
    
    # Unpack all parameters
    ωr, ωz = trap_frequencies(cfg.atom_params, cfg.trap_params);
    Δ0, δ0 = cfg.detuning_params;
    τ = cfg.tspan[end] / 2.0;

    tspan_noise = [0.0:cfg.tspan[end]/1000:cfg.tspan[end];];
    nodes = (tspan_noise, );
    first_laser_phase_amplitudes  = cfg.error_options["laser_noise"] ? cfg.first_laser_phase_amplitudes  : zero(cfg.first_laser_phase_amplitudes);
    second_laser_phase_amplitudes = cfg.error_options["laser_noise"] ? cfg.second_laser_phase_amplitudes : zero(cfg.second_laser_phase_amplitudes);
    ϕ_sec = t -> 0.0;
    ϕ_first = t -> t < τ ? 0.0 : cfg.ξ;

    Γ0, Γ1, Γl   = cfg.error_options["spontaneous_decay_intermediate"] ? cfg.decay_params[1:3] : zeros(3)
    Γr           = cfg.error_options["spontaneous_decay_rydberg"]      ? cfg.decay_params[4]   :  0.0
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
                cfg.atom_centers[1], cfg.atom_centers[2],
                ωr, ωz,
                cfg.error_options,
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
