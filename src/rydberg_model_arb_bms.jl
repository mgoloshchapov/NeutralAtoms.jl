function simulation_one_shift(shift, cfg_t::RydbergConfig; ode_kwargs...    )
    samples = samples_generate(
        cfg_t.trap_params,
        cfg_t.atom_params,
        cfg_t.n_samples;
        harmonic=true
        )[1]

    ωr, ωz = trap_frequencies(cfg_t.atom_params, cfg_t.trap_params);
    Δ0, δ0 = cfg_t.detuning_params;

    tspan_noise = [0.0:cfg_t.tspan[end]/1000:cfg_t.tspan[end];];
    nodes = (tspan_noise, );
    red_laser_phase_amplitudes  = cfg_t.red_laser_phase_amplitudes;
    blue_laser_phase_amplitudes = cfg_t.blue_laser_phase_amplitudes;

    Γ0, Γ1, Γl   = cfg_t.spontaneous_decay_intermediate ? cfg_t.decay_params[1:3] : zeros(3)
    Γr           = cfg_t.spontaneous_decay_rydberg      ? cfg_t.decay_params[4]   :  0.0
    decay_params = [Γ0, Γ1, Γl, Γr]
    J = JumpOperators(decay_params)

    ρ0 = cfg_t.ψ0 ⊗ dagger(cfg_t.ψ0);
    #Density matrix averaged over realizations of laser noise and atom dynamics.
    ρ  = [zero(ρ0) for _ ∈ 1:length(cfg_t.tspan)];
    ρt  = [zero(ρ0) for _ ∈ 1:length(cfg_t.tspan)];
    #Second moment for error estimation of level populations. 
    ρ2 = [zero(ρ0) for _ ∈ 1:length(cfg_t.tspan)];

    function __simulation(sample)
        H = GenerateHamiltonian(
            sample .+ vcat(shift, [0.,0.,0.]), 
            ωr, ωz,
            cfg_t.free_motion,
            cfg_t.atom_motion,
            cfg_t.laser_noise,
        
            tspan_noise,
            cfg_t.f,
            red_laser_phase_amplitudes,
            blue_laser_phase_amplitudes,
            nodes,
        
            cfg_t.red_laser_params,
            cfg_t.blue_laser_params,
        
            Δ0, 
            δ0
            )

        ρt .= timeevolution.master_dynamic(cfg_t.tspan, ρ0, H, J; ode_kwargs...)[2];
        ρ  .+= ρt
        ρ2 .+= ρt .^ 2
    end

    for sample in ProgressBars.ProgressBar(samples)
        __simulation(sample)
    end;

    return ρ ./ cfg_t.n_samples, ρ2 ./ cfg_t.n_samples
end;
