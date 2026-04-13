"""
    Ω_red(laser_params)

Return the on-axis red-laser coupling used by the legacy arbitrary-beam model.
"""
function Ω_red(laser_params)
    Ω0, w0, z0 = laser_params;
    return Ω0 
end;

"""
    Ω_blue(x, y, z, laser_params)

Return the blue-laser coupling for the legacy arbitrary-beam model.

The beam shape is selected from `laser_params` and dispatched to either
`simple_flattopLG_field` or `simple_flattopHG_field`.
"""
function Ω_blue(x,y,z, laser_params)
    Ω0, w0, z0, beam_type, n, m = laser_params;
    if beam_type == "simp_flattop_LG"
        return simple_flattopLG_field(x,y,z,laser_params)
    elseif beam_type == "simp_flattop_HG"
        return simple_flattopHG_field(x,y,z,laser_params)
    end;
end;

"""
    simulation_blue_intens(tspan, ψ0, atom_params, trap_params, samples,
        red_laser_params, blue_laser_params, detuning_params, decay_params;
        atom_motion=true, free_motion=true, spontaneous_decay=true,
        parallel=false)

Simulate a single-atom excitation sequence with an explicitly shaped blue beam.

This is an older, specialized path for beam-profile studies. For the main
single-atom workflow, prefer `simulation` with `RydbergConfig`.
"""
function simulation_blue_intens(
    tspan, ψ0,  
    atom_params,    trap_params,
    samples,
    red_laser_params,    blue_laser_params,
    detuning_params,    decay_params;
    atom_motion = true,    free_motion = true,   
    spontaneous_decay = true, parallel = false
    )

    N = length(samples);
    ωr, ωz = trap_frequencies(atom_params, trap_params);
    Δ0, δ0 = detuning_params;

    if spontaneous_decay
        decay_params_temp = decay_params;
    else
        decay_params_temp = [0.0, 0.0];
    end;

    Γg, Γgt = decay_params_temp;
    J, Jdagger = [sqrt(Γg)*σgp, sqrt(Γgt)*σgtp], [sqrt(Γg)*σpg, sqrt(Γgt)*σpgt];

    ρ0 = ψ0 ⊗ dagger(ψ0);

    #Density matrix averaged over realizations of laser noise and atom dynamics.
    ρ_mean = [zero(ψ0 ⊗ dagger(ψ0)) for _ ∈ 1:length(tspan)];
    ρ_temp = [zero(ψ0 ⊗ dagger(ψ0)) for _ ∈ 1:length(tspan)];

    #Second moment for error estimation of level populations. 
    ρ2_mean = [zero(ψ0 ⊗ dagger(ψ0)) for _ ∈ 1:length(tspan)];
    
    for i ∈ 1:N
        if atom_motion
            #Atom initial conditions
            xi, yi, zi, vxi, vyi, vzi = samples[i];
        else
            xi, yi, zi, vxi, vyi, vzi = zeros(6);
        end;
        
        #Atom trajectories
        X = t -> R(t, xi, vxi, ωr; free=free_motion);
        Y = t -> R(t, yi, vyi, ωr; free=free_motion);
        Z = t -> R(t, zi, vzi, ωz; free=free_motion);
        Vz = t -> V(t, zi, vzi, ωz; free=free_motion);

        #Hamiltonian params trajectories
        Ht = TimeDependentSum(
        [
            t -> -Δ(Vz(t), red_laser_params) - Δ0;
            t -> -δ(Vz(t), red_laser_params, blue_laser_params[1:3]; parallel=parallel) - δ0;
            t ->  Ω_red(red_laser_params) / 2.0;
            t ->  conj(Ω_red(red_laser_params) / 2.0);
            t -> Ω_blue(X(t), Y(t), Z(t), blue_laser_params) / 2.0;
            t -> conj(Ω_blue(X(t), Y(t), Z(t), blue_laser_params) / 2.0);
        ],
        operators
        );

        # #Returns hamiltonian and jump operators in a form required by timeevolution.master_dynamic
        function super_operator(t, rho)
            return Ht, J, Jdagger;
        end;
        
        _, ρ_temp = timeevolution.master_dynamic(tspan, ρ0, super_operator);

        ρ_mean = ρ_mean + ρ_temp;
        ρ2_mean = ρ2_mean + ρ_temp .^ 2;
    end;

    return ρ_mean/N, ρ2_mean/N
end;
