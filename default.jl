using NeutralAtoms
using QuantumOptics
using Serialization


function get_default_configs()
       # Atom params
       m = 86.9091835;     
       T = 100.0;
       atom_params = [m, T];

       # Trap params
       U0 = 1000.0;
       w0 = 1.0;
       λ0 = 0.813;
       M2 = 1.0;
       z0 = w0_to_z0(w0, λ0, M2);
       trap_params = [U0, w0, z0];

       #Params for laser phase noise
       h0 = 13.0 * 1e-6;    #MHz^2/MHz
       hg1 = 25.0 * 1e-6;   #MHz^2/MHz
       hg2 = 10.0e3 * 1e-6; #MHz^2/MHz
       fg1 = 130.0 * 1e-3;  #MHz
       fg2 = 234.0 * 1e-3;  #MHz
       σg1 = 18.0 * 1e-3;   #MHz
       σg2 = 1.5 * 1e-3;    #MHz


       red_laser_phase_params  = [h0, [hg1, hg2], [σg1, σg2], [fg1, fg2]];
       blue_laser_phase_params = [h0, [hg1, hg2], [σg1, σg2], [fg1, fg2]];
       f = [0.01:0.0025:1.0;];
       red_laser_phase_amplitudes = ϕ_amplitudes(f, red_laser_phase_params);
       blue_laser_phase_amplitudes = ϕ_amplitudes(f, blue_laser_phase_params);


       ### Excitation beam parameters
       λr = 0.795;
       λb = 0.475;
       wr = 50.0;
       wb = 10.0;
       zr = w0_to_z0(wr, λr);
       zb = w0_to_z0(wb, λb);

       Δ0 = 2.0*π * 1000.0;
       Ω = 2π * 2.0;
       Ωr = sqrt(2* Δ0 * Ω);
       Ωb = sqrt(2* Δ0 * Ω);
       # Orientation and flat-top param, n=1 - gauss
       # θr = π/2;
       θr = π/2;
       θb = 0.0;
       nr = 1;
       nb = 1;

       red_laser_params = [Ωr, wr, zr, θr, nr];
       blue_laser_params = [Ωb, wb, zb, θb, nb];

       detuning_params = [Δ0, δ_twophoton(Ωr, Ωb, Δ0)];
       Γ = 2.0*π * 5.75;
       Γ0, Γ1, Γl = Γ/4, Γ/4, 2*Γ/4;
       # Quasiclassical calculations of BBR-induced depopulation rates and effective lifetimes
       # of Rydberg nS, nP, and nD alkali-metal atoms with n ≤ 80. T = 300, n=60, S_1/2, Rb87
       τr = 111.74; 
       Γr = 1/τr;
       decay_params = [Γ0, Γ1, Γl, Γr];

       # Simulation params
       T0 = T_twophoton(Ωr, Ωb, Δ0)
       # tspan = [0.0:T0/10:2*T0;];
       tspan = [0.0:T0/20:5*T0;];
       ψ0 = ket_1;
       n_samples = 20;

       atom_motion = true;
       free_motion = true;
       laser_noise = false;
       spontaneous_decay_intermediate = true;
       spontaneous_decay_rydberg = true;
       # spontaneous_decay_intermediate = false;
       # spontaneous_decay_rydberg = false;

       cfg = NeutralAtoms.RydbergConfig(
              tspan,
              ψ0,

              atom_params,
              trap_params,
              n_samples,

              f,
              red_laser_phase_amplitudes,
              blue_laser_phase_amplitudes,

              red_laser_params,
              blue_laser_params,

              detuning_params,
              decay_params,

              atom_motion,
              free_motion,
              laser_noise,
              spontaneous_decay_intermediate,
              spontaneous_decay_rydberg
              );

       d = 2.0;
       atom_centers = [[-d/2, 0.0, 0.0], [d/2, 0.0, 0.0]]
       c6 = 2π * 135298
       ΔtoΩ = 0.377371
       Ωτ = 4.29268
       ξ = 3.90242
       ket_pos = (ket_0 + ket_1)/sqrt(2)
       ψ0_cz = ket_pos ⊗ ket_pos
       Ω_twophoton = (2π/T0)
       τ = 2π / (Ω_twophoton * sqrt(ΔtoΩ^2 + 2.0))
       ϕ2 = 2*τ * ΔtoΩ * Ω_twophoton;
       ϕ1 = (ϕ2 - π)/2
       # tspan_cz = [0.0:τ/10:2*τ;];
       tspan_cz = [0.0, 2*τ];


       cfg_czlp = NeutralAtoms.CZLPConfig(
              tspan_cz,
              ψ0_cz,

              atom_params,
              trap_params,
              n_samples,

              f,
              red_laser_phase_amplitudes,
              blue_laser_phase_amplitudes,

              red_laser_params,
              blue_laser_params,

              detuning_params,
              decay_params,

              atom_motion,
              free_motion,
              laser_noise,
              spontaneous_decay_intermediate,
              spontaneous_decay_rydberg,

              atom_centers,
              c6,
              ΔtoΩ,
              Ωτ,
              ξ,
              ϕ1
       )

       return cfg, cfg_czlp
end