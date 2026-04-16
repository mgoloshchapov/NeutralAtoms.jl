using NeutralAtoms
using QuantumOptics
using Serialization
using DataFrames
using CSV
using DelimitedFiles 

function _get_6P_config()
       m = 86.9091835;            # Atom params
       T = 1.0 # 70.0;
       atom_params = [m, T];

       U0 = 800 * 0.65; #797.7 changed 26.02.2026
       w0 = 1.42; # Trap params
       λ0 = 0.813;
       M2 = 1.0;
       z0 = w0_to_z0(w0, λ0, M2);
       trap_params = [U0, w0, z0];

       h0 = 13.0 * 1e-6;    #MHz^2/MHz #Params for laser phase noise
       hg1 = 25.0 * 1e-6;   #MHz^2/MHz
       hg2 = 10.0e3 * 1e-6; #MHz^2/MHz
       fg1 = 130.0 * 1e-3;  #MHz
       fg2 = 234.0 * 1e-3;  #MHz
       σg1 = 18.0 * 1e-3;   #MHz
       σg2 = 1.5 * 1e-3;    #MHz
       first_laser_phase_params  = [h0, [hg1, hg2], [σg1, σg2], [fg1, fg2]];
       second_laser_phase_params = [h0, [hg1, hg2], [σg1, σg2], [fg1, fg2]];
       f = [0.01:0.0025:1.0;];
       first_laser_phase_amplitudes = ϕ_amplitudes(f, first_laser_phase_params);
       second_laser_phase_amplitudes = ϕ_amplitudes(f, second_laser_phase_params);

       λr = 1.012;        ### Excitation beam parameters
       λb = 0.42;
       wr = 20.0;
       wb = 8.0 #10.0;
       zr = w0_to_z0(wr, λr);
       zb = w0_to_z0(wb, λb);

       Δ0 = 2.0*π * 2000 #870 #1600.0 #1000.0 #* 1600
       Ω = 2π * 2.3 # 2.0;
       a = 1 #2 #sqrt(1.864) #       #a = 2.0 
       Ωr = a * sqrt(2* Δ0 * Ω);
       Ωb = 1/a * sqrt(2* Δ0 * Ω);
       #Ωr =  2.0*π * 268 #88.4077 #268 #232.6 #88.4077
       #Ωb =  2.0*π * 36 #36 #18.02 #48.0929
       # Orientation and flat-top param, n=1 - gauss
       # θr = π/2;
       θr = π;
       θb = 0.0;
       nr = 1;
       nb = 1;

       first_laser_params = Dict("Ω" => Ωr,"w0" => wr,"z0" => zr,
              "θ" => θr,"n_sg" => nr,"type" => "gauss") #first_laser_params = [Ωr, wr, zr, θr] #, nr];
       second_laser_params = Dict("Ω" => Ωb,"w0" => wb,"z0" => zb,
              "θ" => θb,"n_sg" => nb,"type" => "gauss")      #second_laser_params = [Ωb, wb, zb, θb] #, nb];

       detuning_params = [Δ0, -δ_twophoton(Ωr, Ωb, Δ0)];
       Γ = 2.0*π * 5.75;
       Γ0, Γ1, Γl = Γ/4, Γ/4, 2*Γ/4;
       # Quasiclassical calculations of BBR-induced depopulation rates and effective lifetimes
       # of Rydberg nS, nP, and nD alkali-metal atoms with n ≤ 80. T = 300, n=60, S_1/2, Rb87
       τr = 111.74; 
       Γr = 1/τr;
       decay_params = [Γ0, Γ1, Γl, Γr];

       T0 = T_twophoton(Ωr, Ωb, Δ0) # Simulation params
       # tspan = [0.0:T0/10:2*T0;];
       #tspan = [0.0:T0/20:5*T0;];
       tspan = [0.0:T0/20:3.0;]; #15.0
       ψ0 = ket_1;
       n_samples = 20;
       shift = [0.0,0.0,0.0]

       atom_motion = true;
       free_motion = true;
       laser_noise = false;
       spontaneous_decay_intermediate = true #false;
       spontaneous_decay_rydberg = true #false;
       # spontaneous_decay_intermediate = false;
       # spontaneous_decay_rydberg = false;

       cfg = NeutralAtoms.RydbergConfig(
              tspan,
              ψ0,

              atom_params,
              trap_params,
              n_samples,

              f,
              first_laser_phase_amplitudes,
              second_laser_phase_amplitudes,

              first_laser_params,
              second_laser_params,
              shift,

              detuning_params,
              decay_params,

              atom_motion,
              free_motion,
              laser_noise,
              spontaneous_decay_intermediate,
              spontaneous_decay_rydberg
              );
 
       d = 2.7 #2.0;
       atom_centers = [[-d/2, 0.0, 0.0], [d/2, 0.0, 0.0]]
       #atom_centers = [[0.0,-d/2, 0.0], [0.0, d/2, 0.0]]
       #atom_centers = [[0.0, 0.0,-d/2], [0.0, 0.0,d/2]]
       c6 = 2π * 135298
       
       """Ωτ = 4.278785545408966  # с учетом блокады
       ΔtoΩ = 0.38378019520864026 
       ξ = 3.9162081717218746"""
       ΔtoΩ = 0.377371 #идеальные
       Ωτ = 4.29268
       ξ = 3.90242
       
       ket_pos = (ket_0 + ket_1)/sqrt(2)
       ψ0_cz = ket_pos ⊗ ket_pos

       Ω_twophoton = (2π/T0)
       #τ = 2π / (Ω_twophoton * sqrt(ΔtoΩ^2 + 2.0))
       τ = Ωτ/Ω_twophoton
       VV = c6 / d^6 #3.4^6 #       Ω_twophoton / 2 / V
       #ΔtoΩ = ΔtoΩ + Ω_twophoton/(2*VV)
       detuning_params = [Δ0, ΔtoΩ * Ω_twophoton - δ_twophoton(Ωr, Ωb, Δ0)];
       
       ϕ2 = 2*τ * ΔtoΩ * Ω_twophoton;
       ϕ1 = (ϕ2 - π)/2  
       # tspan_cz = [0.0:τ/10:2*τ;];
       tspan_cz = [0.0, 2*τ];
       #detuning_params = [Δ0, ΔtoΩ*Ω_twophoton + δ_twophoton(Ωr, Ωb, Δ0)];
       
       cfg_czlp = NeutralAtoms.CZLPConfig(
              tspan_cz,
              ψ0_cz,

              atom_params,
              trap_params,
              n_samples,

              f,
              first_laser_phase_amplitudes,
              second_laser_phase_amplitudes,

              first_laser_params,
              second_laser_params,

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

function get_6P_config()
       _, cfg_CZ_6P = _get_6P_config()

       a = cfg_CZ_6P.second_laser_params #["z"] = z;
       cfg_CZ_6P.second_laser_params = cfg_CZ_6P.first_laser_params
       cfg_CZ_6P.first_laser_params = a

       Γ =  8.475 #2.0*π * 5.75;
       Γ0, Γ1, Γl = Γ/4, Γ/4, 2*Γ/4;
       # Quasiclassical calculations of BBR-induced depopulation rates and effective lifetimes
       # of Rydberg nS, nP, and nD alkali-metal atoms with n ≤ 80. T = 300, n=60, S_1/2, Rb87
       τr = 445.74; 
       Γr = 1/τr;
       cfg_CZ_6P.decay_params = [Γ0, Γ1, Γl, Γr];

       return cfg_CZ_6P
end;

function get_5P_config()
       _, cfg_CZ_5P = get_default_config()
       ket_pos = (ket_0 + ket_1) / sqrt(2)
       cfg_CZ_5P.ψ0 = ket_pos ⊗ ket_pos 

       w = 2. ;
       bl_lsr_prms = cfg_CZ_5P.second_laser_params
       Ω0, w0, z0 = bl_lsr_prms["Ω"], bl_lsr_prms["w0"], bl_lsr_prms["z0"];
       z = z0*(w/w0)^2        #cfg_CZ_6P.second_laser_params["type"] = "gauss" 
       cfg_CZ_5P.second_laser_params["w0"] = w ; 
       cfg_CZ_5P.second_laser_params["z"] = z;

       Γ = 2.0*π * 5.75;
       Γ0, Γ1, Γl = Γ/4, Γ/4, 2*Γ/4;
       # Quasiclassical calculations of BBR-induced depopulation rates and effective lifetimes
       # of Rydberg nS, nP, and nD alkali-metal atoms with n ≤ 80. T = 300, n=60, S_1/2, Rb87
       τr = 111.74; 
       Γr = 1/τr;
       cfg_CZ_5P.decay_params = [Γ0, Γ1, Γl, Γr];

       return cfg_CZ_5P
end;