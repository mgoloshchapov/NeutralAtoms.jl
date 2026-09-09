using NeutralAtoms
using QuantumOptics
using Serialization

function cfg_grad(ΔtoΩ, ξ)
       _, cfg_ = _get_6P_config(ΔtoΩ, ξ)
       cfg_.error_options = Dict("laser_noise" => false,"spontaneous_decay_intermediate" => false,"spontaneous_decay_rydberg" => false,
       "atom_motion" => false,"free_motion" => false,"xy_motion" => false,"z_motion" => false,"Doppler" => false)
       cfg_.atom_params[2]=0.1
       return cfg_
end

function get_6P_config()
       ΔtoΩ = 0.37737
       ξ = 3.90242   
       ξ = 3.982 #3.955 #3.6
       return _get_6P_config(ΔtoΩ, ξ)
end

function _get_6P_config(ΔtoΩ, ξ)
       m = 86.9091835;            # Atom params
       T = 70.0 # 70.0;
       atom_params = [m, T];
       # Trap params
       U0 = 800.0*0.65; #max * эффективность использования
       w0 = 1.42; 
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

       ### Excitation beam parameters
       λ1 = 0.42;
       λ2 = 1.012;
       w1 = 15.0 #8.0;
       w2 = 20.0;
       z1 = w0_to_z0(w1, λ1);
       z2 = w0_to_z0(w2, λ2);

       Δ0 = 2.0*π * 2000.0  # 2000 
       Ω = 2π * 2.36  # 2.38
       a = 1 #*sqrt(7) # ξ = 3.9564
       Ω1 = a * sqrt(2* Δ0 * Ω);
       Ω2 = 1/a * sqrt(2* Δ0 * Ω);

       # Orientation and flat-top param, n=1 - gauss
       θ1 = π/2 #0.0;
       θ2 = -π/2 #π;
       n1 = 1;
       n2 = 1;

       first_laser_params = Dict("Ω" => Ω1,"w0" => w1,"z0" => z1,
              "θ" => θ1,"n_sg" => n1,"type" => "gauss") #first_laser_params = [Ωr, wr, zr, θr] #, nr];
       second_laser_params = Dict("Ω" => Ω2,"w0" => w2,"z0" => z2,
              "θ" => θ2,"n_sg" => n2,"type" => "gauss")      #second_laser_params = [Ωb, wb, zb, θb] #, nb];

       detuning_params = [Δ0, -δ_twophoton(Ω1, Ω2, Δ0)];
       Γ = 8.475;
       Γ0, Γ1, Γl = Γ/4, Γ/4, 2*Γ/4; #уточнить
       # Quasiclassical calculations of BBR-induced depopulation rates and effective lifetimes
       # of Rydberg nS, nP, and nD alkali-metal atoms with n ≤ 80. T = 300, n=60, S_1/2, Rb87
       τr = 445.0; #from arc getStateLifetime: n3=74, l3=0, j3=0.5 
       Γr = 1/τr;
       decay_params = [Γ0, Γ1, Γl, Γr];
       
       # Simulation params
       T0 = 2.5 #T_twophoton(Ω1, Ω2, Δ0)
       tspan = [0.0:T0/200:T0;];
       #tspan = [0.0:T0/20:5.5*T0;];
       ψ0 = ket_1;
       n_samples = 20;
       shift = [0.0,0.0,0.0]
       error_options = Dict("laser_noise" => false,
                        "blockade"=>false,
                        "spontaneous_decay_intermediate" => true,
                        "spontaneous_decay_rydberg" => true,
                        "atom_motion" => true,
                        "free_motion" => true,
                        "xy_motion" => true,
                        "z_motion" => true,
                        "Doppler" => true
                        )

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

              error_options
              );
 
       d = 3.6;

       #atom_centers = [[0.0, 0.0, -0.5*d], [0.0, 0.0, 0.5*d]]
       atom_centers = [[-0.5*d, 0.0, 0.0 ], [0.5*d, 0.0, 0.0]]

       c6 = 2π * 1.6e6  # 135298 # для n=74 1.6e6  #
       #ΔtoΩ = 0.377371
       Ωτ = 4.29268
       #ξ = 3.9564 # 3.90242
       ket_pos = (ket_0 + ket_1)/sqrt(2)
       ψ0_cz = ket_pos ⊗ ket_pos
       
       Ω_twophoton = (2π/T_twophoton(Ω1, Ω2, Δ0))
       τ = 2π / (Ω_twophoton * sqrt(ΔtoΩ^2 + 2.0))
       VV = c6 / d^6 #3.4^6 #       Ω_twophoton / 2 / V
       #ΔtoΩ = ΔtoΩ + Ω_twophoton/(2*VV)
       detuning_params = [Δ0, ΔtoΩ * Ω_twophoton - δ_twophoton(Ω1, Ω2, Δ0)];
       
       ϕ2 = 2*τ * ΔtoΩ * Ω_twophoton;
       ϕ1 = (ϕ2 - π)/2  
       ϕ_RZ = 1.9053
       tspan_cz = [0.0, 2*τ];
       
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

              error_options,

              atom_centers,
              c6,
              ΔtoΩ,
              Ωτ,
              ξ,
              ϕ_RZ
       )

       return cfg, cfg_czlp
end