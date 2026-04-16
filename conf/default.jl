using NeutralAtoms
using QuantumOptics
using Serialization
using DataFrames
using CSV
using DelimitedFiles 


function get_default_configs()
       # Atom params
       m = 86.9091835;     
       T = 100.0 # 70.0;
       atom_params = [m, T];

       # Trap params
       #U0 = 1000.0; 
       U0 = 800 * 0.65; #797.7 changed 26.02.2026
       #w0 = 1.0;
       w0 = 1.42; 
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


       first_laser_phase_params  = [h0, [hg1, hg2], [σg1, σg2], [fg1, fg2]];
       second_laser_phase_params = [h0, [hg1, hg2], [σg1, σg2], [fg1, fg2]];
       f = [0.01:0.0025:1.0;];
       first_laser_phase_amplitudes = ϕ_amplitudes(f, first_laser_phase_params);
       second_laser_phase_amplitudes = ϕ_amplitudes(f, second_laser_phase_params);


       ### Excitation beam parameters
       λr = 0.795;
       λb = 0.475;
       wr = 50.0;
       wb = 1.2 #2.0 #10.0;
       zr = w0_to_z0(wr, λr);
       zb = w0_to_z0(wb, λb);

       Δ0 = 2.0*π * 1711 #870 #1600.0 #1000.0 #* 1600
       Ω = 2π * 2.41 # 2.0;
       a = 2 #sqrt(1.864) #       #a = 2.0 
       Ωr = a * sqrt(2* Δ0 * Ω);
       Ωb = 1/a * sqrt(2* Δ0 * Ω);
       #Ωr =  2.0*π * 268 #88.4077 #268 #232.6 #88.4077
       #Ωb =  2.0*π * 36 #36 #18.02 #48.0929
       # Orientation and flat-top param, n=1 - gauss
       # θr = π/2;
       θr = π/2;
       θb = 0.0;
       nr = 1;
       nb = 1;

       #first_laser_params = [Ωr, wr, zr, θr] #, nr];
       #second_laser_params = [Ωb, wb, zb, θb] #, nb];
       first_laser_params = Dict("Ω" => Ωr,"w0" => wr,"z0" => zr,
              "θ" => θr,"n_sg" => nr,"type" => "gauss")
       second_laser_params = Dict("Ω" => Ωb,"w0" => wb,"z0" => zb,
              "θ" => θb,"n_sg" => nb,"type" => "gauss")

       detuning_params = [Δ0, -δ_twophoton(Ωr, Ωb, Δ0)];
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
       #tspan = [0.0:T0/20:5*T0;];
       tspan = [0.0:T0/20:3.0;]; #15.0
       ψ0 = ket_1;
       n_samples = 20;

       atom_motion = true;
       free_motion = true;
       laser_noise = false;
       spontaneous_decay_intermediate = true #false;
       spontaneous_decay_rydberg = true #false;
       # spontaneous_decay_intermediate = false;
       # spontaneous_decay_rydberg = false;
       err_optns = Dict("laser_noise" => false,
                        "spontaneous_decay_intermediate" => true,
                        "spontaneous_decay_rydberg" => true,
                        "atom_motion" => true,
                        "free_motion" => true,
                        "xy_osc" => true,
                        "z_osc" => true,
                        "Doppler" => true
                        )
       shift = [0.0,0.0,0.0]

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
 
       d = 3.4 #2.0;
       #atom_centers = [[-d/2, 0.0, 0.0], [d/2, 0.0, 0.0]]
       atom_centers = [[0.0,-d/2, 0.0], [0.0, d/2, 0.0]]
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

function get_default_cxy()
       df = DataFrame(readdlm("flattop_data.dat"), :auto);
       corrX = -92.1
       corrY = -88.7
       x = df.x2 .+ corrX
       y = df.x4 .+ corrY
       z = 1 .- df.x6;

       xx = round.(3.8*x, digits=2);
       yy = round.(3.8*y, digits=2); #yy[1]
       x_span = maximum(xx)-minimum(xx);
       y_span = maximum(yy)-minimum(yy);
       z0 = sqrt.(z) .- 0.168 #0.166
       z_z = (z0) ./ maximum(z0)
       zz = z_z * 1.05 ; #1.02;
       len_x, len_y = 0,0
       if (xx[1] == xx[31])
       len_x = 30
       len_y = 33
       end;
       n_max = 20
       m_max = n_max
       dx = x_span/len_x
       dy = y_span/len_y
       x_0 = [-x_span/2 : x_span/(len_x-1) : x_span/2;];
       y_0 = [-y_span/2 : y_span/(len_y-1) : y_span/2;];

       w = 2. ;
       cxy = NeutralAtoms.decomposition_2d(xx, yy, zz, w,dx,dy);
       c_xy = 1.05 * cxy; # renormalization
       return c_xy
end

function get_ideal_cxy()
       N = 18
       M = 2
       c_xy = 0.0 .* zeros(M+1,N+1) #NeutralAtoms.HG_coefficients(M+1,N+1)
       n=trunc(Int,N/2)
       m=trunc(Int,M/2)
       for j in 0:n
              for i in 0:m
                     c_xy[2*i+1,2*j+1] = NeutralAtoms.HG_coeff_big(j,n) * NeutralAtoms.HG_coeff_big(i,m) #for i in 1:(N+1)
              end
       end
       return c_xy 
end

function get_default_cfg2gauss(w)
       function corr_om(ϕ, d,ww,zz)
              E = gauss_field(0.0,0.0,0.0, ww, zz) + exp(1.0im*ϕ) * gauss_field(0.0,d,0.0, ww, zz)
              return 1/abs(E)
       end;
       _, cfg_CZ_2gauss = get_default_configs();
       ϕ0 =  0.0 #π/2+
       bl_lsr_prms = cfg_CZ_2gauss.second_laser_params
       Ω, w0, z0 = bl_lsr_prms["Ω"], bl_lsr_prms["w0"], bl_lsr_prms["z0"];
       z = z0*(w0/w)^2
       cfg_CZ_2gauss.second_laser_params["type"] = "2gauss" 
       cfg_CZ_2gauss.second_laser_params["rel_phase"] = ϕ0
       cfg_CZ_2gauss.second_laser_params["Ω"] = Ω*corr_om(ϕ0,1.7+1.7,w,z)
       cfg_CZ_2gauss.second_laser_params["Ω1"] = Ω*corr_om(ϕ0,-1.7-1.7,w,z)
       cfg_CZ_2gauss.second_laser_params["w1"] = cfg_CZ_2gauss.second_laser_params["w0"] = w ; 
       cfg_CZ_2gauss.second_laser_params["z1"] = cfg_CZ_2gauss.second_laser_params["z"] = z
       cfg_CZ_2gauss.second_laser_params["beams_centers"] = [[0.0,1.7,0.0],[0.0,-1.7,0.0]] ;
              
       T = 100.1 #100.0;
       cfg_CZ_2gauss.atom_params = [86.9091835, T]
       cfg_CZ_2gauss.spontaneous_decay_intermediate = false #true; #
       cfg_CZ_2gauss.spontaneous_decay_rydberg = false #true; #
       ket_pos = (ket_0 + ket_1) / sqrt(2)
       cfg_CZ_2gauss.ψ0 = ket_pos ⊗ ket_pos

       return cfg_CZ_2gauss
end