#Basis states
const basis = NLevelBasis(5);
const ket_0 = nlevelstate(basis, 1);
const ket_1 = nlevelstate(basis, 2);
const ket_r = nlevelstate(basis, 3);
const ket_p = nlevelstate(basis, 4);
const ket_l = nlevelstate(basis, 5);

#Operators
const Id  = dense(identityoperator(basis));
const σ0p = ket_0 ⊗ dagger(ket_p);
const σp0 = ket_p ⊗ dagger(ket_0);
const σ1p = ket_1 ⊗ dagger(ket_p);
const σp1 = ket_p ⊗ dagger(ket_1);
const σpr = ket_p ⊗ dagger(ket_r);
const σrp = ket_r ⊗ dagger(ket_p);
const np  = ket_p ⊗ dagger(ket_p);
const nr  = ket_r ⊗ dagger(ket_r);
const σlp = ket_l ⊗ dagger(ket_p);
const σpl = ket_p ⊗ dagger(ket_l);

const σlr = ket_l ⊗ dagger(ket_r);
const σrl = ket_r ⊗ dagger(ket_l);

const operators = [np, nr, σ1p, σp1, σpr, σrp];


@inline function get_atom_trajectories(
    sample::Vector{Float64}, 
    ωr::Float64,
    ωz::Float64,
    atom_motion::Bool,
    free_motion::Bool
    )
    x, y, z, vx, vy, vz = atom_motion ? sample : zeros(Float64, 6);
    X  = t -> R(t, x, vx, ωr; free=free_motion);
    Y  = t -> R(t, y, vy, ωr; free=free_motion);
    Z  = t -> R(t, z, vz, ωz; free=free_motion);
    Vx = t -> V(t, x, vx, ωr; free=free_motion);
    Vy = t -> V(t, y, vy, ωr; free=free_motion);
    Vz = t -> V(t, z, vz, ωz; free=free_motion);
    return X, Y, Z, Vx, Vy, Vz;
end

@inline function gauss_field(x, y, z, w0, z0; n0=1, θ0=0)
    return A(x, y, z, w0, z0; n=n0, θ=θ0) .* A_phase(x, y, z, w0, z0; θ=θ0)
end;

@inline function Ω(x, y, z, laser_params) 
    if (laser_params["type"] == "gauss")
        w0 = laser_params["w0"]
        z0 = laser_params["z0"]
        θ = laser_params["θ"]
        n_sg = laser_params["n_sg"]
        return laser_params["Ω"] .* gauss_field(x, y, z, w0, z0; n0=n_sg, θ0=θ) #A(x, y, z, w0, z0; n=n_sg, θ=θ) .* A_phase(x, y, z, w0, z0; θ=θ);
     
    elseif (laser_params["type"] == "flattop")
        θ = laser_params["θ"]
        xx, zz = x*cos(θ) - z*sin(θ), x*sin(θ) + z*cos(θ)
        Ω_w_z = [laser_params["Ω"], laser_params["w0"],laser_params["z0"]] 
        return reconstruct_HG_field_2d(xx, y, zz, Ω_w_z, laser_params["coeffs_xy"]) #, laser_params["coeffs_xy"])
    
    elseif (laser_params["type"] == "2gauss")
        θ = laser_params["θ"]

        x1, y1, z1 = laser_params["beams_centers"][1]
        E1 = laser_params["Ω"] .* gauss_field(x.-x1, y.-y1, z.-z1, laser_params["w0"], laser_params["z0"]; θ0=θ)
        
        x2, y2, z2 = laser_params["beams_centers"][2]
        E2 = laser_params["Ω1"] .* gauss_field(x.-x2, y.-y2, z.-z2, laser_params["w1"], laser_params["z1"]; θ0=θ)

        return E1 .+ exp(1.0im * laser_params["rel_phase"]) .* E2
    else
        throw(error("Unsupported type of beam"))
    end;
        #elseif (laser_params["type"] == "flattop")
        
end;

#Due to Doppler shift for red laser
@inline function Δ(vx, vz, laser_params)
    w0, z0, θ =  laser_params["w0"], laser_params["z0"], laser_params["θ"] #laser_params[2:4]
    k = 2.0 * z0/w0^2;

    Δx = k * sin(θ) * vx
    Δz = k * cos(θ) * vz
    return Δx + Δz
end;

#Due to Doppler shifts for red and blue lasers
@inline function δ(vx, vz, red_laser_params, blue_laser_params)
    wr0, zr0, θr = red_laser_params["w0"], red_laser_params["z0"], red_laser_params["θ"]  #red_laser_params[2:4];
    wb0, zb0, θb = blue_laser_params["w0"], blue_laser_params["z0"], blue_laser_params["θ"]  #blue_laser_params[2:4];
    
    kr = 2.0 * zr0/wr0^2;
    kb = 2.0 * zb0/wb0^2;

    δx = (kr*sin(θr) + kb*sin(θb))*vx
    δz = (kr*cos(θr) + kb*cos(θb))*vz 

    return δx + δz
end;

#Jump operators for master equation 
@inline function JumpOperators(decay_params)
    Γ0, Γ1, Γl, Γr = decay_params;
    decay_operators = [sqrt(Γ0)*σ0p, sqrt(Γ1)*σ1p, sqrt(Γl)*σlp, sqrt(Γr)*σlr]
    return decay_operators
end;

@inline function GenerateHamiltonian(
    sample, 
    ωr, ωz,
    free_motion,
    atom_motion,
    laser_noise,

    tspan_noise,
    f,
    red_laser_phase_amplitudes,
    blue_laser_phase_amplitudes,
    nodes,

    red_laser_params,
    blue_laser_params,
    
    Δ0, 
    δ0
    )
    # Trajectories
    X, Y, Z, Vx, Vy, Vz = get_atom_trajectories(
        sample, 
        ωr, ωz, 
        atom_motion, 
        free_motion);

    

    # Interpolate phase noise traces to pass to hamiltonian
    if laser_noise
        # Generate phase noise traces for red and blue lasers
        ϕ_red_res  = ϕ(tspan_noise, f, red_laser_phase_amplitudes);
        ϕ_blue_res = ϕ(tspan_noise, f, blue_laser_phase_amplitudes);

        ϕ_red  = interpolate(nodes, ϕ_red_res, Gridded(Linear()));
        ϕ_blue = interpolate(nodes, ϕ_blue_res, Gridded(Linear()));
    else
        ϕ_red  = t -> 0.0;
        ϕ_blue = t -> 0.0;
    end

    # Hamiltonian params trajectories
    Ωr = t -> exp(1.0im * ϕ_red(t)) * Ω(X(t), Y(t), Z(t), red_laser_params);
    Ωb = t -> exp(1.0im * ϕ_blue(t)) * Ω(X(t), Y(t), Z(t), blue_laser_params);

    H = TimeDependentSum(
        [
            t -> -Δ(Vx(t), Vz(t), red_laser_params) - Δ0,
            t -> -δ(Vx(t), Vz(t), red_laser_params, blue_laser_params) - δ0,
            t -> Ωr(t) / 2.0,
            t -> conj(Ωr(t)) / 2.0,
            t -> Ωb(t) / 2.0,
            t -> conj(Ωb(t)) / 2.0,
        ],
        operators
        );

    return H
end;


function simulation(
    cfg::RydbergConfig; 
    temperature_calibrate=false, 
    ode_kwargs...
    )

    if temperature_calibrate
        cfg_t = calibrate_two_photon(cfg)
    else
        cfg_t = deepcopy(cfg)
    end

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
            sample, 
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

function Ω_twophoton(Ωr, Ωb, Δ)
    return abs(Ωb * Ωr / (2.0 * Δ))
end;

function T_twophoton(Ωr, Ωb, Δ)
    return 2.0*π / Ω_twophoton(Ωr, Ωb, Δ)
end;

function δ_twophoton(Ωr, Ωb, Δ)
    return (abs(Ωr)^2 - abs(Ωb)^2)/(4.0 * Δ)
end;

function Ωr_required(Ω, Ωb, Δ)
    return 2.0 * Δ * Ω / abs(Ωb)
end;


function calibrate_two_photon(cfg::RydbergConfig, n_samples=1000)
    cfg_calibrated = deepcopy(cfg)
    Ωr = cfg.red_laser_params[1]
    Ωb = cfg.blue_laser_params[1]

    samples = samples_generate(cfg.trap_params, cfg.atom_params, n_samples)[1]

    # Correct single-photon Rabi frequencies to match temperature averaged two-photon Rabi frequency
    samples_Ω2 = [Ω_twophoton(Ω(x, y, z, cfg_calibrated.red_laser_params), Ω(x, y, z, cfg_calibrated.blue_laser_params), Δ(vx, vz, cfg_calibrated.red_laser_params)) for (x, y, z, vx, _, vz) in samples]
    factor_Ω2 = sqrt((sum(samples_Ω2) / length(samples)) / Ω_twophoton(Ωr, Ωb, cfg.detuning_params[1]))
    Ωr_cor, Ωb_cor = Ωr / factor_Ω2, Ωb / factor_Ω2
    cfg_calibrated.red_laser_params[1]  = Ωr_cor
    cfg_calibrated.blue_laser_params[1] = Ωb_cor

    # Correct resonance detuning to match temperature averaged AC Stark shifts
    samples_δ = [δ_twophoton(Ω(x, y, z, cfg_calibrated.red_laser_params), Ω(x, y, z, cfg_calibrated.blue_laser_params), Δ(vx, vz, cfg_calibrated.red_laser_params)) for (x, y, z, vx, _, vz) in samples]
    δ_cor = sum(samples_δ) / length(samples_δ)
    δ_ideal = δ_twophoton(Ωr, Ωb, cfg.detuning_params[1])
    cfg_calibrated.detuning_params[2] += δ_cor - δ_ideal

    return cfg_calibrated
end