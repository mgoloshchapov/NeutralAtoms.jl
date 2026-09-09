#Basis states
const basis = NLevelBasis(5);
const ket_0 = nlevelstate(basis, 1);
const ket_1 = nlevelstate(basis, 2);
const ket_r = nlevelstate(basis, 3);
const ket_p = nlevelstate(basis, 4);
const ket_l = nlevelstate(basis, 5);

@doc "Computational ground state `|0⟩` in the package's five-level basis." ket_0
@doc "Computational ground state `|1⟩` in the package's five-level basis." ket_1
@doc "Target Rydberg state `|r⟩` in the package's five-level basis." ket_r
@doc "Intermediate excited state `|p⟩` used by the two-photon model." ket_p
@doc "Loss state `|l⟩` used to collect population outside the ideal manifold." ket_l

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

"""
    get_atom_trajectories(sample, ωr, ωz, err_optns)

    get_atom_trajectories(sample, center, ωr, ωz, err_optns)

Construct position and velocity trajectories for one sampled atom.

The returned closures are used internally when assembling time-dependent
Hamiltonians for single- and two-atom simulations.
"""
@inline function get_atom_trajectories(
    sample::AbstractVector{<:Real},
    ωr::Float64,
    ωz::Float64,
    err_optns::AbstractDict{String, <:Any})
    return get_atom_trajectories(sample, zeros(Float64, 3), ωr, ωz, err_optns)
end

@inline function get_atom_trajectories(
    sample::AbstractVector{<:Real},
    center::AbstractVector{<:Real},
    ωr::Float64,
    ωz::Float64,
    err_optns::AbstractDict{String, <:Any})
    
    free_motion = err_optns["free_motion"]
    atom_motion = err_optns["atom_motion"]
    cx, cy, cz = center

    x, y, z, vx, vy, vz = atom_motion ? sample : zeros(Float64, 6);

    if err_optns["Doppler"]   
        Vx = t -> V(t, x, vx, ωr; free=free_motion);
        Vy = t -> V(t, y, vy, ωr; free=free_motion);
        Vz = t -> V(t, z, vz, ωz; free=free_motion);
    else
        Vx = t -> 0.0 
        Vy = t -> 0.0 
        Vz = t -> 0.0
    end

    if err_optns["xy_motion"]
        X  = t -> cx + R(t, x, vx, ωr; free=free_motion);
        Y  = t -> cy + R(t, y, vy, ωr; free=free_motion);
    else
        X = t -> cx
        Y = t -> cy
    end

    if err_optns["z_motion"]
        Z  = t -> cz + R(t, z, vz, ωz; free=free_motion);
    else
        Z  = t -> cz
    end

    return X, Y, Z, Vx, Vy, Vz;
end

"""
    gauss_field(x, y, z, w0, z0; n0=1, θ0=0)

Return the normalized complex Gaussian field envelope used by the refactored
dict-based beam parameter API.

This is the Gaussian-beam specialization that `Ω` uses when
`laser_params["type"] == "gauss"`.
"""
@inline function gauss_field(x, y, z, w0, z0; n0=1, θ0=0)
    return A(x, y, z, w0, z0; n=n0, θ=θ0) .* A_phase(x, y, z, w0, z0; θ=θ0)
end;

"""
    Ω(x, y, z, laser_params)

Return the position-dependent complex Rabi coupling for one excitation beam.

Depending on `laser_params`, this helper evaluates either a Gaussian beam or one
of the simple flat-top beam parameterizations used by the package.
"""
@inline function Ω(x, y, z, laser_params) 
    if (laser_params["type"] == "gauss")
        w0 = laser_params["w0"]
        z0 = laser_params["z0"]
        θ = laser_params["θ"]
        n_sg = laser_params["n_sg"]
        return laser_params["Ω"] .* gauss_field(x, y, z, w0, z0; n0=n_sg, θ0=θ) #A(x, y, z, w0, z0; n=n_sg, θ=θ) .* A_phase(x, y, z, w0, z0; θ=θ);
     
    elseif (laser_params["type"] == "flattop_HG")
        θ = laser_params["θ"]
        xx, zz = x*cos(θ) - z*sin(θ), x*sin(θ) + z*cos(θ)
        Ω_w_z = [laser_params["Ω"], laser_params["w0"],laser_params["z0"]] 
        return reconstruct_HG_field_2d(xx, y, zz, Ω_w_z, laser_params["coeffs_xy"]) #, laser_params["coeffs_xy"])

    elseif (laser_params["type"] == "flattop_LG")
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
        
end;

"""
    Δ(vx, vz, laser_params)

Return the single-photon Doppler shift induced by atomic motion for the first
laser.
"""
@inline function Δ(vx, vz, laser_params)
    w0, z0, θ =  laser_params["w0"], laser_params["z0"], laser_params["θ"] #laser_params[2:4]
    k = 2.0 * z0/w0^2;
    #Δx = k * sin(θ) * vx #    Δz = k * cos(θ) * vz
    return k * (sin(θ) * vx + cos(θ) * vz) #Δx + Δz
end;

"""
    δ(vx, vz, red_laser_params, blue_laser_params)

Return the two-photon Doppler detuning induced by atomic motion for the first and
second excitation beams.
"""
@inline function δ(vx, vz, first_laser_params, second_laser_params)
    wr0, zr0, θr = first_laser_params["w0"], first_laser_params["z0"], first_laser_params["θ"]  #first_laser_params[2:4];
    wb0, zb0, θb = second_laser_params["w0"], second_laser_params["z0"], second_laser_params["θ"]  #second_laser_params[2:4];
    
    kr = 2.0 * zr0/wr0^2;
    kb = 2.0 * zb0/wb0^2;

    δx = (kr*sin(θr) + kb*sin(θb))*vx
    δz = (kr*cos(θr) + kb*cos(θb))*vz 

    return δx + δz
end;

"""
    JumpOperators(decay_params)

Construct the Lindblad jump operators for the single-atom model.

`decay_params` is ordered as `[Γ0, Γ1, Γl, Γr]`.
"""
@inline function JumpOperators(decay_params)
    Γ0, Γ1, Γl, Γr = decay_params;
    decay_operators = [sqrt(Γ0)*σ0p, sqrt(Γ1)*σ1p, sqrt(Γl)*σlp, sqrt(Γr)*σlr]
    return decay_operators
end;

"""
    GenerateHamiltonian(sample, center, ωr, ωz, error_options,
        tspan_noise, f, red_laser_phase_amplitudes, blue_laser_phase_amplitudes,
        nodes, red_laser_params, blue_laser_params, Δ0, δ0)

Assemble the full time-dependent single-atom Hamiltonian for one Monte Carlo
sample.

This combines thermal motion, Doppler shifts, beam inhomogeneity, and optional
laser phase noise into the effective two-photon model.
"""
@inline function GenerateHamiltonian(
    sample,
    center,
    ωr, ωz,
    error_options,
    tspan_noise,    f,
    first_laser_phase_amplitudes,
    second_laser_phase_amplitudes,
    nodes,
    first_laser_params,
    second_laser_params,
    Δ0,     δ0    )

    # Trajectories
    X, Y, Z, Vx, Vy, Vz = get_atom_trajectories(sample, center, ωr, ωz, error_options);

    sigma_coeffs = [
        t -> Δ(Vx(t), Vz(t), first_laser_params) - Δ0,
        t -> δ(Vx(t), Vz(t), first_laser_params, second_laser_params) - δ0,]

    # Interpolate phase noise traces to pass to hamiltonian
    if error_options["laser_noise"]
        # Generate phase noise traces for first and second lasers
        ϕ_first_res  = ϕ(tspan_noise, f, first_laser_phase_amplitudes);
        ϕ_second_res = ϕ(tspan_noise, f, second_laser_phase_amplitudes);

        ϕ_first  = interpolate(nodes, ϕ_first_res, Gridded(Linear()));
        ϕ_second = interpolate(nodes, ϕ_second_res, Gridded(Linear()));
    else
        ϕ_first  = t -> 0.0;
        ϕ_second = t -> 0.0;
    end

    # Hamiltonian params trajectories 
    Ωr = t -> exp(1.0im * ϕ_first(t)) * Ω(X(t), Y(t), Z(t), first_laser_params);
    Ωb = t -> exp(1.0im * ϕ_second(t)) * Ω(X(t), Y(t), Z(t), second_laser_params);

    H = TimeDependentSum(
        [sigma_coeffs;[
            t -> Ωr(t) / 2.0,
            t -> conj(Ωr(t)) / 2.0,
            t -> Ωb(t) / 2.0,
            t -> conj(Ωb(t)) / 2.0,
        ]],
        operators
        );

    return H
end;

"""
    simulation(cfg::RydbergConfig; temperature_calibrate=false, ode_kwargs...)

Simulate single-atom two-photon Rydberg excitation and average over Monte Carlo
realizations.

The model combines the mechanisms emphasized in
[arXiv:1802.10424](https://arxiv.org/abs/1802.10424): finite-temperature atomic
motion and Doppler dephasing, spontaneous emission from the intermediate and
Rydberg states, and stochastic laser phase noise.

# Arguments
- `cfg::RydbergConfig`: simulation configuration.

# Keywords
- `temperature_calibrate = false`: if `true`, first adjust the ideal laser
  parameters with `calibrate_two_photon`.
- `ode_kwargs...`: keyword arguments forwarded to
  `timeevolution.master_dynamic`.

# Returns
- `(ρ, ρ2)`, the first and second moments of the density-matrix trajectory.
"""
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

    Γ0, Γ1, Γl   = cfg_t.error_options["spontaneous_decay_intermediate"] ? cfg_t.decay_params[1:3] : zeros(3) 
    Γr           = cfg_t.error_options["spontaneous_decay_rydberg"]      ? cfg_t.decay_params[4]   :  0.0
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
            cfg.shift,
            ωr, ωz,
            cfg_t.error_options,
        
            tspan_noise,
            cfg_t.f,
            cfg_t.first_laser_phase_amplitudes,
            cfg_t.second_laser_phase_amplitudes,
            nodes,
        
            cfg_t.first_laser_params,
            cfg_t.second_laser_params,
        
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


"""
    Ω_twophoton(Ωr, Ωb, Δ)

Return the effective two-photon Rabi frequency

```math
\\Omega_{\\mathrm{2ph}} = \\left| \\frac{\\Omega_r \\Omega_b}{2 \\Delta} \\right|.
```
"""
function Ω_twophoton(Ω1, Ω2, Δ)
    return abs(Ω1 * Ω2 / (2.0 * Δ))
end;

"""
    T_twophoton(Ωr, Ωb, Δ)

Return the effective two-photon Rabi period `2π / Ω_twophoton(Ωr, Ωb, Δ)`.
"""
function T_twophoton(Ω1, Ω2, Δ)
    return 2.0*π / Ω_twophoton(Ω1, Ω2, Δ)
end;

"""
    δ_twophoton(Ωr, Ωb, Δ)

Return the differential AC Stark shift of the effective two-photon transition:

```math
\\delta_{\\mathrm{2ph}} = \\frac{|\\Omega_r|^2 - |\\Omega_b|^2}{4 \\Delta}.
```
"""
function δ_twophoton(Ω1, Ω2, Δ)
    return (abs(Ω1)^2 - abs(Ω2)^2)/(4.0 * Δ)
end;

"""
    Ωr_required(Ω, Ωb, Δ)

Return the red-laser single-photon Rabi frequency needed to realize an effective
two-photon coupling `Ω` with blue-laser coupling `Ωb` and detuning `Δ`.
"""
function Ωr_required(Ω, Ωb, Δ)
    return 2.0 * Δ * Ω / abs(Ωb)
end;

"""
    calibrate_two_photon(cfg::RydbergConfig, n_samples=1000)

Temperature-average the effective two-photon coupling and AC Stark shift, then
return a corrected copy of `cfg`.

This helper is useful when matching an ideal pulse design to a finite-temperature
ensemble before running `simulation`.
"""
function calibrate_two_photon(cfg::RydbergConfig, n_samples=1000)
    cfg_calibrated = deepcopy(cfg)
    Ωr = cfg.first_laser_params["Ω"]
    Ωb = cfg.second_laser_params["Ω"]

    samples = samples_generate(cfg.trap_params, cfg.atom_params, n_samples)[1]

    # Correct single-photon Rabi frequencies to match temperature averaged two-photon Rabi frequency
    samples_Ω2 = [Ω_twophoton(Ω(x, y, z, cfg_calibrated.first_laser_params), Ω(x, y, z, cfg_calibrated.second_laser_params), Δ(vx, vz, cfg_calibrated.first_laser_params)) for (x, y, z, vx, _, vz) in samples]
    factor_Ω2 = sqrt((sum(samples_Ω2) / length(samples)) / Ω_twophoton(Ωr, Ωb, cfg.detuning_params[1]))
    Ωr_cor, Ωb_cor = Ωr / factor_Ω2, Ωb / factor_Ω2
    cfg_calibrated.first_laser_params["Ω"]  = Ωr_cor
    cfg_calibrated.second_laser_params["Ω"] = Ωb_cor

    # Correct resonance detuning to match temperature averaged AC Stark shifts
    samples_δ = [δ_twophoton(Ω(x, y, z, cfg_calibrated.first_laser_params), Ω(x, y, z, cfg_calibrated.second_laser_params), Δ(vx, vz, cfg_calibrated.first_laser_params)) for (x, y, z, vx, _, vz) in samples]
    δ_cor = sum(samples_δ) / length(samples_δ)
    δ_ideal = δ_twophoton(Ωr, Ωb, cfg.detuning_params[1])
    cfg_calibrated.detuning_params[2] += δ_cor - δ_ideal

    return cfg_calibrated
end
