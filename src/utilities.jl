#Beam waist radius
function w(z, w0, z0)
    return w0 .* sqrt.(1.0 .+ (z ./z0) .^2);
end;

"""
    w0_to_z0(w0, λ, M2=1.0)

Return Rayleigh length given beam waist radius `w0`, wavelength `λ` and `M2` parameter.

### Input

- `w0` -- beam waist radius in ``μm``
- `λ` -- beam wavelength in ``μm``
- `M2` -- (optional, default: `1.0`) M2 parameter of beam. For ideal Gaussian beam M2=1.0

### Output

Rayleigh length in ``μm``
"""
function w0_to_z0(w0, λ, M2=1.0)
    return π*w0^2/λ / M2;
end;


#Amplitude of gaussian beam with |E0|=1
function A(x, y, z, w0, z0; n=1, θ=0.0)
    xt, yt, zt = sqrt(x^2 + y^2), 0.0, z 
    xt, zt = xt*cos(θ)-zt*sin(θ), xt*sin(θ) + zt*cos(θ)
    return (w0 ./ w(zt, w0, z0)) .* exp.(- ((xt .^2 .+ yt .^2) ./ (w(zt, w0, z0) .^2)) .^ n)
end;


#Intensity of gaussian beam with |E0|=1
function I(x, y, z, w0, z0; n=1, θ=0.0)
    xt, yt, zt = sqrt(x^2 + y^2), 0.0, z 
    xt, zt = xt*cos(θ) - zt*sin(θ), xt*sin(θ) + zt*cos(θ)
    return ((w0 ./ w(zt, w0, z0)) .* exp.(-((xt .^2 .+ yt .^2) ./ (w(zt, w0, z0) .^2)).^n)) .^2
end;


#Phase of gaussian beam
function A_phase(x, y, z, w0, z0; θ=0.0)
    xt, yt, zt = sqrt(x^2 + y^2), 0.0, z 
    xt, zt = xt*cos(θ)-zt*sin(θ), xt*sin(θ) + zt*cos(θ)
    k = 2.0 * z0 / w0^2;
    return exp.(-1.0im * k * zt .* (0.5*(xt .^2 + yt .^ 2) ./ (zt .^2  + z0 .^2)) + 1.0im * atan.(zt ./ z0));
end;



#Complex amplitude of gaussian beam with |E0|=1
function E(x, y, z, w0, z0;n=1, θ=0.0)
    return A(x,y,z,w0,z0;n=n, θ=θ) .* A_phase(x,y,z,w0,z0; θ=θ)
end;


"""
    trap_frequencies(atom_params, trap_params)

Calculates trap frequencies from atom and trap parameters 

### Input

- `atom_params` -- vector [atom mass in a.u., atom temperature in ``\\mu K``]
- `trap_params` -- vector [trap depth ``U_{0}`` in ``\\mu K``, beam waist radius in ``\\mu m``, beam Rayleigh length in ``\\mu m``]

### Output

Vector of radial and longitudinal trap frequencies [`ωr`, `ωz`]

"""
function trap_frequencies(atom_params, trap_params)
    m, T = atom_params;
    U0, w0, z0 = trap_params;
    ωr = vconst/w0 * sqrt(4 * U0/m);
    ωz = vconst/z0 * sqrt(2 * U0/m);
    
    return ωr, ωz;
end;


function get_rydberg_probs(ρ, ρ2, eps=1e-12)
    probs_dict = OrderedCollections.OrderedDict{String, Vector{Float64}}();

    names = ["0", "1", "r", "p", "l"];
    states = [ket_0, ket_1, ket_r, ket_p, ket_l];
    for i in 1:5
        P = real(expect(states[i] ⊗ dagger(states[i]), ρ))
        P2 = real(expect(states[i] ⊗ dagger(states[i]), ρ2))
        S = @. sqrt(P2 - P^2 .+ eps) / length(ρ)
        probs_dict["P"*names[i]] = P
        probs_dict["S"*names[i]] = S 
    end 

    return probs_dict
end

function get_two_qubit_probs(ρ, ρ2, eps=1e-12)
    probs_dict = OrderedCollections.OrderedDict{String, Vector{Float64}}();
    names = ["00", "01", "10", "11"];

    states = [
        ket_0 ⊗ ket_0, 
        ket_0 ⊗ ket_1, 
        ket_1 ⊗ ket_0, 
        ket_1 ⊗ ket_1
        ];
    for i in 1:4
        P = real(expect(states[i] ⊗ dagger(states[i]), ρ))
        P2 = real(expect(states[i] ⊗ dagger(states[i]), ρ2))
        S = @. sqrt(P2 - P^2 .+ eps) / length(ρ)
        probs_dict["P"*names[i]] = P
        probs_dict["S"*names[i]] = S 
    end 

    return probs_dict
end

function plot_rydberg_probs(tspan, probs_dict)
    names = ["0", "1", "r", "p", "l"];
    colors = ["lightblue", "blue", "red", "orange", "green"];

    plt = Plots.plot()
    for i in 1:5
        P = probs_dict["P"*names[i]]
        S = probs_dict["S"*names[i]]
        plot!(
            tspan, [P P], fillrange=[P+S P-S], 
            ylim=(0.0, 1.0), xlim=(minimum(tspan), maximum(tspan)), 
            fillalpha=0.25, c=colors[i], 
            label=[nothing "P" * names[i]], linewidth=3
            )
    end
    xlabel!("Time, μs")
    ylabel!("Probability")
    title!("Rydberg Rabi oscillations")

    display(plt)
end

function plot_two_qubit_probs(tspan, probs_dict)
    names = ["00", "01", "10", "11"];
    colors = Plots.cgrad(:bam, 4, categorical = true)

    plt = Plots.plot()
    for i in 1:4
        P = probs_dict["P"*names[i]]
        S = probs_dict["S"*names[i]]
        plot!(
            tspan, [P P], fillrange=[P+S P-S], 
            ylim=(0.0, 1.0), xlim=(minimum(tspan), maximum(tspan)), 
            fillalpha=0.25, c=colors[i], 
            label=[nothing "P" * names[i]], linewidth=3
            )
    end
    xlabel!("Time, μs")
    ylabel!("Probability")
    title!("Two-qubit probabilities")

    display(plt)
end



# function calibrate_rabi(trap_params, atom_params, laser_params; n_samples=1000)
#     Ω, w, z, θ, n = laser_params
#     samples = samples_generate(trap_params, atom_params, n_samples)[1]
#     Ω_samples = [A(sample[1], sample[2], sample[3], w, z; n=n, θ=θ) for sample in samples];
#     Ω2_samples = [A(sample[1], sample[2], sample[3], w, z; n=n, θ=θ)^2 for sample in samples];
#     factor = sum(Ω_samples) / length(Ω_samples)
#     factor2 = sqrt(sum(Ω_samples^2) / length(Ω_samples))
#     return factor, factor2
# end



mutable struct RydbergConfig
    tspan::Vector{Float64}
    ψ0::Ket{NLevelBasis{Int64}, Vector{ComplexF64}}

    atom_params::Vector{Float64}
    trap_params::Vector{Float64}
    n_samples::Int64

    f::Vector{Float64}
    red_laser_phase_amplitudes::Vector{Float64}
    blue_laser_phase_amplitudes::Vector{Float64}
    
    red_laser_params::Vector{Float64}
    blue_laser_params::Vector{Float64}
    
    detuning_params::Vector{Float64}
    decay_params::Vector{Float64}

    atom_motion::Bool
    free_motion::Bool
    laser_noise::Bool
    spontaneous_decay_intermediate::Bool
    spontaneous_decay_rydberg::Bool
end


mutable struct CZLPConfig
    tspan::Vector{Float64}
    ψ0

    atom_params::Vector{Float64}
    trap_params::Vector{Float64}
    n_samples::Int64

    f::Vector{Float64}
    red_laser_phase_amplitudes::Vector{Float64}
    blue_laser_phase_amplitudes::Vector{Float64}
    
    red_laser_params::Vector{Float64}
    blue_laser_params::Vector{Float64}
    
    detuning_params::Vector{Float64}
    decay_params::Vector{Float64}

    atom_motion::Bool
    free_motion::Bool
    laser_noise::Bool
    spontaneous_decay_intermediate::Bool
    spontaneous_decay_rydberg::Bool

    atom_centers::Vector{Vector{Float64}}
    c6::Float64
    ΔtoΩ::Float64
    Ωτ::Float64
    ξ::Float64
    ϕ_RZ::Float64
end