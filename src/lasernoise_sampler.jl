@doc doc"""
    Sϕ(f, laser_phase_params)

Return the laser phase-noise spectral density sampled on the frequency grid `f`.

`laser_phase_params` is expected to be `[h0, hg, σg, fg]`: a white-noise floor
`h0` plus optional Gaussian servo bumps with amplitudes `hg`, widths `σg`, and
center frequencies `fg`. The returned spectrum is divided by `f^2`, matching the
conversion from frequency-noise to phase-noise spectral density used in the
single-atom Rydberg model.
"""
function Sϕ(f, laser_phase_params)
    h0, hg, σg, fg = laser_phase_params;
    res = 2.0 * h0 * ones(length(f));
    
    if length(hg) > 0
        for i ∈ [1:length(hg);]
            res = res .+ 2*hg[i] .* exp.(-(f .- fg[i]).^2 ./ (2 * σg[i]^2));
        end;
    end;
        
    return res ./ (f .^ 2)
end;

@doc doc"""
    ϕ_amplitudes(f, laser_phase_params)

Convert a phase-noise spectral-density model into cosine-series amplitudes.

The amplitudes are suitable for `ϕ`, which samples random initial phases and
constructs a time-domain phase-noise realization for one laser trajectory.
"""
function ϕ_amplitudes(f, laser_phase_params)
    h0, hg, σg, fg = laser_phase_params;
    Δf = f[2]-f[1];
    res = 2.0 * h0 * ones(length(f));
    
    if length(hg) > 0
        for i ∈ [1:length(hg);]
            res = res .+ 2*hg[i] .* exp.(-(f .- fg[i]).^2 ./ (2 * σg[i]^2));
        end;
    end;
    
    return 2.0 .* sqrt.(Δf * res) ./ f;
end;

"""
    ϕ(tspan, f, amplitudes)

Sample one stochastic laser phase-noise trajectory on `tspan`.

# Arguments
- `tspan`: time grid in `μs`.
- `f`: frequency grid in `MHz`.
- `amplitudes`: Fourier amplitudes, typically produced by `ϕ_amplitudes`.

# Returns
- Vector of phase samples in radians.
"""
function ϕ(tspan, f, amplitudes)
    N = length(f);
    ϕf = rand(Uniform(0.0, 2.0*π), N); #generate random phases for components
    res = vec(sum(amplitudes .* cos.(2*π * f .* tspan' .+ ϕf), dims=1));

    return res
end;
