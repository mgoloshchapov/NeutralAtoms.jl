"""
    is_zero(x)

Return `true` when `x == 0`.

This small helper is used when locating the first failed recapture event in
`release_evolve`.
"""
function is_zero(x)
    return x == 0
end;

"""
    release_evolve(tspan, cord, atom_params, trap_params; eps=1e-3)

Propagate one sampled atom during the release window and record recapture.

# Arguments
- `tspan`: release times in `μs`.
- `cord`: sampled initial condition `[x, y, z, vx, vy, vz]`.
- `atom_params`: `[m, T]`, with mass in atomic mass units and temperature in
  `μK`.
- `trap_params`: `[U0, w0, z0]`, with `U0` in `μK` and lengths in `μm`.

# Keywords
- `eps = 1e-3`: recapture cutoff. Atoms with total energy above `U0 * (1 - eps)`
  are treated as lost.

# Returns
- Boolean-like vector whose `i`th entry is `1` if the atom remains recaptured at
  `tspan[i]` and `0` otherwise.
"""
function release_evolve(tspan, cord, atom_params, trap_params; eps=1e-3)
    xi, yi, zi, vxi, vyi, vzi = cord;
    m, T = atom_params;
    U0, w0, z0 = trap_params;
    
    x = xi .+ vxi * tspan;
    # y = yi .+ vyi * tspan - g0 * tspan .^2;
    y = yi .+ vyi * tspan;
    z = zi .+ vzi * tspan;
    
    kinetic = K(cord, trap_params, m);
    potential = U0 .* (1.0 .- A.(x, y, z, w0, z0) .^2);
    recap = (kinetic .+ potential) .< U0 * (1.0-eps);
    
    idx = findfirst(is_zero, recap);
    
    #Changed != nothing to !isnothing
    if !isnothing(idx)
        recap[idx:end] .= 0;
    end;
    
    return recap
end; 


"""
    release_recapture(tspan, trap_params, atom_params, N; freq=10, skip=1000, eps=1e-3, harmonic=true)

Estimate a release-and-recapture curve by averaging over `N` sampled atoms.

This is the package's temperature-estimation entry point for optical-tweezer
experiments: sample thermal initial conditions, evolve each atom during the
release window, and average the recapture outcomes.

# Arguments
- `tspan`: release times in `μs`.
- `trap_params`: `[U0, w0, z0]`, with `U0` in `μK` and lengths in `μm`.
- `atom_params`: `[m, T]`, with mass in atomic mass units and temperature in
  `μK`.
- `N`: number of Monte Carlo samples.

# Keywords
- `freq = 10`: number of Metropolis updates skipped between retained samples.
- `skip = 1000`: burn-in length for the Metropolis sampler.
- `eps = 1e-3`: recapture cutoff passed to `release_evolve`.
- `harmonic = true`: when `true`, use the harmonic approximation for fast
  Gaussian sampling. When `false`, use the Metropolis sampler in the full trap
  potential.

# Returns
- `(recapture, acc_rate)`, where `recapture` is the recapture probability at each
  time in `tspan` and `acc_rate` is the Metropolis acceptance rate. For
  `harmonic = true`, the acceptance rate is `1`.
"""
function release_recapture(tspan, trap_params, atom_params, N; freq=10, skip=1000, eps=1e-3, harmonic=true)
    samples, acc_rate = samples_generate(trap_params, atom_params, N; freq=freq, skip=skip, harmonic=harmonic);
    recapture = zeros(length(tspan));
    
    for i ∈ 1:N
        recapture .+= release_evolve(tspan, samples[i], atom_params, trap_params; eps=eps);
    end;
    
    return recapture ./ N, acc_rate
end;
