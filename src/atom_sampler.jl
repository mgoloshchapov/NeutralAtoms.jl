#Constants and scales
#------------------------------------
const c = ustrip(u"m/s", c_0);  #Speed of light
const kB = ustrip(u"J/K", k_B)  #Boltzmann constant
const mu = ustrip(u"kg", m_u);  #Unit of atomic mass

const E0 = kB * 1e-6;       #Characteristic energy in μK
const g0 = 9.81 * 1e-6;     #Gravity free fall acceleration
const vconst = sqrt(E0/mu); #Useful constant for kinetic energy
const r0 = 1e-6;            #Characteristic distance in m
#------------------------------------

"""
    Π(cord, trap_params)

Return the Gaussian-tweezer potential energy at phase-space point `cord`.
"""
function Π(cord, trap_params)
    U0, w0, z0 = trap_params;
    x, y, z, vx, vy, vz = cord;
    return U0 .* (1.0 .- A(x, y, z, w0, z0) .^2);
end;


"""
    Π_Harmonic(cord, trap_params)

Return the harmonic approximation to the tweezer potential used for fast thermal
sampling near the trap center.
"""
function Π_Harmonic(cord, trap_params)
    U0, w0, z0 = trap_params;
    x, y, z, vx, vy, vz = cord;
    
    r2 = x .^2 .+ y .^2;
    return U0 .* (2.0*r2 ./w0^2 + (z ./z0).^2);
end;



"""
    K(cord, trap_params, m)

Return the kinetic energy associated with `cord` for an atom of mass `m`.
"""
function K(cord, trap_params, m)
    U0, w0, z0 = trap_params;
    x, y, z, vx, vy, vz = cord;
    return m/vconst^2 *(vx .^2 + vy .^2 + vz .^2) / 2.0
end;



"""
    H(cord, trap_params, m; harmonic=false)

Return the total energy of one sampled atom in the tweezer model.

Set `harmonic = true` to use `Π_Harmonic`; otherwise the full Gaussian
potential `Π` is used.
"""
function H(cord, trap_params, m; harmonic=false)
    if harmonic
        return Π_Harmonic(cord, trap_params) .+ K(cord, trap_params, m)
    else
        return Π(cord, trap_params) .+ K(cord, trap_params, m)
    end;
end;


"""
    prob_boltzmann(cord, trap_params, atom_params; harmonic=false)

Return the unnormalized Boltzmann weight used by the Metropolis sampler in
`samples_generate`.
"""
function prob_boltzmann(cord, trap_params, atom_params; harmonic=false)    
    m, T = atom_params;
    return exp.(- H(cord, trap_params, m; harmonic) ./ T)
end;


"""
    samples_generate(trap_params, atom_params, N; freq=10, skip=1000, harmonic=true)

Generate Monte Carlo samples of initial atom coordinates and velocities.

# Arguments
- `trap_params`: `[U0, w0, z0]`, with `U0` in `μK` and lengths in `μm`.
- `atom_params`: `[m, T]`, with mass in atomic mass units and temperature in
  `μK`.
- `N`: number of Monte Carlo samples.

# Keywords
- `freq = 10`: number of Metropolis updates skipped between retained samples.
- `skip = 1000`: burn-in length for the Metropolis sampler.
- `harmonic = true`: use the harmonic approximation instead of the full
  Metropolis sampler.
- `eps = 1e-2`: reject samples with total energy above `U0 * (1 - eps)`.

# Returns
- `(samples, acc_rate)`, where each sample is `[x, y, z, vx, vy, vz]`. For
  `harmonic = true`, the acceptance rate is reported as `1`.
"""
function samples_generate(trap_params, atom_params, N; freq=10, skip=1000, harmonic=true, eps=1e-2)
    U0, w0, z0 = trap_params;
    m, T = atom_params;

    if harmonic
        mean = zeros(6);
        cov = Diagonal(([T*w0^2/(4.0*U0),T*w0^2/(4.0*U0),T*z0^2/(2.0*U0),vconst^2*T/m,vconst^2*T/m,vconst^2*T/m]));
        d = MvNormal(mean, cov);
        
        samples = Vector{Vector{Float64}}();
        while length(samples) < N
            cord = rand(d)  
            if H(cord, trap_params, m) < U0 * (1-eps)
                push!(samples, cord)
            end
        end
        return samples, 1
    else
        mean = zeros(6);
        vstep = vconst*sqrt(T/m);
        rstep = sqrt(T/U0)/2;
        cov = Diagonal(([w0*rstep, w0*rstep, z0*sqrt(2)*rstep, vstep, vstep, vstep]) .^ 2);
        d = MvNormal(mean, cov);
        samples = [[0.0, 0.0, 0.0, vstep/sqrt(3), vstep/sqrt(3), vstep/sqrt(3)]];
        u_acc = rand(Uniform(0.0, 1.0), N*freq + skip);
        acc_rate = 0;
        for i ∈ 1:N*freq + skip - 1
            cord_last = samples[end];
            cord_new = cord_last + rand(d);
            p_acc = prob_boltzmann(cord_new, trap_params, atom_params; harmonic)/prob_boltzmann(cord_last, trap_params, atom_params; harmonic);
            
            if p_acc > u_acc[i] && H(cord_new, trap_params, m; harmonic) < U0 * (1-eps)
                push!(samples, cord_new);
                acc_rate += 1; 
            else
                push!(samples, cord_last);
            end;
        end;

        return samples[1+skip:freq:end], acc_rate/(N*freq + skip)
    end;
end;


"""
    R(t, ri, vi, ω; free=false)

Return the position of one Cartesian coordinate at time `t`.

For `free = false`, the coordinate follows harmonic motion with frequency `ω`;
for `free = true`, it follows ballistic motion.
"""
function R(t, ri, vi, ω; free=false)
    # if free
    #     return ri + vi * t;
    # else
    #     return ri * cos(ω * t) + vi/ω * sin(ω * t);
    # end;
    return free ? ri + vi * t : ri * cos(ω * t) + vi/ω * sin(ω * t);
end;    


"""
    V(t, ri, vi, ω; free=false)

Return the velocity of one Cartesian coordinate at time `t`.

This is the velocity companion to `R`.
"""
function V(t, ri, vi, ω; free=false)
    # if free
    #     return vi;
    # else
    #     return vi * cos(ω * t) - ri * ω * sin(ω * t);
    # end;
    return free ? vi : vi * cos(ω * t) - ri * ω * sin(ω * t);
end;       


@doc doc"""
    get_trap_params(ωr, ωz, U0, λ; m = 87.0, dif=true)

Infer Gaussian-tweezer parameters from target trap frequencies.

The trap is assumed to be formed by a single Gaussian beam with equal radial
frequencies in `x` and `y`.

# Arguments
- `ωr`: radial oscillation frequency in `MHz`.
- `ωz`: longitudinal oscillation frequency in `MHz`.
- `U0`: trap depth in `μK`.
- `λ`: trap wavelength in `μm`.

# Keywords
- `m = 87.0`: atom mass in atomic mass units.
- `dif = true`: when `true`, infer the trap geometry from the frequency
  difference and return the implied depth.

# Returns
- `(w0, z0, U)`: waist, Rayleigh length, and effective depth in package units.
"""
function get_trap_params(ωr, ωz, U0, λ; m = 87.0, dif=true)
    if dif
        Δω = ωr - ωz;
        w0 = λ/(1.0 - Δω/ωr) /(sqrt(2.0) * π);
        z0 = π*w0^2 / λ;
        Utemp = m * (ωr / (2.0 * vconst))^2;
    else
        w0 = vconst * sqrt(4.0 * U0 / m) / ωr ;
        z0 = vconst * sqrt(2.0 * U0 / m) / ωz ;
        Utemp = U0;
    end;
    return w0, z0, Utemp
end;
