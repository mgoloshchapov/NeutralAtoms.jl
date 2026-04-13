#_genlaguerre(n, α, x) = binomial(n+α,n) * HypergeometricFunctions.M(-n, α+1, x)
#Ω0, w0, z0, beam_type, n = laser_params; []

"""
    Gouy_phase((l, m, z, z_r))

Return the Gouy phase for a mode with transverse indices `l` and `m`.
"""
function Gouy_phase(lmzzr_gouy)
    l, m, z, z_r = lmzzr_gouy
    return ((l+m+1) * atan(z/z_r))
end

_genlaguerre(n, α, x) = binomial(n+α,n) * HypergeometricFunctions.M(-n, α+1, x)
"""
    LG(r2, w, k)

Return the radial Laguerre-Gauss factor used by the simple flat-top beam model.
"""
function LG(r2, w, k)
    #L_k = _genlaguerre(k,0, 2 .* r.^2 ./ w^2)
    #norm = 2^(n/2) * factorial(n)^0.5 * (pi/2)^(1/4)
    return exp.(-r2 ./ w^2) .* _genlaguerre(k,0, 2 .* r2 ./ w^2) #Hn.(2^0.5 .* x ./ w ) #./ norm
end

"""
    HG(x, w, n)

Return the Hermite-Gauss factor used by the simple flat-top beam model.
"""
function HG(x, w, n)
    her_arg = zeros(n+1)
    her_arg[n+1] = 1
    Hn = Hermite(her_arg)
    #norm = 2^(n/2) * factorial(n)^0.5 * (pi/2)^(1/4)
    return exp.(-x.^2 ./ w^2) * Hn.(2^0.5 .* x ./ w ) #./ norm
end

"""
    LG_coeff(k, n)

Return the coefficient multiplying the `k`th Laguerre-Gauss term in the simple
flat-top expansion of order `n`.
"""
function LG_coeff(k,n)
    res = 0.0
    for i in k:n
        res +=  factorial(i) / (2^(i) * factorial(i - k))  
    end;
    return (-1.0)^k * res / factorial(k)
end

"""
    HG_coeff(k, n)

Return the coefficient multiplying the `k`th Hermite-Gauss term in the simple
flat-top expansion of order `n`.
"""
function HG_coeff(k,n)
    res = 0.0
    for i in k:n
        res += factorial(2*i) / (factorial(i) * 2^(3*i) * factorial(i - k) * factorial(2*k))
    end;
    return res
end

"""
    HG_coefficients(p, q)

Return the coefficient matrix used by the separable Hermite-Gauss flat-top
construction.
"""
function HG_coefficients(p,q)
    c = zeros(p,q)
    for i in 1:p
        for j in 1:q
            c[i,j] =  HG_coeff(i, p-1) * HG_coeff(j, q-1)
        end;
    end;
    return c
end

"""
    simple_flattopHG_field(x, y, z, laser_params)

Return the complex field amplitude of the package's simple Hermite-Gauss
flat-top beam model.

This helper is used by `Ω` when the blue excitation beam is described by a
Hermite-Gauss superposition instead of a single Gaussian mode.
"""
function simple_flattopHG_field(x,y,z,laser_params)
    Ω, w0, zr, beam_type, n, m, sqz = laser_params;
    w = w0 .* sqrt.(1.0 .+ (z ./zr) .^2);
    k = 2 * zr / w0^2
    phase0 = k .* z .+ k ./ 2 .* z .* (x.^2 .+ y.^2) ./ (z.^2 .+ zr.^2) 

    E_sum = 0.0 + 0.0im
    #c = HG_coefficients(trunc(Int, n/2)+1, trunc(Int, m/2)+1)
    for i in 0:2:n
        for j in 0:2:m
            gouy = (i+j+1) .* atan.(z ./ zr)
            cij = HG_coeff(trunc(Int,i/2),trunc(Int,n/2)) * HG_coeff(trunc(Int,j/2),trunc(Int,m/2))
            E_sum += cij * HG.(x, w, i) * HG.(y, w * sqz, j) * exp.(-1.0im * gouy) 
            #E_sum += c[trunc(Int, i/2) + 1, trunc(Int, j/2) + 1] * HG.(x, w, i) * HG.(y, w, j) * exp.(-1.0im * gouy) 
        end;
    end;
    return Ω .* w0 .* E_sum ./ w .* exp.(-1.0im * phase0) 
end

"""
    simple_flattopLG_field(x, y, z, laser_params)

Return the complex field amplitude of the package's simple Laguerre-Gauss
flat-top beam model.
"""
function simple_flattopLG_field(x,y,z,laser_params)
    Ω, w0, zr, beam_type, l, m= laser_params;
    w = w0 .* sqrt.(1.0 .+ (z ./zr) .^2);
    k = 2 * zr / w0^2
    phase0 = k .* z .+ k ./ 2 .* z .* (x.^2 .+ y.^2) ./ (z.^2 .+ zr.^2) 

    E_sum = 0.0 + 0.0im
    Nd2 = trunc(Int, l/2)
    for i in 0:Nd2
        E_sum += LG_coeff(i, Nd2) * LG.(x.^2 .+ y.^2, w, i) * exp.(-1.0im * (i+1) .* atan.(z ./ zr)) 
    end;
    return Ω .* w0 ./ w .* exp.(-1.0im * phase0) * E_sum
end
