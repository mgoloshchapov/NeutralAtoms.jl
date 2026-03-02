#_genlaguerre(n, α, x) = binomial(n+α,n) * HypergeometricFunctions.M(-n, α+1, x)
#Ω0, w0, z0, beam_type, n = laser_params; []

function Gouy_phase(lmzzr_gouy)
    l, m, z, z_r = lmzzr_gouy
    return ((l+m+1) * atan(z/z_r))
end

_genlaguerre(n, α, x) = binomial(n+α,n) * HypergeometricFunctions.M(-n, α+1, x)
function LG(r2, w, k)
    #L_k = _genlaguerre(k,0, 2 .* r.^2 ./ w^2)
    #norm = 2^(n/2) * factorial(n)^0.5 * (pi/2)^(1/4)
    return exp.(-r2 ./ w^2) .* _genlaguerre(k,0, 2 .* r2 ./ w^2) #Hn.(2^0.5 .* x ./ w ) #./ norm
end

@inline function HG(x, w, n)
    her_arg = zeros(n+1)
    her_arg[n+1] = 1
    Hn = Hermite(her_arg)
    #norm = 2^(n/2) * factorial(n)^0.5 * (pi/2)^(1/4)
    #return exp(-x^2 / w^2) * Hn(2^0.5 * x / w ) #./ norm
    return exp.(-x.^2 ./ w.^2) .* Hn.(2^0.5 .* x ./ w ) #./ norm
end

function LG_coeff(k,n)
    res = 0.0
    for i in k:n
        res +=  factorial(i) / (2^(i) * factorial(i - k))  
    end;
    return (-1.0)^k * res / factorial(k)
end

function HG_coeff(k,n)
    res = 0.0
    for i in k:n
        res += factorial(2*i) / (factorial(i) * 2^(3*i) * factorial(i - k) * factorial(2*k))
    end;
    return res
end

function HG_coefficients(p,q)
    c = zeros(p,q)
    for i in 1:p
        for j in 1:q
            c[i,j] =  HG_coeff(i, p-1) * HG_coeff(j, q-1)
        end;
    end;
    return c
end

function simple_flattopHG_field(x,y,z,laser_params)
    Ω, w0, zr, θ, n, m, sqz = laser_params; 
    #x, z = x0*cos(θ) - z0*sin(θ), x0*sin(θ) + z0*cos(θ)
    w = w0 .* sqrt.(1.0 .+ (z ./zr) .^2);
    k = 2 * zr / w0^2
    phase0 = k .* z .+ k ./ 2 .* z .* (x.^2 .+ y.^2) ./ (z.^2 .+ zr.^2) 

    E_sum = 0.0 + 0.0im
    #c = HG_coefficients(trunc(Int, n/2)+1, trunc(Int, m/2)+1)
    for i in 0:2:n
        for j in 0:2:m
            gouy = (i+j+1) .* atan.(z ./ zr)
            cij = HG_coeff(trunc(Int,i/2),trunc(Int,n/2)) * HG_coeff(trunc(Int,j/2),trunc(Int,m/2))
            E_sum += cij * HG.(x, w, trunc(Int,i)) * HG.(y, w * sqz, trunc(Int,j)) * exp.(-1.0im * gouy) 
            #E_sum += c[trunc(Int, i/2) + 1, trunc(Int, j/2) + 1] * HG.(x, w, i) * HG.(y, w, j) * exp.(-1.0im * gouy) 
        end;
    end;
    return Ω .* w0 .* E_sum ./ w .* exp.(-1.0im * phase0) 
end;

function simple_flattopLG_field(x,y,z,laser_params)
    Ω, w0, zr, θ, l = laser_params; #Ω, w0, zr, θ, l, m = laser_params;
    w = w0 .* sqrt.(1.0 .+ (z ./zr) .^2);
    k = 2 * zr / w0^2
    phase0 = k .* z .+ k ./ 2 .* z .* (x.^2 .+ y.^2) ./ (z.^2 .+ zr.^2) 

    E_sum = 0.0 + 0.0im
    Nd2 = trunc(Int, l/2)
    for i in 0:Nd2
        E_sum += LG_coeff(trunc(Int,i), Nd2) * LG.(x.^2 .+ y.^2, w, trunc(Int,i)) * exp.(-1.0im * (i+1) .* atan.(z ./ zr)) 
    end;
    return Ω .* w0 ./ w .* exp.(-1.0im * phase0) * E_sum
end

function decomposition_2d(x::Vector{Float64}, y::Vector{Float64}, F::Vector{Float64},
                         w::Float64, dx::Float64, dy::Float64; n_max=20, m_max=20)
    c = zeros(n_max+1, m_max+1); #m_max
    for n in 0:n_max
        for m in 0:m_max
            #norm = 1/(2^n * sqrt(2*π)*factorial(n))
            c[n+1, m+1] = sum(F .* HG(x, w, n) .* HG(y, w, m)) * dx * dy  / (w^2 * 2^(n+m) * (π/2) * factorial(n)*factorial(m))
        end;
    end;
    return c
end

function reconstruct_HG_field_2d(x::Float64,y::Float64,z::Float64, Ω_w_z::Vector{Float64}, c_xy::Matrix{Float64}) 
    #laser_params::Dict{String, Any}) ## cy::Vector{Float64}
    #c_xy = laser_params["coeffs_xy"]
    size_c = size(c_xy)
    Ω, w0, zr = Ω_w_z #laser_params["Ω"], laser_params["w0"], laser_params["z0"]; # wx, wy = ws 
    z0 = z / zr
    w = w0 * sqrt(1.0 + z0*z0);
    k = 2 * zr / w0^2
    phase0 = k * z + (k / 2) * z * (x^2 + y^2) / (z^2 + zr^2) 
    
    eiz = exp(-1.0im * atan(z0))
    HG_x = [HG(x, w, i-1) for i in 1:(size_c[1])]
    HG_y = [HG(y, w, j-1) for j in 1:(size_c[2])]

    E_sum = (0.0 + 0.0im) #* x
    #zr_arr = (0.0 + 0.0im) #* x
    ph_i = 1.0 + 0.0im #.+ zr_arr
    for i in 1:(size_c[1])
        ph_j = 1.0 + 0.0im # .+ zr_arr

        for j in 1:(size_c[2])
            #E_sum .+= c_xy[i,j] * HG_x[i] * HG_y[j] * exp.(-1.0im * (i + j - 1) * atan.(z ./ zr))
            E_sum += c_xy[i,j] * HG_x[i] * HG_y[j] * ph_j * ph_i
            ph_j *= eiz
        end;
        ph_i *= eiz
    end;
    return Ω * w0 / w * E_sum * exp(-1.0im * phase0) * eiz
end;