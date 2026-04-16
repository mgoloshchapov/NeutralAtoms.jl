module NeutralAtoms
using Distributed

using Distributions, Random
using PhysicalConstants.CODATA2018: c_0, k_B, m_u
using Unitful
using LinearAlgebra
using QuantumOptics
using SplitApplyCombine
using Interpolations
using Polynomials
using SpecialPolynomials
using HypergeometricFunctions
using OrderedCollections
using Plots
using ProgressBars
using MPI
using OrdinaryDiffEq
using Colors
using Statistics

export 

    w0_to_z0, trap_frequencies, E, I,
    release_recapture,
    samples_generate, R, V, get_trap_params, H, samples_visualise,
    Sϕ, ϕ_amplitudes, ϕ,
    Ω_twophoton, T_twophoton, δ_twophoton, Ωr_required, 
    ket_0, ket_1, ket_r, ket_p, ket_l,
    
    simple_flattopHG_field, simple_flattopLG_field,
    HG_coeff, gauss_field, HG_coefficients, 
    decomposition_HG_2d, reconstruct_HG_field_2d,

    simulation, RydbergConfig, get_rydberg_probs, plot_rydberg_probs,
    simulation_czlp, CZLPConfig, get_two_qubit_probs, plot_two_qubit_probs,

    get_gate, project_on_qubit, get_parity_fidelity, 
    get_parity, get_p_prob_and_fid, get_fidelity_phi,
    get_fidelity_osc, get_parity_osc
        
include("utilities.jl")
include("basic_experiments.jl")
include("lasernoise_sampler.jl")
include("atom_sampler.jl")
include("rydberg_model.jl")
include("arbitrary_beams.jl")
include("cz_model.jl")
include("fidelity.jl")
include("gates.jl")

end