"""
    NeutralAtoms

Simulation tools for neutral-atom experiments in optical tweezers.

The package focuses on the workflow used in two-photon Rydberg experiments:
trap characterization, thermal sampling of atom motion, stochastic laser phase
noise, Lindblad dynamics for single-atom excitation, and blockade-mediated
two-qubit phase-gate simulations.
"""
module NeutralAtoms
using JLD2
using Distributed
using Markdown
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
using OrdinaryDiffEq
using Colors
using Statistics

export 

    w0_to_z0, trap_frequencies, E, I,
    release_recapture,
    samples_generate, R, V, get_trap_params, H,
    Sϕ, ϕ_amplitudes, ϕ,
    Ω_twophoton, T_twophoton, δ_twophoton, Ωr_required, 
    ket_0, ket_1, ket_r, ket_p, ket_l,
    
    simple_flattopHG_field, simple_flattopLG_field,
    HG_coeff, gauss_field, HG_coefficients, 
    decomposition_HG_2d, reconstruct_HG_field_2d,
    
    pure_simulation_czlp,
    simulation, RydbergConfig, get_rydberg_probs, plot_rydberg_probs,
    simulation_czlp, CZLPConfig, get_two_qubit_probs, plot_two_qubit_probs,

    get_gate, project_on_qubit, get_parity_osc,
    get_fidelity_with_rz_phi, CZ_caliration, PhiPlus_fidelity_osc,
    get_rydberg_fidelity_configs, get_rydberg_error_budget, get_cz_error_budget,

    save_with_JLD2, load_with_JLD2
        
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
