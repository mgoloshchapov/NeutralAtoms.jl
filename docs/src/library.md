# API Reference

```@meta
CurrentModule = NeutralAtoms
```

This reference is organized by concepts rather than source-file order. Start
with the tutorials in [Examples](examples.md) if you want the recommended
workflow first.

## Package

```@docs
NeutralAtoms
```

## Beam And Trap Optics

```@docs
w
w0_to_z0
A
I
A_phase
E
trap_frequencies
get_trap_params
```

## Thermal Sampling And Trajectories

```@docs
Π
Π_Harmonic
K
H
prob_boltzmann
samples_generate
R
V
is_zero
release_evolve
release_recapture
```

## Laser Phase Noise

```@docs
Sϕ
ϕ_amplitudes
ϕ
```

## Single-Atom Rydberg Simulation

```@docs
RydbergConfig
ket_0
ket_1
ket_r
ket_p
ket_l
get_atom_trajectories
Ω_old
Ω
Δ
δ
Hamiltonian
JumpOperators
GenerateHamiltonian
simulation
Ω_twophoton
T_twophoton
δ_twophoton
Ωr_required
calibrate_two_photon
simulation_mpi
get_rydberg_probs
plot_rydberg_probs
```

## Arbitrary Beam Shaping

```@docs
Gouy_phase
LG
HG
LG_coeff
HG_coeff
HG_coefficients
simple_flattopHG_field
simple_flattopLG_field
Ω_red
Ω_blue
simulation_blue_intens
```

## CZ And Two-Qubit Simulation

```@docs
CZLPConfig
JumpOperatorsTwo
get_V
GenerateHamiltonianTwo
get_blockade_stark_shift_factor
simulation_czlp
get_two_qubit_probs
plot_two_qubit_probs
```

## Gate And Fidelity Utilities

```@docs
get_gate
project_on_qubit
get_rydberg_fidelity_configs
get_rydberg_infidelity
plot_rydberg_infidelity
get_parity_fidelity
get_parity
get_parity_fidelity_temp
get_cz_infidelity
plot_cz_infidelity
```
