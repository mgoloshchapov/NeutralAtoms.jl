"""
    get_gate(U)

Embed a `2 × 2` qubit unitary `U` into the five-level single-atom basis.

Only the computational subspace spanned by `ket_0` and `ket_1` is modified; the
other levels are left unchanged.
"""
function get_gate(U)
    op = dense(identityoperator(basis))
    op.data[1:2, 1:2] = U 
    return op
end


"""
    project_on_qubit(ρ, n=1)

Project a one- or two-atom density operator onto the computational qubit
subspace.

# Arguments
- `ρ`: density operator in the five-level basis.
- `n = 1`: number of atoms to project.

# Returns
- Density operator restricted to the `|0⟩, |1⟩` subspace for each atom.
"""
function project_on_qubit(ρ, n=1)
    basis_sub = SubspaceBasis(basis, [ket_0, ket_1])
    P = reduce(⊗, [projector(basis_sub, basis) for _ in 1:n])
    return P * ρ * dagger(P)
end;


Hadamard = get_gate(Matrix{ComplexF64}([1 1; 1 -1]/sqrt(2)))
X = get_gate(Matrix{ComplexF64}([
        0 1; 
        1 0
        ]))
Y = get_gate(Matrix{ComplexF64}([
    0     -1.0im; 
    1.0im      0
    ]));
Z = get_gate(Matrix{ComplexF64}([
    1  0; 
    0 -1
    ]));
RX = θ -> get_gate(Matrix{ComplexF64}([cos(θ/2) -1.0im*sin(θ/2); -1.0im*sin(θ/2) cos(θ/2)]));
RY = θ -> get_gate(Matrix{ComplexF64}([cos(θ/2) -sin(θ/2); sin(θ/2) cos(θ/2)]));
RZ = θ -> get_gate(Matrix{ComplexF64}([
    exp(-1.0im*θ/2)             0; 
    0               exp(1.0im*θ/2)
    ]));

# CNOT = Matrix{ComplexF64}([
#     1 0 0 0; 
#     0 1 0 0;
#     0 0 0 1;
#     0 0 1 0]);
# CZ = Matrix{ComplexF64}([
#     1 0 0 0;
#     0 1 0 0;
#     0 0 1 0;
#     0 0 0 -1;
# ])
