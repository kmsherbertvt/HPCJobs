# Legend

### Hydrogen Chains and Cycles
Matrices with filename "{c?}H{\d+}{E|S}{J|B|P}{m?}.npy" correspond to hydrogen toy models.

The label "H{N}" indicates a line of N equally-spaced hydrogens,
    while  "cH{N}" indicates a ring.

The next letter indicates whether the geometry is near equilibrium (E) or stretched (S).
For stretched geometries, the distance between nearest neighbors is 3.0Å.
For equilibrium geometries, the geomtries are NOT optimized, so they are only _close_ to equilibrium.
The exact seaparations are given in the section below.

The next letter indicates the qubit mapping used:
    Jordan-Wigner (J), Bravyi-Kitaev (B) or parity (P).
The presence of an m after the mapping indicates whether a two-qubit tapering was applied
    to account for conservation of particles in each spin sector.

#### Equilibrium Geometries
- H2: 0.75Å
- H4: 0.90Å
- H6: 0.95Å
- H8: 0.95Å
- cH4: 1.25Å
- cH6: 1.00Å
- cH8: 1.05Å

## Hubbard Chains and Cycles
Matrices with filename "{c?}L{\d+}{t|U}{A|C|M}{J|B|P}{m?}.npy"
    correspond to 1d Hubbard models with particle-hole symmetry.

The letter "L" stands for "Lattice", since "H" for "Hubbard" was already taken for hydrogen.
The prefix "c" stands for "cyclic", ie. periodic boundary conditions.

The next letter is either a little "t" or a big "U",
    indicating which factor dominates in the model.
The actual value of t is always 1.0; the variable is the dimensionless parameter `u = U/4t`.
In the so-called "t" model, u=0.1, and in the so-called "U" model, u=1.0.

The next letter, "A", "C", or "M", represents the orbital basis.
- "A" is for "atomic". This is the one for which two-body terms are local, ie. the usual for Hubbard.
- "C" is for "core", ie. diagonalizing the hopping terms. Thus, the reference is the ground-state for u=0.
- "M" is for "molecular", ie. true Hartree-Fock.

The next letter indicates the qubit mapping used:
    Jordan-Wigner (J), Bravyi-Kitaev (B) or parity (P).
The presence of an m after the mapping indicates whether a two-qubit tapering was applied
    to account for conservation of particles in each spin sector.


## Special Matrices
- "lih30": LiH with a bond separation of 3.0Å.
    Frozen core and Qiskit parity mapping with 2-qubit tapering. Total of 4 qubits.
- "H215": H2 with a bond separation of 1.5Å.
    OpenFermion parity mapping with 2-qubit tapering. Total of 4 qubits.