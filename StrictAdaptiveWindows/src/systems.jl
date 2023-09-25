#=

Call me crazy but I don't think we should try to standardize the process.
For Ayush's matrices, there is no way I'm getting N and S matrices,
    and anyway I don't think that's a good thing to standardize. Matrices are bad.

So let us let this file simply implement the key parts of the system:
    H       (needed to generate observable)
    ψ_REF   (needed to generate reference state)
    ψ_FCI   (needed for accuracy comparisons)

    REF
    FCI
    FES

It may be worthwhile folding this struct into the main module but
    the process for building each one probably isn't, so let's leave it here for now.

NOTE: farsight...for systems with a degenerate ground state,
    we *pry* want FES to give the first *non-degenerate* energy.
But that might not always be well-defined.

=#

import CtrlVQE

using Serialization: serialize, deserialize
using LinearAlgebra: eigen, Hermitian, diag
using Parameters: @with_kw
using NPZ: npzread

import ..MATRIXPATH, ..SYSTEMPATH

import ..as_dict

function load_system(
    code::String;
    matrixpath=MATRIXPATH, systempath=SYSTEMPATH,
)
    systemfile = "$systempath/$code"
    # LOAD SYSTEM IF POSSIBLE
    if isfile(systemfile)
        system = MolecularSystem(; deserialize(systemfile)...)
    # OTHERWISE, CREATE IT AND THEN SAVE IT
    else
        system = make_system(code; matrixpath=matrixpath)
        serialize(systemfile, as_dict(system))
    end
    return system
end


const Hchain_pattern = r"^(c?H)(\d+)(.)(P|J|B)(m|n)?$"
function make_system(code::String; matrixpath=matrixpath=MATRIXPATH)
    # TEST H chain PATTERN
    matched = match(Hchain_pattern, code)
    if !isnothing(matched)
        cfg, num, geo, map, tap = matched

        isnothing(tap) && (tap = "")
        N = parse(Int, num)

        return MolecularSystem(
            "$matrixpath/$code.npy",
            "$matrixpath/$(cfg)$(num)_$(map)$(tap)_N.npy",
            "$matrixpath/$(cfg)$(num)_$(map)$(tap)_S.npy",
            N, 0,
        )
    end

    return MolecularSystem("$matrixpath/$code.npy")
end







@with_kw struct MolecularSystem{F}
    H::Matrix{Complex{F}}
    i_REF::Int                      # INDEX OF REFERENCE BASIS STATE
    n::Int                          # NUMBER OF QUBITS
    # FREQUENTLY USED STATEVECTORS
    ψ_REF::Vector{Complex{F}}
    ψ_FCI::Vector{Complex{F}}
    # ENERGIES
    REF::F
    FCI::F
    FES::F
end

function MolecularSystem(H::AbstractMatrix, i_REF::Int)
    F = real(eltype(H))

    N = size(H,1)
    n = round(Int, log2(size(H,1)))

    ψ_REF = CtrlVQE.LinearAlgebraTools.basisvector(N, i_REF)
    Λ, U = eigen(Hermitian(H))
    return MolecularSystem(
        convert(Matrix{Complex{F}}, H),
        i_REF, n,
        convert(Vector{Complex{F}}, ψ_REF),                 # ψ_REF
        convert(Vector{Complex{F}}, U[:,1]),                # ψ_FCI
        real(H[i_REF, i_REF]),                              # REF
        Λ[1],                                               # FCI
        Λ[2],                                               # FES
    )
end

function MolecularSystem(H_file::String)
    H = npzread(H_file)
    i_REF = argmin(real.(diag(H)))
    return MolecularSystem(H, i_REF)
end

function MolecularSystem(H_file::String, N_file::String, S_file::String, N::Int, S::Int)
    H = npzread(H_file)
    Nd = real.(diag(npzread(N_file)))
    Sd = real.(diag(npzread(S_file)))

    σ = sortperm(real.(diag(H)))
    sector = (N .≈ Nd) .&& (S .≈ Sd)
    i_REF = σ[findfirst(i -> sector[i], σ)]

    isnothing(i_REF) && error("No basis vector matches given fock sector.")
    return MolecularSystem(H, i_REF)
end