#=
# New type

In this file I will define a new `struct` to store all the relevant quantities of the problem.

Similarly to what I have done with the internal Bosonic Josephson junction, I will define a parent type `Control` and then I will define a subtype for each of the different Hamiltonians I am going to encounter.

The child types will be:
- `ControlFull`: for the full Hamiltonian
- `ControlSTA`: for the intermediate Hamiltonian

Regardless of the Hamiltonian, the `Control` will contain the following fields:
- `N`: number of particles
- `T`: total time of the evolution
- `Ωf`: final value of the control function
- `U`: interaction strength

Having everything set up in  dimensionless fashion, there is no need for the initial value of the control function $ \Omega_0 $, because it is always equal to 1.
=#

abstract type Control end
"""
ControlFull(N::Int64, J0::Float64, Jf::Float64, U::Float64, T::Float64)
Control Parameter for the eSTA protocol, it also contains the parameters that characterize the protocol like the number of corrections we want as well as the number of maximum states we are interested in
"""
struct ControlFull <: Control
    N::Int64
    Ωf::Float64
    U::Float64
    T::Float64
    nλ::Int64
    states::AbstractVector{Int64}
end

"""
ControlSTA(N::Int64, T::Float64, J0::Float64, Jf::Float64, U::Float64)
Derived type of ControlParameter for the STA protocol
"""
struct ControlSTA <: Control
    N::Int64
    Ωf::Float64
    U::Float64
    T::Float64
end

#=
# Notation about Λ

We have that $ \Lambda = UN/2 $
=#

"""
    Corrs(c::Control,n::Int64, kn::AbstractVector{ComplexF64}, gn::ComplexF64)
This struct contains the value of the corrections kn and gn for a given control object `c`, calculated for a given state `n`.
"""
struct Corrs
    c::Control
    n::Int64
    kn::AbstractVector{ComplexF64}
    gn::ComplexF64
end
Hess(c::Corrs) = c.kn * c.kn'
v(c::Corrs) = conj(c.gn) * c.kn |> real
Hess(cs::AbstractArray{Corrs, 1}) = [Hess(c) for c in cs] 
v(cs::AbstractArray{Corrs, 1}) = [v(c) for c in cs]
