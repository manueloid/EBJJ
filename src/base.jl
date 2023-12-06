#=
# New type
In this file I will define a new `struct` to store all the relevant quantities of the problem.

Similarly to what I have done with the internal Bosonic Josephson junction, I will define a parent type `ControlParameter` and then I will define a subtype for each of the different Hamiltonians I am going to encounter.

The child types will be:
- `ControlFull`: for the full Hamiltonian
- `ControlSTA`: for the intermediate Hamiltonian

Regardless of the Hamiltonian, the `ControlParameter` will contain the following fields:
- `N`: number of particles
- `T`: total time of the evolution
- `J0`: initial value of the control function
- `Jf`: final value of the control function
- `U`: interaction strength


=#

abstract type Control end
"""
ControlFull(N::Int64, J0::Float64, Jf::Float64, U::Float64, T::Float64)
Control Parameter for the eSTA protocol, it also contains the parameters that characterize the protocol like the number of corrections we want as well as the number of maximum states we are interested in
"""
struct ControlFull <: Control
    N::Int64
    J0::Float64
    Jf::Float64
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
    J0::Float64
    Jf::Float64
    U::Float64
    T::Float64
end

#=
## Type constructor

Here I will implement the type constructor for the `ControlParameter` type.

I will start from a constructor that does not take any input and returns a Control type with some default values.
The default values will be:

- $ N = 30 $
- $ T = 0.05\pi $
- $ J_0 = 0.1 $
- $ J_f = 0.05 $
- $ U = 0.02 $

After that - since most of the time I am only interested in changing the total time and the number of particle - I will define a constructor that takes only those two parameters as input and then it will set the other parameters to default values.

From the paper, the constant values are:

I will use multiple dispatch in such a way that the order of the parameters does not matter.
=#

ControlFull() = ControlFull(30, 0.1, 0.05, 0.02, 0.05 * pi, 1, 2:2)
ControlFull(c::ControlFull, states::AbstractVector{Int64}) = ControlFull(c.N, c.J0, c.Jf, c.U, c.T, c.nλ, states)
ControlSTA() = ControlSTA(30, 0.1, 0.05, 0.02, 0.05 * pi)
ControlSTA(c::ControlFull) = ControlSTA(c.N, c.J0, c.Jf, c.U, c.T)
ControlFull(N::Int64, T::Float64) = ControlFull(N, 0.1, 0.05, 0.02, T, 1, 2:2)
ControlFull(T::Float64, N::Int64) = ControlFull(N, 0.1, 0.05, 0.02, T, 1, 2:2)
ControlSTA(N::Int64, T::Float64) = ControlSTA(N, 0.1, 0.05, 0.02, T)
ControlSTA(T::Float64, N::Int64) = ControlSTA(N, 0.1, 0.05, 0.02, T)
c_time(c::ControlFull, t::Float64) = ControlFull(c.N, c.J0, c.Jf, c.U, t, c.nλ, c.states)
c_time(c::ControlSTA, t::Float64) = ControlSTA(c.N, c.J0, c.Jf, c.U, t)

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
