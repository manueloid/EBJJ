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
Derived type of ControlParameter for the full Hamiltonian
"""
struct ControlFull <: Control
    N::Int64
    J0::Float64
    Jf::Float64
    U::Float64
    T::Float64
end

"""
ControlSTA(N::Int64, T::Float64, J0::Float64, Jf::Float64, U::Float64)
Derived type of ControlParameter for the intermediate Hamiltonian
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

ControlFull() = ControlFull(30, 0.1, 0.05, 0.02, 0.05 * pi)
ControlSTA() = ControlSTA(30, 0.1, 0.05, 0.02, 0.05 * pi)
ControlFull(N::Int64, T::Float64) = ControlFull(N, 0.1, 0.05, 0.02, T)
ControlFull(T::Float64, N::Int64) = ControlFull(N, 0.1, 0.05, 0.02, T)
ControlSTA(N::Int64, T::Float64) = ControlSTA(N, 0.1, 0.05, 0.02, T)
ControlSTA(T::Float64, N::Int64) = ControlSTA(N, 0.1, 0.05, 0.02, T)
