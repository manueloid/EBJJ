#=
# Constant quantities for the simulation

## 1. Introduction

In this code I will implement the ConstantQuantity type, to store all those quantities that will be used in the simulation of the quantum time evolution.

In general, the plan is to simulate a system of $ N $ bosonic particles, where the system is driven by a control function $ J(t) $.
In particular the initial and final value of the control function $ J(t) $ are given and will be fixed for different final times $ t_f $.

This means the initial and final state of the system are fixed, the only thing that changes is the shape of the control function $ J(t) $, and the respective time evolution operator.
Another fixed quantities are the operator $ J_x $ and $ J_z $, which are square matrics of size $ N $.

In the implementation of the types, I will follow the same steps I used in the other eSTA calculations.

The only difference comes from the fact that I am going to use the set up a new package instead of using just a collection of functions.
This hopefully will help speed up the loading time of the package, and will make it easier to use and to test.
=#
#=
## 2. Implementation

Here I will define the new types and the methods to create them.
Ideally I would have the construction method based on a `Control` type, from the `EBJJ` package.

The Hamiltonian of the system is given by 
$$
    H_{BH} = -2J(t) \hat{J}_x + U \hat{J}_z^2 
\tag{1}
$$

I am going to define also some constructors for this type.
I will start from a generic one that takes the parameters of the system as a input and returns the corresponding `ConstantQuantity` type.
This first constructor will take only number types as input, I will specialize this in the remainder of the code to take also `Control` types as input.
=#

struct ConstantQuantity
    Jx::Operator
    Jy::Operator
    Jz::Operator
    ψ0::Ket
    ψf::Ket
    css::Ket
end
function ConstantQuantity(N::Int64, Ωf::Float64, U::Float64)
    Jx = sigmax(SpinBasis(N / 2)) / 2 |> dense       # Jx operator
    Jy = sigmay(SpinBasis(N / 2)) / 2 |> dense       # Jy operator
    Jz = sigmaz(SpinBasis(N / 2)) / 2 |> dense       # Jz is a diagonal matrix 
    ψ0 = eigenstates(- 2 * Jx + U * Jz^2)[2][1]      # ground state of the system at time 0 
    ψf = eigenstates(- 2 * Ωf * Jx + U * Jz^2)[2][1] # ground state of the system at time T
    css = eigenstates(- 2 * Jx)[2][1] # CSS initial state
    return ConstantQuantity(Jx, Jy, Jz, ψ0, ψf, css)
end
ConstantQuantity(c::Control) = ConstantQuantity(c.N, c.Ωf, c.U)
#=
I want to define a new struct that will be the combination of the `Control` type and the `ConstantQuantity` type.

By doing that I will be able to use less variables when calling functions. This will make the functions more readable.

I will define a new `ControlQuantity` type.
=#

"""
    ControlQuantity(c::Control, q::ConstantQuantity)
New type that contains both the control parameters of the sysytem as well as the constant quantities like the operators and the ground states.
This is going to be used to reduce the number of variables passed to the functions.
"""
struct ControlQuantity
    c::Control
    q::ConstantQuantity
end
ControlQuantity(c::Control) = ControlQuantity(c, ConstantQuantity(c))
"""
    Error(time_err::Float64, mod_err::Float64)
New type where I store the errors that I am going to use when calculating the robustness
Little and useless breakdown of the fields:
- `time_err`: error that will be used in the calculation of the time noise robustness
- `mod_err`:  error that will be used in the calculation of the modulation noise robustness
"""
struct Error
    time_err::Float64
    mod_err::Float64
end
abstract type AbstractError end
"""
    TimeError(err::Float64)
New type that contains only the error that will be used in the calculation of the time noise robustness
It has only one field `err`
"""
struct TimeError <: AbstractError
    err::Float64
end
"""
    ModError(err::Float64)
New type that contains only the error that will be used in the calculation of the modulation noise robustness
It has only one field `err`
"""
struct ModError <: AbstractError
    err::Float64
end

# Now I am going to overload the `-` function to be able to have the inverse of the `TimeError` and `ModError` types
import Base: -
-(e::TimeError) = TimeError(-e.err)
-(e::ModError) = ModError(-e.err)

Error(δ::TimeError, ε::ModError) = Error(δ.err, ε.err)
Error(δ::TimeError) = Error(δ.err, 0.0)
Error(ε::ModError) = Error(0.0, ε.err)
