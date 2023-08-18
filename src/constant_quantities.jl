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

using QuantumOptics
using EBJJ

struct ConstantQuantity
    Jz::Operator
    Jx::Operator
    ψ0::Ket
    ψf::Ket
end
function ConstantQuantity(N::Int64, J0::Float64, Jf::Float64, U::Float64)
    Jz = sigmaz(SpinBasis(N / 2)) / 2 |> dense       # Jz is a diagonal matrix 
    Jx = sigmax(SpinBasis(N / 2)) / 2 |> dense       # Jx operator
    ψ0 = eigenstates(-2.0 * J0 * Jx + U * Jz^2)[2][1] # Initial state
    ψf = eigenstates(-2.0 * Jf * Jx + U * Jz^2)[2][1] # Final state
    return ConstantQuantity(Jz, Jx, ψ0, ψf)
end
ConstantQuantity(c::Control) = ConstantQuantity(c.N, c.J0, c.Jf, c.U)
