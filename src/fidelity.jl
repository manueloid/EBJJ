#=
# Fidelity calculation
=#
#=
## 1. Introduction

In this code I will implement all the functions to calculate the fidelity given the eSTA corrections.
Those functions will take as input the `ConstantQuantity` type amidst other, and will return the fidelity.

First I will start by defining a function that calculates the fidelity for a single final time and as a function of time.
This is not the final function but rather a first implementation of the code that will also help testing it out.

I don't know if it makes sense to define a function that takes the primitive types as arguments, I would rather define a function that simulate the time evolution of the system assuming all the other parameters have been set.

Following this train of thought, I do believe the function should take as input the `ConstantQuantity` type, and the `Control` type.

Then I should use some kind of filter:
1. If the time array that is passed to the function start with a 0.0, then I should calculate the time evolution of the fidelity. 
This should be true also if I pass an integer number and a final time, meaning I want to check the evolution of the fidelity as a time progresses.

2. If the time array does not start with a 0.0, then the function should calculate the fidelity at the final time, without caring about the evolution of the fidelity.
This should be true also when I pass an integer number and the two final times, meaning I want the fidelity at final time for each final time in the array `t0, ..., tf`.

=#
#=
## 2. Implementation

I think the best way to implement the whole thing is to first define a function that given all the operators and the parameters, returns the Hamiltonian of the system.
The computational cost for the definition of the control function $ J(t) $ is negligible, so I will define it inside the function.

Then I will define a function that takes the Hamiltonian and the initial state, and returns the time evolution operator.
=#
#=
### 2.1. Single final time

The system is mainly controlled by the control function $ J(t) $.
The other parameters that influences the system are the number of particles $ N $, the initial and final value of the control function $ J_0 $ and $ J_f $, and the interaction strength $ U $, as well as the final time $ t_f$.
_Coincidentally_ these are all the values in the `Control` type.

The main source of memory consumption is, I think, the definition of the operators and the eigenstates.
A lot of memory can be saved by defining the operators and the eigenstates only once, and then passing them to the function.

In theory though, this is not a big issue when I want to calculate the evolution of the fidelity as a function of time, but I will do that regardless as it gives me some flexibility.
=#
#=
# Fidelity calculation
=#
#=
## 1. Introduction

In this code I will implement all the functions to calculate the fidelity given the eSTA corrections.
Those functions will take as input the `ConstantQuantity` type amidst other, and will return the fidelity.

First I will start by defining a function that calculates the fidelity for a single final time and as a function of time.
This is not the final function but rather a first implementation of the code that will also help testing it out.

I don't know if it makes sense to define a function that takes the primitive types as arguments, I would rather define a function that simulate the time evolution of the system assuming all the other parameters have been set.

Following this train of thought, I do believe the function should take as input the `ConstantQuantity` type, and the `Control` type.

Then I should use some kind of filter:
1. If the time array that is passed to the function start with a 0.0, then I should calculate the time evolution of the fidelity. 
This should be true also if I pass an integer number and a final time, meaning I want to check the evolution of the fidelity as a time progresses.

2. If the time array does not start with a 0.0, then the function should calculate the fidelity at the final time, without caring about the evolution of the fidelity.
This should be true also when I pass an integer number and the two final times, meaning I want the fidelity at final time for each final time in the array `t0, ..., tf`.

=#
#=
## 2. Implementation

I think the best way to implement the whole thing is to first define a function that given all the operators and the parameters, returns the Hamiltonian of the system.
The computational cost for the definition of the control function $ J(t) $ is negligible, so I will define it inside the function.

Then I will define a function that takes the Hamiltonian and the initial state, and returns the time evolution operator.
=#
#=
### 2.1. Single final time

The system is mainly controlled by the control function $ J(t) $.
The other parameters that influences the system are the number of particles $ N $, the initial and final value of the control function $ J_0 $ and $ J_f $, and the interaction strength $ U $, as well as the final time $ t_f$.
_Coincidentally_ these are all the values in the `Control` type.

The main source of memory consumption is, I think, the definition of the operators and the eigenstates.
A lot of memory can be saved by defining the operators and the eigenstates only once, and then passing them to the function.

In theory though, this is not a big issue when I want to calculate the evolution of the fidelity as a function of time, but I will do that regardless as it gives me some flexibility.
=#

#=
## 3 Robustness

I defined a new type called `Error` so that I did not need to define any functions that calculate the robustness, I am going to pass the errors to the fidelity function an then it will calculate the fidelity given the corresponding error.

In this way I do not have to define multiple fidelities functions, one is enough
=#

"""
    fidelity(q::ConstantQuantity, c::Control, e::Error)
Returns the fidelity of the protocol given the parameters of the system in `ConstantQuantity` and `Control` types, and the error `e` in the control function.
"""
function fidelity(q::ConstantQuantity, c::ControlSTA, e::Error=Error(0.0, 0.0))
    ε, δ = e.mod_err, e.time_err
    Ω(t) = control_function(t + δ, c) * (1 + ε)
    H(t, psi) = -2 * Ω(t) * q.Jx + c.U * q.Jz^2
    fid(t, psi) = abs2.(dagger(q.ψf) * psi) # Function that calculates the fidelity of the system
    return timeevolution.schroedinger_dynamic([0.0, c.T], q.ψ0, H; fout=fid)[2][end]
end
"""
    fidelity(q::ConstantQuantity, c::Control, corrs = Vector{Float64}, e::Error = Error(0.0, 0.0))
Returns the fidelity of the eSTA approach when the corrections have already been calculated, given the parameters of the system and an error `e`
"""
function fidelity(q::ConstantQuantity, c::ControlFull, corrs::Vector{Float64}, e::Error=Error(0.0, 0.0))
    ε, δ = e.mod_err, e.time_err
    Ω(t) = control_functionX(t + δ, c, corrs) * (1 + ε)
    H(t, psi) = -2 * Ω(t) * q.Jx + c.U * q.Jz^2
    fid(t, psi) = abs2.(dagger(q.ψf) * psi) # Function that calculates the fidelity of the system
    return timeevolution.schroedinger_dynamic([0.0, c.T], q.ψ0, H; fout=fid)[2][end]
end
fidelity(q::ConstantQuantity, c::ControlFull, corrs::AbstractVector{Corrs}, e::Error=Error(0.0, 0.0)) = fidelity(q, c, corrections(corrs), e)
fidelity(q::ConstantQuantity, c::ControlFull, e::Error=Error(0.0, 0.0)) = fidelity(q, c, corrections(c), e)
robustness(q::ConstantQuantity, c::Control, δ::TimeError) = (fidelity(q, c, Error(0.0, δ.err)) - fidelity(q, c, Error(0.0, -δ.err))) / (2 * δ.err)
robustness(q::ConstantQuantity, c::Control, ε::ModError) = (fidelity(q, c, Error(ε.err, 0.0)) - fidelity(q, c, Error(-ε.err, 0.0))) / (2 * ε.err)
