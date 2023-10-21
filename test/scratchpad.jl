#=
# Non numerical tests

In this file I will implement a series of tests that will basically check that the plots of the quantities I am interested in are correct.

This file is not a proper set of test, but rather a qualitative check of the results.

I will check the auxiliary functions, the control function, the Hamiltonian of the system and finally the fidelity.
=#
#=
# 1 - Auxiliary functions

Here I need to check that the auxiliary functions work as expected.
I will plot a series of auxiliary functions for different final times to see how they look.
=#

using Plots
tfs = range(0.01, 0.05, length=10) |> collect # Final times
time_evo(final_time::Float64) = range(0.0, final_time, length=1000) |> collect # return a range of points in the interval [0, tf]
plot()
for tf in tfs
    c = ControlFull(10, 50.0, 40.0, 0.49, tf)
    b(t) = auxiliary(t, c)
    plot!(time_evo(tfs[end]), b.(time_evo(tfs[end])))
end
plot!()

#=
# 2 - Control function

Now I am going to check how the control function behaves.

I will follow the same steps I used earlier.
=#

plot()
for tf in tfs
    c = ControlFull(10, 50.0, 40.0, 0.49, tf)
    J(t) = control_function(t, c)
    plot!(time_evo(tfs[end]), J.(time_evo(tfs[end])))
end
plot!()

#=
## 2.1 Control Function with corrections

I have also implmentend a function that return the value of the corrected control function given a list of corrections.
Here I am going to check if it works as expected.

The maximum values of the control function depend on the final time, so I will define the corrections as something of the form `1/tf * rand(5)`, hopefully the magnitude of the corrections will be similar to the one of the control function $J$.
=#

plot()
for tf in tfs
    c = ControlFull(10, 50.0, 40.0, 0.49, tf)
    J(t) = control_function(t, c)
    corr = 1 / tf * rand(5)
    J_corr(t) = control_function(t, c, corr)
    plot!(time_evo(tfs[end]), J_corr.(time_evo(tfs[end])))
    plot!(time_evo(tfs[end]), J.(time_evo(tfs[end])))
end
plot!()


#=
# 3 - Hamiltonian 

I have to figure out how to check if the Hamiltonian is correct.
I think I am on the right track because the fidelity seems to consistently go to 1 for longer final times and drop for shorter times.
=#

#=
# 4 - Fidelity

I will define a new function that evaluates both the corrections and the corresponding fidelity internally.
I will lay out everything inside the body of the function, hoping that I can get some nice result
=#


"""
    fidelity_tf(tfs::Vector{Float64}, c::Control)
Return the final fidelity of the system for each final time in the array `tfs`, given the system parameters in the `Control` type `c`.
It internally evaluates the corrections and the corresponding fidelity, as well as giving the fidelity of the STA protocol.
The Hamiltonian of the system is defined internally in the for loop.
"""
function fidelity(tfs::Vector{Float64}, c::Control)
    # Definition of the quantities that are not depending on the final time
    q = ConstantQuantity(c)         # Give the initial and final state as well as the operators
    ψ0, ψf = q.ψ0, q.ψf             # Initial and final states
    Jx, Jz = q.Jx, q.Jz             # Spin operators
    fid_sta = zeros(Float64, length(tfs)) # Array that will contain the fidelity of the STA protocol
    fid = zeros(Float64, length(tfs))     # Array that will contain the fidelity of the eSTA protocol
    Threads.@threads for index in eachindex(tfs)
        tf = tfs[index]
        cl = ControlFull(c.N, c.J0, c.Jf, c.U, tf) # local Control parameter
        corrs = corrections([2, 4, 6, 8], cl, 1)         # Corrections
        println("calculating fidelity for $tf")
        fid_sta[index] = EBJJ.fidelity(q, cl, tf)
        fid[index] = EBJJ.fidelity(q, cl, tf, corrs)
    end
    return fid_sta, fid
end
using QuantumOptics

J(γ::Int64, N::Int64, U) = U * N / (2.0 * γ)
N = 10
U = 0.40
J0 = J(10, N, U);
Jf = J(40, N, U);
tf = 0.5
tfs = range(0.06, tf, length=10) |> collect
np = 10:10:20 |> collect
c = ControlFull(N, J0, Jf, U, tf);
sta, esta = fidelity(tfs, c)
using Plots
plot(tfs, sta)
plot!(tfs, esta)

"""
    fidelity(q::ConstantQuantity, c::Control, e::Error)
Returns the fidelity of the protocol given the parameters of the system in `ConstantQuantity` and `Control` types, and the error `e` in the control function.
"""
function fidelity(q::ConstantQuantity, c::ControlSTA, e::Error=Error(0.0, 0.0))
    ε, δ = e.mod_err, e.time_err
    J(t) = control_function(t + δ, c) * (1 + ε)
    H(t, psi) = -2.0 * J(t) * q.Jx + c.U * q.Jz^2
    fid(t, psi) = abs2.(dagger(q.ψf) * psi) # Function that calculates the fidelity of the system
    return timeevolution.schroedinger_dynamic([0.0, c.T], q.ψ0, H; fout=fid)[2][end]
end
"""
    fidelity(q::ConstantQuantity, c::Control, corrs = Vector{Float64}, e::Error = Error(0.0, 0.0))
Returns the fidelity of the eSTA approach when the corrections have already been calculated, given the parameters of the system and an error `e`
"""
function fidelity2(q::ConstantQuantity, c::ControlFull, e::Error=Error(0.0, 0.0))
    ε, δ = e.mod_err, e.time_err
    corrs = corrections(2, c, 1)
    J(t) = control_function(t + δ, c, corrs) * (1 + ε)
    H(t, psi) = -2.0 * J(t) * q.Jx + c.U * q.Jz^2
    fid(t, psi) = abs2.(dagger(q.ψf) * psi) # Function that calculates the fidelity of the system
    return timeevolution.schroedinger_dynamic([0.0, c.T], q.ψ0, H; fout=fid)[2][end]
end

function test(q::ConstantQuantity, c::ControlFull)
    corrs = corrections(2, c, 1)
    return fidelity(q, c, corrs)
end
