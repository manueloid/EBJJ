
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


using EBJJ 
"""
    fidelities(c::Control, tfs::AbstractVector{Float64})
Calculate the fidelity of the control object `c` at the final times `tfs` for different values of the final time
"""
function fidelities(c::Control, tfs::AbstractVector{Float64})
    q = ConstantQuantity(c)
    fid = zeros(length(tfs))
    Threads.@threads for index in eachindex(tfs)
        cl = EBJJ.c_time(c, tfs[index])
        fid[index] = EBJJ.fidelity(q, cl)
    end
    return fid
end
"""
    robustnesses(c::Control, tfs::AbstractVector{Float64}, ε::AbstractError)
    Calculate the robustness of the control object `c` at the final times `tfs` for different values of the final time in the array `tfs`
"""
function robustnesses(c::Control, tfs::AbstractVector{Float64}, ε::AbstractError)
    q = ConstantQuantity(c)
    rob = zeros(length(tfs))
    Threads.@threads for index in eachindex(tfs)
        cl = EBJJ.c_time(c, tfs[index])
        rob[index] = EBJJ.robustness(q, cl, ε)
    end
    return rob
end

max_state = 12
nλ = 5
U = 0.06
N = 10
J0 = 0.245
Jf = 0.1225
t0, tf = 0.1, 0.9
tfs = range(t0, tf, length=19) 
c = ControlFull(N, J0, Jf, U, tf, nλ, 2:2:max_state);
cs = ControlSTA(c);
δ, ε = TimeError(1e-7), ModError(1e-7)
esta = fidelities(c, tfs)
sta = fidelities(cs, tfs)
using Plots
plot(tfs, esta, label="esta")
plot!(tfs, sta, label="sta")
