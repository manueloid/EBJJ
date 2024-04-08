
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

using Plots, EBJJ

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

using EBJJ, ProgressMeter
"""
    fidelities(c::Control, tfs::AbstractVector{Float64})
Calculate the fidelity of the control object `c` at the final times `tfs` for different values of the final time
"""
function fidelities(c::Control, tfs::AbstractVector{Float64})
    q = ConstantQuantity(c)
    fid = zeros(length(tfs))
    p = Progress(length(tfs), 1, "Calculating fidelities")
    Threads.@threads for index in eachindex(tfs)
        cl = EBJJ.c_time(c, tfs[index])
        fid[index] = EBJJ.fidelity(q, cl)
        next!(p)
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
    p = Progress(length(tfs), 1, "Calculating robustnesses")
    Threads.@threads for index in eachindex(tfs)
        cl = EBJJ.c_time(c, tfs[index])
        rob[index] = EBJJ.robustness(q, cl, ε) |> abs
        next!(p)
    end
    return rob
end
"""
    kngn(c::Control, tfs::AbstractVector{Float64})
Return the value of both the kn and the gn for a control object `c` at the final times `tfs` for different values of the final time in the array `tfs`
"""
function kngn(c::Control, tfs::AbstractVector{Float64})
    corrs = []
    Threads.@threads for index in eachindex(tfs)
        cl = EBJJ.c_time(c, tfs[index])
        push!(corrs, EBJJ.corrections(cl))
    end
    return corrs
end

using EBJJ, Plots, DelimitedFiles
max_state = 2
nλ = 1
U = .1
N = 200
t0, tf = 0.005,  1.
Ωf = 0.1
tfs = range(t0, tf, length=100) 
c = ControlFull(N, Ωf, U, tf, nλ, 2:2:max_state);
cs = ControlSTA(c);
q = ConstantQuantity(c)

fid_esta = fidelities(c, tfs)
fid_sta = fidelities(cs, tfs)

plot()
plot(tfs, fid_esta, label="eSTA", xlim = (0.02, 0.3) )
plot!(tfs, fid_sta, label="STA" )

some = readdlm("/home/manueloid/Repos/ExternalBJJ/data/fidelity/fid200esta01.dat")
some2 = readdlm("/home/manueloid/Repos/ExternalBJJ/data/fidelity/fid200sta01.dat")

# plot(tfs, fid_esta, label="eSTA")
plot!(some[:,1], some[:,2], label="eSTA data", style = :dash)
plot!(some2[:,1], some2[:,2], label="STA data", style = :dash)
 savefig("/home/manueloid/Desktop/compare200_01.pdf")

rob_esta = robustnesses(c, tfs, ModError(1.e-7))
rob_sta = robustnesses(cs, tfs, ModError(1.e-7))
plot()
plot(tfs, rob_esta, label="eSTA")
plot!(tfs, rob_sta, label="STA")
