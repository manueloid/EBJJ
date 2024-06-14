
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

using EBJJ, ProgressMeter, QuantumOptics
"""
    oat(q::ConstantQuantity, c::ControlFull, corrs::Vector{Float64}, e::Error = Error(0.0, 0.0))
Return the final state of the evolution starting from the initial CSS state and evolving it using the Hamiltonian of the system under study.
"""
function oat(q::ConstantQuantity, c::ControlFull, corrs::Vector{Float64}, e::Error=Error(0.0, 0.0))
    ε, δ = e.mod_err, e.time_err
    css_state = q.css
    Ω(t) = control_functionX(t + δ, c, corrs) * (1 + ε)
    H(t, psi) = -2 * Ω(t) * q.Jx + c.U * q.Jz^2
    return timeevolution.schroedinger_dynamic([0.0, c.T], css_state, H)[end][end] # Returns the final state 
end
function oat(q::ConstantQuantity, c::ControlSTA, e::Error=Error(0.0, 0.0))
    ε, δ = e.mod_err, e.time_err
    css_state = eigenstates(-2 * q.Jx)[end][1]
    Ω(t) = control_functionX(t + δ, c) * (1 + ε)
    H(t, psi) = -2 * Ω(t) * q.Jx + c.U * q.Jz^2
    return timeevolution.schroedinger_dynamic([0.0, c.T], css_state, H)[end][end] # Returns the final state 
end
oat(q::ConstantQuantity,c::ControlFull, e::Error=Error(0.0, 0.0)) = fidelity_css(q, c, corrections(corrections(c)), e)
"""
    our(q::ConstantQuantity, c::ControlFull, corrs::Vector{Float64}, e::Error = Error(0.0, 0.0))
Start the simulation from the ground state of _our_ Hamiltonian and evolve it using the Hamiltonian of the system under study.
It returns the state at the final time.
"""
function our(q::ConstantQuantity, c::ControlFull, corrs::Vector{Float64}, e::Error=Error(0.0, 0.0))
    ε, δ = e.mod_err, e.time_err
    in = q.ψ0
    Ω(t) = control_functionX(t + δ, c, corrs) * (1 + ε)
    H(t, psi) = -2 * Ω(t) * q.Jx + c.U * q.Jz^2
    return timeevolution.schroedinger_dynamic([0.0, c.T], in, H)[end][end]
end
function our(q::ConstantQuantity, c::ControlSTA, e::Error=Error(0.0, 0.0))
    ε, δ = e.mod_err, e.time_err
    in = q.ψ0
    Ω(t) = control_functionX(t + δ, c) * (1 + ε)
    H(t, psi) = -2 * Ω(t) * q.Jx + c.U * q.Jz^2
    return timeevolution.schroedinger_dynamic([0.0, c.T], in, H)[end][end]
end


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
function fidelities_css(c::Control, tfs::AbstractVector{Float64})
    q = ConstantQuantity(c)
    fid = zeros(length(tfs))
    p = Progress(length(tfs), 1, "Calculating fidelities")
    Threads.@threads for index in eachindex(tfs)
        cl = EBJJ.c_time(c, tfs[index])
        fid[index] = fidelity_css(q, cl)
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
nλ = 2
U = .05
N = 50
t0, tf = 0.002, 0.2
Ωf = 0.2
tfs = range(t0, tf, length=10) 
c = ControlFull(N, Ωf, U, t0, nλ, 2:2:max_state);
cs = ControlSTA(c);
q = ConstantQuantity(c)

squeezing(state::Ket, q::ConstantQuantity, α::Float64) = dagger(state) * ( cos(α) * q.Jz + sin(α) * q.Jy)^2 * state - (dagger(state) * ( cos(α) * q.Jz + sin(α) * q.Jy) * state)^2 |> real |> s -> s / (c.N / 4)
ξN(psi, q) = (dagger(psi) * (q.Jz)^2 * psi - (dagger(psi) * (q.Jz) * psi)^2) / (c.N / 4) |> real
sta_state = our(q,cs)
corrs = corrections(corrections(c))
esta_state = oat(q, c, corrs)


αrange = 0.0:0.01:π
sq = [squeezing(final, q, α) for α in αrange]
sq_sta, sq_esta = ξN.([sta_state, esta_state], Ref(q))
min = findmin(sq)
αrange[min[2]]
plot(αrange, sq, label = "Squeezing", xlabel = "α", ylabel = "Squeezing")
hline!([sq_sta], label = "sta")
hline!([sq_esta], label = "esta")




# fid_esta = fidelities(c, tfs)
fid_css_sta = fidelities_css(cs, tfs)
fid_css_esta = fidelities_css(c, tfs)

fid_sta = fidelities(cs, tfs)
fid_esta = fidelities(c, tfs)

plot()
plot(tfs , fid_css_esta, label="eSTA", title = "Fidelity starting from the CSS state")
plot!(tfs, fid_css_sta, label="STA")
# savefig("/home/manueloid/Desktop/fid_css.png")
# Save data to files 

plot()
plot(tfs, fid_esta, label="eSTA", title = "Fidelity starting from ψ0")
plot!(tfs, fid_sta, label="STA")
# savefig("/home/manueloid/Desktop/fid_standard.png")

# Control function
j_sta(t) = control_function(t, cs)
corrs = corrections(corrections(c))
j_esta(t) = control_function(t,c, corrs)
trange = range(0.0, c.T, length = 1000)

plot()
plot(trange, j_esta.(trange), label="eSTA", title = "Control function, tf = $(c.T)")
plot!(trange, j_sta.(trange), label="STA")

savefig("/home/manueloid/Desktop/control_function.png")

writedlm("/home/manueloid/Desktop/test_esta.dat", hcat(tfs, fid_esta))
writedlm("/home/manueloid/Desktop/test_sta.dat", hcat(tfs, fid_sta))

some = readdlm("/home/manueloid/Repos/ExternalBJJ/data/fidelity/fid50esta04.dat")
some2 = readdlm("/home/manueloid/Repos/ExternalBJJ/data/fidelity/fid50sta04.dat")
# plot(tfs, fid_esta, label="eSTA")
plot!(some[:,1], some[:,2], label="eSTA data", style = :dash)
plot!(some2[:,1], some2[:,2], label="STA data", style = :dash)

rob_esta = robustnesses(c, tfs, ModError(1.e-7))
rob_sta = robustnesses(cs, tfs, ModError(1.e-7))
plot()
plot(tfs, rob_esta, label="eSTA")
plot!(tfs, rob_sta, label="STA")

using ForwardDiff, EBJJ
ts = range(0.0, c.T, length = 1000)
f(t) = control_functionX(t,c, [-1.0, 1.0])
g(t) = control_functionX(t,c)
plot(ts, f.(ts))
plot!(ts, g.(ts))

