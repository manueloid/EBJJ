
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

function ΔJ_oat(ψ::Ket, q::ConstantQuantity)
    Jz, Jy = q.Jz, q.Jy
    num = ψ' * (Jz * Jy + Jy * Jz) * ψ
    den = ψ' * Jz^2 * ψ - ψ' * Jy^2 * ψ
    α = 1 / 2 * atan(num / den) |> real
    return cos(α)^2 * ψ' * Jz^2 * ψ + sin(α)^2 * ψ' * Jy^2 * ψ + cos(α) * sin(α) * ψ' * (Jz * Jy + Jy * Jz) * ψ
end
function css_oat(q::ConstantQuantity, c::Control)
    in = q.css # Initial state
    H(t, psi) = c.U * q.Jz^2
    ΔJ(t, psi) = ΔJ_oat(psi, q) / (c.N / 4) |> real
    # return timeevolution.schroedinger_dynamic([0.0, c.T], in, H)[end][end] # Returns the final state 
    return timeevolution.schroedinger_dynamic([0.0, c.T], in, H; fout=ΔJ)[2][end]
end
function gs_oat(q::ConstantQuantity, c::Control)
    in = q.ψ0
    H(t, psi) = c.U * q.Jz^2
    ΔJ(t, psi) = ΔJ_oat(psi, q) / (c.N / 4) |> real
    # return timeevolution.schroedinger_dynamic([0.0, c.T], in, H)[end][end] # Returns the final state 
    return timeevolution.schroedinger_dynamic([0.0, c.T], in, H; fout=ΔJ)[2][end]
end
function gs_ibjj(q::ConstantQuantity, c::ControlFull, corrs::Vector{Float64})
    in = q.ψ0
    Ω(t) = control_functionX(t, c, corrs)
    H(t, psi) = -2 * Ω(t) * q.Jx + c.U * q.Jz^2
    ΔJ(t, psi) = (dagger(psi) * (q.Jz)^2 * psi - (dagger(psi) * (q.Jz) * psi)^2) / (c.N / 4) |> real
    return timeevolution.schroedinger_dynamic([0.0, c.T], in, H; fout=ΔJ)[2][end]
end
function gs_ibjj(q::ConstantQuantity, c::ControlSTA)
    in = q.ψ0
    Ω(t) = control_functionX(t, c)
    H(t, psi) = -2 * Ω(t) * q.Jx + c.U * q.Jz^2
    ΔJ(t, psi) = (dagger(psi) * (q.Jz)^2 * psi - (dagger(psi) * (q.Jz) * psi)^2) / (c.N / 4) |> real
    return timeevolution.schroedinger_dynamic([0.0, c.T], in, H; fout=ΔJ)[2][end]
end

using EBJJ, Plots, DelimitedFiles
max_state = 2
nλ = 2
U = 0.05
N = 50
t0, tf = 0.002, 0.6
Ωf = 0.1
tfs = range(t0, tf, length=10)
c = ControlFull(N, Ωf, U, tf, nλ, 2:2:max_state);
cs = ControlSTA(c);
q = ConstantQuantity(c)

function squeezing_css(c::Control, tfs=AbstractVector{Float64})
    q = ConstantQuantity(c)
    ξs = zeros(4, length(tfs))
    p = Progress(length(tfs), 1, "Calculating ξN")
    Threads.@threads for index in eachindex(tfs)
        cl = EBJJ.c_time(c, tfs[index])
        corrs = corrections(corrections(cl))
        cs = ControlSTA(cl)
        ξs[1, index] = css_oat(q, cl)
        ξs[2,index] = gs_oat(q, cl)
        ξs[3, index] = gs_ibjj(q, cs)
        ξs[4, index] = gs_ibjj(q, cl, corrs)
        next!(p)
    end
    return ξs
end
all = squeezing_css(c, tfs)

plot(tfs, all)

function squeezing_ξ(c::Control, tfs::AbstractVector{Float64})
    q = ConstantQuantity(c)
    ξN = zeros(length(tfs))
    p = Progress(length(tfs), 1, "Calculating ξN")
    Threads.@threads for index in eachindex(tfs)
        cl = EBJJ.c_time(c, tfs[index])
        ξN[index] = squeezing_ξ(q, cl)
        next!(p)
    end
    return ξN
end

ξN_css = squeezing_css(cs, tfs)
ξN_esta = squeezing_ξ(c, tfs)
ξN_sta = squeezing_ξ(cs, tfs)

plot()
plot(tfs, ξN_css, label="eSTA", title="Squeezing parameter")
plot!(tfs, ξN_esta, label="eSTA")
plot!(tfs, ξN_sta, label="STA")
# savefig("/home/manueloid/Desktop/ξN.png")


# fid_esta = fidelities(c, tfs)
fid_css_sta = fidelities_css(cs, tfs)
fid_css_esta = fidelities_css(c, tfs)

fid_sta = fidelities(cs, tfs)
fid_esta = fidelities(c, tfs)

plot()
plot(tfs, fid_css_esta, label="eSTA", title="Fidelity starting from the CSS state")
plot!(tfs, fid_css_sta, label="STA")
# savefig("/home/manueloid/Desktop/fid_css.png")
# Save data to files 

plot()
plot(tfs, fid_esta, label="eSTA", title="Fidelity starting from ψ0")
plot!(tfs, fid_sta, label="STA")
# savefig("/home/manueloid/Desktop/fid_standard.png")

# Control function
j_sta(t) = control_function(t, cs)
corrs = corrections(corrections(c))
j_esta(t) = control_function(t, c, corrs)
trange = range(0.0, c.T, length=1000)

plot()
plot(trange, j_esta.(trange), label="eSTA", title="Control function, tf = $(c.T)")
plot!(trange, j_sta.(trange), label="STA")

savefig("/home/manueloid/Desktop/control_function.png")

writedlm("/home/manueloid/Desktop/test_esta.dat", hcat(tfs, fid_esta))
writedlm("/home/manueloid/Desktop/test_sta.dat", hcat(tfs, fid_sta))

some = readdlm("/home/manueloid/Repos/ExternalBJJ/data/fidelity/fid50esta04.dat")
some2 = readdlm("/home/manueloid/Repos/ExternalBJJ/data/fidelity/fid50sta04.dat")
# plot(tfs, fid_esta, label="eSTA")
plot!(some[:, 1], some[:, 2], label="eSTA data", style=:dash)
plot!(some2[:, 1], some2[:, 2], label="STA data", style=:dash)

rob_esta = robustnesses(c, tfs, ModError(1.e-7))
rob_sta = robustnesses(cs, tfs, ModError(1.e-7))
plot()
plot(tfs, rob_esta, label="eSTA")
plot!(tfs, rob_sta, label="STA")

using ForwardDiff, EBJJ
ts = range(0.0, c.T, length=1000)
f(t) = control_functionX(t, c, [-1.0, 1.0])
g(t) = control_functionX(t, c)
plot(ts, f.(ts))
plot!(ts, g.(ts))

