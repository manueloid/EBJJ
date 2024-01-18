
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
"""
    corrections(c::Control, tfs::AbstractVector{Float64})
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


using EBJJ, Plots
max_state = 10
nλ = 8
U = 0.5
N = 10
t0, tf = 0.01, 4.0
Ωf = 0.5
tfs = range(t0, tf, length=49) 
c = ControlFull(N, Ωf, U, tf, nλ, 2:2:max_state);
cs = ControlSTA(c);
fid_esta = fidelities(c, tfs)
fid_sta = fidelities(cs, tfs)
plot(tfs, fid_esta, label="eSTA")
plot!(tfs, fid_sta, label="STA")

scatter()
for  i in eachindex(tfs)
    hess, vec = EBJJ.Hess(corr[i]), EBJJ.v(corr[i]) # Value of the hessian and the vector at the first final time, for different states n
    push!(hess, sum(hess))
    hessian_plot = scatter!(
        ones(length(hess))*tfs[i],
        [h[1]|>real for h in hess],
        legend=false,
   )
end
plot!()
scatter()
final_vec = []
final_hess = []
for  i in eachindex(tfs)
    hess, vec = EBJJ.Hess(corr[i]), EBJJ.v(corr[i]) # Value of the hessian and the vector at the first final time, for different states n
    push!(final_vec, sum(vec))
    push!(final_hess, sum(hess))
    push!(vec, sum(vec))
    hessian_plot = scatter!(
        ones(length(vec))*tfs[i],
        [h[1]|>real for h in vec],
        legend=false,
   )
end
plot!()

using LinearAlgebra
for  i in eachindex(tfs)
    hess, vec = EBJJ.Hess(corr[i]), EBJJ.v(corr[i]) # Value of the hessian and the vector at the first final time, for different states n
    num = 
    hessian_plot = scatter!(
        ones(length(vec))*tfs[i],
        [h[1]|>real for h in vec],
        legend=false,
   )
end
plot!()

plot(tfs, reduce(vcat, final_vec), label="vec")
plot!(tfs, reduce(vcat, final_hess|>real), label="hess")

#=
# Fidelity for the eSTA approach

There is a problem with the fidelity in the sense that the STA seems to be working great even in a regime where it should not be working.
What I can do is simulate the system in a range where STA should not work that well and calculate the overlap beetween the numerical simulation and the analitical solution of the STA.

I will take the `position` function from another file and I will also define a new function to calculate the time evolution of the system given the parameters of the system and the control function.
The fidelity function will be taken from anothe file as well.
=#

using SpecialPolynomials, QuadGK, EBJJ, Test, HCubature, QuantumOptics
he(n::Int64, x::Float64) = SpecialPolynomials.basis(Hermite, n)(x)
"""
    f2(t::Float64, c::Control)
Return the value of the argument of the Gaussian term of the wave function for a certain time `t`, given the parameters of the system in the `Control` type.
"""
function f2(t::Float64, c::Control)
    h, Λ = 2.0 / c.N, EBJJ.Λ(c)
    # Definition of the auxiliary functions
    b(t) = auxiliary(t, c)
    db(t) = EBJJ.auxiliary_1d(t, c)
    return sqrt(1 / (2Λ)) * 1 / (h * b(t)^2) - im * db(t) / (2Λ * h * b(t))
end
"""
    φ(t::Float64, c::Control)
Function that returns the value of the phase factor for the given control parameters.
"""
function φ(t::Float64, c::Control)
    h, Λ = 2.0 / c.N, EBJJ.Λ(c)
    imag_phase_integrand(t::Float64) = 2Λ * h * real(f2(t, c)) # Integrand of the phase factor
    return quadgk(τ -> imag_phase_integrand(τ), 0.0, t, atol=1e-7)[1]
end
"""
    position(n::Int64, z::Float64, α::ComplexF64, h::Float64)
Fourier transform of the wave function in momentum representation, for a general complex number `α`.
In this case I needed to introduce the parameter `h` coming from the Fourier transform.
This function does not have the imaginaray phase term.
"""
function position(n::Int64, z::Float64, α::ComplexF64, h::Float64)
    r = sqrt(real(α))
    return (r^2 / pi)^(1 / 4) * # Normalisation term
           (2^n * factorial(n))^(-1 / 2) * # Hermite polynomial normalisation term
           im^n / sqrt(h) * 1 / sqrt(α) * (conj(α) / α)^(n / 2) * # Fourier transform term
           exp(-z^2 / (2 * α)) * # Gaussian term
           he(n, z * r / abs(α)) # Hermite polynomial
end
"""
    position(n::Int64, z::Float64, t::Float64, c::Control)
Function that returns the Fourier transform of the wave function in position representation for a certain time `t`, given the parameters of the system in the `Control` type.
"""
function position(n::Int64, z::Float64, t::Float64, c::Control)
    h = 2.0 / c.N
    α = f2(t, c)
    return position(n, z, α, h) * exp(-im * (n + 1 / 2) * φ(t, c)) # Phase factor
end
max_state = 2
nλ = 2
U = 0.1
N = 50
t0, tf = 0.01, 0.75
Ωf = 0.5
tfs = range(t0, tf, length=49) 
c = ControlFull(N, Ωf, U, tf, nλ, 2:2:max_state);
cs = ControlSTA(c);
q = ConstantQuantity(cs) 
"""
    fidelity(q::ConstantQuantity, c::ControlSTA)
Returns the time evolution 
"""
function fidelity(q::ConstantQuantity, c::ControlSTA; nsteps=100)
    time = range(0.0, c.T, length=nsteps)
    Ω(t) = control_function(t, c)
    H(t, psi) = - Ω(t) * q.Jx + c.U * q.Jz^2
    # fid(t, psi) = abs2.(dagger(q.ψf) * psi) # Function that calculates the fidelity of the system
    time, solutions =  timeevolution.schroedinger_dynamic(time, q.ψ0, H)
end
sol =fidelity(q,cs)
zrange = range(-cs.N/2, cs.N/2, length=cs.N + 1)
analytic(t) = [position(0, z, t, cs) for z in zrange] .* sqrt(2 / cs.N)
times = sol[1]
solutions = [ s.data for s in sol[2] ]
analytic_time = analytic.(times)
using Plots
plot(times, [s .* conj(a) |> sum |> abs for (s, a) in zip(solutions, analytic_time)])
EBJJ.fidelity(q, c)

#=
These are the most relevant values that I need to start the comparison:
1. The final time
2. the number of steps I want to calculate the fidelity for
3. 
=#
