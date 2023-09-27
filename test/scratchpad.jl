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



J(γ::Int64, N::Int64, U) = U * N / (2.0 * γ)
N = 10
J0 = J(10, N, 49.0);
Jf = J(40, N, 49.0);
tf = 0.1 # in seconds
U = 0.49;
fid = zeros(10)
fid_corr = zeros(10)
tfs = range(0.04, tf, length=10) |> collect
c = ControlFull(N, J0, Jf, U, tf);
q = ConstantQuantity(c);

for (index, tf) in enumerate(tfs)
    c = ControlFull(N, J0, Jf, U, tf)
    J(t::Float64) = control_function(t, c)
    corr = corrections([2, 4, 6], c)
    fid[index] = fidelity(q, c)
    fid_corr[index] = fidelity(q, c, +corr)
end

using Plots
plot(tfs, fid)
plot!(tfs, fid_corr)


fid = fidelity(q, c)
corr = corrections([2], c)
fidelity(q, c, corr)

b(t) = auxiliary(t, c)
db(t) = EBJJ.auxiliary_1d(t, c)
d2b(t) = EBJJ.auxiliary_2d(t, c)
J1(t) = control_function(t, c)
tarr = range(0.0, c.T, length=100) |> collect
plot(tarr, b.(tarr))
plot!(tarr, db.(tarr))
plot!(tarr, d2b.(tarr))
plot!(tarr, J1.(tarr))
