#=
### 1 Look at this fidelity

I will first try to see where the problem is from.
I will implement some check to see if increasing the number of particles produces an increase fidelity given a fixed final time.

On the other hand I will try to see if fixing the number of particle and increasing the final time produces an increase in the fidelity.
=#

function fidelities(N::Vector{Int64}, tf::Float64)
    global fid = Vector{Float64}(undef, length(N))
    for (index, np) in enumerate(N)
        c = ControlFull(np, tf)
        q = ConstantQuantity(c)
        # corrs = corrections([2], c) # For the moment I do not care about the corrections
        fid[index] = fidelity(q, c)
    end
    return fid
end
function fidelities(N::Int64, tf::Vector{Float64})
    global fid = Vector{Float64}(undef, length(tf))
    c = ControlFull(N, tf[1])
    q = ConstantQuantity(c)
    for (index, t) in enumerate(tf)
        c = ControlFull(N, t)
        # corrs = corrections([2], c) # For the moment I do not care about the corrections
        fid[index] = fidelity(q, c)
    end
    return fid
end
n = 10:2:50 |> collect
tfs = range(0.01, 0.5, length=100) |> collect
fidelities(N::Vector{Int64}, tfs::Vector{Float64}) = [fidelities(n, tfs) for n in N]

fid = fidelities(n, tfs)

fid = hvcat(size(fid, 1), fid...)
using Plots
heatmap(n, tfs, fid)


heatmap(tfs, fid, label="Fidelity", xlabel="Final time", ylabel="Fidelity", title="Fidelity as a function of the final time")

using EBJJ
c = ControlFull(40, 0.01, 0.1, 0.01, 1.0);
q = ConstantQuantity(c);
corrs = corrections([2], c);
fidelity(q, c, corrs)
fidelity(q, c)


#=
#### 1.1. Conclusion of the analysis

It looks like the fidelity does not scale with the number of particles and the final time.
I will try to see if the problem is in the control function or in the Hamiltonian.
=#

#=
### 2 Testing the control function and its corrections

Here I am going to see if the control function $ J(t) $ is working as expected and I will also check how the corrections are changing the shape of it.

I will also check the auxiliary function, hopefully everything works out as expected.

I would like to get the same result they have in the paper and in order to do that I am going to use their same values, even though it is not as straightforward as they define everything in a weird way.

Nevertheless, the first quantity I need to define is the Rabi time $ t_{rabi} $, which is defined as $ t_{rabi} = \frac{\pi}{J_0} $, where $ J_0 $ is the value of the control function for $ t = 0 $.
Then they go on to define this value $ \gamma_0 = \frac{NU}{2J_0} $, which is the value of the control function at $ t = 0 $, and similarly for $ \gamma_f$.

Later on, they use some experimental values, so I am going to use them.
In particular, we have 
- $ N = 100 $ 
- $ U = 0.49 Hz $
- $ t_f $ = 20 ms $
- $ \gamma_0 = 10 $
- $ \gamma_f = 100 $
=#

J(γ::Int64, N::Int64, U) = U * N / (2.0 * γ)
N = 10
J0 = J(10, N, 0.49);
Jf = J(100, N, 0.49);
tf = 0.05 # in seconds
U = 0.49;
c = ControlFull(N, J0, Jf, U, tf);
J(t::Float64) = control_function(c)(t);
tarr(tf::Float64) = range(0.0, tf, length=1000);
plot(tarr(tf) ./ tf, J.(tarr(tf)))
q = ConstantQuantity(c);
corr = corrections([2], c)

fid = fidelity(q, c)
fid_corr = fidelity(q, c, N^2 * corr)

εs = range(-500.0, 500.0, length=100) |> collect
fidε = [fidelity(q, c, ε * corr) for ε in εs]
plot(εs, fidε)


fid = fidelity(q, c)
fidelity(q, c, corr)
