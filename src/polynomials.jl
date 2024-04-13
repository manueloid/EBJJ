#=
# Polynomials used in the calculations

In this file I am going to define all the polynomials that go into the definition of the control function $ J(t) $.
In particular I am going to define the auxiliary function $ b(t) $ and its derivatives, and the control function $ J(t) $.

The fac that I am going to define the derivatives explicitely, allows me to not use the package `ForwardDiff` in the main file.
=#
#=
### 1. Auxiliary function

The auxiliary function is the one that will go in the Ermakov equation, and it is defined as
$$
    b(t) = 1 + 
       6 (\gamma - 1) \left( \frac{t}{T} \right)^3 
    - 15 (\gamma - 1) \left( \frac{t}{T} \right)^4 
    + 10 (\gamma - 1) \left( \frac{t}{T} \right)^5
    \tag{1}
$$
where $ T $ is the total time of the process and $ \gamma $ in this case is given by $ \gamma = \left( \frac{J_0}{J_f} \right)^{1/4} $.

#### 1.1 Auxiliary function derivatives

I am going to explicitely define the derivatives of the auxiliary function $ b(t) $.
For the first derivative, the result is
$$
    \dot{b}(t) = \frac{1}{T} \left( 
        30 (\gamma - 1) \left( \frac{t}{T} \right)^2 
      - 60 (\gamma - 1) \left( \frac{t}{T} \right)^3
      + 30 (\gamma - 1) \left( \frac{t}{T} \right)^4 
      \right)
    \tag{2}
$$
and for the second derivative, we have
$$
    \ddot{b}(t) = \frac{1}{T^2} \left( 
        60 (\gamma - 1) \left( \frac{t}{T} \right)
      - 180 (\gamma - 1) \left( \frac{t}{T} \right)^2
      + 120 (\gamma - 1) \left( \frac{t}{T} \right)^3 
      \right).
    \tag{3}
$$
Moreover, all of these functions are piecewise defined, so that they are zero outside the interval $ [0, T] $.
=#
#=
### 2 Control function

The control function is defined as
$$
    J(t) = \frac{J_0}{b(t)^4} - \frac{\ddot{b}(t)}{2 b(t) U N}
    \tag{4}
$$
where $ U $ is the energy of the system and $ N $ is the number of particles.
=#

"""
    piecewise(f::Function, lim::AbstractVector{Number})
    Take a function and returns its piecewise version given the interval limits `t₀` and `tₑ`
"""
function piecewise(f::Function,t, t₀::Number, tₑ::Number)
    if t < t₀
        return f(t₀)
    elseif t > tₑ
        return f(tₑ)
    else return f(t)
    end
end
piecewise(f::Function, t, interval::AbstractVector{Number}) = piecewise(f, t, interval[1], interval[end])
"""
    auxiliary(t, tf::Float64, σ::Float64) 
Return the value of the auxiliary function `b(t)` at time `t` given the final time `tf` and the parameter `σ`.
"""
function auxiliary(t::Float64, tf::Float64, σ::Float64)
    b(t) = 1 +
           10 * (σ - 1) * (t / tf)^3 -
           15 * (σ - 1) * (t / tf)^4 +
           6 * (σ - 1) * (t / tf)^5
    return piecewise(b, t, 0.0, tf)
end
"""
    auxiliary(t::Float64, c::Control)
Return the value of the auxiliary function at time `t` given the control parameter `c`.
It internally defines the parameter `σ` as `σ = (J0 / Jf)^0.25 givent the control parameter `c`.
"""
function auxiliary(t::Float64, c::Control)
    σ = (1/c.Ωf)^0.25 # calculate the value of σ, given the boundary conditions
    return auxiliary(t, c.T, σ) # return the value of the auxiliary function
end
"""
    auxiliary_1d(t::Float64, c::Control)
Return the value of the first derivative of the auxiliary function at time `t` given the control parameter `c`.
"""
function auxiliary_1d(t::Float64, c::Control)
    σ = (1/c.Ωf)^0.25 # calculate the value of σ, given the boundary conditions
    f(t) = 1 / c.T * (
        30 * (σ - 1) * (t / c.T)^2 -
        60 * (σ - 1) * (t / c.T)^3 +
        30 * (σ - 1) * (t / c.T)^4
    )
    return piecewise(f, t, 0.0, c.T)
end
"""
    auxiliary_2d(t::Float64, c::Control)
Return the value of the second derivative of the auxiliary function at time `t` given the control parameter `c`.
"""
function auxiliary_2d(t::Float64, c::Control)
    σ = (1/c.Ωf)^0.25 # calculate the value of σ, given the boundary conditions
    f(t) = 1 / c.T^2 * (
        60 * (σ - 1) * (t / c.T) -
        180 * (σ - 1) * (t / c.T)^2 +
        120 * (σ - 1) * (t / c.T)^3
    )
    return piecewise(f, t, 0.0, c.T)
end
"""
    control_function(t::Float64, c::Control)
Return the control function `J(t)` at time `t` given the control parameter `c`.
The control function is defined as
    J(t) = J₀ / b(t)⁴ - b̈(t) / (2 b(t) U N)
where b(t) is the auxiliary function.
"""
function control_function(t::Float64, c::Control)
    f(t) =  1 / auxiliary(t, c)^4 - auxiliary_2d(t, c) / (2 * auxiliary(t, c) * c.N * c.U)
    return piecewise(f, t, 0.0, c.T)
end
"""
    correction_polyin(c::Control, correction_vector::Array{Float64,1})
Given an array of values `[ y₁, ..., yₙ]`, return the Lagrange polynomial that interpolates the points
`[(0,0) , (t₁, y₁), ..., (tₙ, yₙ), (T, 0)]` at time `t`.
The function is defined in the interval [0, c.final_time]
"""
function correction_polyin(t::Float64, c::Control, correction_vector::Array{Float64,1})
    ys = [0.0; correction_vector; 0.0] # faster than vcat([]...)
    xs = range(0.0, c.T, length=length(ys)) |> collect
    return Lagrange(xs, ys)(t)
end
"""
    control_function(t::Float64, c::Control, correction_vector::Array{Float64,1})
Return the value of the corrected control function, given a list of coefficients `correction_vector`.
The vector is used to obtain the Lagrange polynomial and then add it to the control function.
"""
function control_function(t::Float64, c::Control, correction_vector::Array{Float64,1})
    J_corr(t::Float64) = control_function(t, c) + correction_polyin(t, c, correction_vector)
    return piecewise(J_corr, t, 0.0, c.T)
end
control_function(t::Float64, c::ControlFull) = control_function(t, c, corrections(corrections(c))
