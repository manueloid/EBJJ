"""
`piecewise(t, cp::Control, f::Function)` take a function `f` and generates a piecewise function of the form 
`f(t) = f(0) if t < 0
		f(t) if 0 <= t <= tf
		f(tf) if t > tf`
"""
function piecewise(t, tf::Float64, f::Function)
    if t < 0
        return f(0.0)
    elseif t > tf
        return f(tf)
    else
        return f(t)
    end
end
"""
`auxiliary(t, tf::Float64, σ::Float64)` 
Return the value of the auxiliary function at time `t`, for the boundary conditions σ.
"""
function auxiliary(t, tf::Float64, σ::Float64)
    b(t) = 1 +
           10 * (σ - 1) * (t / tf)^3 -
           15 * (σ - 1) * (t / tf)^4 +
           6 * (σ - 1) * (t / tf)^5
    return piecewise(t, tf, b)
end
"""
`auxiliary(t, tf::Float64, J0::Float64, Jf::Float64)`
	Return the value of the auxiliary function at time `t` given the boundary conditions `J0` and `Jf`.
"""
function auxiliary(t, tf::Float64, J0::Float64, Jf::Float64)
    σ = (J0 / Jf)^0.25 # calculate the value of σ, given the boundary conditions
    return auxiliary(t, tf, σ) # return the value of the auxiliary function
end
"""
`auxiliary(t, c::Control)`
	Return the auxiliary function at time `t`, given the control parameter `c`.
"""
function auxiliary(t, c::Control)
    return auxiliary(t, c.T, c.J0, c.Jf) # return the value of the auxiliary function
end
"""
`control_function(t, tf::Float64,  J0::Float64, Jf::Float64; N::Int64=100, U::Float64=0.49)`
Given the boundary conditions of the control function, the total time of the process and the number of particles and the interaction strength, it returns the value of the control function at time `t`.
"""
function control_function(t, tf::Float64, J0::Float64, Jf::Float64, N::Int64=30, U::Float64=0.02)
    σ = (J0 / Jf)^0.25 # calculate the value of σ, given the boundary conditions
    b(t) = auxiliary(t, tf, σ) # define the auxiliary function
    db(t) = ForwardDiff.derivative(b, t) # define the first derivative of the auxiliary function
    ddb(t) = ForwardDiff.derivative(db, t) # define the second derivative of the auxiliary function
    J(t) = J0 / b(t)^4 - ddb(t) / (2 * b(t) * U * N) # define the control function
    return piecewise(t, tf, J) # return the value of the control function as a piecewise function
end
"""
`control_function(t, c::Control)`
Given the control parameter `c` returns the value of the polynomials at time `t`.
This is the multiple dispatch version of the previous function `control_function`.
"""
control_function(t, c::Control) = control_function(t, c.T, c.J0, c.Jf, c.N, c.U)
