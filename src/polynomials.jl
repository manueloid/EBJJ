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
auxiliary(c::Control) = t -> auxiliary(t, c) # auxiliary function at the final time
"""
    control_function(auxiliary::Function, ddauxiliary::Function, c::Control)
Given the auxiliary function and its second derivative, it returns the control function as an anonymous function.
"""
function control_function(auxiliary::Function, ddauxiliary::Function, c::Control)
    J(t::Float64) = c.J0 / auxiliary(t)^4 - ddauxiliary(t) / (2 * auxiliary(t) * c.U * c.N) # define the control function
    return J # return the control function
end
