using EBJJ, QuantumOptics
todecibel(x) = 10 * log10(x)
l(t, c::ControlSTA) = 1 + ( c.Ωf  - 1) * t / c.T
"""
    ξN(q::ConstantQuantity, c::ControlSTA)
Returns the evolution of the squeezing parameter `ξN = ΔJz² / (N/4)` from time 0 to final time `c.T`.
"""
function squeezing_ξ(q::ConstantQuantity, c::ControlSTA; linear::Bool = false)
    Ω(t) = linear ? l(t, c) : control_function(t, c)
    H(t, psi) = -2 * Ω(t) * q.Jx + c.U * q.Jz^2
    ΔJ(t, psi) = (dagger(psi) * (q.Jz)^2 * psi - (dagger(psi) * (q.Jz) * psi)^2) / (c.N / 4) |> real
    return timeevolution.schroedinger_dynamic([0.0:c.T/1000:c.T;], q.ψ0, H; fout=ΔJ)[2]
end
"""
    squeezing_ξX(q::ConstantQuantity, c::ControlSTA)
Returns the evolution of the squeezing parameter `ξN = ΔJz² / (N/4)` from time 0 to final time `c.T`, for the X protocol of the STA control.
"""
function squeezing_ξX(q::ConstantQuantity, c::ControlSTA)
    Ω(t) = control_functionX(t, c)
    H(t, psi) = -2 * Ω(t) * q.Jx + c.U * q.Jz^2
    ΔJ(t, psi) = (dagger(psi) * (q.Jz)^2 * psi - (dagger(psi) * (q.Jz) * psi)^2) / (c.N / 4) |> real
    return timeevolution.schroedinger_dynamic([0.0:c.T/1000:c.T;], q.ψ0, H; fout=ΔJ)[2]
end
"""
    ξN(q::ConstantQuantity, c::ControlFull)
Returns the evolution of the squeezing parameter `ξN = ΔJz² / (N/4)` from time 0 to final time `c.T`.
"""
function squeezing_ξ(q::ConstantQuantity, c::ControlFull)
    corrs = corrections(corrections(c))
    Ω(t) = control_function(t, c, corrs)
    H(t, psi) = -2 * Ω(t) * q.Jx + c.U * q.Jz^2
    ΔJ(t, psi) = (dagger(psi) * (q.Jz)^2 * psi - (dagger(psi) * (q.Jz) * psi)^2) / (c.N / 4) |> real
    return timeevolution.schroedinger_dynamic([0.0:c.T/1000:c.T;], q.ψ0, H; fout=ΔJ)[2]
end
"""
    α(q::ConstantQuantity, c::ControlSTA)
Returns the evolution of the squeezing parameter `α = ⟨ψ| 2Jₓ/N|ψ⟩` from time 0 to final time `c.T`.
"""
function squeezing_α(q::ConstantQuantity, c::ControlSTA; linear::Bool = false)
    Ω(t) = linear ? l(t, c) : control_function(t, c)
    H(t, psi) = -2 * Ω(t) * q.Jx + c.U * q.Jz^2
    squeezing(t, psi) = 2 / N * dagger(psi) * q.Jx * psi |> real
    return timeevolution.schroedinger_dynamic([0.0:c.T/1000:c.T;], q.ψ0, H; fout=squeezing)[2]
end
"""
    squeezing_αX(q::ConstantQuantity, c::ControlSTA)
Returns the evolution of the squeezing parameter `α = ⟨ψ| 2Jₓ/N|ψ⟩` from time 0 to final time `c.T` and for the X version of the STA protocol
"""
function squeezing_αX(q::ConstantQuantity, c::ControlSTA)
    Ω(t) = control_functionX(t, c)
    H(t, psi) = -2 * Ω(t) * q.Jx + c.U * q.Jz^2
    squeezing(t, psi) = 2 / N * dagger(psi) * q.Jx * psi |> real
    return timeevolution.schroedinger_dynamic([0.0:c.T/1000:c.T;], q.ψ0, H; fout=squeezing)[2]
end
"""
    α(q::ConstantQuantity, c::ControlFull)
Returns the evolution of the squeezing parameter `α = ⟨ψ| 2Jₓ/N|ψ⟩` from time 0 to final time `c.T`.
"""
function squeezing_α(q::ConstantQuantity, c::ControlFull)
    corrs = corrections(corrections(c))
    Ω(t) = control_function(t, c, corrs)
    H(t, psi) = -2 * Ω(t) * q.Jx + c.U * q.Jz^2
    squeezing(t, psi) = 2 / N * dagger(psi) * q.Jx * psi |> real
    return timeevolution.schroedinger_dynamic([0.0:c.T/1000:c.T;], q.ψ0, H; fout=squeezing)[2]
end
"""
    fidelity_evo(q::ConstantQuantity, c::ControlSTA)
Returns the evolution of the fidelity from time 0 to final time `c.T`.
"""
function fidelity_evo(q::ConstantQuantity, c::ControlSTA; linear::Bool = false)
    Ω(t) = linear ? l(t, c) : control_function(t, c)
    H(t, psi) = -2 * Ω(t) * q.Jx + c.U * q.Jz^2
    fid(t, psi) = abs2.(dagger(q.ψf) * psi) # Function that calculates the fidelity of the system
    return timeevolution.schroedinger_dynamic([0.0:c.T/1000:c.T;], q.ψ0, H; fout = fid)[2]
end
"""
    fidelity_evoX(q::ConstantQuantity, c::ControlSTA)
Returns the evolution of the fidelity from time 0 to final time `c.T` and for the X protocol of the STA control.
"""
function fidelity_evoX(q::ConstantQuantity, c::ControlSTA)
    Ω(t) = control_functionX(t, c)
    H(t, psi) = -2 * Ω(t) * q.Jx + c.U * q.Jz^2
    fid(t, psi) = abs2.(dagger(q.ψf) * psi) # Function that calculates the fidelity of the system
    return timeevolution.schroedinger_dynamic([0.0:c.T/1000:c.T;], q.ψ0, H; fout = fid)[2]
end
"""
    fidelity_evo(q::ConstantQuantity, c::ControlFull)
Returns the evolution of the fidelity from time 0 to final time `c.T`.
"""
function fidelity_evo(q::ConstantQuantity, c::ControlFull)
    corrs = corrections(corrections(c))
    Ω(t) = control_functionX(t, c, corrs)
    H(t, psi) = -2 * Ω(t) * q.Jx + c.U * q.Jz^2
    fid(t, psi) = abs2.(dagger(q.ψf) * psi) # Function that calculates the fidelity of the system
    return timeevolution.schroedinger_dynamic([0.0:c.T/1000:c.T;], q.ψ0, H; fout = fid)[2]
end

# Data generation 

max_state = 2
nλ = 2
U = 0.4
N = 50
t0, tf = 0.005, 0.35
Ωf = 0.1
tfs = range(t0, tf, length=100) 
c = ControlFull(N, Ωf, U, tf, nλ, 2:2:max_state);
cs = ControlSTA(c);
q = ConstantQuantity(c)

tf = c.T * U
ts = 0.0:tf/1000:tf
ξN_sta = squeezing_ξ(q, cs)
ξN_esta = squeezing_ξ(q, c)
ξN_ad = squeezing_ξ(q, cs, linear = true)
ξN_staX = squeezing_ξX(q, cs)
α_sta = squeezing_α(q, cs)
α_esta = squeezing_α(q, c)
α_ad = squeezing_α(q, cs, linear = true)
α_staX = squeezing_α(q, cs, linear = true)
fid_sta = fidelity_evo(q, cs)
fid_esta = fidelity_evo(q, c)
fid_ad = fidelity_evo(q, cs, linear = true)
fid_staX = fidelity_evoX(q, cs)
ξs_esta = todecibel.(ξN_esta .^ 2 ./ (α_esta .^ 2))
ξs_sta = todecibel.(ξN_sta .^ 2 ./ (α_sta .^ 2))
ξs_ad = todecibel.(ξN_ad .^ 2 ./ (α_ad .^ 2))
ξs_staX = todecibel.(ξN_staX .^ 2 ./ (α_staX .^ 2))
begin 
    corrs = corrections(corrections(c))
    Ω(t) = control_functionX(t, c, corrs)
    cf_esta = Ω.(0.0:c.T/1000:c.T)
end
cf_sta = control_function.(0.0:c.T/1000:c.T, Ref(cs))
cf_ad = l.(0.0:c.T/1000:c.T, Ref(cs))
cf_staX = control_functionX.(0.0:c.T/1000:c.T, Ref(cs))
# Style and plotting
using PGFPlotsX, Colors
colors = (
    black=colorant"#000000", # STA	
    red=colorant"#FF0000", # eSTA Full Hamiltonian with Hessian 
    blue=colorant"#0000FF", # eSTA Intermediate Hamiltonian with Hessian
    yellow=colorant"#FFa000", # eSTA Full Hamiltonian with original version
    green=colorant"#008000"# eSTA Intermediate Hamiltonian original version
)
styles = (
    solid="solid",                         # eSTA Intermediate Hamiltonian with Hessian
    dot_dash="dash pattern={on 4pt off 1pt on 1pt off 1pt}",#STA
    ldash="dash pattern={on 2pt off 2pt}",# eSTA Full Hamiltonian with Hessian 
    dash=" dashed",                        # eSTA Full Hamiltonian with original version
    dot="dash pattern={on 1pt off 1pt}",  # eSTA Intermediate Hamiltonian original version
)
common_style = @pgf {vertical_sep = "0pt", xticklabels_at = "edge bottom", xlabels_at = "edge bottom"}
esta_opt = @pgf {color = colors.red, line_width = 1, style = styles.solid}
sta_opt = @pgf{color = colors.black, line_width = 1, style = styles.dash}
ad_opt = @pgf {color = colors.green, line_width = 1, style = styles.dot_dash}
extra_opt = @pgf {color = colors.yellow, line_width = 1, style = styles.ldash}
styles_plot = [esta_opt, sta_opt, ad_opt, extra_opt]
PGFPlotsX.enable_interactive(false)
gr = @pgf GroupPlot(
    {
        group_style = { 
            group_size = " 1 by 3",
            vertical_sep = "10pt",
            xticklabels_at = "edge bottom",
            xlabels_at = "edge bottom"
            },
        enlarge_y_limits = "0.01",
        enlarge_x_limits = "false",
        ticklabel_style = "/pgf/number format/fixed",
        max_space_between_ticks = "40pt",
        try_min_ticks = 3,
        xtick_distance = "$(tf/3)",
        xticklabel_style = {
            "scaled ticks=false",
            "/pgf/number format/fixed", 
            "/pgf/number format/precision=3",
            },
        width = "\\textwidth",
        height = "0.5\\textwidth",
        # ylabel_style = "at ={(rel axis cs: -0.08,0.4)}",
        clip = false,
    },
    {
        ylabel = "\$\\Omega\$(t)",
        raw"extra description/.code={\node[below left,inner sep=0pt] at (rel axis cs: -0.08,1.0) {(a)};}"
    },
    # Control function
    Plot(esta_opt, Table(ts, cf_esta)),
    Plot(sta_opt, Table(ts, cf_sta)),
    Plot(ad_opt, Table(ts, cf_ad)),
    Plot(extra_opt, Table(ts, cf_staX)),
    # α Squeezing
    {
        ylabel = "\$\\xi_s^2\$ ",
        raw"extra description/.code={\node[below left,inner sep=0pt] at (rel axis cs: -0.08,1.0) {(b)};}"
    },
    Plot(esta_opt, Table(ts, ξs_esta)),
    Plot(sta_opt, Table(ts, ξs_sta)),
    Plot(ad_opt, Table(ts, ξs_ad)),
    Plot(extra_opt, Table(ts, ξs_staX)),
    # Fidelity
    {
        ylabel = "\$F\$",
        xlabel = "\$\\chi t\$",
        ymin = min(fid_esta...),  ymax = 1.0,
        raw"extra description/.code={\node[below left,inner sep=0pt] at (rel axis cs: -0.08,1.0) {(c)};}"
    },
    Plot(esta_opt, Table(ts, fid_esta)),
    Plot(sta_opt, Table(ts, fid_sta)),
    Plot(ad_opt, Table(ts, fid_ad)),
    Plot(extra_opt, Table(ts, fid_staX)),
    )
display(homedir() * "/Repos/ExternalBJJ/Documents/Paper/Fig_2_evolution_negative.pdf", gr)
