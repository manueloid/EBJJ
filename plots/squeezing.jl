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
function css_ibjj(q::ConstantQuantity, c::ControlSTA)
    in = q.css
    Ω(t) = control_function(t, c)
    H(t, psi) = -2 * Ω(t) * q.Jx + c.U * q.Jz^2
    ΔJ(t, psi) = ΔJ_oat(psi, q) / (c.N / 4) |> real
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
todecibel(x) = 10 * log10(x)

using EBJJ, Plots, DelimitedFiles
t0,tf = 0.05, 1.0
tfs = range(t0, tf, length=102)
N = 50
U, Ωf = 0.4, 0.1
# plot_name = "/home/manueloid/Desktop/sqN_" * string(U) * "_" * string(U * N / (2*Ωf)) * ".png"

function squeezing(c::Control, tfs=AbstractVector{Float64})
    q = ConstantQuantity(c)
    ξs = zeros(length(tfs), 5)
    p = Progress(length(tfs), 1, "Calculating ξN")
    Threads.@threads for index in eachindex(tfs)
        cl = EBJJ.c_time(c, tfs[index])
        corrs = corrections(corrections(cl))
        cs = ControlSTA(cl)
        ξs[index, 1] = css_oat(q, cl)
        ξs[index, 2] = gs_oat(q, cl)
        ξs[index, 3] = css_ibjj(q, cs)
        ξs[index, 4] = gs_ibjj(q, cs)
        ξs[index, 5] = gs_ibjj(q, cl, corrs)
        next!(p)
    end
    return ξs
end
function squeezing(U::Float64, Ωf::Float64, tfs=AbstractVector{Float64})
     N, nλ, max_state = 50, 2, 2
     c = ControlFull(N, Ωf, U, tfs[end], nλ, 2:2:max_state)
     return squeezing(c, tfs)
end
every = squeezing(U, Ωf, tfs)

every_log = todecibel.(every)
tfs = tfs * U

using PGFPlotsX, Colors
colors = (
    black=colorant"#000000",
    red=colorant"#FF0000",
    blue=colorant"#0000FF",
    yellow=colorant"#FFa000",
    green=colorant"#018000"
)
styles = (
    solid="solid",                         
    dot_dash="dash pattern={on 4pt off 1pt on 1pt off 1pt}",
    ldash="dash pattern={on 2pt off 2pt}",
    dash=" dashed",                      
    dot="dash pattern={on 1pt off 1pt}",
)
common_style = @pgf {vertical_sep = "0pt", xticklabels_at = "edge bottom", xlabels_at = "edge bottom"}
esta_opt = @pgf {color = colors.red, line_width = 1, style = styles.solid}
sta_opt = @pgf{color = colors.black, line_width = 1, style = styles.dash}
css_oat_style = @pgf {color = colors.green, line_width = 1, style = styles.dot_dash}
css_ibjj_style = @pgf {color = colors.yellow, line_width = 1, style = styles.ldash}
gs_oat_style = @pgf {color = colors.blue, line_width = 1, style = styles.dot}
PGFPlotsX.enable_interactive(false)
canvas = @pgf Axis({
        name = "canvas",
        width = "8cm",
        height = "8cm",
        xlabel = raw"$\chi t_f$",
        ylabel = raw"$\xi _N^2$[dB]",
        # raw"extra description/.code={\node[below left,inner sep=0pt] at (rel axis cs: -0.08,1.0) {(a)};}",
        xticklabel_style = {
            "scaled ticks=false",
            "/pgf/number format/fixed", 
            "/pgf/number format/precision=3",
            },
        try_min_ticks = 3,
        max_space_between_ticks = "40pt",
        xmax = 0.05,
        # ymax = 0.8,
        # xtick_distance = 0.006 /4,
        enlarge_x_limits = "false",
        enlarge_y_limits = "0.01",
    },
    {}, Plot(esta_opt, Table(tfs, every_log[:, 5])),
    {}, Plot(sta_opt,Table(tfs, every_log[:, 4])),
    {}, Plot(css_oat_style, Table(tfs, every_log[:, 1])),
    {}, Plot(gs_oat_style, Table(tfs, every_log[:, 2])),
    {}, Plot(css_ibjj_style, Table(tfs, every_log[:, 3])),
)
display("/tmp/some.pdf", canvas)
