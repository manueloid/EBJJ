using EBJJ, ProgressMeter, QuantumOptics, DelimitedFiles
function ΔJ_oat(ψ::Ket, q::ConstantQuantity)
    Jz, Jy = q.Jz, q.Jy
    num = ψ' * (Jz * Jy + Jy * Jz) * ψ
    den = ψ' * Jz^2 * ψ - ψ' * Jy^2 * ψ
    α = 1 / 2 * atan(num / den) |> real
    J(α) = cos(α)^2 * ψ' * Jz^2 * ψ + sin(α)^2 * ψ' * Jy^2 * ψ + cos(α) * sin(α) * ψ' * (Jz * Jy + Jy * Jz) * ψ |> real
    return min(J(α), J(α + π / 2))
end
function css_oat(q::ConstantQuantity, c::Control)
    in = q.css # Initial state
    H(t, psi) = c.U * q.Jz^2
    ΔJ(t, psi) = ΔJ_oat(psi, q) |> real
    ξN = timeevolution.schroedinger_dynamic([0.0, c.T], in, H; fout=ΔJ)[2][end]
    return ξN
end
function gs_oat(q::ConstantQuantity, c::Control)
    in = q.ψ0
    H(t, psi) = c.U * q.Jz^2
    ΔJ(t, psi) = ΔJ_oat(psi, q) |> real
    ξN = timeevolution.schroedinger_dynamic([0.0, c.T], in, H; fout=ΔJ)[2][end]
    return ξN
end
function gs_ibjj(q::ConstantQuantity, c::ControlSTA)
    in = q.ψ0
    Ω(t) = control_function(t, c)
    H(t, psi) = -2 * Ω(t) * q.Jx + c.U * q.Jz^2
    ΔJ(t, psi) = ΔJ_oat(psi, q) |> real
    ξN = timeevolution.schroedinger_dynamic([0.0, c.T], in, H; fout=ΔJ)[2][end]
    return ξN
end
function gs_ibjjX(q::ConstantQuantity, c::ControlFull, corrs::Vector{Float64})
    in = q.ψ0
    Ω(t) = control_functionX(t, c, corrs)
    H(t, psi) = -2 * Ω(t) * q.Jx + c.U * q.Jz^2
    ΔJ(t, psi) = ΔJ_oat(psi, q) |> real
    ξN = timeevolution.schroedinger_dynamic([0.0, c.T], in, H; fout=ΔJ)[2][end]
    return ξN
end
function gs_ibjjX(q::ConstantQuantity, c::ControlSTA)
    in = q.ψ0
    Ω(t) = control_functionX(t, c)
    H(t, psi) = -2 * Ω(t) * q.Jx + c.U * q.Jz^2
    ΔJ(t, psi) = ΔJ_oat(psi, q) |> real
    ξN = timeevolution.schroedinger_dynamic([0.0, c.T], in, H; fout=ΔJ)[2][end]
    return ξN
end
todecibel(x) = 10 * log10(x)
# do not need to calculate the CSS with the IBJJ hamiltonian
# just add the calculation for the non X version of STA

N = 400
Λ0 = 10
U, Ωf = 2Λ0 / N, 0.1
t0, tf = 0.0001 / U, 0.005 / U
tfs = range(t0, tf, length=202)
data_file = "/home/manueloid/Desktop/Ueda_" * string(round(U, digits=2)) * "_" * string(U * N / (2)) * ".dat"

function squeezing(c::Control, tfs)
    q = ConstantQuantity(c)
    ξN = zeros(length(tfs), 5)
    p = Progress(length(tfs), 1, "Calculating ξN")
    Threads.@threads for index in eachindex(tfs)
        cl = EBJJ.c_time(c, tfs[index])
        corrs = corrections(corrections(cl))
        cs = ControlSTA(cl)
        ξN[index, 1] = css_oat(q, cl)
        ξN[index, 2] = gs_oat(q, cl)
        ξN[index, 3] = gs_ibjj(q, cs)
        ξN[index, 4] = gs_ibjjX(q, cs)
        ξN[index, 5] = gs_ibjjX(q, cl, corrs)
        next!(p)
    end
    return ξN
end
function squeezing(U::Float64, Ωf::Float64,N::Int, tfs=AbstractVector{Float64})
    nλ, max_state =  2, 2
    c = ControlFull(N, Ωf, U, tfs[end], nλ, 2:2:max_state)
    return squeezing(c, tfs) 
end
ξN = squeezing(U, Ωf,N, tfs) ./ (N / 4)

writedlm(data_file, hcat(tfs,ξN))

ξN_data = hcat(tfs,ξN)
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
common_style = @pgf {vertical_sep = "0pt", xticklcbels_at = "edge bottom", xlabels_at = "edge bottom"}
esta_opt = @pgf {color = colors.red, line_width = 1, style = styles.solid}
sta_opt = @pgf{color = colors.black, line_width = 1, style = styles.dash}
css_oat_style = @pgf {color = colors.green, line_width = 1, style = styles.dot_dash}
staX_opt = @pgf {color = colors.yellow, line_width = 1, style = styles.ldash}
gs_oat_style = @pgf {color = colors.blue, line_width = 1, style = styles.dot}
#=
Breakdown of colors
- black: STA index 4
- red: eSTA, index 5
- blue: ψ0 with OAT, index 2
- green: CSS with OAT, index 1
- yellow: STA non X, index 3
=#
PGFPlotsX.enable_interactive(false)
canvas = @pgf Axis({
        name = "canvas",
        width = "8cm",
        height = "8cm",
        xlabel = raw"$\chi t_f$",
        ylabel = raw"$\xi _{N,max}^2$",
        # raw"extra description/.code={\node[below left,inner sep=0pt] at (rel axis cs: -0.08,1.0) {(a)};}",
        ticklabel_style = {
            "scaled ticks=false",
            "/pgf/number format/fixed",
            "/pgf/number format/precision=3",
        },
        # try_min_ticks = 3,
        max_space_between_ticks = "40pt",
        # xmin = 0.01, 
        ymax = 0.6,
        xmin = 0.0,
        xtick_distance = 0.005 /5,
        enlarge_x_limits = "false",
        enlarge_y_limits = "0.01",
    },
    {}, Plot(esta_opt, Table(tfs * U, ξN_data[:, 6])),
    {}, Plot(staX_opt, Table(tfs * U, ξN_data[:, 5])),
    {}, Plot(css_oat_style, Table(tfs * U, ξN_data[:, 2])),
    {}, Plot(gs_oat_style, Table(tfs * U, ξN_data[:, 3])),
    {}, Plot(sta_opt, Table(tfs * U, ξN_data[:, 4])),
    VLine({dashed, black}, tfs[54]* U),
    "\\node at (rel axis cs:.8,0.75) {\$ \\Lambda_0 = $( N * U / 2) \$};",
)
plot_name = "/home/manueloid/Desktop/Ueda_" * string(round(U, digits=2)) * "_" * string(U * N / (2)) * ".pdf"
# plot_name = "/tmp/fig.pdf"
display(plot_name, canvas)

sqs = todecibel.(ξN) ./ Jx
plot_name = "/home/manueloid/Desktop/Wineland_" * string(U) * "_" * string(U * N / (2 * Ωf)) * ".pdf"
canvas = @pgf Axis({
        name = "canvas",
        width = "8cm",
        height = "8cm",
        xlabel = raw"$\chi t_f$",
        ylabel = raw"$\xi _s^2$[dB]",
        # raw"extra description/.code={\node[below left,inner sep=0pt] at (rel axis cs: -0.08,1.0) {(a)};}",
        ticklabel_style = {
            "scaled ticks=false",
            "/pgf/number format/fixed",
            "/pgf/number format/precision=3",
        },
        try_min_ticks = 3,
        max_space_between_ticks = "40pt",
        # xmax = 0.05,
        # ymax = 0.8,
        # xtick_distance = 0.006 /4,
        enlarge_x_limits = "false",
        enlarge_y_limits = "0.01",
    },
    {}, Plot(esta_opt, Table(tfs, sqs[:, 5])),
    {}, Plot(sta_opt, Table(tfs, sqs[:, 4])),
    {}, Plot(css_oat_style, Table(tfs, sqs[:, 1])),
    {}, Plot(gs_oat_style, Table(tfs, sqs[:, 2])),
    {}, Plot(css_ibjj_style, Table(tfs, sqs[:, 3])),
)
display(plot_name, canvas)



