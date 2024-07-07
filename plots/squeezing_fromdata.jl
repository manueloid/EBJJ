using DelimitedFiles, PGFPlotsX
N = 400
Λ0 = 2.5
U, Ωf = 2Λ0 / N, 0.1
data_file = "/home/manueloid/Desktop/Ueda_" * string(round(U, digits=2)) * "_" * string(U * N / (2)) * ".dat"
ξN_data1 = readdlm(data_file)
ξN_data1[:, 1] = ξN_data1[:, 1] * U
Λ0 = 10.0
U, Ωf = 2Λ0 / N, 0.1
data_file = "/home/manueloid/Desktop/Ueda_" * string(round(U, digits=2)) * "_" * string(U * N / (2)) * ".dat"
ξN_data2 = readdlm(data_file)
ξN_data2[:, 1] = ξN_data2[:, 1] * U

using PGFPlotsX, Colors
PGFPlotsX.enable_interactive(false)
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

function plots(ξN_data)
    tfs = ξN_data[:, 1]
    pl = @pgf [
        Plot(esta_opt, Table(tfs, ξN_data[:, 6])),
        Plot(staX_opt, Table(tfs, ξN_data[:, 5])),
        Plot(css_oat_style, Table(tfs, ξN_data[:, 2])),
        Plot(gs_oat_style, Table(tfs, ξN_data[:, 3])),
        Plot(sta_opt, Table(tfs, ξN_data[:, 4])),
    ]
    return pl
end
canvas1 = @pgf Axis({
        name = "canvas",
        width = "8cm",
        height = "8cm",
        # xlabel = raw"$\chi t_f$",
        ylabel = raw"$\xi _{N,max}^2$",
        ticklabel_style = {
            "scaled ticks=false",
            "/pgf/number format/fixed",
            "/pgf/number format/precision=3",
        },
        # try_min_ticks = 3,
        max_space_between_ticks = "40pt",
        # xmin = 0.01, 
        ymax = 0.60,
        xmin = 0.0,
        xtick_distance = 0.005 / 5,
        enlarge_x_limits = "false",
        enlarge_y_limits = "0.01",
        raw"extra description/.code={\node[below left,inner sep=0pt] at (rel axis cs: -0.14,1.) {(a)};}"
    },
    plots(ξN_data1),
    "\\node at (rel axis cs:.8,0.9) {\$ \\Lambda_0 = $( N * U / 2) \$};",
)
canvas2 = @pgf Axis({
        name = "canvas",
        width = "8cm",
        height = "8cm",
        # xlabel = raw"$\chi t_f$",
        ylabel = raw"$\xi _{N,max}^2$",
        ticklabel_style = {
            "scaled ticks=false",
            "/pgf/number format/fixed",
            "/pgf/number format/precision=3",
        },
        # try_min_ticks = 3,
        max_space_between_ticks = "40pt",
        # xmin = 0.01, 
        ymax = 0.33,
        xmin = 0.0,
        xtick_distance = 0.005 / 5,
        enlarge_x_limits = "false",
        enlarge_y_limits = "0.01",
        raw"extra description/.code={\node[below left,inner sep=0pt] at (rel axis cs: -0.14,1.) {(b)};}"
    },
    plots(ξN_data2),
    "\\node at (rel axis cs:.8,0.9) {\$ \\Lambda_0 = $( N * U / 2) \$};",
)
plot_name = "/home/manueloid/Desktop/TwistTurnIBJJ_Fig_New_" * string(Λ0) * ".pdf"
# display(plot_name, canvas2)
gp = @pgf GroupPlot(
    {
        group_style = {
            group_size = " 1 by 2",
            vertical_sep = "15pt",
            # xticklabels_at = "edge bottom",
            xlabels_at = "edge bottom"
        },
        ylabel = raw"$\xi _{N,max}^2$",
        xlabel = raw"$\chi t_f$",
        enlarge_y_limits = "0.01",
        enlarge_x_limits = "false",
        ticklabel_style = "/pgf/number format/fixed",
        max_space_between_ticks = "40pt",
        try_min_ticks = 3,
        # xtick_distance="$(tf/4)",
        xticklabel_style = {
            "scaled ticks=false",
            "/pgf/number format/fixed",
            "/pgf/number format/precision=3",
        },
        width = "8cm",
        height = "8cm",
        ylabel_style = "at ={(rel axis cs: -0.14,0.4)}",
        clip = true,
    },
    {
        ymax = 0.60,
        raw"extra description/.code={\node[below left,inner sep=0pt] at (rel axis cs: -0.14,1.) {(a)};}"
    },
    plots(ξN_data1),
    "\\node at (rel axis cs:.8,0.9) {\$ \\Lambda_0 = 2.5\$};",
    {
        ymax = 0.33,
        raw"extra description/.code={\node[below left,inner sep=0pt] at (rel axis cs: -0.14,1.) {(b)};}"
    },
    plots(ξN_data2),
    "\\node at (rel axis cs:.8,0.9) {\$ \\Lambda_0 = $( N * U / 2) \$};",
)
plot_name = "/home/manueloid/Desktop/TwistTurnIBJJ_Fig_New_combined.pdf"
display(plot_name, gp)
