using DelimitedFiles, PGFPlotsX
U,N = 0.1, 50
data_file = "/home/manueloid/Desktop/Ueda_" * string(round(U, digits=2)) * "_" * string(U * N / (2)) * ".dat"
ξN = readdlm(data_file)
tfs = range(0.001, 0.04, length = size(ξN, 1))
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
canvas = @pgf Axis({
        name = "canvas",
        width = "8cm",
        height = "8cm",
        xlabel = raw"$\chi t_f$",
        ylabel = raw"$\xi _N^2$",
        # raw"extra description/.code={\node[below left,inner sep=0pt] at (rel axis cs: -0.08,1.0) {(a)};}",
        ticklabel_style = {
            "scaled ticks=false",
            "/pgf/number format/fixed",
            "/pgf/number format/precision=3",
        },
        try_min_ticks = 3,
        max_space_between_ticks = "40pt",
        # xmin = 0.01, 
        ymax = 0.6,
        xtick_distance = 0.04 /2,
        enlarge_x_limits = "false",
        enlarge_y_limits = "0.01",
    },
    {}, Plot(esta_opt, Table(tfs, ξN[:, 5])),
    {}, Plot(staX_opt, Table(tfs, ξN[:, 4])),
    {}, Plot(css_oat_style, Table(tfs, ξN[:, 1])),
    {}, Plot(gs_oat_style, Table(tfs, ξN[:, 2])),
    {}, Plot(sta_opt, Table(tfs, ξN[:, 3])),
    "\\node at (rel axis cs:.8,0.75) {\$ \\Lambda_0 = $( N * U / 2) \$};"
)
plot_name = "/home/manueloid/Desktop/Ueda_" * string(round(U, digits=2)) * "_" * string(U * N / (2)) * ".pdf"
display(plot_name, canvas)
