using DelimitedFiles, PGFPlotsX, Colors
path = homedir() * "/Repos/ExternalBJJ/data/fidelity/"
filenames = readdir(path)
files = path .* filenames
pattern = r"nfid(\d+)([a-zA-Z]+)(\d+)\.dat"
match_result = match(pattern, files[2])
match_result.captures
fileloc(nn, type, u) = path * "nfid" * string(nn) * type * replace(string(u), "." => "") * ".dat"
types = ["estaX", "sta", "staX", "refl"]

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
styles_plot = [esta_opt, sta_opt, extra_opt, ad_opt]

xmax(val) = @pgf {xmax = val}
PGFPlotsX.enable_interactive(false)
function tabulars(nn, uu)
    tabs = []
    for type in types
        file = fileloc(nn, type, uu)
        data = readdlm(file)
        push!(tabs, data)
    end
    return tabs
end
function plottering(nn, uu)
    plots, ymin = [], 0.0
    for (type, style) in zip(types, styles_plot)
        file = fileloc(nn, type, uu)
        data = readdlm(file)
        plot = @pgf Plot(style, Table(data[1: end, 1], data[1:end, 2]))
        push!(plots, plot)
        ymin = data[:, 2][1]
    end
    println(ymin)
    # push!(plots, ["\\node at (rel axis cs:.8,0.85) {\$ N = $nn \$};"])
    push!(plots, ["\\node at (rel axis cs:.8,0.25) {\$ \\Lambda_0 = $( nn * uu / 2) \$};"])
    ax = @pgf Axis({ymin = ymin}, plots)
    return ax
end

xlabel = @pgf {xlabel = raw"$\chi t_f$"}
ylabel = @pgf {ylabel = "F"}
title = @pgf {title = raw"$\chi t_f$"}
xm50 = @pgf {xmin = 0.0, xmax = 0.04, xtick_distance = "$(0.04 / 2)"}
xm200 = @pgf {xmin = 0.0, xmax = 0.01, xtick_distance = "$(0.01 / 2)"}
xm400 = @pgf {xmin = 0.0, xmax = 0.005, xtick_distance = "$(0.005 / 2)"}
ym10 = @pgf {ymin = 0.85, ymax = 1.0}
ym5 = @pgf {ymin = 0.85, ymax = 1.0}
ym2 = @pgf {ymin = 0.85, ymax = 1.0}
gr = @pgf GroupPlot(
    {
        group_style = {
            group_size = " 3 by 3",
            vertical_sep = "10pt",
            horizontal_sep = "20pt",
            xticklabels_at = "edge bottom",
            xlabels_at = "edge bottom",
            yticklabels_at = "edge left",
            ylabels_at = "edge left",
        },
        xticklabel_style = {
            "scaled ticks=false",
            "/pgf/number format/fixed", 
            "/pgf/number format/precision=4",
            },
        minor_tick_num = 2,
        enlarge_y_limits = "0.01",
        enlarge_x_limits = "false",
        ticklabel_style = "/pgf/number format/fixed",
        max_space_between_ticks = "60pt",
        try_min_ticks = 4,
        ymax = 1.0,
        width = "7cm",
        height = "7cm",
        # xtick_distance = "$(c.T/2)",
        # width = "0.5\\textwidth",
        # height = ".5\\textwidth",
        # ylabel_style = "at ={(rel axis cs: -0.18,0.5)}",
        # clip = false,
    },
    merge!(plottering(50, 0.4), ylabel, xm50, ym10,  {title = "N = 50"}),
    merge!(plottering(200, 0.1), xm200, ym10, {title = "N = 200"}),
    merge!(plottering(400, 0.05), xm400, ym10, {title = "N = 400"}),
    # Second row
    merge!(plottering(50, 0.2), ylabel, xm50, ym5),
    merge!(plottering(200, 0.05), xm200, ym5),
    merge!(plottering(400, 0.025), xm400, ym5),
    # Third row
    merge!(plottering(50, 0.1), ylabel, xlabel, xm50, ym2 ),
    merge!(plottering(200, 0.025), xlabel, xm200, ym2),
    merge!(plottering(400, 0.0125), xlabel, xm400, ym2),
)
# display("/tmp/fig.pdf", gr)
display(homedir() * "/Repos/ExternalBJJ/Documents/Paper/TwistTurnIBJJ_Fig4.pdf", gr)
