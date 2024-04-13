using DelimitedFiles, PGFPlotsX, Colors
path = homedir() * "/Repos/ExternalBJJ/data/sensitivity/"
filenames = readdir(path)
files = path .* filenames
errtypes = ["time", "amp"]
amp_files = filter(file -> occursin("amp", file), files)
time_files = filter(file -> occursin("time", file), files)

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
styles_plot = [esta_opt, sta_opt, extra_opt]

amp_plots, time_plots = [], []
for (amp, time, style) in zip(amp_files, time_files, styles_plot)
    data = readdlm(amp)
    plot = @pgf Plot(style, Table(data[:, 1], data[:, 2]))
    push!(amp_plots, plot)
    data = readdlm(time)
    plot = @pgf Plot(style, Table(data[:, 1], data[:, 2]))
    push!(time_plots, plot)
end

gr = @pgf GroupPlot(
    {
        group_style = {
            group_size = " 1 by 2",
            vertical_sep = "10pt",
            horizontal_sep = "10pt",
            xticklabels_at = "edge bottom",
            xlabels_at = "edge bottom",
            yticklabels_at = "edge left",
            ylabels_at = "edge left",
        },
        enlarge_y_limits = "0.01",
        enlarge_x_limits = "false",
        ticklabel_style = "/pgf/number format/fixed",
        max_space_between_ticks = "84pt",
        try_min_ticks = 4,
        # xtick_distance = "$(c.T/2)",
        width = "\\textwidth",
        height = "0.5\\textwidth",
        # ylabel_style = "at ={(rel axis cs: -0.18,0.5)}",
        # clip = false,
    },
    {ymax = 1, xmin = 0.01,xmax = 0.15, ylabel = raw"$S_m$",
    raw"extra description/.code={\node[below left,inner sep=0pt] at (rel axis cs: -0.08,1.0) {(a)};}"
    },
    amp_plots,
    {ymax = 1, xmin = 0.01,xmax = 0.15, ylabel = raw"$S_t$", xlabel = raw"$\chi t$",
    raw"extra description/.code={\node[below left,inner sep=0pt] at (rel axis cs: -0.08,1.0) {(b)};}"
    },
    time_plots
)
display(homedir() * "/Repos/ExternalBJJ/Documents/Paper/gfx/sensitivity.pdf", gr)
