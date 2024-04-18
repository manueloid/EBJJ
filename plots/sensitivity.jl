using DelimitedFiles, PGFPlotsX, Colors
path = homedir() * "/Repos/ExternalBJJ/data/sensitivity/"
filenames = readdir(path)
files = path .* filenames
errtypes = ["time", "amp"]
nsens = filter(file -> occursin("nsens", file), files)
esta = readdlm(nsens[1])
sta = readdlm(nsens[2])
staX = readdlm(nsens[3])

#=
Column 1 time
Column 2 fidelity
Column 3 sensitivity amplitude
Column 4 sensitivity time
Column 5 eta
=#

colors = (
    black=colorant"#000000", # STA	
    red=colorant"#FF0000", # eSTA Full Hamiltonian with Hessian 
    blue=colorant"#0000FF", # eSTA Intermediate Hamiltonian with Hessian
    yellow=colorant"#FFa000", # eSTA Full Hamiltonian with original version
    green=colorant"#018000"# eSTA Intermediate Hamiltonian original version
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
st_plots5, sm_plots5, eta_plots5 = [], [], []
for (values, style) in zip(nsens[2:2:end], styles_plot)
    data = readdlm(values)
    time, sm, st, eta = data[5:end, 1], data[5:end, 2], data[5:end, 3], data[5:end, 5]
    push!(st_plots5, Plot(style, Table(time, st)))
    push!(sm_plots5, Plot(style, Table(time, sm)))
    push!(eta_plots5, Plot(style, Table(time, eta)))
end
st_plots20, sm_plots20, eta_plots20 = [], [], []
for (values, style) in zip(nsens[1:2:end], styles_plot)
    data = readdlm(values)
    time, sm, st, eta = data[5:end, 1], data[5:end, 2], data[5:end, 3], data[5:end, 5]
    push!(st_plots20, Plot(style, Table(time, st)))
    push!(sm_plots20, Plot(style, Table(time, sm)))
    push!(eta_plots20, Plot(style, Table(time, eta)))
end
PGFPlotsX.enable_interactive(false)
# The best plan for this moment is to set up all the axis and only in a second time I am going to set up some functions to read the data and make the plots
# This is the big one 
function plotting(eta, sm, st)
    canvas = @pgf Axis({
        name = "canvas",
        width = "8cm",
        height = "8cm",
        xlabel = raw"$\chi t_f$",
        ylabel = raw"$\eta$",
        # raw"extra description/.code={\node[below left,inner sep=0pt] at (rel axis cs: -0.08,1.0) {(a)};}",
        xticklabel_style = {
            "scaled ticks=false",
            "/pgf/number format/fixed", 
            "/pgf/number format/precision=3",
            },
        try_min_ticks = 3,
        max_space_between_ticks = "40pt",
        xmax = 0.006,
        ymax = 0.8,
        xtick_distance = 0.006 /4,
        enlarge_x_limits = "false",
        enlarge_y_limits = "0.01",
    },
    eta,
    )
    inset1 = @pgf Axis({
        name = "inset",
        at = "{(canvas.north east)}",
        width = "4.5cm",
        height = "3cm",
        xshift = "-13pt",
        yshift = "-11pt",
        anchor = "north east",
        xtick = "\\empty",
        xlabel = "{}",
        ylabel = raw"$S_t$",
        xmax = 0.006,
        ymax = 0.8,
        try_min_ticks = 2, 
        enlarge_x_limits = "false",
        enlarge_y_limits = "0.01",
    },
    st, 
    )
    # The following two are the insets.
    inset2 = @pgf Axis({
        at = "{(inset.south east)}",
        width = "4.5cm",
        height = "3cm",
        anchor = "north east",
        xticklabel_style = {
            "scaled ticks=false",
            "/pgf/number format/fixed", 
            "/pgf/number format/precision=3",
            },
        xlabel= raw"$\chi t_f$",
        # xtick = "\\empty",
        # xlabel = "{}",
        ylabel = raw"$S_m$",
        xmax = 0.006,
        ymax = 0.8,
        try_min_ticks = 2, 
        xtick_distance = 0.006/2, # xticklabel_style
        enlarge_x_limits = "false",
        enlarge_y_limits = "0.01",
    },
    sm)
    return [canvas, inset1, inset2]
end
canvas5, inset1_5, inset2_5 = plotting(eta_plots5, sm_plots5, st_plots5)
canvas20, inset1_20, inset2_20 = plotting(eta_plots20, sm_plots20, st_plots20)
canvas5["xlabel"] = "\\empty"
tkz = @pgf TikzPicture({
    },
    merge!(canvas20, {raw"extra description/.code={\node[below left,inner sep=0pt] at (rel axis cs: -0.08,0.1) {(b)};}"}), 
    inset1_20,
    inset2_20,
    merge!(canvas5, {at = "{(canvas.north west)}", yshift = "15pt", raw"extra description/.code={\node[below left,inner sep=0pt] at (rel axis cs: -0.08,0.1) {(a)};}"}),
    inset1_5,
    inset2_5,
)
# display("/tmp/plot.pdf", tkz)
display(homedir() * "/Repos/ExternalBJJ/Documents/Paper/Fig_4_sensitivity_full.pdf", tkz)
