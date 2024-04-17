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
st_plots, sm_plots, eta_plots = [], [], []
for (values, style) in zip(nsens[1:2:end], styles_plot)
    data = readdlm(values)
    time, sm, st, eta = data[5:end, 1], data[5:end, 2], data[5:end, 3], data[5:end, 5]
    push!(st_plots, Plot(style, Table(time, st)))
    push!(sm_plots, Plot(style, Table(time, sm)))
    push!(eta_plots, Plot(style, Table(time, eta)))
end
    
    

amp_files = filter(file -> occursin("amp", file), files)
time_files = filter(file -> occursin("time", file), files)

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

amp_plots, time_plots = [], []
for (amp, time, style) in zip(amp_files, time_files, styles_plot)
    data = readdlm(amp)
    plot = @pgf Plot(style, Table(data[:, 1] .* 0.05, data[:, 2]))
    push!(amp_plots, plot)
    data = readdlm(time)
    plot = @pgf Plot(style, Table(data[:, 1] .* 0.05, data[:, 2]))
    push!(time_plots, plot)
end
PGFPlotsX.enable_interactive(false)

# The best plan for this moment is to set up all the axis and only in a second time I am going to set up some functions to read the data and make the plots
# This is the big one 
canvas = @pgf Axis({
        name = "canvas",
        width = "\\textwidth",
        xlabel = raw"$\chi t_f$",
        ylabel = raw"$\eta$",
        xticklabel_style = {
            "scaled ticks=false",
            "/pgf/number format/fixed", 
            "/pgf/number format/precision=3",
            },
        try_min_ticks = 4,
        max_space_between_ticks = "40pt",
        xmax = 0.003,
        ymax = 1.0,
        xtick_distance = 0.003 /4,
        enlarge_x_limits = "false",
        enlarge_y_limits = "0.01",
    },
    eta_plots,
)
inset1 = @pgf Axis({
        name = "inset",
        at = "{(canvas.north east)}",
        width = "0.5\\textwidth",
        height = "0.30\\textwidth",
        xshift = "-10pt",
        yshift = "-10pt",
        anchor = "north east",
        xtick = "\\empty",
        xlabel = "{}",
        ylabel = raw"$S_t$",
        xmax = 0.003,
        ymax = 1.0,
        try_min_ticks = 3, 
        enlarge_x_limits = "false",
        enlarge_y_limits = "0.01",
    },
    st_plots
)
# The following two are the insets.
inset2 = @pgf Axis({
        at = "{(inset.south east)}",
        width = "5cm", width = "0.5\\textwidth",
        height = "0.30\\textwidth",
        anchor = "north east",
        xticklabel_style = {
            "scaled ticks=false",
            "/pgf/number format/fixed", 
            "/pgf/number format/precision=3",
            },
        xlabel= raw"$\chi t_f$",
        ylabel = raw"$S_m$",
        xmax = 0.003,
        ymax = 0.5,
        xtick_distance = 0.003/2, # xticklabel_style
        try_min_ticks = 3, 
        enlarge_x_limits = "false",
        enlarge_y_limits = "0.01",
    },
    sm_plots)
time, eta = esta[:, 1], esta[:, 5]
tkz = @pgf TikzPicture({},
    canvas,
    inset1,
    inset2
)
# display("/tmp/plot.pdf", tkz)
display(homedir() * "/Repos/ExternalBJJ/Documents/Paper/Fig_3_sensitivity400_5.pdf", tkz)
