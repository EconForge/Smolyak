using Smolyak
import PyPlot


# function to plot a 2d or 3d grid
function plot(sg::SmolyakGrid)
    g = sg.grid
    mu  = sg.mu
    d = sg.d
    mu_str = isa(mu, Array) ? replace(strip(string(mu)),'\n'," \\times ") : mu
    if size(g, 2) == 2
        PyPlot.scatter(g[:, 1], g[:, 2])
        PyPlot.title("Smolyak grid: \$d=$(d), \\; \\mu=$(mu_str)\$")
    elseif size(g, 2) == 3
        PyPlot.scatter3D(g[:, 1], g[:, 2], g[:, 3], "b.")
        PyPlot.title("Smolyak grid: \$d=$(d), \\; \\mu=$(mu_str)\$")
    else
        error("ERROR: can only plot 2d or 3d grids")
    end
end
