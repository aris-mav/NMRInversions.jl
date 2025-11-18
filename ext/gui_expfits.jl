## Plots for expfits
"""
    plot(res::expfit_struct...; kwargs...)
Plot the results of an `expfit` call. \n

This function can take any number of `expfit_struct` structures as input. \n
e.g. `plot(data1)` plots the data and the corresponding fit,
but `plot(data1, data2)` or `plot(data1, data2, data3)` work as well,
and they plot all of the results on the same plot.

Keyword (optional) arguments:
- `markersize` : The size of the markers (default is 7).
- `normeq` : Whether to plot the normalised form of the equation or not (default is `true`).
- `yscale` : The scale of the y-axis (default is `identity`).

"""
function Makie.plot(res::expfit_struct...; kwargs...)

    f = Figure(size=(800, 400))
    plot!(f, res...; kwargs...)

    return f

end

function Makie.plot(res_mat::AbstractVecOrMat{NMRInversions.expfit_struct}; kwargs...)

    f = Figure(size=(800, 400))
    plot!(f, res_mat...; kwargs...)

    return f

end

"""
    plot!(fig, res...; names, markersize, normeq)

Plots the results of an `expfit` call on a figure or a grid position object.

The arguments are:

- `fig` : The figure or grid position object.
- `res` : One or more `expfit_struct` structures containing the fit results.

Keyword (optional) arguments:
- `markersize` : The size of the markers (default is 7).
- `normeq` : Whether to plot the normalised form of the equation or not (default is `true`).
- `yscale` : The scale of the y-axis (default is `identity`).

Note that the res inputs are not a vector, but individual `expfit_struct` structures,
like: plot!(fig, data1, data2, data3).
If you want to use a vector of `expfit_struct` structures, make sure to 
splat it by using `...` in the function call (e.g. `plot!(fig, [data1, data2, data3]...)`).
"""
function Makie.plot!(fig::Union{Makie.Figure,Makie.GridPosition}, res::expfit_struct...;
     markersize=7, normeq=true, yscale = identity)

    # Make axes
    if res[1].seq in [NMRInversions.IR]
        ax = Axis(fig[1, 1], xlabel="time (s)", ylabel="Signal (a.u.)",yscale=yscale)

    elseif res[1].seq in [NMRInversions.CPMG]
        ax = Axis(fig[1, 1], xlabel="time (s)", ylabel="Signal (a.u.)",yscale=yscale)

    elseif res[1].seq in [NMRInversions.PFG]
        ax = Axis(fig[1, 1], xlabel="b factor (s/m² e-9)", ylabel="Signal (a.u.)",yscale=yscale)
    end

    for r in res
        if isempty(r.title)
            scatter!(ax, r.x, r.y, markersize=markersize)
        else
            scatter!(ax, r.x, r.y, markersize=markersize, label= r.title )
        end
        lines!(ax, r.xfit, r.yfit, label=getfield(r, normeq == true ? :eqn : :eq))
    end

    axislegend(ax, position=(res[1].seq == IR ? :rb : :rt), nbanks=2)

end
