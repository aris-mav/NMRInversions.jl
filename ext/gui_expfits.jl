## Plots for expfits
using GLMakie, NMRInversions
"""
    plot(res::expfit_struct...; kwargs...)
Plot the results of an `expfit` call. \n

This function can take any number of `expfit_struct` structures as input. \n
e.g. `plot(data1)` plots the data and the corresponding fit,
but `plot(data1, data2)` or `plot(data1, data2, data3)` work as well,
and they plot all of the results on the same plot.

Keyword (optional) arguments:
- `markersize` : The size of the markers (default is 7).
- `normeq` : Whether to plot the normalised form of the equation or not (default is `false`).
- `yscale` : The scale of the y-axis (default is `identity`).

"""
function Makie.plot(res::ExpfitData...; kwargs...)

    f = Figure(size=(800, 400))
    plot!(f, res...; kwargs...)

    return f

end

function Makie.plot(res_mat::AbstractVecOrMat{NMRInversions.ExpfitData}; kwargs...)

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
function Makie.plot!(fig::Union{Makie.Figure,Makie.GridPosition}, res::ExpfitData...;
     markersize=7, normeq=false, yscale = identity)

    xlbl = if res[1].input.axes[1] isa PFG
        "b factor (s/m² e-9)"
    else
        "time (s)"
    end

    ax = Axis(
        fig[1, 1], 
        xlabel=xlbl, 
        ylabel="Signal (a.u.)",
        yscale=yscale,
    )

    for r in res

        x = r.input.axes[1]
        y = r.input.data

        if isempty(r.title)
            scatter!(ax, x, y, markersize=markersize)
        else
            scatter!(ax, x, y, markersize=markersize, label= r.title )
        end

        xfit = NMRInversions.logrange( x[1], x[end], 2^10)
        yfit = NMRInversions.mexp(r.u, typeof(r.input.axes[1])(xfit))

        lines!(ax, xfit, yfit, label=getfield(r, normeq == true ? :eqn : :eq))
    end

    fig[2, 1] = Legend(fig, ax, orientation = :vertical,
                       tellwidth = false, tellheight = true,
                       nbanks=(res[1].title == "" ? 1 : 2)
                       )
end
