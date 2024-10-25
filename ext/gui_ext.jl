module gui_ext

using NMRInversions, GLMakie, PolygonOps, LinearAlgebra, NativeFileDialog


## Plots for expfits
"""
    plot(res::expfit_struct...; kwargs...)
Plot the results of an `expfit` call. \n

This function can take any number of `expfit_struct` structures as input. \n
e.g. `plot(data1)` plots the data and the corresponding fit,
but `plot(data1, data2)` or `plot(data1, data2, data3)` work as well,
and they plot all of the results on the same plot.

This function creates a figure, 
and calls `plot!(fig, res)` for each of the results in `res`,
so please refer to the `plot!` documentation for more information.
"""
function Makie.plot(res::expfit_struct...; kwargs...)

    f = Figure(size=(500, 500))
    plot!(f, res...; kwargs...)

    return f

end

"""
    plot!(fig, res...; names, markersize, normeq)

Plots the results of an `expfit` call on a figure or a grid position object.

The arguments are:

- `fig` : The figure or grid position object.
- `res` : One or more `expfit_struct` structures containing the fit results.
- `names` : Vector with the names of the data (default is a vector of "Data" for each result).
- `markersize` : The size of the markers (default is 7).
- `normeq` : Whether to plot the normalised form of the equation or not (default is `true`).

Note that the res inputs are not a vector, but individual `expfit_struct` structures,
like: plot!(fig, data1, data2, data3).
If you want to use a vector of `expfit_struct` structures, make sure to 
splat it by using `...` in the function call (e.g. `plot!(fig, [data1, data2, data3]...)`).
"""
function Makie.plot!(fig::Union{Makie.Figure,Makie.GridPosition}, res::expfit_struct...;
    names=["Data" for _ in res], markersize=7, normeq=true)

    # Make axes
    if res[1].seq in [NMRInversions.IR]
        ax = Axis(fig[1, 1], xlabel="time (s)", ylabel="Signal (a.u.)")

    elseif res[1].seq in [NMRInversions.CPMG]
        ax = Axis(fig[1, 1], xlabel="time (s)", ylabel="Signal (a.u.)")

    elseif res[1].seq in [NMRInversions.PFG]
        ax = Axis(fig[1, 1], xlabel="b factor (s/m² e-9)", ylabel="Signal (a.u.)")
    end

    for (i, r) in enumerate(res)
        xfit = collect(range(0, 1.1 * maximum(r.x), 500))
        scatter!(ax, r.x, r.y, markersize=markersize, label=names[i])
        lines!(ax, xfit, mexp(r.seq, r.u, xfit), label=getfield(r, normeq == true ? :eqn : :eq))
    end

    axislegend(ax, position=(res[1].seq == IR ? :rb : :rt), nbanks=2)

end


## Plots for 1D inversions

"""
    plot(res::NMRInversions.inv_out_1D...; kwargs...)

Plot the results contained in a `inv_out_1D` structure.
This function can take any number of `inv_out_1D` structures as input. 

"""
function Makie.plot(res::NMRInversions.inv_out_1D...; kwargs...)

    f = Figure(size=(500, 500))
    plot!(f, res...; kwargs...)

    return f

end


"""
    plot!(fig, res...; title)

Plot the results contained in a `inv_out_1D` structure on a figure or a grid position object.
This function can take any number of `inv_out_1D` structures as input.

The arguments are:

- `fig` : The figure or grid position object.
- `res` : One or more `inv_out_1D` structures containing the fit results.

Keyword (optional) arguments:
- `title` : Title of the plot.
"""
function Makie.plot!(fig::Union{Makie.Figure,Makie.GridPosition}, res::NMRInversions.inv_out_1D...;
    title="")

    # Make axes
    if res[1].seq in [NMRInversions.IR]
        ax1 = Axis(fig[:, 1], xlabel="time (s)", ylabel="Signal (a.u.)")
        ax2 = Axis(fig[:, 2], xlabel="T₁ (s)", xscale=log10)

    elseif res[1].seq in [NMRInversions.CPMG]
        ax1 = Axis(fig[:, 1], xlabel="time (s)", ylabel="Signal (a.u.)")
        ax2 = Axis(fig[:, 2], xlabel="T₂ (s)", xscale=log10)

    elseif res[1].seq in [NMRInversions.PFG]
        ax1 = Axis(fig[:, 1], xlabel="b factor (s/m² e-9)", ylabel="Signal (a.u.)")
        ax2 = Axis(fig[:, 2], xlabel="D (m²/s)", xscale=log10)
    end

    for r in res
        draw_on_axes(ax1, ax2, r)
    end

end

function draw_on_axes(ax1, ax2, res::NMRInversions.inv_out_1D)

    scatter!(ax1, res.x, res.y)
    lines!(ax1, res.xfit, res.yfit)
    lines!(ax2, res.X, res.f)

end


## Plots for 2D inversions

"""
    plot(results::AbstractMatrix{inv_out_2D}; kwargs...)

Plot the results contained in a matrix of `inv_out_2D` structures.
The arrangement is determined by the shape of the matrix.

The arguments are:

- `results` : The `inv_out_2D` structure containing the fit results.
kwargs are passed onto `plot!`.
"""
function Makie.plot(
    res_mat::AbstractMatrix{NMRInversions.inv_out_2D},
    kwargs...)

    f = Figure(size=(400, 400) .* reverse(size(res_mat)))

    for (ind, res) in pairs(res_mat)
        plot!(f[Tuple(ind)...], res; title=res.title, kwargs...)
    end

    return f

end


"""
    plot!(fig, res::NMRInversions.inv_out_2D; title, clmap)

Plot the results contained in a `inv_out_2D` structure on a figure or a grid position object.

The arguments are:

- `fig` : The figure or grid position object.
- `res` : The `inv_out_2D` structure containing the fit results.

Keyword (optional) arguments:
- `title` : Title of the plot (default: "").
- `clmap` : Color map of the plot (default: :tempo).
"""
function Makie.plot!(fig::Union{Makie.Figure,Makie.GridPosition}, res::NMRInversions.inv_out_2D;
    title="", clmap=:tempo)

    x = res.X_indir
    y = res.X_dir

    xplot = collect(range(0, 1, length(x)))
    yplot = collect(range(0, 1, length(y)))

    # Make axes

    axmain = Axis(fig[3:10, 1:8])
    axtop = Axis(fig[1:2, 1:8], limits=((xplot[1], xplot[end]), (0, nothing)), title=title, titlesize=17)
    axright = Axis(fig[3:10, 9:10], limits=((0, nothing), (yplot[1], yplot[end])))

    hidedecorations!(axtop)
    hidedecorations!(axright)
    hidedecorations!(axmain)

    axmain_values = Axis(fig[3:10, 1:8],
        xscale=log10, yscale=log10,
        xlabel=L"T_1 \, \textrm{(s)}", ylabel=L"T_2 \,\textrm{(s)}",
        xlabelsize=23, ylabelsize=23,
        limits=(x[1], x[end], y[1], y[end])
    )

    Makie.deactivate_interaction!(axmain_values, :rectanglezoom) # Disable zoom
    Makie.deactivate_interaction!(axmain, :rectanglezoom) # Disable zoom
    Makie.deactivate_interaction!(axright, :rectanglezoom) # Disable zoom
    Makie.deactivate_interaction!(axtop, :rectanglezoom) # Disable zoom

    linkxaxes!(axmain, axtop)
    linkyaxes!(axmain, axright)

    draw_on_axes(axmain, axtop, axright, res, clmap)

end

function draw_on_axes(axmain, axtop, axright, res, clmap)

    empty!(axmain)
    empty!(axtop)
    empty!(axright)

    z = res.F' .* res.filter'
    x = range(0, 1, size(z, 1))
    y = range(0, 1, size(z, 2))

    # Plots
    contourf!(axmain, x, y, z, colormap=clmap, levels=50)
    lines!(axtop, x, vec(sum(z, dims=2)), colormap=clmap, colorrange=(1, 10), color=5, alpha=0.5)
    lines!(axright, vec(sum(z, dims=1)), y, colormap=clmap, colorrange=(1, 10), color=5, alpha=0.5)

    # Plot diagonal line
    lines!(axmain, [(0, 0), (1, 1)], color=:black, linewidth=1)

    #Create a matrix for all the discrete points in the space
    points = [[i, j] for i in x, j in y]
    mask = zeros(size(points))

    for (i, polygon) in enumerate(res.selections)

        # Selected points only
        mask .= [PolygonOps.inpolygon(p, polygon; in=1, on=1, out=0) for p in points]
        spo = mask .* z

        indir_dist = vec(sum(spo, dims=2))
        dir_dist = vec(sum(spo, dims=1))

        #draw current polygon
        lines!(axmain, Point2f.(polygon), linestyle=:dash, colormap=:tab10, colorrange=(1, 10), color=i, alpha=0.9)

        T1 = indir_dist ⋅ res.X_indir / sum(spo)
        T2 = dir_dist ⋅ res.X_dir / sum(spo)

        # Plot peak label
        xc = vec(sum(spo, dims=2)) ⋅ x / sum(spo)
        yc = vec(sum(spo, dims=1)) ⋅ y / sum(spo)

        sc = scatter!(axmain, xc, yc, markersize=15,
            marker=collect('a':'z')[i],
            colormap=:tab10, colorrange=(1, 10),
            color=i,
            glowcolor=:white, glowwidth=4,
            label=" : T₁/T₂ = $(round(T1/T2 , digits=1)) \n    Volume = $(round(sum(spo)/sum(z) *100 ,digits = 1))%")


        text!(axmain, 2.5 * (1 / 100), (99 - 13 * (i - 1)) * (1 / 100),
            text=
            collect('a':'z')[i] *
            " : T₁/T₂ = $(round(T1/T2 , digits=1)) \n     Volume = $(round(sum(spo)/sum(z) *100 ,digits = 1))%",
            align=(:left, :top), colormap=:tab10, colorrange=(1, 10), color=i)

        # lines!(axtop, indir_dist[], linestyle=:dash, colormap=:tab10, colorrange=(1, 10), color=i, alpha=0.5)
        # lines!(axright, dir_dist[], 1:length(res.X_dir), linestyle=:dash, colormap=:tab10, colorrange=(1, 10), color=i, alpha=0.5)
        # axislegend( axmain,  framevisible=false)

    end

end

function dynamic_plots(axmain, axtop, axright, selection, polygon)

    scatter!(axmain, selection, colormap=:tab10, colorrange=(1, 10), color=8)
    lines!(axmain, polygon, colormap=:tab10, colorrange=(1, 10), color=8)

end

function draw_3d(ax3d, res, sp)

    surface!(ax3d, res.F' .* res.filter', colormap=:tempo)

    surface!(ax3d, spo, colormap=:blues,
        colorrange=(minimum(filter(x -> x != 0, vec(z))), maximum(z)),
        lowclip=:transparent)

end


"""
    Makie.plot(res::inv_out_2D)
Run the GUI to plot the results and select peaks you want to label.
"""
function Makie.plot(res::NMRInversions.inv_out_2D)
    # begin
    gui = Figure(size=(900, 500))
    Makie.plot!(gui[2:10, 1:9], res)

    axmain = gui.content[1]
    axtop = gui.content[2]
    axright = gui.content[3]
    # Legend(gui[5, 10:19], axmain)
    
    
    # Buttons
    labelb = Button(gui[1, 10:15]; label="Label current selection")
    clearb = Button(gui[2, 10:15]; label="Clear current selection")
    resetb = Button(gui[1, 15:19]; label="Reset everything")
    filterb = Button(gui[2, 15:19]; label="Filter-out unselected")
    saveb = Button(gui[3, 10:19]; label="Save and exit")

    # Title textbox
    tb = Textbox(gui[1, 2:7], 
                 placeholder="Insert a title for the plot, then press enter.", 
                 width=300, 
                 reset_on_defocus=true,
                 )

    z = res.F' .* res.filter'
    x = collect(range(0, 1, size(z, 1)))
    y = collect(range(0, 1, size(z, 2)))

    #Create a matrix for all the discrete points in the space
    points = [[i, j] for i in x, j in y] #important, used by inpolygon later
    mask = Observable(zeros(size(points)))

    # Selected points only
    spo = @lift(($mask .* $z))

    # Selected points Direct dimension distribution
    dir_dist = @lift(vec(sum($spo, dims=1)))
    indir_dist = @lift(vec(sum($spo, dims=2)))

    # Initialise Selection vectors
    selection = Observable(Point{2,Float32}[])
    polygon = Observable(Point{2,Float32}[])

    # Dynamic plots
    dynamic_plots(axmain, axtop, axright, selection, polygon)


    begin ## SELECTING POINTS

        # Pick points by clicking the plot
        spoint = select_point(gui.content[1].scene, marker=:circle)
        # Update selection Observable by pushing the selected point to it
        on(spoint) do _
            push!(selection[], spoint[])
            selection[] = selection[] # Update explicitly so that the plots are updated automatically

            if size(selection[], 1) > 2
                polygon[] = vcat(selection[], [selection[][1]])
                mask[] = [PolygonOps.inpolygon(p, polygon[]; in=1, on=1, out=0) for p in points]
            end
        end

    end ## SELECTING POINTS


    begin ## BUTTON CLICKS

        on(labelb.clicks) do _

            if size(selection[], 1) < 3
                @warn("You need to make a selection.")
            else

                push!(res.selections, polygon[])

                draw_on_axes(axmain, axtop, axright, res, :tempo)
                dynamic_plots(axmain, axtop, axright, selection, polygon)

                # Clear for next selection
                selection[] = Point{2,Float32}[]
                polygon[] = Point{2,Float32}[]
                mask[] = zeros(size(points))

                reset_limits!(axmain)
            end

        end

        on(clearb.clicks) do _

            selection[] = Point{2,Float32}[]
            polygon[] = Point{2,Float32}[]
            mask[] = zeros(size(points))

        end

        on(filterb.clicks) do _

            if size(selection[], 1) < 3
                @warn("You need to make a selection.")
            else

                res.filter .= res.filter .* [PolygonOps.inpolygon(p, polygon[]; in=1, on=1, out=0) for p in points]'
            end

            draw_on_axes(axmain, axtop, axright, res, :tempo)
            dynamic_plots(axmain, axtop, axright, selection, polygon)

            selection[] = Point{2,Float32}[]
            polygon[] = Point{2,Float32}[]
            mask[] = zeros(size(points))

        end


        on(resetb.clicks) do _
            res.filter .= ones(size(res.F))
            empty!(res.selections)

            draw_on_axes(axmain, axtop, axright, res, :tempo)
            dynamic_plots(axmain, axtop, axright, selection, polygon)

            selection[] = Point{2,Float32}[]
            polygon[] = Point{2,Float32}[]
            mask[] = zeros(size(points))
        end

        on(saveb.clicks) do _

            if !isnothing(tb.stored_string[])
                res.title = tb.stored_string[]
            end

            f = plot(permutedims([res])) # permutedims so that it's a matrix
            savedir = NMRInversions.save_file(res.title, filterlist = "png")

            if savedir == ""
                display("Please enter a name for your file on the file dialog.")
            else
                save(savedir, f)
            end

        end

    end ## BUTTON CLICKS

    gui
end

end


