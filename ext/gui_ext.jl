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

    for (i, r) in enumerate(res)
        scatter!(ax, r.x, r.y, markersize=markersize, label= r.title )
        lines!(ax, r.xfit, r.yfit, label=getfield(r, normeq == true ? :eqn : :eq))
    end

    axislegend(ax, position=(res[1].seq == IR ? :rb : :rt), nbanks=2)

end


## Plots for 1D inversions

"""
    plot(res_mat::VecOrMat{inv_out_1D};kwargs...)

Plot the results contained in a vector of `inv_out_1D` structures.

The keyword (optional) arguments are:

- `selections` : Whether to highlight the selections in the plots (default is `false`).

"""
function Makie.plot(res_mat::AbstractVecOrMat{NMRInversions.inv_out_1D}; selections = false)

    fig = Figure(size=(700, 600))

    res = res_mat[1]
    # Make axes
    if res.seq in [NMRInversions.IR]
        ax1 = Axis(fig[1:5, 1:5], xlabel="time (s)", ylabel="Signal (a.u.)")
        ax2 = Axis(fig[1:5, 6:10], xlabel="T₁ (s)", xscale=log10)
        ax3 = Axis(fig[6:10,1:5], xlabel= "time (s)", ylabel="Residuals (a.u.)")

    elseif res.seq in [NMRInversions.CPMG]
        ax1 = Axis(fig[1:5, 1:5], xlabel="time (s)", ylabel="Signal (a.u.)")
        ax2 = Axis(fig[1:5, 6:10], xlabel="T₂ (s)", xscale=log10)
        ax3 = Axis(fig[6:10,1:5], xlabel= "time (s)", ylabel="Residuals (a.u.)")

    elseif res.seq in [NMRInversions.PFG]
        ax1 = Axis(fig[1:5, 1:5], xlabel="b factor (s/m² e-9)", ylabel="Signal (a.u.)")
        ax2 = Axis(fig[1:5, 6:10], xlabel="D (m²/s)", xscale=log10)
        ax3 = Axis(fig[6:10,1:5], xlabel= "b factor (s/m² e-9)", ylabel="Residuals (a.u.)")
    end

    linkxaxes!(ax1, ax3)

    for r in res_mat
        draw_on_axes(ax1, ax2, ax3, r, selections = selections)
    end

    return fig
end


"""
    plot(res::inv_out_1D)

Run the GUI to plot the 1D inversion results and select peaks you want to label.
"""
function Makie.plot(res::NMRInversions.inv_out_1D)

    fig = plot([res], selections = true)

    slider = IntervalSlider(fig[6,6:10], range = [1:length(res.X)...], startvalues = (1, length(res.X)))

    int_low= lift(slider.interval) do i
        return res.X[i[1]]
    end

    int_high= lift(slider.interval) do i
        return res.X[i[2]]
    end

    weighted_average = lift(slider.interval) do i
        return res.f[i[1]:i[2]]' * res.X[i[1]:i[2]] / sum(res.f[i[1]:i[2]]) 
    end

    labeltext = lift(slider.interval) do i
        return "Selected range from $(round(res.X[i[1]], sigdigits=2)) to $(round(res.X[i[2]], sigdigits=2))
        \n Weighted average = $(round(weighted_average[], sigdigits=2))"
    end

    Label(fig[7,6:10], labeltext, tellwidth=false)
    vlines!(fig.content[2], int_low, color=:red)
    vlines!(fig.content[2], int_high, color=:red)

    button_label = Button(fig[8,6:10], label = "Label current selection")
    button_filter = Button(fig[9,6:10], label = "Filter-out current selection")
    button_reset = Button(fig[10,6:8], label = "Reset selections")
    button_save = Button(fig[10,8:10], label = "Save and exit")

    on(button_label.clicks) do _
        push!(res.selections, slider.interval[])
        empty!(fig.content[1])
        empty!(fig.content[2])
        empty!(fig.content[3])
        draw_on_axes(fig.content[1], fig.content[2], fig.content[3], res, selections = true)
        vlines!(fig.content[2], int_low, color=:red)
        vlines!(fig.content[2], int_high, color=:red)
    end

    on(button_filter.clicks) do _
        empty!(fig.content[1])
        empty!(fig.content[2])
        empty!(fig.content[3])
        filter_range!(res, slider.interval[])
        draw_on_axes(fig.content[1], fig.content[2], fig.content[3], res, selections = true)
        vlines!(fig.content[2], int_low, color=:red)
        vlines!(fig.content[2], int_high, color=:red)
    end

    on(button_reset.clicks) do _
        empty!(res.selections)
        res.filter = ones(length(res.X))
        empty!(fig.content[1])
        empty!(fig.content[2])
        empty!(fig.content[3])
        draw_on_axes(fig.content[1], fig.content[2], fig.content[3], res)
        vlines!(fig.content[2], int_low, color=:red)
        vlines!(fig.content[2], int_high, color=:red)
    end

    on(button_save.clicks) do _
        savedir = NMRInversions.save_file(res.title, filterlist = "png")
        f = plot([res], selections = true)
        save(savedir, f)
    end


    return fig
end


function draw_on_axes(ax1, ax2, ax3, res::NMRInversions.inv_out_1D; selections = false)

    # apply filter to f and update fitted curve and residuals
    f_prime = res.filter .* res.f
    yfit_prime = create_kernel(
        res.seq, res.xfit, (res.seq == PFG ? res.X .* 1e9 : res.X)) * f_prime
    r_prime = create_kernel(
        res.seq, res.x, (res.seq == PFG ? res.X .* 1e9 : res.X)) * f_prime - res.y

    c = length(ax2.scene.plots)
    scatter!(ax1, res.x, real.(res.y), colormap=:tab10, colorrange=(1, 10), color=c)
    lines!(ax1, res.xfit, yfit_prime, colormap=:tab10, colorrange=(1, 10), color=c)
    lines!(ax2, res.X, f_prime, colormap=:tab10, colorrange=(1, 10), color=c)
    lines!(ax3, res.x, r_prime, colormap=:tab10, colorrange=(1, 10), color=c)
    scatter!(ax3, res.x, r_prime, colormap=:tab10, colorrange=(1, 10), color=c)


    if selections

        if res.seq in [NMRInversions.IR]
            Xlabel = ["<T₁> = ", " (s)"]
        elseif res.seq in [NMRInversions.CPMG]
            Xlabel = ["<T₂> = ", " (s)"]
        elseif res.seq in [NMRInversions.PFG]
            Xlabel = ["<D> = ", " (m²/s)"]
        end

        wa = weighted_averages(res)

        for (i,s) in enumerate(res.selections)
            clr = only(
                Makie.numbers_to_colors(
                    [i], 
                    Makie.to_colormap(:tab10), 
                    identity, 
                    Vec2(1, 10), 
                    Makie.automatic, 
                    Makie.automatic, 
                    RGBAf(0, 0, 0, 0)
                )
            )

            band!(
                ax2, 
                res.X[s[1]:s[2]],
                zeros(length(res.f[s[1]:s[2]])),
                (res.f .* res.filter)[s[1]:s[2]],
                color = clr, alpha = 0.5
            )

            height = (1 - (i-1) * 0.1) * maximum(res.f .* res.filter) 
            text!(
                ax2, res.X[2], height, 
                text = Xlabel[1]*"$(round(wa[i], sigdigits = 2))"*Xlabel[2], 
                align = (:left, :top), color = clr
            )

            #=vlines!(ax2, res.X[s[1]], color = clr, alpha = 0.8)=#
            #=vlines!(ax2, res.X[s[2]],color = clr, alpha = 0.8)=#
        end
    end
end


## Plots for 2D inversions

"""
    plot(results::AbstractVecOrMat{inv_out_2D}; kwargs...)

Plot the results contained in a matrix of `inv_out_2D` structures.
The arrangement is determined by the shape of the matrix.

The arguments are:

- `results` : The `inv_out_2D` matrix or vector containing the fit results.

Keyword (optional) arguments:
- `dims` : Dimensions of each plot (default: (400, 400)).
- `title` : Title of the plot (default is the title in the inv_out_2D).
- `colormap` : Color map of the plot (default: :viridis).
- `contf` : Whether to use a filled contour plot (default: false).
- `levels` : Number of contour levels (default: 40).
- `labelsizes` : Size of the axis labels (default: (23, 23)).
- `ticksizes` : Size of the axis ticks (default: (14, 14)).
- `titlesize` : Size of the title (default: 17).
- `titlefont` : Font of the title (default: :bold).
"""
function Makie.plot(
    res_mat::AbstractVecOrMat{NMRInversions.inv_out_2D}; dims = (400,400),
    kwargs...)

    if !isa(res_mat, AbstractMatrix)
        res_mat = permutedims(res_mat)
    end

    f = Figure(size=dims .* reverse(size(res_mat)))

    for (ind, res) in pairs(res_mat)
        plot!(f[Tuple(ind)...], res;  kwargs...)
    end

    return f

end


"""
    plot!(fig, res::NMRInversions.inv_out_2D; kwargs...)

Plot the results contained in a `inv_out_2D` structure on a figure or a grid position object.

The arguments are:

- `fig` : The figure or grid position object.
- `res` : The `inv_out_2D` structure containing the fit results.

Keyword (optional) arguments:
- `title` : Title of the plot (default: "").
- `colormap` : Color map of the plot (default: :viridis).
- `contf` : Whether to use a filled contour plot (default: false).
- `levels` : Number of contour levels (default: 40).
- `labelsizes` : Size of the axis labels (default: (23, 23)).
- `ticksizes` : Size of the axis ticks (default: (14, 14)).
- `titlesize` : Size of the title (default: 17).
- `titlefont` : Font of the title (default: :bold).
"""
function Makie.plot!(fig::Union{Makie.Figure,Makie.GridPosition}, res::NMRInversions.inv_out_2D;
                     title=res.title, colormap=:viridis, contf=false, levels = 40, 
                     labelsizes = (23, 23), ticksizes = (15, 15), titlesize = 17, titlefont = :bold 
                     )

    x = res.X_indir
    y = res.X_dir

    xplot = collect(range(0, 1, length(x)))
    yplot = collect(range(0, 1, length(y)))

    # Make axes

    axmain = Axis(fig[3:10, 1:8])
    axtop = Axis(fig[1:2, 1:8], limits=((xplot[1], xplot[end]), (0, nothing)), 
                 title=title, titlesize=titlesize, titlefont=titlefont)
    axright = Axis(fig[3:10, 9:10], limits=((0, nothing), (yplot[1], yplot[end])))

    hidedecorations!(axtop)
    hidedecorations!(axright)
    hidedecorations!(axmain)

    axmain_values = Axis(fig[3:10, 1:8],
                         xscale=log10, yscale=log10,
                         xlabel=L"T_1 \, \textrm{(s)}", ylabel=L"T_2 \,\textrm{(s)}",
                         xlabelsize=labelsizes[1], ylabelsize=labelsizes[2],
                         limits=(x[1], x[end], y[1], y[end]),
                         xticklabelsize=ticksizes[1], yticklabelsize=ticksizes[2],
                         xtickalign = 1.0, ytickalign = 1.0
                         )

    Makie.deactivate_interaction!(axmain_values, :rectanglezoom) # Disable zoom
    Makie.deactivate_interaction!(axmain, :rectanglezoom) # Disable zoom
    Makie.deactivate_interaction!(axright, :rectanglezoom) # Disable zoom
    Makie.deactivate_interaction!(axtop, :rectanglezoom) # Disable zoom

    linkxaxes!(axmain, axtop)
    linkyaxes!(axmain, axright)

    draw_on_axes(axmain, axtop, axright, res, colormap,contf,levels)

end

function draw_on_axes(axmain, axtop, axright, res, clmap,contf,levels)

    empty!(axmain)
    empty!(axtop)
    empty!(axright)

    z = res.F' .* res.filter'
    x = range(0, 1, size(z, 1))
    y = range(0, 1, size(z, 2))

    # Plots
    if contf == true
        contourf!(axmain, x, y, z, colormap=clmap, levels=levels)
    else
        contour!(axmain, x, y, z, colormap=clmap, levels=levels)
    end
    lines!(axtop, x, vec(sum(z, dims=2)), colormap=clmap, colorrange=(1, 10), color=5, alpha=0.8)
    lines!(axright, vec(sum(z, dims=1)), y, colormap=clmap, colorrange=(1, 10), color=5, alpha=0.8)

    #=band!(axtop, x, zeros(length(x)), vec(sum(z, dims=2)) ,colormap=clmap, colorrange=(1, 10), color = 5 ,alpha=0.5)=#
    #=band!(axright, Point2f.(zeros(length(y)), y), Point2f.(vec(sum(z, dims=1)), y) ,colormap=clmap, colorrange=(1, 10), alpha=0.5)=#

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
    
    # Title textbox
    tb = Textbox(gui[1, 2:7], 
                 width=300, 
                 reset_on_defocus=true,
                 stored_string= res.title == "" ? " " : res.title
                 )
    
    # Buttons
    labelb = Button(gui[1, 10:15]; label="Label current selection")
    clearb = Button(gui[2, 10:15]; label="Cancel current selection")
    resetb = Button(gui[1, 15:19]; label="Reset everything")
    filterb = Button(gui[2, 15:19]; label="Filter-out unselected")
    saveb = Button(gui[3, 10:19]; label="Save and exit")

    Label(gui[5, 10:13], "Colormap:", halign=:right)
    colormenu = Menu(gui[5, 14:18], 
                     options = ["viridis","tempo","heat","curl","balance","solar"],default="viridis")


    Label(gui[6, 10:13], "Levels:", halign=:right)
    levelsbox = Textbox(gui[6, 14:16], width=100, stored_string="40", validator=Int, reset_on_defocus=true)
    levels = Observable(40)

    Label(gui[7, 10:13], "Fill contour:", halign=:right)
    fillcheck = Checkbox(gui[7, 14], checked=false)


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

                draw_on_axes(axmain, axtop, axright, res, colormenu.selection[], fillcheck.checked[],levels[])
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

            draw_on_axes(axmain, axtop, axright, res, colormenu.selection[], fillcheck.checked[],levels[])
            dynamic_plots(axmain, axtop, axright, selection, polygon)

            selection[] = Point{2,Float32}[]
            polygon[] = Point{2,Float32}[]
            mask[] = zeros(size(points))

        end


        on(resetb.clicks) do _
            res.filter .= ones(size(res.F))
            empty!(res.selections)

            draw_on_axes(axmain, axtop, axright, res, colormenu.selection[], fillcheck.checked[], levels[])
            dynamic_plots(axmain, axtop, axright, selection, polygon)

            selection[] = Point{2,Float32}[]
            polygon[] = Point{2,Float32}[]
            mask[] = zeros(size(points))
        end

        on(saveb.clicks) do _

            # permutedims so that res is a matrix
            f = plot(permutedims([res]), levels = levels[], colormap = colormenu.selection[], contf = fillcheck.checked[]) 
            savedir = NMRInversions.save_file(res.title, filterlist = "png")

            if savedir == ""
                display("Please enter a name for your file on the file dialog.")
            else
                save(savedir, f)
            end
        end

        on(tb.stored_string) do str

            if !isnothing(str)
                res.title = str
            end
        end

        on(colormenu.selection) do _
            draw_on_axes(axmain, axtop, axright, res, colormenu.selection[], fillcheck.checked[], levels[])
            dynamic_plots(axmain, axtop, axright, selection, polygon)
        end

        on(fillcheck.checked) do _
            draw_on_axes(axmain, axtop, axright, res, colormenu.selection[], fillcheck.checked[], levels[])
            dynamic_plots(axmain, axtop, axright, selection, polygon)
        end

        on(levelsbox.stored_string) do s
            levels[] = parse(Int, s)
            draw_on_axes(axmain, axtop, axright, res, colormenu.selection[], fillcheck.checked[], levels[])
            dynamic_plots(axmain, axtop, axright, selection, polygon)
        end

    end ## BUTTON CLICKS
    
    
    gui
    

end

end


