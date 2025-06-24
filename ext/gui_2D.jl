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

    seq = res.seq
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

    xlbl = if seq == IRCPMG
        L"T_1 \, \textrm{(s)}"
    elseif seq == PFGCPMG
        L"D \, \textrm{(m^2/s)}"
    end
    ylbl = L"T_2 \,\textrm{(s)}"

    axmain_values = Axis(fig[3:10, 1:8],
                         xscale=log10, yscale=log10,
                         xlabel=xlbl, ylabel=ylbl,
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

    # Plot diagonal line
    if res.seq == IRCPMG
        diag_low = max(x[1],y[1]) 
        diag_high = min(x[end],y[end])  
        lines!(
            axmain_values, 
            [
                (diag_low, diag_low), 
                (diag_high,diag_high)
            ], 
            color=:black, linewidth=1
        )
    end
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


    #Create a matrix for all the discrete points in the space
    points = [[i, j] for i in x, j in y]
    mask = zeros(size(points))

    wa_indir, wa_dir, volumes = weighted_averages(res, silent = true)

    for (i, polygon) in enumerate(res.selections)

        # Selected points only
        mask .= [PolygonOps.inpolygon(p, polygon; in=1, on=1, out=0) for p in points]
        spo = mask .* z

        #draw current polygon
        lines!(axmain, Point2f.(polygon), linestyle=:dash, colormap=:tab10, colorrange=(1, 10), color=i, alpha=0.9)

        # Plot peak label
        xc = vec(sum(spo, dims=2)) ⋅ x / sum(spo)
        yc = vec(sum(spo, dims=1)) ⋅ y / sum(spo)

        lbl = if res.seq == IRCPMG
            " : T₁/T₂ = $(round(wa_indir[i]/wa_dir[i] , digits=1)) \n     Volume = $(round(volumes[i] *100 ,digits = 1))%"
        elseif res.seq == PFGCPMG
            " : D = $(round(wa_indir[i], sigdigits=3)), T₂ = $(round(wa_dir[i], sigdigits=3))\n     Volume = $(round(volumes[i] *100 ,digits = 1))%"

        end

        sc = scatter!(axmain, xc, yc, markersize=15,
            marker=collect('a':'z')[i],
            colormap=:tab10, colorrange=(1, 10),
            color=i,
            glowcolor=:white, glowwidth=4,
            label=lbl)

        text!(axmain, 2.5 * (1 / 100), (99 - 13 * (i - 1)) * (1 / 100),
              text=
              collect('a':'z')[i] * lbl,
              align=(:left, :top), colormap=:tab10, colorrange=(1, 10), color=i,
              #=glowcolor=:white, glowwidth=1,=#
              )

        # lines!(axtop, indir_dist[], linestyle=:dash, colormap=:tab10, colorrange=(1, 10), color=i, alpha=0.5)
        # lines!(axright, dir_dist[], 1:length(res.X_dir), linestyle=:dash, colormap=:tab10, colorrange=(1, 10), color=i, alpha=0.5)
        # axislegend( axmain,  framevisible=false)

        linkxaxes!(axmain, axtop)
        linkyaxes!(axmain, axright)
    end

end

function dynamic_plots(axmain, axtop, axright, selection, polygon)

    scatter!(axmain, selection, colormap=:tab10, colorrange=(1, 10), color=8)
    lines!(axmain, polygon, colormap=:tab10, colorrange=(1, 10), color=8)

    linkxaxes!(axmain, axtop)
    linkyaxes!(axmain, axright)
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
    plot!(gui[2:10, 1:9], res)

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
                save(savedir, f, px_per_unit=2)
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
