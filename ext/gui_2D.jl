## Plots for 2D inversions

using GLMakie: coordinates

"""
    plot(results::AbstractVecOrMat{inv_out_2D}; kwargs...)

Plot the results contained in a matrix of `inv_out_2D` structures.
The arrangement is determined by the shape of the matrix.

The arguments are:

- `results` : The `inv_out_2D` matrix or vector containing the fit results.

Keyword (optional) arguments:
- `dims` : Dimensions of each plot (default: `(400, 400)`).
- `title` : Title of the plot (default: `""`).
- `colormap` : Color map of the plot (default: `:viridis`).
- `contf` : Whether to use a filled contour plot (default : `false`).
- `levels` : Number of contour levels (default: `40`).
- `labelsizes` : Size of the axis labels (default: `(23, 23)`).
- `ticksizes` : Size of the axis ticks (default: `(14, 14)`).
- `titlesize` : Size of the title (default: `17`).
- `titlefont` : Font of the title (default: `:bold`).
- `legendlabelsize` : Size of the legend label (default: `12`)
- `gap` : The gap between the contour plot and the distribution plots (default: `0`)
- `extend_grid` : Whether to extend the gridlines to the plots on the sides (default: `false`)
- `legendposition` : Where to put the legend (default: `:lt` )
- `bandplots` : Whether to use band plots to highlight the selections in the top and right distributions (default: `true`)
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
- `title` : Title of the plot (default: `""`).
- `colormap` : Color map of the plot (default: `:viridis`).
- `contf` : Whether to use a filled contour plot (default : `false`).
- `levels` : Number of contour levels (default: `40`).
- `labelsizes` : Size of the axis labels (default: `(23, 23)`).
- `ticksizes` : Size of the axis ticks (default: `(14, 14)`).
- `titlesize` : Size of the title (default: `17`).
- `titlefont` : Font of the title (default: `:bold`).
- `legendlabelsize` : Size of the legend label (default: `12`)
- `gap` : The gap between the contour plot and the distribution plots (default: `0`)
- `extend_grid` : Whether to extend the gridlines to the plots on the sides (default: `false`)
- `legendposition` : Where to put the legend (default: `:lt` )
- `bandplots` : Whether to use band plots to highlight the selections in the top and right distributions (default: `true`)
"""
function Makie.plot!(fig::Union{Makie.Figure,Makie.GridPosition}, res::NMRInversions.inv_out_2D;
                     title=res.title, colormap=:viridis, contf=false, levels = 40, 
                     labelsizes = (23, 23), ticksizes = (15, 15), titlesize = 17, titlefont = :bold ,
                     legendlabelsize = 12, gap = 0, extend_grid = false, legendposition = :lt,
                     bandplots = true,
                     )

    seq = res.seq
    x = res.X_indir
    y = res.X_dir

    xlbl = if seq == IRCPMG
        L"T_1 \, \textrm{(s)}"
    elseif seq == PFGCPMG
        L"D \, \textrm{(m^2/s)}"
    end

    ylbl = if seq in [IRCPMG, PFGCPMG]
        L"T_2 \,\textrm{(s)}"
    end

    x_low = seq == IRCPMG ? exp10(floor(log10(min(x[1],y[1])))) : exp10(floor(log10(x[1])))
    x_high = seq == IRCPMG ? exp10(ceil(log10(max(x[end],y[end])))) : exp10(ceil(log10(x[end])))
    y_low = seq == IRCPMG ? exp10(floor(log10(min(x[1],y[1])))) : exp10(floor(log10(y[1])))
    y_high = seq == IRCPMG ? exp10(ceil(log10(max(x[end],y[end])))) : exp10(ceil(log10(y[end])))

    gr = fig[1:10, 1:10] = GridLayout()

    # order is important, dont change
    axmain = Axis(
        gr[3:10, 1:8],
        xlabel=xlbl, ylabel=ylbl,
        xlabelsize=labelsizes[1], ylabelsize=labelsizes[2],
        xticklabelsize=ticksizes[1], yticklabelsize=ticksizes[2],
        xtickalign = 0, ytickalign = 0,
        xscale = log10, yscale = log10,
        limits = (x_low, x_high, y_low, y_high)
    )

    axtop = Axis(
        gr[1:2, 1:8], xscale = log10,
        title=title, titlesize=titlesize, titlefont=titlefont,
        ygridvisible=false, 
        limits = (x_low, x_high ,0 , nothing) 
    )
    axright = Axis(
        gr[3:10, 9:10], yscale = log10,
        xgridvisible=false,
        limits = (
            0, nothing, y_low , y_high,
        )
    )

    if abs(log10(x_low)) + abs(log10(x_high)) <= 7
        axmain.xticks = LogTicks(floor(log10(x_low)):ceil(log10(x_high)))
        axtop.xticks = LogTicks(floor(log10(x_low)):ceil(log10(x_high)))
    else
        axmain.xticks = LogTicks(floor(log10(x_low)):2:ceil(log10(x_high)))
        axtop.xticks = LogTicks(floor(log10(x_low)):2:ceil(log10(x_high)))
    end
    if abs(log10(y_low)) + abs(log10(y_high)) <= 7
        axmain.yticks = LogTicks(floor(log10(y_low)):ceil(log10(y_high)))
        axright.yticks = LogTicks(floor(log10(y_low)):ceil(log10(y_high)))
    else
        axmain.yticks = LogTicks(floor(log10(y_low)):2:ceil(log10(y_high)))
        axright.yticks = LogTicks(floor(log10(y_low)):2:ceil(log10(y_high)))
    end

    linkxaxes!(axmain, axtop)
    linkyaxes!(axright, axmain)

    hidedecorations!(axtop, grid= !extend_grid)
    hidedecorations!(axright, grid= !extend_grid)

    colgap!(gr , gap)
    rowgap!(gr , gap)
    if gap <= 3
        hidespines!(axtop, :b)
        hidespines!(axright, :l)
    end

    Makie.deactivate_interaction!(axmain, :rectanglezoom) # Disable zoom
    Makie.deactivate_interaction!(axright, :rectanglezoom) # Disable zoom
    Makie.deactivate_interaction!(axtop, :rectanglezoom) # Disable zoom

    draw_on_axes(axmain, axtop, axright, res, colormap, contf, levels, legendlabelsize, legendposition, bandplots)
end


function draw_on_axes(
    axmain, axtop, axright, res, 
    clmap, contf, levels, legendlabelsize, 
    legendposition, bandplots)

    empty!(axmain)
    empty!(axtop)
    empty!(axright)

    z = res.F' .* res.filter'
    x = res.X_indir
    y = res.X_dir

    # Plots
    if contf == true
        contourf!(axmain, x, y, z, colormap=clmap, levels=levels)
    else
        contour!(axmain, x, y, z, colormap=clmap, levels=levels)
        
    end

    linecolor = :black

    lines!(axtop, x, vec(sum(z, dims=2)), color = linecolor, alpha=0.8, linewidth =1)
    lines!(axright, vec(sum(z, dims=1)), y, color = linecolor, alpha=0.8, linewidth =1)

    #Create a matrix for all the discrete points in the space
    points = [[i, j] for i in x, j in y]
    mask = zeros(size(points))

    wa_indir, wa_dir, volumes = weighted_averages(res, silent = true)

    markers = collect('a':'z')
    colors = cgrad(:tab10)

    for (i, polygon) in enumerate(res.selections)

        # Selected points only
        mask .= [PolygonOps.inpolygon(p, polygon; in=1, on=1, out=0) for p in points]
        spo = mask .* z

        xdist = vec(sum(spo, dims=2))
        ydist = vec(sum(spo, dims=1))

        # Plot peak label
        xc = xdist ⋅ x / sum(spo)
        yc = ydist ⋅ y / sum(spo)

        lbl = if res.seq == IRCPMG
            "T₁/T₂ = $(round(wa_indir[i]/wa_dir[i] , digits=1)) \nVolume = $(round(volumes[i] *100 ,digits = 1))%"
        elseif res.seq == PFGCPMG
            "D = $(round(wa_indir[i], sigdigits=3)), T₂ = $(round(wa_dir[i], sigdigits=3))\nVolume = $(round(volumes[i] *100 ,digits = 1))%"

        end

        alphas = 0.83

        if bandplots == true
            band!(axtop, x, zeros(length(x)), xdist, color = colors[i], alpha = alphas-0.2)
            band!(axright, y, zeros(length(y)), ydist, color = colors[i], alpha = alphas-0.2,
                  direction = :y)
        end 

        lines!(axmain, Point2f.(polygon), linestyle=:dash, color = colors[i], alpha=alphas)

        scatter!(
            axmain, xc, yc, markersize=15,
            marker=markers[i], color = colors[i],
            label=Makie.rich(lbl, color = colors[i]),
            glowcolor=:white, glowwidth=4,
        )
    end

    if res.seq == IRCPMG
        plot_diagonal(axmain,x,y)
    end

    if !isempty(res.selections)
        axislegend(axmain, position=legendposition,
                   framevisible = false,
                   labelsize = legendlabelsize)
    end
end

function plot_diagonal(ax,x,y)

    low = exp10(floor(log10(min(x[1],y[1]))))
    high = exp10(ceil(log10(max(x[end],y[end]))))
    lines!(
        ax, 
        [ (low, low), (high,high) ], 
        color=:black, linewidth=1
    )
end

function dynamic_plots(res, axmain, axtop, axright, selection, poly, xcoord, ycoord)

    color = cgrad(:tab10)[length(res.selections)+1]
    scatter!(axmain, selection, color= color)
    lines!(axmain, poly, color= color)

    vlines!(axmain, xcoord, color = :black, linewidth = 0.5, alpha = 0.5)
    vlines!(axtop, xcoord, color = :black, linewidth = 0.5, alpha = 0.5)
    hlines!(axmain, ycoord, color = :black, linewidth = 0.5, alpha = 0.5)
    hlines!(axright, ycoord, color = :black, linewidth = 0.5, alpha = 0.5)

end

#=function draw_3d(ax3d, res, sp)=#
#==#
#=    surface!(ax3d, res.F' .* res.filter', colormap=:tempo)=#
#==#
#=    surface!(ax3d, spo, colormap=:blues,=#
#=        colorrange=(minimum(filter(x -> x != 0, vec(z))), maximum(z)),=#
#=        lowclip=:transparent)=#
#==#
#=end=#


"""
    Makie.plot(res::inv_out_2D)
Run the GUI to plot the results and select peaks you want to label.
"""
function Makie.plot(res::NMRInversions.inv_out_2D)

    GLMakie.activate!(; title= "2D inversion GUI")

    gui = Figure(size=(900, 500))

    plot!(gui[2:10, 1:9], res, title = "")

    axmain = gui.content[1]
    axtop = gui.content[2]
    axright = gui.content[3]
    
    # Title textbox
    tb = Textbox(gui[1, 1:7], 
                 width=300, 
                 reset_on_defocus=true,
                 stored_string= res.title == "" ? " " : res.title
                 )
    
    # Buttons
    labelb = Button(gui[1, 10:15]; label="Label current selection")
    filterb = Button(gui[2, 10:15]; label="Filter-out selection")
    clearb = Button(gui[3, 10:15]; label="Cancel current selection")

    reset_sel_b = Button(gui[1, 15:19]; label="Reset selections")
    reset_filter_b = Button(gui[2, 15:19]; label="Reset filter")
    saveb = Button(gui[3, 15:19]; label="Save and exit")

    Label(gui[5, 10:13], "Colormap:", halign=:right)
    colormenu = Menu(gui[5, 14:18], 
                     options = ["viridis","tempo","heat","curl","balance","solar"],default="viridis")


    Label(gui[6, 10:13], "Levels:", halign=:right)
    levelsbox = Textbox(gui[6, 14:16], width=100, stored_string="40", validator=Int, reset_on_defocus=true)
    levels = Observable(40)

    Label(gui[7, 10:13], "Fill contour:", halign=:right)
    fillcheck = Checkbox(gui[7, 14], checked=false)

    Label(gui[8, 10:13], "Band plots:", halign=:right)
    bandcheck = Checkbox(gui[8, 14], checked=true)

    Label(gui[9,10:19],"α = $(round(res.alpha, sigdigits=3)) , SNR = $(round(res.SNR, digits=1))")

    xcoord = Observable(0.0)
    ycoord = Observable(0.0)
    coord_label = Observable("")

    xlbl = if res.seq == IRCPMG
        "T₁"
    elseif res.seq == PFGCPMG
        "D"
    end

    ylbl = if res.seq in [IRCPMG, PFGCPMG]
        "T₂"
    end

    on(events(gui).mouseposition) do _
        if is_mouseinside(axmain)
            xcoord[] = (exp10.(mouseposition(axmain)[1]))
            ycoord[] = (exp10.(mouseposition(axmain)[2]))
            coord_label[] = "$xlbl = $(round(xcoord[],sigdigits=2)) , $ylbl = $(round(ycoord[],sigdigits=2))"
        elseif is_mouseinside(axtop)
            xcoord[] = (exp10.(mouseposition(axtop)[1]))
            ycoord[] = 0.0
            coord_label[] = "$xlbl = $(round(xcoord[],sigdigits=2))"
        elseif is_mouseinside(axright)
            xcoord[] = 0.0
            ycoord[] = (exp10.(mouseposition(axright)[2]))
            coord_label[] = "$ylbl = $(round(ycoord[],sigdigits=2))"
        else
            xcoord[] = 0.0
            ycoord[] = 0.0
            coord_label[] = "Hover mouse over plot to view coordinates."
        end
    end

    Label(gui[10,10:19], coord_label)

    z = res.F' .* res.filter'
    x = res.X_indir
    y = res.X_dir

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
    visual_polygon = Observable(Point{2,Float32}[])

    # Dynamic plots
    dynamic_plots(res,axmain, axtop, axright, selection, visual_polygon, xcoord, ycoord)

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

        on(events(gui).mouseposition) do _
            if is_mouseinside(axmain)
                if size(selection[], 1) > 2
                    visual_polygon[] = vcat(selection[], [Point{2,Float32}(xcoord[],ycoord[])], [selection[][1]])
                end
            end
            if !is_mouseinside(axmain)
                if size(selection[], 1) > 2
                    visual_polygon[] =  polygon[] 
                end
            end
        end

    end ## SELECTING POINTS


    begin ## BUTTON CLICKS

        on(labelb.clicks) do _

            if size(selection[], 1) < 3
                @warn("You need to make a selection.")
            else

                delete!.(gui.content[findall(x -> x isa Legend, gui.content)])
                push!(res.selections, polygon[])

                draw_on_axes(axmain, axtop, axright, res, Symbol(colormenu.selection[]), fillcheck.checked[],levels[], 12, :lt, bandcheck.checked[])
                dynamic_plots(res, axmain, axtop, axright, selection, visual_polygon, xcoord, ycoord)

                # Clear for next selection
                selection[] = Point{2,Float32}[]
                polygon[] = Point{2,Float32}[]
                visual_polygon[] = Point{2,Float32}[]
                mask[] = zeros(size(points))

                reset_limits!(axmain)
            end

        end

        on(clearb.clicks) do _

            selection[] = Point{2,Float32}[]
            polygon[] = Point{2,Float32}[]
            visual_polygon[] = Point{2,Float32}[]
            mask[] = zeros(size(points))

        end

        on(filterb.clicks) do _

            if size(selection[], 1) < 3
                @warn("You need to make a selection.")
            else

                res.filter .= res.filter .* [PolygonOps.inpolygon(p, polygon[]; in=0, on=0, out=1) for p in points]'
            end

            delete!.(gui.content[findall(x -> x isa Legend, gui.content)])
            draw_on_axes(axmain, axtop, axright, res, colormenu.selection[], fillcheck.checked[],levels[],12,:lt, bandcheck.checked[])
            dynamic_plots(res, axmain, axtop, axright, selection, visual_polygon, xcoord, ycoord)

            selection[] = Point{2,Float32}[]
            polygon[] = Point{2,Float32}[]
            visual_polygon[] = Point{2,Float32}[]
            mask[] = zeros(size(points))

        end

        on(reset_sel_b.clicks) do _
            empty!(res.selections)

            delete!.(gui.content[findall(x -> x isa Legend, gui.content)])
            draw_on_axes(axmain, axtop, axright, res, colormenu.selection[], fillcheck.checked[], levels[],12,:lt, bandcheck.checked[])
            dynamic_plots(res, axmain, axtop, axright, selection, visual_polygon, xcoord, ycoord)

            selection[] = Point{2,Float32}[]
            polygon[] = Point{2,Float32}[]
            visual_polygon[] = Point{2,Float32}[]
            mask[] = zeros(size(points))
        end

        on(reset_filter_b.clicks) do _
            res.filter .= ones(size(res.F))

            delete!.(gui.content[findall(x -> x isa Legend, gui.content)])
            draw_on_axes(axmain, axtop, axright, res, colormenu.selection[], fillcheck.checked[], levels[],12,:lt, bandcheck.checked[])
            dynamic_plots(res, axmain, axtop, axright, selection, visual_polygon, xcoord, ycoord)

            selection[] = Point{2,Float32}[]
            polygon[] = Point{2,Float32}[]
            visual_polygon[] = Point{2,Float32}[]
            mask[] = zeros(size(points))
        end

        on(saveb.clicks) do _

            # permutedims so that res is a matrix
            f = plot(permutedims([res]), levels = levels[], colormap = Symbol(colormenu.selection[]), contf = fillcheck.checked[], bandplots=bandcheck) 
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
            delete!.(gui.content[findall(x -> x isa Legend, gui.content)])
            draw_on_axes(axmain, axtop, axright, res, colormenu.selection[], fillcheck.checked[], levels[],12,:lt, bandcheck.checked[])
            dynamic_plots(res, axmain, axtop, axright, selection, visual_polygon, xcoord, ycoord)
        end

        on(fillcheck.checked)  do _
            delete!.(gui.content[findall(x -> x isa Legend, gui.content)])
            draw_on_axes(axmain, axtop, axright, res, colormenu.selection[], fillcheck.checked[], levels[],12,:lt, bandcheck.checked[])
            dynamic_plots(res, axmain, axtop, axright, selection, visual_polygon, xcoord, ycoord)
        end

        on(bandcheck.checked)  do _
            delete!.(gui.content[findall(x -> x isa Legend, gui.content)])
            draw_on_axes(axmain, axtop, axright, res, colormenu.selection[], fillcheck.checked[], levels[],12,:lt, bandcheck.checked[])
            dynamic_plots(res, axmain, axtop, axright, selection, visual_polygon, xcoord, ycoord)
        end

        on(levelsbox.stored_string) do s
            levels[] = parse(Int, s)
            draw_on_axes(axmain, axtop, axright, res, colormenu.selection[], fillcheck.checked[], levels[],12,:lt, bandcheck.checked[])
            dynamic_plots(res, axmain, axtop, axright, selection, visual_polygon, xcoord, ycoord)
        end

    end ## BUTTON CLICKS


    gui
    
end


## override function to prevent logscale bug
# https://discourse.julialang.org/t/glmakie-selecting-points-in-log-log-plot-sync-dissimilar-axes/130371/9
function select_point(scene; blocking = false, priority=2, space = :data, kwargs...)
    key = Mouse.left
    waspressed = Observable(false)
    point = Observable([Point2f(0,0)])
    point_ret = Observable(Point2f(0,0))
    # Create an initially hidden  arrow
    plotted_point = scatter!(
        scene, point; space = :data, transformation = Makie.Transformation(scene; transform_func = identity), visible = false, marker = Circle, markersize = 20px,
        color = RGBAf(0.1, 0.1, 0.8, 0.5), kwargs...,
    )

    onany(events(scene).mousebutton, Makie.transform_func_obs(scene), priority=priority) do event, transform_func
        if event.button == key && is_mouseinside(scene)
            mp = mouseposition(scene)
            if event.action == Mouse.press
                waspressed[] = true
                plotted_point[:visible] = true  # start displaying
                point[][1] = mp
                point[] = point[]
                return Consume(blocking)
            end
        end
        if !(event.button == key && event.action == Mouse.press)
            if waspressed[] # User has selected the rectangle
                waspressed[] = false
                final_point = Makie.project(scene, :data, space, point[][1])
                if space == :data
                    final_point = Makie.apply_transform(Makie.inverse_transform(transform_func), final_point)
                end
                point_ret[] = final_point
            end
            plotted_point[:visible] = false
            return Consume(blocking)
        end
        return Consume(false)
    end
    on(events(scene).mouseposition, priority=priority) do event
        if waspressed[]
            mp = mouseposition(scene)
            point[][1] = mp
            point[] = point[] # actually update observable
            return Consume(blocking)
        end
        return Consume(false)
    end

    return point_ret
end

