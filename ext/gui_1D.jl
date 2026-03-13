## Plots for 1D inversions

"""
    plot(res_mat::VecOrMat{inv_out_1D};kwargs...)

Plot the results contained in a vector of `inv_out_1D` structures.

The keyword (optional) arguments are:

- `selections` : Whether to highlight the selections in each of the plots (default is `true`).
- `stacked` : Whether to plot different curves as a stack, instead of superimposed. (default is `true`).

"""
function Makie.plot(res_mat::AbstractVecOrMat{NMRInversions.inv_out_1D}; 
                    selections::Bool = true,
                    stacked::Bool = true, 
                    yscale = identity, 
                    xscale = identity,
                    )

    fig = Figure(size=(700, 600))

    res = res_mat[1]

    position1 = fig[1:5, 1:5]
    position2 = (stacked && length(res_mat) > 1 ) ? fig[1:10, 6:10] : fig[1:5, 6:10]
    position3 = fig[6:10,1:5]

    # Make axes
    if res.seq in [NMRInversions.IR, NMRInversions.SR]
        ax1 = Axis(position1, xlabel="time (s)", ylabel="Signal (a.u.)", yscale = yscale, xscale = xscale)
        ax2 = Axis(position2, xlabel=L"T_1 \, \textrm{(s)}", xscale=log10)
        ax3 = Axis(position3, xlabel= "time (s)", ylabel="Residuals (a.u.)", xscale = xscale)

    elseif res.seq in [NMRInversions.CPMG]
        ax1 = Axis(position1, xlabel="time (s)", ylabel="Signal (a.u.)")
        ax2 = Axis(position2, xlabel=L"T_2 \, \textrm{(s)}", xscale=log10)
        ax3 = Axis(position3, xlabel= "time (s)", ylabel="Residuals (a.u.)")

    elseif res.seq in [NMRInversions.PFG]
        ax1 = Axis(position1, xlabel="b factor (s/m² e-9)", ylabel="Signal (a.u.)")
        ax2 = Axis(position2, xlabel=L"D \, \textrm{(m^2/s)}" , xscale=log10)
        ax3 = Axis(position3, xlabel= "b factor (s/m² e-9)", ylabel="Residuals (a.u.)")
    end

    ax2.xlabelsize = 18

    linkxaxes!(ax1, ax3)

    for (i,r) in enumerate(res_mat)
        draw_on_axes(ax1, ax2, ax3, r, 
                     selections = selections, 
                     offset = stacked ? i * 1.1 : 0,
                     n_plots = length(res_mat),
                     i = i,
                     )
    end

    ax1.yscale = yscale
    ax1.xscale = xscale

    if length(res_mat) > 1
        if stacked
        else
            Legend(fig[6:10,6:10], ax2)
        end
    end

    ax2.limits = (
        minimum(minimum(getfield.(res_mat,:X))),
        maximum(maximum(getfield.(res_mat,:X))),
        stacked ? nothing : minimum(map(res -> -0.025 * maximum(res.f .* res.filter), res_mat)),
        stacked ? nothing : maximum(map(res ->  1.025 * maximum(res.f .* res.filter), res_mat)),
    )


    if stacked 
        for a in (ax1, ax2, ax3)
            a.yticksvisible = false
            a.yticklabelsvisible = false
            a.yticks = 1.1:1.1:length(res_mat) * 1.1
        end
    end

    return fig
end


function redraw(fig::Figure,res::inv_out_1D, int_low, int_high)
    empty!(fig.content[1])
    empty!(fig.content[2])
    empty!(fig.content[3])
    draw_on_axes(fig.content[1], fig.content[2], fig.content[3], res, selections = true)
    vlines!(fig.content[2], int_low, color=:red)
    vlines!(fig.content[2], int_high, color=:red)
end

"""
    plot(res::inv_out_1D)

Run the GUI to plot the 1D inversion results and select peaks you want to label.
"""
function Makie.plot(res::NMRInversions.inv_out_1D)

    GLMakie.activate!(;title= "1D inversion GUI")
    fig = plot([res], selections = true, stacked= false)

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
        \n Weighted average = $(round(weighted_average[], sigdigits=3))"
    end

    Label(fig[7,6:10], labeltext, tellwidth=false)
    vlines!(fig.content[2], int_low, color=:red)
    vlines!(fig.content[2], int_high, color=:red)

    button_label = Button(fig[8,6:8], label = "Label selection")
    button_filter = Button(fig[8,8:10], label = "Filter-out selection")
    button_reset_s = Button(fig[9,6:8], label = "Reset selections")
    button_reset_f = Button(fig[9,8:10], label = "Reset filter")
    button_save = Button(fig[10,8:10], label = "Save and exit")

    log_y = false
    log_x = false
    if !any(real.(res.y) .<= 0)
        button_yscale = Button(fig[10,6:8], label = "Change y scale")

        on(button_yscale.clicks) do _
            if log_y
                fig.content[1].yscale = identity
                log_y = false
            else
                fig.content[1].yscale = log10
                log_y = true
            end

        end
    else 
        button_xscale = Button(fig[10,6:8], label = "Change x scale")

        on(button_xscale.clicks) do _
            if log_x
                fig.content[1].limits = (exp10(floor(log10(res.x[1]))), nothing, nothing, nothing) #this prevents an error
                fig.content[1].xscale = identity
                fig.content[3].xscale = identity
                fig.content[1].limits = (nothing, nothing, nothing, nothing) #undo the above
                log_x = false
            else
                fig.content[1].xscale = log10
                fig.content[3].xscale = log10
                log_x = true
            end

        end
    end

    on(button_label.clicks) do _
        push!(res.selections, slider.interval[])
        redraw(fig, res, int_low, int_high)
    end

    on(button_filter.clicks) do _
        filter_range!(res, slider.interval[])
        redraw(fig, res, int_low, int_high)
    end

    on(button_reset_s.clicks) do _
        empty!(res.selections)
        redraw(fig, res, int_low, int_high)
    end

    on(button_reset_f.clicks) do _
        res.filter = ones(length(res.X))
        redraw(fig, res, int_low, int_high)
    end

    on(button_save.clicks) do _
        savedir = NMRInversions.save_file(res.title, filterlist = "png")
        f = plot([res], selections = true, 
                 yscale = log_y ? log10 : identity,
                 xscale = log_x ? log10 : identity
                 )
        save(savedir, f, px_per_unit = 2.5)
    end


    return fig
end


"""
- `offset` normalises the data to 1 and plots them as `data .+ offset`, used to stack many graphs.
"""
function draw_on_axes(ax1, ax2, ax3, res::NMRInversions.inv_out_1D; 
                      selections = false, offset::Real = 0, n_plots::Int =1, i = 1)

    # apply filter to f and update fitted curve and residuals
    f_prime = res.filter .* res.f

    yfit_prime = create_kernel(
        res.seq, res.xfit, (res.seq == PFG ? res.X .* 1e9 : res.X), y = res.y
    ) * f_prime

    r_prime = create_kernel(
        res.seq, res.x, (res.seq == PFG ? res.X .* 1e9 : res.X), y = res.y
    ) * f_prime - real.(res.y)

    scatter!(ax1, res.x, 
             offset == 0 ? real.(res.y) : (scale_to_one!(real.(res.y)) .+ offset), 
             colormap=:tab10, colorrange=(1, 10), color=i)
    lines!(ax1, res.xfit, 
           offset == 0 ? real.(yfit_prime) : (scale_to_one!(real.(yfit_prime)) .+ offset), 
           colormap=:tab10, colorrange=(1, 10), color=i)
    !selections && lines!(ax2, res.X, 
                          offset == 0 ? real.(f_prime) : (scale_to_one!(real.(f_prime)) .+ offset), 
                          colormap=:tab10, colorrange=(1, 10), color=i, label = res.title)
    lines!(ax3, res.x, 
           offset == 0 ? real.(r_prime) : (scale_to_one!(real.(r_prime)) .* 0.5 .+ offset ), 
           colormap=:tab10, colorrange=(1, 10), color=i)
    scatter!(ax3, res.x, 
             offset == 0 ? real.(r_prime) : (scale_to_one!(real.(r_prime)) .* 0.5 .+ offset ), 
             colormap=:tab10, colorrange=(1, 10), color=i)


    if selections

        if res.seq in [NMRInversions.IR]
            Xlabel = ["<T₁> = ", " s, Area= "," %"]
        elseif res.seq in [NMRInversions.CPMG]
            Xlabel = ["<T₂> = ", " s, Area= "," %"]
        elseif res.seq in [NMRInversions.PFG]
            Xlabel = ["<D> = ", " m²/s, Area= "," %"]
        end

        wa, areas = weighted_averages(res, silent = true)

        for (i,s) in enumerate(res.selections)
            b_low = zeros(length(res.f[s[1]:s[2]])) .+ offset
            b_high = if (offset == 0)
                f_prime[s[1]:s[2]]
            else
                scale_to_one!(f_prime)[s[1]:s[2]] .+ offset
            end

            band!(
                ax2, res.X[s[1]:s[2]], b_low, b_high,
                colormap=:tab10, colorrange=(1, 10), color=i,alpha = 0.35,
            )
        end

        lines!(ax2, res.X, 
               offset == 0 ? real.(f_prime) : (scale_to_one!(real.(f_prime)) .+ offset), 
                          colormap=:tab10, colorrange=(1, 10), color= i, label = res.title)

        for (i,s) in enumerate(res.selections)

            height = (1 - (i-1) * (0.1 + 0.07*log(n_plots))) * maximum(f_prime) + offset 

            sel = collect('a':'z')[i]
            h = f_prime[argmin(abs.(res.X .- wa[i]))]/2 + offset

            scatter!(ax2, 
                     wa[i], 
                     h,
                     markersize=15,
                     marker=sel,
                     colormap=:tab10, colorrange=(1, 10),
                     color=i,
                     glowcolor=:white, glowwidth=4)

            lbl = 
                " " * sel *" : "* Xlabel[1]* "$(round(wa[i], sigdigits = 2))"* 
                Xlabel[2]* "$(round(areas[i]*100, sigdigits = 2))"*Xlabel[3]

            text!(
                ax2, res.X[2], height, 
                text = lbl, 
                align = (:left, :top), 
                colormap=:tab10, colorrange=(1, 10), color=i,
                glowcolor=:white, glowwidth=2
            )
        end
    end
    ax2.limits = (minimum(res.X), maximum(res.X),  -0.025 * maximum(res.f .* res.filter), 1.025 * maximum(res.f .* res.filter))
end
