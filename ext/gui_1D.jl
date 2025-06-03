## Plots for 1D inversions

"""
    plot(res_mat::VecOrMat{inv_out_1D};kwargs...)

Plot the results contained in a vector of `inv_out_1D` structures.

The keyword (optional) arguments are:

- `selections` : Whether to highlight the selections in the plots (default is `false`).

"""
function Makie.plot(res_mat::AbstractVecOrMat{NMRInversions.inv_out_1D}; selections = false , yscale = identity)

    fig = Figure(size=(700, 600))

    res = res_mat[1]
    # Make axes
    if res.seq in [NMRInversions.IR]
        ax1 = Axis(fig[1:5, 1:5], xlabel="time (s)", ylabel="Signal (a.u.)", yscale = yscale)
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

    ax1.yscale = yscale

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
        \n Weighted average = $(round(weighted_average[], sigdigits=3))"
    end

    Label(fig[7,6:10], labeltext, tellwidth=false)
    vlines!(fig.content[2], int_low, color=:red)
    vlines!(fig.content[2], int_high, color=:red)

    button_label = Button(fig[8,6:8], label = "Label selection")
    button_filter = Button(fig[8,8:10], label = "Filter-out selection")

    button_reset = Button(fig[9,6:10], label = "Reset selections")

    log_y = false
    button_yscale = Button(fig[10,6:8], label = "Change y scale")
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
        f = plot([res], selections = true, yscale = log_y ? log10 : identity)
        save(savedir, f)
    end

    on(button_yscale.clicks) do _
        if log_y
            fig.content[1].yscale = identity
            log_y = false
        else
            fig.content[1].yscale = log10
            log_y = true
        end

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

    i = length(ax2.scene.plots)
    scatter!(ax1, res.x, real.(res.y), colormap=:tab10, colorrange=(1, 10), color=i)
    lines!(ax1, res.xfit, yfit_prime, colormap=:tab10, colorrange=(1, 10), color=i)
    lines!(ax2, res.X, f_prime, colormap=:tab10, colorrange=(1, 10), color=i)
    lines!(ax3, res.x, r_prime, colormap=:tab10, colorrange=(1, 10), color=i)
    scatter!(ax3, res.x, r_prime, colormap=:tab10, colorrange=(1, 10), color=i)


    if selections

        if res.seq in [NMRInversions.IR]
            Xlabel = ["<T₁> = ", " (s), Area= "," %"]
        elseif res.seq in [NMRInversions.CPMG]
            Xlabel = ["<T₂> = ", " (s), Area= "," %"]
        elseif res.seq in [NMRInversions.PFG]
            Xlabel = ["<D> = ", " (m²/s), Area= "," %"]
        end

        wa, areas = weighted_averages(res)

        for (i,s) in enumerate(res.selections)

            band!(
                ax2, 
                res.X[s[1]:s[2]],
                zeros(length(res.f[s[1]:s[2]])),
                (res.f .* res.filter)[s[1]:s[2]],
                colormap=:tab10, colorrange=(1, 10), color=i,alpha = 0.5
            )

            height = (1 - (i-1) * 0.1) * maximum(res.f .* res.filter) 
            text!(
                ax2, res.X[2], height, 
                text = Xlabel[1]*"$(round(wa[i], sigdigits = 2))"*Xlabel[2]*"$(round(areas[i]*100, sigdigits = 2))"*Xlabel[3], 
                align = (:left, :top), 
                colormap=:tab10, colorrange=(1, 10), color=i
            )

            #=vlines!(ax2, res.X[s[1]], color = clr, alpha = 0.8)=#
            #=vlines!(ax2, res.X[s[2]],color = clr, alpha = 0.8)=#
        end
    end
end
