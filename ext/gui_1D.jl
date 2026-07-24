## Plots for 1D inversions
using NMRInversions
"""
    plot(res_mat::VecOrMat{InversionData};kwargs...)

Plot the results contained in a vector of `InversionData` structures.

The keyword (optional) arguments are:

- `selections` : Whether to highlight the selections in each of the plots (default is `true`).
- `stacked` : Whether to plot different curves as a stack, instead of superimposed. (default is `true`).

"""
function Makie.plot(
    res_mat::AbstractVecOrMat{NMRInversions.InversionData{1}};
    selections::Bool=true,
    stacked::Bool=true,
    yscale=identity,
    xscale=identity,
)

    fig = Figure(size=(700, 600))

    res = res_mat[1]

    position1 = fig[1:5, 1:5]
    position2 = (stacked && length(res_mat) > 1) ? fig[1:10, 6:10] : fig[1:5, 6:10]
    position3 = fig[7:11, 1:5]

    ax1 = Axis(
        position1,
        xlabel=x_label(res.input.axes[1]), ylabel="Signal (a.u.)",
        yscale=yscale, xscale=xscale
    )
    ax2 = Axis(
        position2,
        xlabel=X_label(res.axes[1]), xscale=log10
    )
    ax3 = Axis(
        position3, xlabel=x_label(res.input.axes[1]), ylabel="Residuals (a.u.)",
        xscale=xscale
    )

    ax2.xlabelsize = 18

    linkxaxes!(ax1, ax3)

    for (i, r) in enumerate(res_mat)
        static_plots(ax1, ax2, ax3, r,
            selections=selections,
            offset=stacked ? i * 1.1 : 0,
            n_plots=length(res_mat),
            i=i,
        )
    end

    ax1.yscale = yscale
    ax1.xscale = xscale

    if length(res_mat) > 1
        if stacked
        else
            Legend(fig[6:11, 6:10], ax2)
        end
    end

    ax2.limits = (
        minimum(minimum(vec) for r in res_mat for vec in r.axes),
        maximum(maximum(vec) for r in res_mat for vec in r.axes),
        stacked ? nothing : minimum(map(res -> -0.025 * maximum(res.data .* res.filter), res_mat)),
        stacked ? nothing : maximum(map(res -> 1.025 * maximum(res.data .* res.filter), res_mat)),
    )


    if stacked
        for a in (ax1, ax2, ax3)
            a.yticksvisible = false
            a.yticklabelsvisible = false
            a.yticks = 1.1:1.1:length(res_mat)*1.1
        end
    end

    return fig
end


function redraw(fig::Figure, res::InversionData{1}, int_low, int_high)
    empty!(fig.content[1])
    empty!(fig.content[2])
    empty!(fig.content[3])
    static_plots(fig.content[1], fig.content[2], fig.content[3], res, selections=true)
    vlines!(fig.content[2], int_low, color=:red)
    vlines!(fig.content[2], int_high, color=:red)
end

"""
    plot(res::InversionData)

Run the GUI to plot the 1D inversion results and select peaks you want to label.
"""
function Makie.plot(res::NMRInversions.InversionData{1})

    GLMakie.activate!(; title="1D inversion GUI")
    fig = plot([res], selections=true, stacked=false)

    slider = IntervalSlider(
        fig[6, 6:10],
        range=[1:length(res.axes[1])...],
        startvalues=(1, length(res.axes[1]))
    )

    interval_low = lift(slider.interval) do i
        return res.axes[1][i[1]]
    end

    interval_high = lift(slider.interval) do i
        return res.axes[1][i[2]]
    end

    weighted_average = lift(slider.interval) do i
        return res.data[i[1]:i[2]]' * res.axes[1][i[1]:i[2]] / sum(res.data[i[1]:i[2]])
    end

    labeltext = lift(slider.interval) do i
        return """
        Selected range from $(round(res.axes[1][i[1]], sigdigits=2)) \
        to $(round(res.axes[1][i[2]], sigdigits=2)).
        Weighted average = $(round(weighted_average[], sigdigits=3))
        """
    end

    Label(fig[7, 6:10], labeltext, tellwidth=false)
    vlines!(fig.content[2], interval_low, color=:red)
    vlines!(fig.content[2], interval_high, color=:red)

    button_label = Button(fig[8, 6:8], label="Label selection")
    button_filter = Button(fig[8, 8:10], label="Filter-out selection")
    button_reset_s = Button(fig[9, 6:8], label="Reset selections")
    button_reset_f = Button(fig[9, 8:10], label="Reset filter")
    button_xscale = Button(fig[10, 6:8], label="Change x scale")
    button_yscale = Button(fig[10, 8:10], label="Change y scale")
    button_save = Button(fig[11, 6:10], label="Save and exit")

    log_y = false
    log_x = false

    on(button_yscale.clicks) do _
        if !any(real.(res.input.data) .<= 0)
            if log_y
                fig.content[1].yscale = identity
                log_y = false
            else
                fig.content[1].yscale = log10
                log_y = true
            end
        else
            @warn("Negative values present. Can't go to log scale.")
        end
    end

    on(button_xscale.clicks) do _
        if log_x
            fig.content[1].limits = (
                exp10(floor(log10(res.axes[1][1]))),
                nothing, nothing, nothing
            ) # looks weird but it prevents an error
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

    on(button_label.clicks) do _
        low, high = slider.interval[]
        selection = [[res.axes[1][low]], [res.axes[1][high]]]
        push!(res.selections, selection)
        redraw(fig, res, interval_low, interval_high)
    end

    on(button_filter.clicks) do _
        res.filter[slider.interval[][1]+1:slider.interval[][2]] .= 0
        scale_filter!(res)
        redraw(fig, res, interval_low, interval_high)
    end

    on(button_reset_s.clicks) do _
        empty!(res.selections)
        redraw(fig, res, interval_low, interval_high)
    end

    on(button_reset_f.clicks) do _
        res.filter = ones(size(res.data))
        redraw(fig, res, interval_low, interval_high)
    end

    on(button_save.clicks) do _
        savedir = NMRInversions.save_file(res.title, filterlist="png")
        f = plot([res], selections=true,
            yscale=log_y ? log10 : identity,
            xscale=log_x ? log10 : identity
        )
        save(savedir, f, px_per_unit=2.5)
    end


    return fig
end


"""
- `offset` normalises the data to 1 and plots them as `data .+ offset`, used to stack many graphs.
"""
function static_plots(ax1, ax2, ax3, res::NMRInversions.InversionData{1};
    selections=false, offset::Real=0, n_plots::Int=1, i=1)

    W_½ = sqrt(res.input.W[1])
    f = res.data .* res.filter
    x = res.input.axes[1]
    g = res.input.data
    wg = W_½ * g

    xfit = typeof(x)(
        exp10.(range(log10(x[1]), log10(x[end]), 512))
    )

    X = (x isa PFG ? res.axes[1] .* 1e9 : res.axes[1])

    yfit = create_kernel(xfit, X, y=wg) * f

    r = W_½ * create_kernel(x, X, y=wg) * f - real.(wg)

    scatter!(
        ax1, x,
        offset == 0 ? real.(g) : (scale_to_one!(real.(g)) .+ offset),
        colormap=:tab10, colorrange=(1, 10), color=i
    )
    lines!(
        ax1, xfit,
        offset == 0 ? real.(yfit) : (scale_to_one!(real.(yfit)) .+ offset),
        colormap=:tab10, colorrange=(1, 10), color=i
    )
    !selections && lines!(
        ax2, res.axes[1],
        offset == 0 ? real.(f) : (scale_to_one!(real.(f)) .+ offset),
        colormap=:tab10, colorrange=(1, 10), color=i, label=res.title
    )
    lines!(
        ax3, x,
        offset == 0 ? real.(r) : (scale_to_one!(real.(r)) .* 0.5 .+ offset),
        colormap=:tab10, colorrange=(1, 10), color=i
    )
    scatter!(
        ax3, x,
        offset == 0 ? real.(r) : (scale_to_one!(real.(r)) .* 0.5 .+ offset),
        colormap=:tab10, colorrange=(1, 10), color=i
    )

    if selections

        if res.input.axes[1] isa Union{IR,SR}
            text_labels = ["<T₁> = ", " s, Area= ", " %"]
        elseif res.input.axes[1] isa CPMG
            text_labels = ["<T₂> = ", " s, Area= ", " %"]
        elseif res.input.axes[1] isa PFG
            text_labels = ["<D> = ", " m²/s, Area= ", " %"]
        end

        wa, areas = NMRInversions.weighted_averages(res, silent=true)

        for (i, idx) in enumerate(NMRInversions._selection_indices(res))

            b_low = zeros(length(res.data[idx])) .+ offset
            b_high = if (offset == 0)
                f[idx]
            else
                scale_to_one!(f)[idx] .+ offset
            end

            band!(
                ax2, res.axes[1][idx], b_low, b_high,
                colormap=:tab10, colorrange=(1, 10), color=i, alpha=0.35,
            )
        end

        lines!(ax2, res.axes[1],
            offset == 0 ? real.(f) : (scale_to_one!(real.(f)) .+ offset),
            colormap=:tab10, colorrange=(1, 10), color=i, label=res.title)

        for i in eachindex(res.selections)

            height = (1 - (i - 1) * (0.1 + 0.07 * log(n_plots))) * maximum(f) + offset

            sel = collect('a':'z')[i]
            h = f[argmin(abs.(res.axes[1] .- wa[i]))] / 2 + offset

            scatter!(ax2,
                wa[i],
                h,
                markersize=15,
                marker=sel,
                colormap=:tab10, colorrange=(1, 10),
                color=i,
                glowcolor=:white, glowwidth=4)

            txt =
                " " * sel * " : " * text_labels[1] * "$(round(wa[i], sigdigits = 2))" *
                text_labels[2] * "$(round(areas[i]*100, sigdigits = 2))" * text_labels[3]

            text!(
                ax2, res.axes[1][2], height,
                text=txt,
                align=(:left, :top),
                colormap=:tab10, colorrange=(1, 10), color=i,
                glowcolor=:white, glowwidth=2
            )
        end
    end

    ax2.limits = (
        minimum(minimum(vec) for vec in res.axes),
        maximum(maximum(vec) for vec in res.axes),
        -0.025 * maximum(f),
        1.025 * maximum(f),
    )

end
