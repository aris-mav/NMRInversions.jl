"""
    plot(data::ExperimentData)

Quickly visualise the contents of a 1D NMR data structure.
"""
function Makie.plot(data::NMRInversions.ExperimentData{1})

    GLMakie.activate!(; title="Data visualisation")

    fig = Figure()
    ax = Axis(fig[:, :])

    x = data.axes[1]

    if any(isa.(data.data, Complex))

        plot!(ax, x, imag.(data.data),
            color=Cycled(2), label="Complex part")

        plot!(ax, x, real.(data.data),
            color=Cycled(1), label="Real part")
    else
        plot!(ax, x, data.data, label="Data points")
    end

    fig[2, 1] = Legend(fig, ax, orientation=:horizontal,
        tellwidth=false, tellheight=true)

    ax.xlabel = lblx(x)
    ax.ylabel = "Signal (a.u.)"

    return fig
end


"""
    plot(data::ExperimentData)

Interactive GUI to visualise the contents of a 2D NMR data structure.
"""
function Makie.plot(data::NMRInversions.ExperimentData{2})

    GLMakie.activate!(; title="Data visualisation")

    dt = real.(data.data)
    xd = data.axes[1]
    xi = data.axes[2]

    if any(map(x -> length(x) > 16000, [xd, xi]))

        throw("Dataset too big to plot, Makie throws errors. \
              Each dimension should have less than 16000 points.")
    end

    fig = Figure()
    ax_dir = Axis(fig[1:4, 1:5], ylabel="Signal (a.u.)",
        title="Direct dimension ($(nameof(typeof(data.axes[1]))))")
    ax_indir = Axis(fig[6:9, 1:5], ylabel="Signal (a.u.)",
        title="Indirect dimension ($(nameof(typeof(data.axes[2]))))")
    ax_full = Axis3(fig[1:9, 6:10], zticklabelsize=0.1, zlabelvisible=false)

    sld = Slider(fig[5, 1:5], range=1:length(xi), startvalue=1)
    sli = Slider(fig[10, 1:5], range=1:length(xd), startvalue=1)

    indir_slice = @lift(dt[$(sli.value), :])
    dir_slice = @lift(dt[:, $(sld.value)])

    scatter!(ax_dir, xd, dir_slice, color=Cycled(1))
    scatter!(ax_indir, xi, indir_slice, color=Cycled(2))
    surface!(ax_full, xd, xi, dt, alpha=0.5)
    scatter!(
        ax_full,
        xd, @lift(fill(xi[$(sld.value)], length(xd))), dir_slice,
        color=Cycled(1), markersize=10
    )
    scatter!(
        ax_full,
        @lift(fill(xd[$(sli.value)], length(xi))), xi, indir_slice,
        color=Cycled(2), markersize=10
    )
    scatter!(
        ax_dir,
        @lift(xd[$(sli.value)]), @lift(dt[$(sli.value), $(sld.value)]),
        color=Cycled(2)
    )
    scatter!(
        ax_indir,
        @lift(xi[$(sld.value)]), @lift(dt[$(sli.value), $(sld.value)]),
        color=Cycled(1)
    )

    y_low = minimum(dt) * 1.1
    y_high = maximum(dt) * 1.1

    ax_dir.limits = (nothing, nothing, y_low, y_high)
    ax_indir.limits = (nothing, nothing, y_low, y_high)

    ax_indir.xlabel = lblx(xi)
    ax_dir.xlabel = lblx(xd)

    ax_full.xlabel = lblx(xd)
    ax_full.ylabel = lblx(xi)

    return fig
end
