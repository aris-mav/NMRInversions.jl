"""
plot(data::input1D)

Quickly visualise the contents of a 1D NMR data structure.
"""
function Makie.plot(data::NMRInversions.input1D)
    fig = Figure()
    ax = Axis(fig[:,:])

    if any(isa.(data.y, Complex))
        plot!(ax,data.x,imag.(data.y),color = Cycled(2), label= "Complex part")
        plot!(ax,data.x,real.(data.y),color = Cycled(1), label= "Real part")
        axislegend(ax)
    else
        plot!(ax,data.x,data.y)
    end

    ax.ylabel = "Signal (a.u.)"
    if data.seq in (IR, SR)
        ax.xlabel = "time (s)"
    elseif data.seq == CPMG
        ax.xlabel = "time (s)"
    elseif data.seq == PFG
        ax.xlabel = "b factor (s/m² e-9)"
    end

    return fig
end


"""
plot(data::input2D)

Interactive GUI to visualise the contents of a 2D NMR data structure.
"""
function Makie.plot(data::NMRInversions.input2D)

    fig = Figure()
    ax_dir = Axis(fig[1:4,1:5], ylabel= "Signal (a.u.)", 
                  title= "Direct dimension")
    ax_indir = Axis(fig[6:9,1:5], ylabel= "Signal (a.u.)",
                    title = "Indirect dimension")
    ax_full = Axis3(fig[1:9,6:10], zticklabelsize = 0.1, zlabelvisible=false)

    d = real.(data.data)
    xd = data.x_direct
    xi = data.x_indirect

    sld = Slider(fig[5, 1:5], range = 1:length(xi), startvalue = 1)
    sli = Slider(fig[10, 1:5], range = 1:length(xd), startvalue = 1)

    indir_slice = @lift(d[$(sli.value),:])
    dir_slice = @lift(d[:,$(sld.value)])

    scatter!(ax_dir, xd, dir_slice, color = Cycled(1))
    scatter!(ax_indir, xi, indir_slice, color = Cycled(2))
    surface!(ax_full,xd, xi, d, alpha=0.5)
    scatter!(ax_full, xd , @lift(fill(xi[$(sld.value)] ,length(xd)) ), dir_slice, color = Cycled(1), markersize = 10)
    scatter!(ax_full, @lift(fill(xd[$(sli.value)] ,length(xi)) ), xi , indir_slice, color = Cycled(2), markersize = 10)

    y_low = minimum(d) * 1.1
    y_high = maximum(d) * 1.1
    ax_dir.limits = (nothing,nothing,y_low,y_high)
    ax_indir.limits = (nothing,nothing,y_low,y_high)

    ax_indir.ylabel = "Signal (a.u.)"
    ax_dir.xlabel = "time (s)"
    ax_full.xlabel = "time (s)"

    if data.seq == IRCPMG
        ax_indir.xlabel = "time (s)"
        ax_full.ylabel = "time (s)"
    elseif data.seq == PFGCPMG
        ax_indir.xlabel = "b factor (s/m² e-9)"
        ax_full.ylabel = "b factor (s/m² e-9)"
    end

    return fig
end
