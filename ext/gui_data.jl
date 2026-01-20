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

    surface!(ax_full,data.x_direct, data.x_indirect, real.(data.data))
    scatter!(ax_dir, data.x_direct, real.(data.data)[:,1],color = Cycled(1))
    scatter!(ax_indir, data.x_indirect, real.(data.data)[1,:],color = Cycled(2))

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
