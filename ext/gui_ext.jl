module gui_ext

using NMRInversions, GLMakie, PolygonOps, LinearAlgebra, NativeFileDialog

function lblx(ax::DataAxis)
    if ax isa Union{IR, SR, CPMG}
        return "time (s)"
    elseif ax isa PFG
        return "b factor (s/m² e-9)"
    end
end

function lblX(ax::DataAxis)
    if ax isa Union{IR, SR}
        return L"T_1 \, \textrm{(s)}"
    elseif ax isa CPMG
        return L"T_2 \, \textrm{(s)}"
    elseif ax isa PFG
        return "D (m²/s)"
    end
end

include("./gui_data.jl")
include("./gui_1D.jl")
include("./gui_2D.jl")
# include("./gui_expfits.jl")

end
