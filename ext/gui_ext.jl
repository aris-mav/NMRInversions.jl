module gui_ext

using NMRInversions, GLMakie, PolygonOps, LinearAlgebra, NativeFileDialog

function lblx(ax::DataAxis)
    if ax isa Union{IR,SR,CPMG}
        return "time (s)"
    elseif ax isa PFG
        return "b factor (s/m² e-9)"
    else
        return ""
    end
end

function lblX(ax::DataAxis)
    if ax isa Union{IR,SR}
        return L"T_1 \, \textrm{(s)}"
    elseif ax isa CPMG
        return L"T_2 \, \textrm{(s)}"
    elseif ax isa PFG
        return L"D \, \textrm{(m^2/s)}"
    else
        return ""
    end
end

function symb(ax::DataAxis)
    if ax isa Union{IR,SR}
        return "T₁"
    elseif ax isa CPMG
        return "T₂"
    elseif ax isa PFG
        return "D"
    else
        return ""
    end
end

include("./gui_data.jl")
include("./gui_1D.jl")
include("./gui_2D.jl")
include("./gui_expfits.jl")

end
