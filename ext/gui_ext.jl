module gui_ext

using NMRInversions, GLMakie, PolygonOps, LinearAlgebra, NativeFileDialog

x_label(::Union{IR,SR,CPMG}) = "time (s)"
x_label(::PFG) = "b factor (s/m² e-9)"
x_label(::Spectrum) = "δ (ppm)"
x_label(::DataAxis) = ""

X_label(::Union{IR,SR}) = L"T_1 \, \textrm{(s)}"
X_label(::CPMG) = L"T_2 \, \textrm{(s)}"
X_label(::PFG) = L"D \, \textrm{(m^2/s)}"
X_label(::Spectrum) = L"\delta \, \textrm{(ppm)}"
X_label(::DataAxis) = ""

symb(::Union{IR,SR}) = "T₁"
symb(::CPMG) = "T₂"
symb(::PFG) = "D"
symb(::Spectrum) = "δ"
symb(::DataAxis) = ""

include("./gui_data.jl")
include("./gui_1D.jl")
include("./gui_2D.jl")
include("./gui_expfits.jl")

end
