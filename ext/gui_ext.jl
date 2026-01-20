module gui_ext

using NMRInversions, GLMakie, PolygonOps, LinearAlgebra, NativeFileDialog

include("./gui_1D.jl")
include("./gui_2D.jl")
include("./gui_expfits.jl")
include("./gui_data.jl")

end
