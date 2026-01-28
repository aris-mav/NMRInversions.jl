module NMRInversions

using DelimitedFiles
using LinearAlgebra
using SparseArrays
using NativeFileDialog
using PolygonOps
using Optim 

include("types.jl")
include("import_data.jl")
include("kernels.jl")
include("finding_alpha.jl")
include("optim_regularizations.jl")
include("L1_regularization.jl")
include("inversion_functions.jl")
include("exp_fits.jl")
include("misc.jl")

end
