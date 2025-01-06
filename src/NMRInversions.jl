module NMRInversions

# Import dependencies
using DelimitedFiles
using LinearAlgebra
using SparseArrays
using NativeFileDialog
using PolygonOps
using Optim 

"""
to do list:
- introduce faf, flip angle fraction, to the kernel functions. 1 would be a perfect pulse, 0 would be no pulse.
- add precompilation

"""

## Include the package files 
include("types.jl")
include("inversions_io.jl")
include("kernels.jl")
include("finding_alpha.jl")
include("optim_regularizations.jl")
include("L1_regularization.jl")
include("inversion_functions.jl")
include("exp_fits.jl")
include("misc.jl")

# Export useful functions
export invert
export create_kernel
export import_csv
export import_spinsolve, import_geospec

# The following functions are defined here and modified from extension files
function select_peaks end
function pubfig end
export select_peaks, pubfig

end
