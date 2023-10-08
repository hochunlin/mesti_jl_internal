###### Update on 20231007
module MESTI

using SparseArrays
using LinearAlgebra
using Statistics
using Printf
using TensorCast
using LazyGrids
if haskey(ENV, "MUMPS_PREFIX")
    using MUMPS3
    using MPI
end

include("build_transverse_function_1d.jl")
include("mesti_build_fdfd_matrix.jl")
include("mesti_matrix_solver.jl")
include("setup_longitudinal.jl")
include("mesti_main.jl")
include("mesti_build_channels.jl")
include("mesti2s.jl")

end
