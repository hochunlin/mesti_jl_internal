ENV["MUMPS_PREFIX"] = "/home/hclin/MUMPS_5.6.0_par_metis/lib"
# Include MESTI module
include("MESTI.jl")
using .MESTI
using LinearAlgebra, SparseArrays, Test

include("matrix_solver_test.jl")
include("interface_t_r_test.jl")
include("unitary_test.jl")
