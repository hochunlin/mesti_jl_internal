using Pkg

Pkg.add("LazyGrids")
Pkg.add("MPI")
Pkg.add("LinearAlgebra")
Pkg.add("SparseArrays")
Pkg.add("Statistics")
Pkg.add("TensorCast")
Pkg.add("Printf")
Pkg.add(PackageSpec(url="http://github.com/wrs28/MUMPS3.jl.git#5.3.3-update")) #MUMPS3.jl, the Julia interface for MUMPS.

# optional 
Pkg.add("MAT")
Pkg.add("Random")
Pkg.add("PyPlot")
Pkg.add("Interpolations")