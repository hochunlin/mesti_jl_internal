# This Julia script demonstrates the usage of the MUMPS solver for solving sparse linear systems.
# It includes testing the solver with randomly generated sparse matrices.
# This script is taken and modified from https://github.com/wrs28/MUMPS3.jl/blob/5.3.3-update/test/basic_solve.jl

# Set the environment variable "MUMPS_PREFIX" to the path of MUMPS library
ENV["MUMPS_PREFIX"] = "/home/hclin/MUMPS_5.6.0_par_metis/lib"

# Import necessary packages
using MUMPS3
using MPI, LinearAlgebra, SparseArrays, Test

# Check if MPI is initialized, and initialize if not
MPI.Initialized() ? nothing : MPI.Init()

# Define the size of matrix and the number of RHS
N, M = 2000, 20

# Test the functionality in a test set
@testset "Basic solve: " begin
    for i ∈ 1:100
        A = sparse(I,N,N) + sprand(N,N,1/N) # Generate a sparse matrix
        y = sprand(N,M,1/sqrt(N*M)) # Generate a sparse random matrix for the right-hand side

        x = mumps_solve(A,y) # Solve the linear system using MUMPS

        @test norm(A*x-y) ≤ sqrt(eps()) # Test the correctness of the solution
    end
end