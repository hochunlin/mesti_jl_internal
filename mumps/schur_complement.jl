# This Julia script demonstrates the usage of the MUMPS solver for solving sparse linear systems.
# It includes testing the Schur complement with randomly generated sparse matrices.
# This script is taken and modified from https://github.com/wrs28/MUMPS3.jl/blob/5.3.3-update/test/schur_complement.jl

# Set the environment variable "MUMPS_PREFIX" to the path of MUMPS library
ENV["MUMPS_PREFIX"] = "/home/hclin/MUMPS_5.6.0_par_metis/lib"

# Import necessary packages
using MUMPS3
using MPI, LinearAlgebra, SparseArrays, Test

# Check if MPI is initialized, and initialize if not
MPI.Initialized() ? nothing : MPI.Init()

# Define the size of matrix and the number of RHS
N, M = 2000, 20

A = sparse(I,N,N) + sprand(N,N,1/N) # Generate a sparse matrix

invA = inv(Matrix(A)) # Calculate the inverse of matrix A

# Test set for Schur complement calculation
@testset "Schur complement: " begin
    for i ∈ 1:100

        y = sprand(N,M,1/sqrt(N*M)) # Generate a sparse random matrix for the right-hand side

        S = mumps_schur_complement(A,y) # Calculate the Schur complement using the MUMPS solver

        schur_inds = unique!(sort(y.rowval)) # Identify the indices corresponding to the rows and columns of the Schur complement

        @test norm(inv(S) - invA[schur_inds,schur_inds]) ≤ sqrt(eps()) # Test the correctness of the Schur complement calculation
    end
end
