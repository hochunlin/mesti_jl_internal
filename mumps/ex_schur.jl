ENV["MUMPS_PREFIX"] = "/home/hclin/MUMPS_5.6.0/lib";

using MUMPS3, MPI, LinearAlgebra, SparseArrays, Random;
MPI.Initialized() ? nothing : MPI.Init()

m = 100; n = 5; p = .5; T = ComplexF64; rng = MersenneTwister(1234);
A = I + sprand(rng,T,m,m,p);
A⁻¹ = inv(Matrix(A));
B = sprand(rng,T,m,n,p);
C = sprand(rng,T,n,m,p);
D = sprand(rng,T,n,n,p);
M = [A B; C D]; #Construct a (p+q)×(p+q) matrix M.
id = Mumps(M, sym=0, par=1); #Initializes a Mumps object with the same type as M
set_schur_centralized_by_column!(id, 101:105); #Specify q in the matrix M.
set_job!(id,1); #Specify that Mumps will perform an analysis.
set_icntl!(id,7,5); #Specify the Metis ordering.
invoke_mumps!(id); #Perform the analysis.
set_job!(id,2); #Specify that Mumps will perform a factorization.
invoke_mumps!(id); #Perform the factorization.
S = get_schur_complement(id); #Take the Schur Complement.
norm(D - C*A⁻¹*B - S) #Check the error of the Schur Complement.
println("norm(D - C*A⁻¹*B - S)=", D - C*A⁻¹*B - S")