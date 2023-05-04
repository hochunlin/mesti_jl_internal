### 20230504 Summary:

1. Rename opts.use_keep_401_1 as opts.parallel_dependency_graph.

2. Add opts.use_single_precision_MUMPS. This enables to use the single-precision MUMPS. We simply convert the double precision to the single precision for the matrices A and B (or K for the APF case).

3. Extend the element data type of matrices A, B, and C: including SparseMatrixCSC{Complex{Int64}, Int64}

   #### mesti_matrix_solver.jl

* Add use_single_precision_MUMP

* Extend the element data type of matrices A, B, and C: including SparseMatrixCSC{Complex{Int64}, Int64}

* Add use_single_precision_MUMPS

  #### mesti2s.jl

* Add use_single_precision_MUMPS


### 20230503 Summary:

1. Add branches for the computation of 2D TM fields Ex(y,z), where only epsilon_xx(y,z) is used. Always use the dimension of epsilon_xx to determine if  2D TM computation is required.
2. Keep as many common lines as possible.
3. Sessions not implemented:
   1. Self-energy
   2. RGF method
   3. codes related to symmetrize_K in mesti2s
   
   #### matrix_build_fdfd_matrix.jl

* Add extra allowed types for input arguments:
  * epsilon_xx: Matrix.
  * epsilon_yy, epsilon_zz, xBC, xPML: Nothing
* Add 2D TM checking term **use_2D_TM** to determine which branch to go if there are differences between 2D and 3D codes.
  
  * use_2D_TM = true if dims(epsilon_xx) = 2
* Add 2D TM branches; No change at 3D branches
* Add a new function handle for 2D TM fields (with "nothing" specified for users) after the main function handle.

  #### setup_longitudinal.jl

* Add extra allowed types for input arguments:

  - kxdx_all: Nothing
* Add extra allowed types for output arguments:

  - side.kxdx_prop: Nothing
* Add 2D TM checking term **use_2D_TM** to determine which branch to go if there are differences between 2D and 3D codes.

  - use_2D_TM = true if kxdx_all = nothing && kLambda_x = nothing && ind_zero_kx = nothing
* Add 2D TM branches; No change at 3D branches

  #### mesti_build_channels.jl

* Add extra allowed types for input arguments:

  * nx_Ex, nx_Ey, xBC, ny_Ey, n0: Nothing
* Add extra allowed types for output arguments:

  * kxdx_all, kxdx_prop: Nothing
* Add 2D TM checking term **use_2D_TM** to determine which branch to go if there are differences between 2D and 3D codes.
  * use_2D_TM = true if nx_Ex == nothing && nx_Ey == nothing && ny_Ex != nothing && ny_Ey == nothing
  * use_2D_TM = true if dims(syst.epsilon_xx) = 2
* Add 2D TM branches; No change at 3D branches
* Add a new function handle for 2D TM fields (with "nothing" specified for users) after the two main function handles.

  #### mesti_matrix_solver.jl

* Edit line 352 by changing "transpose" to its equivalent "permutedims", since the output of the former one has the form Transpose{SparseMatrixCSC{ComplexF64, Int64}} instead of SparseMatrixCSC{ComplexF64, Int64} and the codes throw an error because matrix C does not support this type.
* Exchange **provide_rhs!** and **set_icntl!** at lines 705-706 and 739-739.

  #### mesti.jl

* Add extra allowed types for input arguments:

  * Syst.epsilon_xx: Matrix
  * Syst.epsilon_yy, Syst.epsilon_zz, Syst.xBC, Syst.kx_B: Nothing

* Add 2D TM checking term **use_2D_TM** to determine which branch to go if there are differences between 2D and 3D codes.

  * use_2D_TM = true if dims(epsilon_xx) = 2

* Add 2D TM branches, mainly at adding PMLs, building B and C; No change at 3D branches.

* Edit one line at building matrix C section: 

  from **C_ii[ii] = [C_ii[ii]; [spzeros(M_ii, nxyz_before) transpose(reshape(data, nxyz_data, M_ii)) spzeros(M_ii, nxyz_after)]]**

  to **C_ii[ii] = permutedims([transpose(C_ii[ii]) [spzeros(nxyz_before, M_ii); reshape(data, nxyz_data, M_ii); spzeros(nxyz_after, M_ii)]],(2,1))**

  Since the later concatenating way (the same as building matrix B) seems to be faster than the former one.

  #### mesti2s.jl

* Add 2D TM checking term **use_2D_TM** to determine which branch to go if there are differences between 2D and 3D codes.

  * use_2D_TM = true if dims(epsilon_xx) = 2
* Add 2D TM branches; No change at 3D branches