## Prerequisites for MUMPS installation on Linux

We need the following tools before compiling MUMPS and its Julia interface on Linux.

### Compilers 

The compilation of MUMPS requires both C and Fortran compilers. Although both C and Fortran compilers are included in GNU Compiler Collection (GCC), we recommend using Intel [C](https://www.intel.com/content/www/us/en/developer/tools/oneapi/dpc-compiler.html#gs.gtmcma), [Fortran](https://www.intel.com/content/www/us/en/developer/tools/oneapi/fortran-compiler.html#gs.gtma2f) compilers, and [MPI](https://www.intel.com/content/www/us/en/developer/tools/oneapi/mpi-library.html#gs.gtmtr3) libraries.

On a Linux cluster, if Intel compilers are not loaded, you can load them on a Lmod module system by typing
```shell
module load intel
module load intel intel-mpi
```
In non-Lmod module system, you may type (for Intel oneAPI case),
```shell
source .../intel/oneapi/compiler/latest/env/vars.sh
source  .../intel/oneapi/mpi/latest/env/vars.sh
```
where `...` is the path to Intel folder.



On a local machine, you can download the Intel installer [here](https://www.intel.com/content/www/us/en/developer/tools/oneapi/hpc-toolkit.html). Make sure to use custom installation and choose the following: 
	(a) Intel Inspector
	(b) Intel MPI Library
	(c) Intel oneAPI DPC++/C++ Compiler \& Intel C++ Compiler Classic
	(d) Intel Fortran Compiler \& Intel Fortran Compiler Classic        

After installing the Intel compiler, you need to set the environment variables specific to these compilers. Run the following commands:

```shell
source opt/intel/oneapi/compiler/latest/env/vars.sh
source opt/intel/oneapi/mpi/latest/env/vars.sh
```

These commands will configure the necessary environment variables for the Intel MPI, and Intel compilers.

### BLAS, LAPACK, and ScaLPACK

MUMPS requires both BLAS, LAPACK, and ScaLPACK libraries, which are standard libraries on Linux cluster. These libraries are also included in many implementations, such as MKL and OpenBLAS. 

In the example below, we use MKL. Different clusters install MKL in different paths. You will need the correct path for the linker to find corresponding BLAS, LAPACK, and ScaLPACK libraries. In the provided `Makefile.inc`, we assume the MKL path has been exported to an environment variable called `MKLROOT`. Here we show how to export the correct `MKLROOT`. The BLAS, LAPACK, and ScaLPACK libraries can be found under `$(MKLROOT)/lib/intel64`. 

If you are using the Lmod module system and MKL is installed, you can use 

```shell
module load intel
module load intel-mkl
echo $MKLROOT
```

`MKLROOT`  is assigned and should be printed out, if MKL is available. The `MKLROOT` should be under the path `.../mkl`, where `...` is the Intel directory path.

If Lmod module system is not used in your cluster, then you can also try to find the path like `.../intel/oneapi/mkl`, where `...` is the Intel directory path. Then you can assign the `MKLROOT` by 

```shell
source .../intel/oneapi/mkl/latest/env/vars.sh
echo $MKLROOT
```

`MKLROOT`  is assigned and should be printed out.

On a local machine, please visit the following link: https://www.intel.com/content/www/us/en/developer/tools/oneapi/onemkl.html to acquire MKL. During the installation process, ensure that we select custom installation and choose only the Intel oneAPI Math Kernel Library. This will ensure the necessary libraries are properly installed for MUMPS to function correctly.

After installing the MKL, you can set the environment variable `MKLROOT`. Run the following commands:

```shell
source /opt/intel/oneapi/mkl/latest/env/vars.sh
echo $MKLROOT
```

`MKLROOT`  is assigned and should be printed out.



Note that MUMPS supports shared memory, multithreaded parallelism through the use of multithreaded
BLAS libraries. We provide `Makefile.inc` on Linux and they activate the OpenMP feature. You can use the environment variable `OMP_NUM_THREADS` to set the number of threads.

### METIS

In 3D system, because METIS ordering is more efficient than AMD ordering, we should install the METIS program for graph partitioning (not to be confused with MESTI), including [GKlib](https://github.com/KarypisLab/GKlib) and [METIS](https://github.com/KarypisLab/METIS) . Install them in the following steps:

Install them in the following steps:

(a) Downloading GKlib

```shell
git clone https://github.com/KarypisLab/GKlib.git
```

(b) Installing GKlib

```shell
cd GKlib; make config openmp=set; make install; cd ..;
```

(c) Downloading METIS

```shell
git clone https://github.com/KarypisLab/METIS.git
```

(d) Setting METIS to double precision

```shell
sed -i "43s/32/64/" METIS/include/metis.h
```

(e) Installing METIS

```shell
cd METIS; make config shared=1; make install; cd ..;
```

Then, by default, the library file, header file, and binaries will be installed in `~/local/lib`, `~/local/include`, and `~/local/bin`.

### Running MUMPS in Julia

In some cases, the cluster or local machine cannot find the libraries by itself when you run Julia interface for MUMPS. To solve those issues, you can append those library paths to `LD_PRELOAD` before running Julia. For example, with the Intel oneAPI installed under /opt, you can type,

```shell
source /opt/intel/oneapi/mkl/latest/env/vars.sh
source /opt/intel/oneapi/mpi/latest/env/vars.sh
source /opt/intel/oneapi/compiler/latest/env/vars.sh
export LD_PRELOAD=$LD_PRELOAD:$MKLROOT/lib/intel64/libmkl_intel_lp64.so
export LD_PRELOAD=$LD_PRELOAD:$MKLROOT/lib/intel64/libmkl_intel_thread.so
export LD_PRELOAD=$LD_PRELOAD:/opt/intel/oneapi/inspector/latest/lib64/libiomp5.so
export LD_PRELOAD=$LD_PRELOAD:$MKLROOT/lib/intel64/libmkl_core.so
export LD_PRELOAD=$LD_PRELOAD:$MKLROOT/lib/intel64/libmkl_blacs_intelmpi_lp64.so
export LD_PRELOAD=$LD_PRELOAD:$MKLROOT/lib/intel64/libmkl_scalapack_lp64.so
```



