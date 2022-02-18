# libParanumal
Forked from https://github.com/paranumal/libparanumal

## ACOUSTICS MODULE
> Implemented by
> * Anders Melander (<s144277@student.dtu.dk>)
> * Emil Strøm (<s144259@student.dtu.dk>)

> * Maintained by Nikolas Borrel-Jensen (<nibor@elektro.dtu.dk>)

> Papers
> * Pind, F., Engsig-Karup, A. P., Jeong, C.-H., Hesthaven, J. S., Mejling, M. S., & Strømann-Andersen, J. (2019). Time domain room acoustic simulations using the spectral element method. The Journal of the Acoustical Society of America, 145(6), 3299–3310. https://doi.org/10.1121/1.5109396
> * Melander. (2021). Massively Parallel Nodal Discontinous Galerkin Finite Element Method Simulator for Room Acoustics. <https://doi.org/10.1177/ToBeAssigned>

### INSTALLATION
1. To access DTUs HPC system, log in using the command (see section about HPC system as well) <br>
    `> ssh username@login1.gbar.dtu.dk`
2. Clone the code from git **in the home user directory** `~/` <br>
    `> git clone https://github.com/dtu-act/libparanumal`
3. OCCA 1.0 is present in the libparanumal folder. **IMPORTANT**: do not use the version from git since this version of libParanumal is out-of-sync (OCCA is a third-party library used for compiling code to various platforms, such as CUDA/GPU)
4. Build OCCA and libParanumal. **IMPORTANT**: to exploit the GPU, you should build in an environment with access to GPUs. On DTUs systems, do (see also section below about DTU HPC) <br>
    `> voltash`<br>
    Enter the acoustics folder <br>
    `> cd ~/libparanumal/solvers/acoustics` <br>
    and then execute the build script <br>
    `> ./build_acoustics.sh` <br>

### RUNNING THE CODE
5. For testing the setup, run the following command from location `libparanumal/solvers/acoustics` (output is written into `libparanumal/solvers/acoustics/data`). **NOTE:** You might run out of resources on DTU HPC system when running without the queue, resulting in "out of memory" messages. <br>
    `> ./RUNGPU.bsub` <br>    
6. Several examples including frequency-dependent and independent cases, can be found inside `libparanumal/solvers/acoustics/tests/`

### GMSH
* use Gmsh (http://gmsh.info) to create meshes (see 'Massively Parallel Nodal Discontinous Galerkin Finite Element Method Simulator for 3D Room Acoustics' thesis for details)
* Tutorial: https://www.youtube.com/watch?v=xL2LmDsDLYw

### ONLY FOR USERS OF THE HPC SYSTEM AT THE TECHNICAL UNIVERSITY OF DENMARK
#### HPC system
1. Login using <br>
    `> ssh username@login1.gbar.dtu.dk`
2. change to GPU capabilities by typing <br> 
    `> voltash`
3. you can run interactively by executing (after having compiled OCCA and libParanumal)<br> 
    `> ./RUNGPU.bsub`
4. when running long-running tasks, use the queue <br> 
    `bsub < <the_script>.bsub`

#### References
* VPN: https://www.inside.dtu.dk/en/medarbejder/it-og-telefoni/wifi-og-fjernadgang/vpn-cisco-anyconnect
* HPC: https://www.hpc.dtu.dk/
* HPC GPU parameters: https://www.hpc.dtu.dk/?page_id=2759
* LSF job submission system: http://www.cc.dtu.dk/?page_id=1416
* HPC examples: https://www.hpc.dtu.dk/?page_id=2021
* Gmsh: http://gmsh.info
* Tutorial: https://www.youtube.com/watch?v=xL2LmDsDLYw
* DGFEM simulator: https://github.com/dtu-act/libparanumal

#### Useful UNIX commands
* `> ssh username@login1.gbar.dtu.dk`   # login
* `> voltash`                           # switch to GPU cluster
* `> bsub < <the_script>.bsub`          # add to the queue system
* `> bstat`                             # job status
* `> bkill <id>`                        # kill job
* `> ./the_scripts > logfile.txt`       # run the script and pipe the output to a log file

### Installation on MacOS X
#### Compiling libparanumal

* Install OCCA (see above). Set `-j ‘numthreads’` accordingly when installing, otherwise installation can stall
* Install fortran: https://github.com/fxcoudert/gfortran-for-macOS/releases [not sure this is needed - try without and update this doc]
* Install libomp: https://stackoverflow.com/questions/43555410/enable-openmp-support-in-clang-in-mac-os-x-sierra-mojave
* Install OpenMPI: https://www.open-mpi.org/faq/?category=building
	- ./configure --prefix=/usr/local CC=gcc CXX=g++ FC=gfortran [not sure FC=gfortran is needed - try without and update this doc]

If the following error is shown

    ---[ Error ]------------------------------------------------
        File     : /Users/nikolasborrel/Documents/PhD/libparanumal/occa/src/io/utils.cpp
        Function : write
        Line     : 372
        Message  : Failed to open [c76425872528f884/compilerVendorTest.cpp]

delete folder `~/.occa`. Can also be fixed with chmod 777 (permission is probably messed up)

#### VS CODE
A few references that might be helpful for setting up VS CODE

* https://code.visualstudio.com/docs/editor/debugging
* https://code.visualstudio.com/docs/cpp/config-clang-mac
* https://code.visualstudio.com/docs/cpp/c-cpp-properties-schema-reference
* https://www2.cs.duke.edu/courses/cps108/doc/makefileinfo/sample.html

---
## General
An experimental set of finite element flow solvers for heterogeneous (GPU/CPU) systems. The initial development of libParanumal was performed by the [Parallel Numerical Algorithms Group at Virginia Tech](http://paranumal.com).   

libParanumal is funded in part by the US Department of Energy as part of the activities of the [Center for Efficient Exscale Discretizations](http://ceed.exascaleproject.org). 

Why is it called libParanumal ?: the high-order finite-element implementations in libParanumal are __spectrally__ accurate and rely heavily on __ghost__ elements for MPI communications.

If you use libParanumal as part of a research project see Section 8 below for papers to reference.

---
### 1. Overview 

Brief summary of major features:

A. Supported elements:
  - Triangles, quadrilaterals, tetrahedra, hexahedra.
  - Lagrange basis functions up to degree 15.
  - Partial support for Bezier-Bernstein basis functions.
  
B. Mesh wrangling:
  - Gmsh format file loaders.
  - Load balanced geometric partitioning using space filling curves (Hilbert or Morton ordering). 
  - Clustered partitioning for multirate time stepping.
  
C. Elliptic solver:
  - Linear Poisson and screened Poisson potential solvers.
  - GPU optimized matrix-vector products.
  - Hybrid p-type multigrid and algebraic multigrid  preconditioned conjugate gradient solver.
  - Sparse matrix or nearly matrix-free algebraic multigrid for coarse levels of multigrid hierarchy.

D. Heterogeneous accelerated flow solvers:
  - Linearized Euler equations.
  - Isothermal compressible Navier-Stokes solver with:
     * Upwind discontinuous Galerkin discretization in space.
     * Dormand-Prince adaptive Runge-Kutta integration in time.
  - Isothermal Galerkin-Boltzmann gas dynamics solver with:
     * Penalty flux DG discretization in space.
     * Adaptive semi-analytic (pointwise exponential) integration in time.
     * Multiaxial quasi-perfectly matched absorbing layer far field boundary condition.
  - Incompressible Navier-Stokes solver with:
     * Choice of continuous FEM or interior penalty DG in space.
     * Extrapolation-BDF integration in time.
     * Sub-cycling (Operator Integration Factor Splitting) for advection.

E. Dependencies:
   - Message Passing Interface (MPI).
      * The libParanumal makefiles assume that mpic++ are installed and visible in your path.     
   - Open Concurrent Compute Abstraction (OCCA) 
      * OCCA must be installed.
      * OCCA will try to detect if any of these execution models are installed OpenMP, CUDA, OpenCL, HIP.
      * If OCCA does not detect any of these it will default to Serial execution.
      * You will need to adjust the libParnumal setup input files to choose the execution model and compute device appropriate for your system.
      * The OCCA github repo is [here](https://github.com/libocca/occa)
      * The OCCA webpage is [here](http://libocca.org)
      

---
### 2. Code block diagram 
<img src="http://www.math.vt.edu/people/tcew/libParanumalNekDiagramFA18-crop.png" width="1024" >

---
### 3. Clone: libParanumal
`git clone https://github.com/paranumal/libparanumal`

---
### 4. OCCA dependency (currently OCCA 1.0 forked by Noel Chalmers) 
`git clone https://github.com/noelchalmers/occa`

#### 4-1. Build OCCA 
`cd occa`    
export OCCA_DIR=\`pwd\`  
`export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$OCCA_DIR/lib`    
`make -j`    
`cd ../  `  

---
### 5. Running the codes: 

The elliptic solver and flow solvers reside in sub-directories of the solver directory. Each sub-directory includes makefile, src directory, data directory (including header files for defining boundary conditions), okl kernel directory, and setups directory. The setups directory includes a number of example input files that specify input parameters for the solver.

#### 5-1. Build libParanumal elliptic example
  
`cd libparanumal/solvers/elliptic`    
`make -j  `  

#### 5-2. Run elliptic example with provided quadrilateral set up file on a single device:
  
`./ellipticMain setups/setupQuad2D.rc`  

#### 5-3. Run the same example with two devices:

`mpiexec -n 2 ./ellipticMain setups/setupQuad2D.rc`  
 
---

### 6. License

The MIT License (MIT)

Copyright (c) 2017 Tim Warburton, Noel Chalmers, Jesse Chan, Ali Karakus

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

### 7. References

Discontinuous Galerkin Boltzmann (bns) solver [arXiv version](https://arxiv.org/abs/1805.02082): `Karakus, A., Chalmers, N., Hesthaven, J.S. and Warburton, T., 2018. Discontinuous Galerkin Discretizations of the Boltzmann Equations in 2D: semi-analytic time stepping and absorbing boundary layers. arXiv preprint arXiv:1805.02082.`

Incompressible Navier-Stokes (discontinuous) Galerkin (ins) solver [arXiv version](https://arxiv.org/abs/1801.00246): `Karakus, A., Chalmers, N., Swirydowicz, K. and Warburton, T., 2017. GPU Acceleration of a High-Order Discontinuous Galerkin Incompressible Flow Solver. arXiv preprint arXiv:1801.00246.`

Optimization of elliptic mat-vec operations for (elliptic) solver on hexes [arXiv version](https://arxiv.org/abs/1711.00903): `Świrydowicz, K., Chalmers, N., Karakus, A. and Warburton, T., 2017. Acceleration of tensor-product operations for high-order finite element methods. arXiv preprint arXiv:1711.00903.`


