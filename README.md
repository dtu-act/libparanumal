# libParanumal - ACOUSTICS MODULE
---
Forked from https://github.com/paranumal/libparanumal

Implemented by [Anders Melander (DTU)](mailto:adame@dtu.dk) and [Emil Strøm (DTU)](mailto:s144259@student.dtu.dk) <br>
Maintained by [Nikolas Borrel-Jensen (DTU)](mailto:nibor@dtu.dk)

Papers
* Pind, F., Engsig-Karup, A. P., Jeong, C.-H., Hesthaven, J. S., Mejling, M. S., & Strømann-Andersen, J. (2019). Time domain room acoustic simulations using the spectral element method. The Journal of the Acoustical Society of America, 145(6), 3299–3310. [https://doi.org/10.1121/1.5109396](https://doi.org/10.1121/1.5109396)
* Melander. (2021). Massively Parallel Nodal Discontinous Galerkin Finite Element Method Simulator for Room Acoustics. [https://doi.org/10.1177/ToBeAssigned](https://doi.org/10.1177/ToBeAssigned)

---

## INSTALLATION
1. Clone the code from git in you home directory <br>
    `> git clone https://github.com/dtu-act/libparanumal`
2. Build OCCA and libParanumal. To exploit the GPU, you should build in an environment with access to GPUs (see section DTU HPC) <br>
    `> voltash` (enabling GPU environment) <br>
    `> cd ~/libparanumal/solvers/acoustics` <br>
    and then execute the build script <br>
    `> ./build_acoustics.sh` <br>
3. To build the tests, run <br>
    `> tests/build_acoustics_tests.sh` <br>

**NOTE**:
* Do not use other versions of OCCA than the one from the root folder, since libParanumal is incompatible with newer versions.
* If you want to run on CPUs, remove the loading of the `CUDA` module inside `build_acoustics.sh` and rebuild.

## RUNNING THE CODE
* Run an example by (output written to `simulationSetups/output`)<br>
    `> simulationSetups/RUN_EXAMPLE.sh` <br>
* Run the tests and make sure all tests pass. Several examples including frequency independent and dependent cases can be found inside the test folder <br>
    `> tests/run_tests.sh` <br>    

**NOTE**: You might run out of resources ("out of memory messages) on DTU HPC if the task is not excecuted using the queue. <br>

## GENERATING DATA FOR TRAINING MACHINE LEARNING MODELS
Please see [instructions](https://github.com/dtu-act/libparanumal/tree/master/solvers/acoustics/simulationSetups/deeponet).

## SETTINGS FILE FORMATS

### Main settings file
The main settings file contains the parameters for simulations, such as source and source positions, material parameters, path to the mesh etc.

The application is executed by <br>
    `> mpirun -n 1 ./acousticsMain path/to/settings/file` <br>

where the argument is the path to the settings file containing the following properties:
```
[SIMULATION_ID]
cube_500hz_p4_5ppw_freq_indep

[MESH FILE]
../../meshes/tests/cube_500hz_p4_5ppw_freq_indep.msh

; [MESH FILE IC] # OPTIONAL: mesh determining source positions indexed with argument when calling main
; ../../meshes/tests/cube_500hz_p4_5ppw_freq_indep_src_pos.msh

[OUTPUT DIRECTORY]
simulationSetups/output

[POLYNOMIAL DEGREE]
4

[CURVILINEAR MESH] # 0: Off, 1: On. Mesh order MUST match [POLYNOMIAL DEGREE]
0

[TIME INTEGRATOR] # LSERK4 | EIRK4 | DOPRI5 (broken)
LSERK4

[FINAL TIME]
0.1

[RECEIVER] # OPTIONAL
simulationSetups/setupdata/receivers2_cube.dat

; [DT] # OPTIONAL: set time step resolution explicitly (CFL will be not be used)
; 9.296099356709925e-07 

[CFL] # Courant condition. Use CFL=1 for freq.-indep BC and CFL=0.2 for freq.-dep. BC. NOTE: overwritten then [DT] is set
1.0

[RHO] # Density of the medium
1.2

[C] # Speed of sound in the medium
343.0

[Z_IND] # Z value for frequency independent boundary condition
7400

[LRVECTFIT] # Generated using vectorfitDriverLR.m
simulationSetups/setupdata/LRDATA14.dat

[ERVECTFIT] # Generated using vectorfitDriverER.m
simulationSetups/setupdata/ERDATA14.dat

[FREQUENCY]
500 # [Hz]

[WRITE_WAVE_FIELD] # NONE | XDMF | H5 | H5Compact
NONE

[TEMPORAL_PPW_OUTPUT] # [WRITE_WAVE_FIELD] != NONE: temporal resolution of the output wave field
2

[MESH_RECTILINEAR_PPW] # [WRITE_WAVE_FIELD] != NONE: Resolution of the mesh for the initial condition
2

[SOURCE_TYPE] # GAUSSIAN | GRF (gaussian random fields)
GAUSSIAN

; [GRF_LENGTH_SCALE] # GRF only: length scale, bigger values lead to smoother functions
; 0.3

; [SXYZ] # OPTIONAL: Width of initial pulse (will be computed automatically if not defined)
; 0.4

[SX] # x coordinate of initial pulse. Ignored when: [GAUSSIAN_SOURCE_POS] = MESH or [SOURCE_TYPE] = GRF
0.5

[SY] # y coordinate of initial pulse. Ignored when: [GAUSSIAN_SOURCE_POS] = MESH or [SOURCE_TYPE] = GRF
0.5

[SZ] # z coordinate of initial pulse. Ignored when: [GAUSSIAN_SOURCE_POS] = MESH or [SOURCE_TYPE] = GRF
0.5
```

#### Overall settings
* `[MESH FILE]`: The mesh (see section "Mesh generation with Gmsh") for the simulation.
* `[RECEIVER]`: The IRs are written to disk corresponding to the source positions pointed to. Can be skipped if only the wave field output is needed.
* `[LRVECTFIT]` / `[ERVECTFIT]`: Miki's model is used for modeling frequency dependent boundaries and the fitted parameters are contained in a separate file with the path given here. The scripts for generating the coefficients are located inside `simulationSetups/vector_fitting_tools/`.

#### Settings for outputting all mesh pressures for all timestep
* `[WRITE_WAVE_FIELD]`: The full pressure field can be exported for each time step and output format is given here. `XDMF` is a compact format and can be visualized using a wide range of applications, such as ParaView. `H5Compact` is collabsing all timestep into a single tag in the H5 file for faster loading, but cannot be visualized directly. `H5` is similar to `H5Compact` but also writes the tetrahedral connection matrix.
* `[TEMPORAL_PPW_OUTPUT]`: The temporal sampling resolution can be set here, skipping some of the oversampled time instances.
* `[MESH FILE IC]`: Source positions are determined by the mesh nodes in this `.msh` created in Gmsh. The source index is determined by the last argument when invoking the application and multiple simulations can be triggered by using a task array: \verb|./acousticsMain setup_file $src_index|, where \verb|$src_index| is a bash variable triggered by a task array.
* If Gaussian Random Fields (GRFs) are to be used as intital source term, set `#define INCLUDE_GRF 1` in `acoustics.h` and link the armadillo library with `-larmadillo` inside the makefile.

### Receiver positions
The receiver file set in `[RECEIVER]` includes the receiver position `x,y,z` locations for each receiver. The first line indicates the number of receivers.
```
2
0.1 0.1 0.1
0.1 0.4 0.3
```

### Frequency dependent boundaries settings file
For frequency dependent boundaries, Miki's model is used and the coefficient can be fitted using the Matlab script `simulationSetups/vector_fitting_tools/`. The generated file has the format

```
14 10 2
0.037628829927
-0.075844970252
...
-3745.107750216643
-17416.897171651282
0.002138957511
-----
sigma = 47700
dmat = 0.050000000000
freqRange = [50,2000]
```

## MESH GENERATION WITH GMSH
The mesh file can be created using [Gmsh](http://gmsh.info). For more complicated geometries, Sketchup or similar can be used to create the geometry and then be imported to Gmsh. The steps for generating the mesh with appropriate element size and boundary conditions in Gmsh is shown below. The geometry is defined inside the `.geo` files and a few examples can be found in `meshes/geo/`. <br>

See also [Gmsh](http://gmsh.info) and [Gmsh video tutorial](https://www.youtube.com/watch?v=xL2LmDsDLYw).

### Geometries and materials
For simple geometries, it is easiest to modify the `.geo` file directly. For more complicated geometries, it can be useful to use Sketchup with the plugin for exporting to Gmsh. 

1. For setting boundary types, open a `.geo` file in Gmsh and choose `Modules->Geometry->Edit Script`. The text file looks like the following, where the dimension can easily be changed:
```
cl__1 = 1.0;
xdim = 4;
ydim = 2.7;
zdim = 3;

Point(1) = {0, 0, 0, cl__1};
Point(2) = {xdim, 0, 0, cl__1};
Point(3) = {xdim, ydim, 0, cl__1};
Point(4) = {0, ydim, 0, cl__1};
Point(5) = {0, 0, zdim, cl__1};
Point(6) = {xdim, 0, zdim, cl__1};
Point(7) = {xdim, ydim, zdim, cl__1};
Point(8) = {0, ydim, zdim, cl__1};
Line(9) = {1, 2};
Line Loop(22) = {16, 13, 14, 15};
Plane Surface(22) = {22};
...
Physical Surface("Frequency Independent",2) = {22, 24, 26, 28, 30, 32};
Physical Volume(10) = {34};
```
2. The field `Physical Surface("Frequency Independent",2)` is defining the type of boundaries for the corresponding plane surfaces, where the number in the tuple is defining the boundary type (the string is just for convenience). The boundary types are
    * 1: Perfectly reflecting (Neumann) boundaries.
    * 2: Frequency independent impedance boundaries.
    * 3: Local-reacting frequency dependent impedance boundaries.
    * 4: Extended-reaction frequency dependent impedance boundaries.

### Element size
The input to libParanumal is a mesh discretized in terms of elements. Depending on the polynomial order chosen in libParanumal, the mesh resolution should be chosen accordingly.

<bf>Example:</bf>
Assume that we have maximun frequency $f_\text{max} = 1000 \text{Hz}$, speed of sound $c = 343$ m/s, points per wavelength $\text{ppw} = 5$ and polynomial order $P = 4$. Then the element resolution $\Delta x_{\text{elem}}$ is calculated as $\Delta x_{\text{elem}} = \frac{c}{f_{\text{max}}\times \text{ppw}} \times P = 0.2744 \text{ m}$.

The element size can be set as follows:

1. In the top menu, select `Tools->Options` and in the popup window choose `Geometry->General` and set `Global model scaling` to 1.
2. in the popup window choose `Mesh->General` and set *both* fields in `Min/max element size` to $\Delta x$.
3. Close the popup window.

### Mesh generation
Having set the boundary materials and element size, the geometry can now be meshed.

1. Choose `Modules->Mesh->3D` and then optimize the mesh by choosing `Optimize 3D (Netgen)`. Note, that if the element resolution is changes, 1D, 2D and 3D needs to be meshed for the changes to take effect.
2. Save the file by choosing `File->Export` and save the file with `.msh` extension. Choose `Version 2 (ASCII)` and click ok.

## GAUSSIAN RANDOM FIELDS
Support for GRFs as initial conditions is not enabled as default. As explained earlier, set the flag `#define INCLUDE_GRF 1` in `acoustics.h` and link the armadillo library with `-larmadillo` inside the makefile. Also, the static library [https://github.com/bigladder/btwxt](https://github.com/bigladder/btwxt) used for interpolating from static grids to Gaussian quadrature points should be linked by assigning `-L$(LIBSDIR) -lbtwxt` to the `LIBS` environment variable . The static library `libbtwxt.a` should be located inside `libparanumal/libs/` and the header files located inside `librapanumal/include/btwxt`. To build the library, do the following:

1. Clone the code from git into e.g. `libparanumal/include/` <br>
    `> git clone https://github.com/bigladder/btwxt`
2. Setup with cmake <br>
    `> mkdir build/` <br>
    `> cd build` <br>
    `> cmake --install ../src ` <br>
3. Build the static library from the generated makefile <br>
    `> cd src/` <br>
    `> make` <br>
4. Move the library to the libs folder <br>
    `> mv libbtwxt.a /path/to/libparanumal/libs/` <br>
5. Move the header files into the root of the `btwxt` folder <br>
    `> cd /path/to/libparanumal/include` <br>
    `> mv btwxt/src/btwxt.h .` <br>
    `> mv btwxt/src/error.h .` <br>
    `> mv btwxt/src/griddeddata.h .` <br>
    `> mv btwxt/src/gridpoint.h .` <br>
 6. The remaining files can be cleaned up. <br>

## USING DTU HPC
* `> ssh username@login1.gbar.dtu.dk`   # login
* `> voltash`                           # switch to GPU cluster
* `> bsub < <the_script>.sh`            # add to the queue system (see `simulationSetups/run_cube.sh`)
* `> bstat`                             # job status
* `> bkill <id>`                        # kill job
* `> ./the_scripts > logfile.txt`       # run the script and pipe the output to a log file

## DEVELOPMENT USING VSCODE - REMOTE IDE
Setting up an environment for effective and easy development is essential. With VS Code, it is possible to work on your local computer, but running the code [remotely](https://code.visualstudio.com/docs/remote/remote-overview) using the [Remote Development](https://marketplace.visualstudio.com/items?itemName=ms-vscode-remote.vscode-remote-extensionpack) extension.

[Here](https://docs.ccv.brown.edu/oscar/connecting-to-oscar/remote-ide) you can find good documentation how to setup such an environment (adjust to your environment).

### DTU remote setup
The `~/.ssh/config` should look like this for connecting to DTU HPC system:

```
Host gbar-remote-ide
User <your username>
IdentityFile ~/.ssh/gbar
Hostname nilogin.gbar.dtu.dk
```

In VSCode, select `Remote-SSH: Connect to Host...` from `View->Command Palette` and select `gbar-remote-ide`. You might need to change the timeout limit when entering the password.

## Some references
* [VPN](https://www.inside.dtu.dk/en/medarbejder/it-og-telefoni/wifi-og-fjernadgang/vpn-cisco-anyconnect)
* [HPC](https://www.hpc.dtu.dk/)
* [HPC GPU parameters](https://www.hpc.dtu.dk/?page_id=2759)
* [LSF job submission system](http://www.cc.dtu.dk/?page_id=1416)
* [HPC examples](https://www.hpc.dtu.dk/?page_id=2021)

---
# GENERAL
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


