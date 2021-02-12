## INSTALLATION
1. If you want to access DTUs HPC system, login using command (see section about HPC system as well) <br>
    `> ssh username@login1.gbar.dtu.dk`
2. Create a folder `libparanumal` **in the home user directory** (i.e. `~/', otherwise you should modify the build script) and clone the code from git into the folder using the command <br>
    `> git clone https://github.com/dtu-act/libparanumal`
3. Copy OCCA 1.0 into libparanumal folder. **IMPORTANT**: do not use version from git, ask the acoustic group at DTU (OCCA is a third-party library used for compiling code to various platforms, such as CUDA/GPU)
4. Build OCCA and libParanumal. **IMPORTANT**: to exploit the GPU, you should build in an environment with access to GPUs. On DTUs systems, do (see also section below about DTU HPC) <br>
    `> voltash`<br>
    Enter the acoustics folder <br>
    `> cd ~/libparanumal/solvers/acoustics` <br>
    and then execute the build script <br>
    `> ./build_acoustics.sh` <br>
    
5. Now libParanumal should be compiled, For testing the setup, run the following command from location `libparanumal/solvers/acoustics` (output is written into `libparanumal/solvers/acoustics/data`) <br>
    `> ./RUNGPU.bsub` <br>    
6. Several examples including frequency dependent and independent cases can be found inside `libparanumal/solvers/acoustics/tests/`

## HPC system
1. Login using <br>
    `> ssh username@login1.gbar.dtu.dk`
2. change to GPU capabilities by typing <br> 
    `> voltash`
3. you can run interactively by executing (after having compiled OCCA and libParanumal)<br> 
    `> ./RUNGPU.bsub`
4. when running long-running tasks, use the queue <br> 
    `bsub < <the_script>.bsub`

## GMSH
* use Gmsh (http://gmsh.info) to create meshes (see 'Massively Parallel Nodal Discontinous Galerkin Finite Element Method Simulator for 3D Room Acoustics' thesis for details)
* Tutorial: https://www.youtube.com/watch?v=xL2LmDsDLYw

## References
* VPN: https://www.inside.dtu.dk/en/medarbejder/it-og-telefoni/wifi-og-fjernadgang/vpn-cisco-anyconnect
* HPC: https://www.hpc.dtu.dk/
* HPC GPU parameters: https://www.hpc.dtu.dk/?page_id=2759
* LSF job submission system: http://www.cc.dtu.dk/?page_id=1416
* HPC examples: https://www.hpc.dtu.dk/?page_id=2021
* Gmsh: http://gmsh.info
* Tutorial: https://www.youtube.com/watch?v=xL2LmDsDLYw
* DGFEM simulator: https://github.com/dtu-act/libparanumal

## OCCA
Download and build OCCA v1.0 (ask the DTU acoustics group for the correct version: nibor@elektro.dtu.dk or hsllo@elektro.dtu.dk)

## Useful UNIX commands
* `> ssh username@login1.gbar.dtu.dk`   # login
* `> voltash`                           # switch to GPU cluster
* `> bsub < <the_script>.bsub`          # add to the queue system
* `> bstat`                             # job status
* `> bkill <id>`                        # kill job
* `> ./the_scripts > logfile.txt`       # run the script and pipe the output to a log file