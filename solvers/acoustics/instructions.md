INSTALLATION
1) git clone https://github.com/dtu-act/libparanumal
2) copy OCCA 1.0 into libparanumal folder (do not use version from git, ask the acoustic group at DTU)
3) build OCCA and acousticsMain by running 'libparanumal/solvers/acoustics/build_acoustics.sh'
    - NOTE: to exploit the GPU, you should build in an environment with access to GPUs (see section about DTU HPC)
4) run RUNGPU.bsub for a test example

HPC system
1) Login using ssh username@login1.gbar.dtu.dk 
2) change to GPU capabilities by typing voltash
3) you can run interactively by executing ./RUNGPU.bsub
4) when running long-running tasks, use the queue
    - bsub < <the_script>.bsub

GMSH
- use Gmsh (http://gmsh.info) to create meshes (see Massively Parallel Nodal Discontinous Galerkin Finite Element Method Simulator for 3D Room Acoustics for details)
- Tutorial: https://www.youtube.com/watch?v=xL2LmDsDLYw


OTHER NOTES

VPN: https://www.inside.dtu.dk/en/medarbejder/it-og-telefoni/wifi-og-fjernadgang/vpn-cisco-anyconnect
HPC: https://www.hpc.dtu.dk/
HPC GPU parameters: https://www.hpc.dtu.dk/?page_id=2759
LSF job submission system: http://www.cc.dtu.dk/?page_id=1416
HPC examples: https://www.hpc.dtu.dk/?page_id=2021

Gmsh: http://gmsh.info
Tutorial: https://www.youtube.com/watch?v=xL2LmDsDLYw

DGFEM simulator: https://github.com/dtu-act/libparanumal

Download and build occa: use version 1.0 of the OCCA lib (ask nibor@elektro.dtu.dk for the correct version of the code)

Terminal commands
> ssh username@login1.gbar.dtu.dk   # login
> bsub < <the_script>.bsub          # run on the queue system
> bstat                             # status
> bkill <id>                        # kill


# RUN INTERACTIVELY
> voltash
> ./the_scripts > logfile.txt