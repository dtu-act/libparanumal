#!/bin/bash

./boltzmannMainQuad3D ../../meshes/cubed_grid_layers_medium.msh 3 1 0 0 1 1
./boltzmannMainQuad3D ../../meshes/cubed_grid_layers_medium.msh 3 0 1 0 1 1
./boltzmannMainQuad3D ../../meshes/cubed_grid_layers_medium.msh 3 1 1 0 1 1
./boltzmannMainQuad3D ../../meshes/cubed_grid_layers_medium.msh 3 1 1 0 1 2
./boltzmannMainQuad3D ../../meshes/cubed_grid_layers_medium.msh 3 1 1 1 1 2
./boltzmannMainQuad3D ../../meshes/cubed_grid_layers_medium.msh 3 1 1 1 2 2