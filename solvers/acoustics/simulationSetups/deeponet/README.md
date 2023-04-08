# DATA GENERATION FOR DEEPONET
This text explains how to generate the data used for training and evaluating the DeepONet from the paper "Sound propagation in realistic interactive 3D scenes with parameterized sources using deep neural operators" by Nikolas Borrel-Jensen, Somdatta Goswami, Allan P. Engsig-Karup, George Em Karniadakis, and Cheol-Ho Jeong (PNAS).

The Discontinuous Gallerkin Finite Element method (DG-FEM) method is used for generating the data and should be installed first by following the [instructions](https://github.com/dtu-act/libparanumal). Information about the setting files can also be found here.

The data used for the results from the paper can be provided by contacting Cheol-Ho Jeong <chje@dtu.dk> or Finn T. Agerkvist <ftag@dtu.dk>.

## OVERVIEW
Four geometries are considered:

* Cube 2m x 2m x 2m
* Furnished room
* L-shape room
* Dome

For each of the geometries, data are generated for training, validation and testing:

* Training data: generated with a spatial resolution of 6 points per wavelength for source positions distributed with 1/5 of a wavelength inside a subdomain of the geometry.
* Validation data: generated with a spatial resolution of 5 points per wavelength for source positions distributed with roughly a full wavelength inside the same subdomain as the traininng data.
* Test data: generated with a spatial resolution of 5 points per wavelength for five source positions inside the same subdomain as the traininng data.

For training and validation data, the source positions are defined by a GMESH mesh, where the vertices are determining the positions.

## TRAINING AND VALIDATION DATA
The setup scripts for generating the training and validation data are located inside the folder `setup_train_val_srcpos_mesh/` and the simulation scripts for running the executable with the appropriate setup scripts are located in the `deeponet` root folder

* run_cube2x2x2_train and run_cube2x2x2_val
* run_furnished_train and run_furnished_val
* run_Lshape_train and run_Lshape_val
* run_dome_train and run_dome_val

These scripts are using a job array passing an index argument from 1 to the number of source position to the executable telling the application with vertix index from the mesh file defined for the parameter `[MESH FILE IC]`.

## TEST DATA
The setup scripts for generating the test data are located inside the folder `setup_test_srspos5/` and the simulation scripts for running the executable with the appropriate setup scripts are located in the `deeponet/` root folder

* run_cube2x2x2_test
* run_furnished_test
* run_Lshape_test
* run_dome_test

These scripts are using a job array passing the path to one of the five setup scripts corresponding to one of the five source positions.

## INPUT FILES
### Geometries
The mesh geometry is defined in the setup file by the tag `[MESH FILE]` and the mesh defining the source positions is defined by the tag `[MESH FILE IC]`. The mesh files used for the simulations can be found in the PNAS Supplementary Information.

## MORE DETAILS

### Frequency independent and dependent materials
Frequency independent impedances is defined by the tag `[Z_IND]` in the setup script, whereas the coefficients for fitting a frequency-dependent material (using Miki's model) are referenced by tag `[LRVECTFIT]` and defined inside `deeponet/freq_dep_lr.dat`. Note that the material type is defined for each surface of the geometry defined in GMESH `.geo` file.

### Resolution
The spatial resolution is determined by the mesh referenced in `[MESH FILE]`. The input to the libParanumal DG-FEM solver is a mesh discretized in terms of elements and hence the actual spatial discretization resolution is given by the element resolution and the polynomial order determined by the tag `[POLYNOMIAL DEGREE]` in the setup file. Hence, depending on the polynomial order chosen in libParanumal, the mesh resolution should be chosen accordingly.

<bf>Example:</bf>
Assume that we have maximun frequency $f_\text{max} = 1000 \text{Hz}$, speed of sound $c = 343$ m/s, points per wavelength $\text{ppw} = 5$ and polynomial order $P = 4$. Then the element resolution $\Delta x_{\text{elem}}$ is calculated as $\Delta x_{\text{elem}} = \frac{c}{f_{\text{max}}\times \text{ppw}} \times P = 0.2744 \text{ m}$.

### Ensuring corresponding temporal resolutions
We need to ensure the same temporal resolutions between the training, validation and test datasets, since we are not interested in temporal interpolation for the trained DeepONet model (i.e. we are overfitting in the temporal dimension). The temporal resolution is determined from the spatial resolution and the Courant condition (CFL) and therefore, special care is required when using different grid resolutions for training and validation/test data.

The temporal resolution $\Delta t$ can be set explicitly in the setup file by the `[DT]` tag and has been set in the validation and testing data to be equal to the $\Delta t$ used for the training data. Stability is ensured, since the training data has higher spatial resolution (6 points per wavelenght) compared to validation and test data (5 points per wavelength).

Moreover, since sampling at the Nyquist limit is enought for reconstructing the discrete impulse response, the tag `[TEMPORAL_PPW_OUTPUT]` has been set to 2 points per wavelength pruning the data written to disk to only include time slices of approximately this resolution.