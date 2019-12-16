# 1D-DG Project

## To compile :
"cmake .." in build then "make"

## To clean :
- clean outputs : "make data-clean" in build

- clean .o and .mod : "make clean" in build

- clean everything : "make cleanall" in build

## To execute a forward run :

"./FORWARD" or "./FORWARD parameter-file.par" from build if there is a parameter file to specify

To execute a fwi run :

"./FWI" or "./FWI parameter-file.par" from build if there is a parameter file to specify

NB : a default parameter file is given in the root of the project and is called "param.par"

## Dependencies :
To compile and run the code you will need :

- gfortran
- cmake
- BLAS and LAPACK libraries

For animations you can have :

- gnuplot

or

- python3 (faster)

- conda to set the python environment properly

To set python properly you have to :
- create the environment from the envornment.yml in the animation_script directory : `conda env create -f environment.yml`
- set the environment : `conda activate 1D_DG_env`
- set gnuplot booleen in param file to false