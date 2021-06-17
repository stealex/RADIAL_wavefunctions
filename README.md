# RADIAL_wavefunctions

This repository contains the famous "RADIAL" package of Salvat et al. Additionally, there is a `C++` interface to this package for the computation of bound and scattering wave-functions.

## Compile the whole thing
There is a list of prequisites:
- `Linux` distro with `gcc` (version >= 7)
- `gfortran` (seomthing recent)
- `GNUmake`

You should be able to build the program using:
```
make
```
in the root directory of the project.
In case of problems contact me.

## Running
The functonality, although advanced, is not completed. You need to create a folder called `Nuclei` at the same level with the root directory of this project. Inside `Nuclei` you must create a folder for each nucleus. The name of this folder must be the "A" of the nucleus followed by it's name (e.q. `76Ge` for Germanium-76). For example:
```
cd project_source
mkdir ../Nuclei/76Ge
```

Having this structure set-up, take a look into `wavefunctions.conf` file. It is fairly well documented inline. Modify it to suit your needs and run:
```
./bin/main wavefunctions.conf
```
This will create a set of files under `../Nucleus/76Ge` that can be imported in `Mathematica` using the `DBDPhaseSpace` package.
