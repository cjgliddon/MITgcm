# HPIMM: the High Pressure Icy Moon Module of MITgcm

HPIMM is a package for MITgcm that enables it to realistically simulate the geophysical fluid dynamics of high-pressure icy ocean worlds. 

## How to Use

The main function of the `buildmitgcm.sh` script is to create an *experiment directory*, which contains a new compiled MITgcm executable (which runs the model) with user-specified code modifications, as well as a subdirectory of input files for defining model parameters and sample scripts for analysis. The basic command line structure for compiling is:
```
$ .\buildmitgcm.sh <expname> <gridsize> <ncores> compilebuild
```

`<expname>` is just the name of the directory to be created. The `<gridsize>` argument specifies a particular size file, `SIZE.h`, which will be used in model compilation. `<ncores>` is the number of cores to be used in running the model. **Important:** The script assumes a standard formatting for the `<gridsize>` argument and size file names. `<gridsize>` should have the format `{string}(nz,ny,nx)`, where `{string}` is any string identifier. The corresponding size file has the name `SIZE.h.{string}x{nx}y{ny}.{ncores}p` -- note that because the tiling of the MITgcm grid depends on the number of processors, we must specify different size files for runs with different numbers of processors (even if the total number of grid cells is the same).

The code modifications which the user wants to use in compiling MITgcm are stored in a directory `code_now`. This directory includes a `packages.conf` file which specifies which packages to include in the compilation; the `hpimm` and `shelfice` packages should be included in this file. Additionally, project-specific input files are copied from the folder `input_now`. The most important file in this folder is the Matlab script `gendata.m`, which when run automatically generates the binary files necessary to initialize a model run for a given set of planetary parameters and which modifies the input namelists accordingly.