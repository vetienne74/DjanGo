
<img src=".//misc/fileForReadme/DjanGoLogo.png"  alt="./misc/fileForReadme/DjanGoLogo.png">

# Welcome to DjanGo

**Version 1.0**

Contact: Vincent Etienne / Email: vetienne@rocketmail.com

# Overview

## Description

DjanGo is a versatile C++ tool for modelling wave propagation
with finite-difference & finite-element numerical methods. Initially developped for geophysical applications (seismic wave propagation), it includes now capabilities towards sound and music modelling. However, it should be possible to adapt DjanGo to model any kind of waves (mechanical, electromagnetic and so on). Anyone who is interested in using DjanGo or contributing to its development is invited to contact me.  

Features of DjanGo
- Wave propagation formulation
  - Acoustic and elastic wave formulation
  - 1st and 2nd order equation formulation
  - 1D and 2D geometries
- Boundary conditions
  - Free surface, rigid surface
  - Absorbing boundaries: CPML and Sponge
- Parallelism
  - Multi threading (OpenMP)
  - Subdomain decomposition (MPI)

FDM characteristics
- Staggered grids
- 2nd to 16th order in space
- Second-order leap frog time integration scheme

FEM characteristics
- Continuous Galerkin Spectral Element method
- Discontinuous Galerkin Spectral Element method with full p-adaptivity
- Segment (1D) and quadrangle (2D) elements
- GLL nodes (arbitrarily number)
- Second-order leap frog time integration scheme

Intensive test cases have been designed to validate accuracy and convergence.

## Versions

Version      | Description | Release date
------------ | ----------- | ------------
v1.0         |  **Initial version** | December, 2023


# Environment set-up

## Basic requirements

To use DjanGo on your computer
* C++ compiler with OpenMP support
* MPI library

:bell: **The Linux environment is my favorite and all scripts in this projects are shell scripts that run on Linux.**

DjanGo runs on CPU, development for GPU acceleration is on the way.

To display figures and hear sounds with the provided material, Matlab is necessary.

## Environment script (mandatory)

In order to compile and run DjanGo, you need to source one of the files in the directory `./env`

`cd ./env`

Example to set up the environment for DjanGo with GCC compiler:

`source ./source setEnvNeptuneGcc.sh`

[Display command output](misc/fileForReadme/setEnvDjanGo.txt)

:bell: **For a new system, you would need to create a file for your system** (take example from one of the existing files)

# Compilation

## Makefile

To compile DjanGo, go to `./build`, and use the command `make`

[Display command output](misc/fileForReadme/make.txt)

Executable can be found in `./bin/django`

:bell: If DjanGo environment has not been set (see [Environment script (mandatory)](#environment-script-mandatory)), compilation will abort.

By default, DjanGo is compiled in single in precision

To compile in double precision: `make precision=double`

Some additional tools are provided and need also to be compiled. Go to `./misc/tool/build`, and use the command `make`

[Display command output](misc/fileForReadme/makeTools.txt)

# Validation

## Validation tests

To check that DhanGo has correctly been built and works fine, go to `./validation` and launch

`sh runValidationTests.sh`

[Display command output](misc/fileForReadme/runValidationTests.txt)

This script runs a set a light test cases and should complete within few minutes (even on a laptop).

You should get in the ouptput report (displayed on the terminal)

* All tests marked as PASSED (437 tests passed)
* No test marked as FAILED

Check the summary at the end of report to have a quick look on this.

:bell: These tests are intended for validation purpose only, they do not allow for performance measurements.

## Validated hardware, operating systems and compilers

DjanGo has been successfully tested on the hardware, operating systems and compilers listed below.

# Execution

## Usage

DjanGo can be launched from a terminal with all configuration parameters within a single xml file.

**To get help on the parameters**

`./bin/django -h`

[Display command output](misc/fileForReadme/commandLineParam.txt)


# Directory description

appli: applications

bin: binaries

build: compilation directory (Makefile is there)

env: environment scripts

include: all include *.h files

misc: miscellaneous dirs, including tools

script: validation scripts

src: all source *.cpp files

validation: various validation tests organised in categories (see 0readme in subfolders)
