# beam-cohesive-fracture

A secant-based solution to beam elements (with 6 DOFs) connected to a specialized continuum element with a crack down the middle in the local x-direction. The crack is modeled with the PPR cohesive zone model. This version of the code allows the user to specify an internal temperature based on some distribution.

## Prerequisites

To build the code, you will need to have the following installed:

- [GNU Make](https://www.gnu.org/software/make/)
- [GNU Fortran](https://gcc.gnu.org/fortran/)

## Building

To build the code, simply run `make` in the root directory of the project. This will create an executable called `secant-base-3.out`.

## Running

To run the code, simply run `./secant-base-3.out input.inp` in the root directory of the project. This will create the following list of files with output in the root directory of the project:

```bash
input.inp_boundary-disp.txt
input.inp_boundary-force.txt
input.inp_coh-stress.txt
input.inp_extmesh.txt
input.inp_internal-disp.txt
input.inp_internal-force.txt
input.inp_intmesh.txt
input.inp_separation.txt
```
