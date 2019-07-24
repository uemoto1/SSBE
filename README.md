# SSBE
Experimental Implementation of SALMON's Semiconductor Bloch Equation

## Structure

### SALMON SBE (SSBE)
The Houston SBE solver for SALMON program:

```
salmon_sbe
├── Makefile                GNU makefile for build
├── Makefile.inc.gnu        (Setting for gfortran+OpenMPI+BLAS/LAPACK)
├── Makefile.inc.intel-knl  (Setting for ifort+IntelMPI+IntelMKL)
├── input_parameter.f90     Input file parser
├── input_parameter.py      (...and its generator script)
├── pulse.f90               Pulse waveform generator
├── sbe_solver.f90          Semiconductor Bloch Solver
├── main.f90                Main routine
│
├── common                  SALMON unified libraries ....
│   ├── pack_unpack.f90
│   └── structures.f90
├── io
│   └── salmon_file.f90
├── math
│   └── salmon_math.f90
├── misc
│   └── unusedvar.f90
├── parallel
│   ├── salmon_communication.f90
│   └── salmon_parallel.f90
│
└── SSBE                    Binary
```

### TDSE SBE (SSBE)
The Bloch SBE solver for TDSE program:
```
tdse_sbe
├── ...
```

### Build

On source code directory (`salmon_sbe/` or `tdse_sbe/`),
rename the platform-dependent files (`.gnu` or `.intel-knl`) to `Makefile.inc`, and execute `make` command: 
```
cp Makefile.inc.YOUR_PLATFORM Makefile.inc
make
```
