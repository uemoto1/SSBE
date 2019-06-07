# SSBE
Example Implementation of Semiconductor Bloch Equation for SALMON program

## Structure

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

