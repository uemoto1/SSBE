# SSBE
Example Implementation of Semiconductor Bloch Equation for SALMON program

## Structure

```
salmon_sbe
├── Makefile                GNU makefile for build
├── Makefile.inc.gnu        (Configulation for gfortran/openmpi)
├── Makefile.inc.intel-knl  (Configulation for gfortran/openmpi)
├── input_parameter.f90     Input file parser
├── input_parameter.py      (...and its generator script)
├── pulse.f90               Pulse waveform generator
├── sbe_solver.f90          Semiconductor Bloch Solver
├── main.f90                Mainroutine
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

