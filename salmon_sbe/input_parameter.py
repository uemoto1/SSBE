#!/usr/bin/env python
import os

input_parameter = {
    'control': {
        'sysname':          ['character(64)',   [],     "untitled"],
        'directory':        ['character(64)',   [],     "./"],
        'read_sbe_gs_bin':  ['character(64)',   [],     "n"]
    },
    'system': {
        'al_vec1':      ['real(8)',         [3],    0.],
        'al_vec2':      ['real(8)',         [3],    0.],
        'al_vec3':      ['real(8)',         [3],    0.],
        'nstate':       ['integer',         [],     0],
        'nelec':        ['integer',         [],     0],
    },
    'kgrid': {
        'num_kgrid_gs':     ['integer',     [3],    0],
        'num_kgrid_sbe':    ['integer',     [3],    0],
    },
    'analysis': {
        'e_min_dielec':     ['real(8)',     [],     0.],
        'e_max_dielec':     ['real(8)',     [],     1.],
        'n_dielec':         ['integer',     [],     1000],
        'gamma_dielec':     ['real(8)',     [],     0.005],
        'out_rt_step' :     ['integer',     [],     10],
        'out_dielec' :      ['character(1)',[],     "y"],
    },
    'tgrid': {
        'dt':               ['real(8)',     [],     0.08],
        'nt':               ['integer',     [],     5000],
    },
    'emfield': {
        'epdir_re1':         ['real(8)',    [3],    [0.0, 0.0, 1.0]],
        'epdir_im1':         ['real(8)',    [3],    [0.0, 0.0, 0.0]],
        'rlaser_int_wcm2_1': ['real(8)',    [],     1e10],
        'omega1':            ['real(8)',    [],     0.05698],
        'pulse_tw1':         ['real(8)',    [],     413.5],
        'phi_cep1':          ['real(8)',    [],     0.00],
    },
}












template = """! This file is automatically created by input_parameter.py
module input_parameter
    use salmon_file, only: open_filehandle
    implicit none

{DEF_VARIABLE}

contains

    subroutine read_input()
        implicit none
        integer :: ret, fh
        character(256) :: tmp

{DEF_NAMELIST}

{VARIABLE_DEFAULT}

        fh = open_filehandle('.namelist.tmp')
        do while (.true.)
            read(*, '(a)', iostat=ret) tmp
            if (ret < 0) exit ! End of file
            tmp = adjustl(tmp)
            if (tmp(1:1) .ne. '!') write(fh, '(a)') trim(tmp)
        end do
{READ_NAMELIST}

        close(fh)

{VAR_DUMP}
    end subroutine read_input
end module input_parameter
"""





def indent(text, depth=0):
    return "\n".join([" " * depth + line for line in text.splitlines()])

def f90dim(dimlist):
    if 0 < len(dimlist):
        return "(" + ", ".join(["1:%d" % i for i in dimlist]) + ")"
    else:
        return ""

def f90val(value):
    return repr(value).replace("[", "(/").replace("]", "/)")

def f90fmt(f90type):
    if "character" in f90type:
        return "a"
    elif "integer" in f90type:
        return "i7"
    elif "real" in f90type:
        return "es12.4e3"
    else:
        raise TypeError


def_variable = ""
read_namelist = ""
def_namelist = ""
default_variable = ""
var_dump = ""

for group in input_parameter:
    def_namelist += "namelist/{GROUP}/ &\n& {VARLIST}\n".format(
        GROUP=group,
        VARLIST=", &\n& ".join(input_parameter[group].keys()) 
    )
    read_namelist += "rewind(fh); read(fh, nml={GROUP}, iostat=ret)\n".format(
        GROUP=group
    )

    for name in input_parameter[group]:
        f90type, dimlist, default = input_parameter[group][name]
        def_variable += "{F90TYPE} :: {NAME}{F90DIM}\n".format(
            F90TYPE=f90type, NAME=name, F90DIM=f90dim(dimlist)    
        )
        default_variable += "{NAME} = {VALUE}\n".format(
            NAME=name, VALUE=f90val(default)
        )
        var_dump += """write(*, '("# {NAME} =",99(1x,{F90FMT}))') {NAME}\n""".format(
            NAME=name, F90FMT=f90fmt(f90type)
        )
        
            

title, ext = os.path.splitext(__file__)
f90file = "%s.f90" % title
inpfile = "%s.inp" % title

with open(f90file, "w") as fh:
    fh.write(template.format(
        DEF_VARIABLE=indent(def_variable, 4),
        DEF_NAMELIST=indent(def_namelist, 8),
        READ_NAMELIST=indent(read_namelist, 8),
        VARIABLE_DEFAULT=indent(default_variable, 8),
        VAR_DUMP=indent(var_dump, 8)
    ))

with open(inpfile, "w") as fh:
    for group in input_parameter:
        fh.write("&%s\n" % group)
        for name in input_parameter[group]:
            f90type, dimlist, default = input_parameter[group][name]
            fh.write("\t%s = %s\n" % (
                name, 
                repr(default).replace("[","").replace("]","")
            ))
        fh.write("/\n\n")

1