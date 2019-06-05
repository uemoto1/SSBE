input_parameter = {
    'control': {
        'sysname': ['character(64)', []],
        'directory': ['character(64)', []],
    },
    'system': {
        'al_vec1': ['real(8)', [3]],
        'al_vec2': ['real(8)', [3]],
        'al_vec3': ['real(8)', [3]],
        'nstate': ['integer', []],
        'nelec': ['integer', []],
    },
    'kgrid': {
        'num_kgrid': ['integer', [3]],
    },
}

template = """
module input_parameter
    implicit none

{DEF_VARIABLE}

contains

    subroutine read_input()
        implicit none
        integer :: ret, fh
        character(256) :: tmp

{DEF_NAMELIST}

        fh = open_filehandle('.namelist.tmp')
        do while (.true.)
            read(*, '(a)', iostat=ret) tmp
            if (ret < 0) exit ! End of File

            tmp = trim(adjustl(tmp))
            if (tmp(1:1) <> '!') write(fh, '(a)') tmp
        end do
{READ_NAMELIST}
        close(fh)
    end subroutine read_input
end module input_parameter
"""

def indent(text, depth=0):
    return "\n".join([" " * depth + line for line in text.splitlines()])

def f90dim(dimlist):
    if 0 < len(dimlist):
        return "(" + ", ".join([str(i) for i in dimlist]) + ")"
    else:
        return ""


def_variable = ""
read_namelist = ""
def_namelist = ""

for group in sorted(input_parameter.keys()):
    def_namelist += "namelist/{GROUP}/ &\n& {VARLIST}\n".format(
        GROUP=group,
        VARLIST=", &\n& ".join(sorted(input_parameter[group].keys()))    
    )
    read_namelist += "rewind(fh); read(fh, nml={GROUP})\n".format(
        GROUP=group
    )
    for name in sorted(input_parameter[group].keys()):
        f90type, dimlist = input_parameter[group][name]
        def_variable += "{F90TYPE} :: {NAME}{F90DIM}\n".format(
            F90TYPE=f90type, NAME=name, F90DIM=f90dim(dimlist)    
        )
            

print(template.format(
    DEF_VARIABLE=indent(def_variable, 4),
    DEF_NAMELIST=indent(def_namelist, 8),
    READ_NAMELIST=indent(read_namelist, 8)
))