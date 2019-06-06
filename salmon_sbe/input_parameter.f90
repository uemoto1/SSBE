! This file is automaticall created by input_parameter.py
module input_parameter
    use salmon_file, only: open_filehandle
    implicit none

    character(64) :: directory
    character(64) :: sysname
    integer :: num_kgrid(3)
    real(8) :: al_vec1(3)
    real(8) :: al_vec2(3)
    real(8) :: al_vec3(3)
    integer :: nelec
    integer :: nstate

contains

    subroutine read_input()
        implicit none
        integer :: ret, fh
        character(256) :: tmp

        namelist/control/ &
        & directory, &
        & sysname
        namelist/kgrid/ &
        & num_kgrid
        namelist/system/ &
        & al_vec1, &
        & al_vec2, &
        & al_vec3, &
        & nelec, &
        & nstate

        fh = open_filehandle('.namelist.tmp')
        do while (.true.)
            read(*, '(a)', iostat=ret) tmp
            if (ret < 0) exit ! End of file
            tmp = adjustl(tmp)
            if (tmp(1:1) .ne. '!') write(fh, '(a)') trim(tmp)
        end do
        rewind(fh); read(fh, nml=control)
        rewind(fh); read(fh, nml=kgrid)
        rewind(fh); read(fh, nml=system)
        close(fh)
    end subroutine read_input
end module input_parameter

