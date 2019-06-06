! This file is automatically created by input_parameter.py
module input_parameter
    use salmon_file, only: open_filehandle
    implicit none

    character(64) :: directory
    character(64) :: sysname
    integer :: num_kgrid(1:3)
    real(8) :: al_vec3(1:3)
    real(8) :: al_vec2(1:3)
    real(8) :: al_vec1(1:3)
    integer :: nstate
    integer :: nelec

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
        & al_vec3, &
        & al_vec2, &
        & al_vec1, &
        & nstate, &
        & nelec

        directory = ''
        sysname = './'
        num_kgrid = 0
        al_vec3 = 0.0
        al_vec2 = 0.0
        al_vec1 = 0.0
        nstate = 0
        nelec = 0

        fh = open_filehandle('.namelist.tmp')
        do while (.true.)
            read(*, '(a)', iostat=ret) tmp
            if (ret < 0) exit ! End of file
            tmp = adjustl(tmp)
            if (tmp(1:1) .ne. '!') write(fh, '(a)') trim(tmp)
        end do
        rewind(fh); read(fh, nml=control, iostat=ret)
        rewind(fh); read(fh, nml=kgrid, iostat=ret)
        rewind(fh); read(fh, nml=system, iostat=ret)

        close(fh)

        write(*, '("# directory =",99(1x,a))') directory
        write(*, '("# sysname =",99(1x,a))') sysname
        write(*, '("# num_kgrid =",99(1x,i7))') num_kgrid
        write(*, '("# al_vec3 =",99(1x,f7.3))') al_vec3
        write(*, '("# al_vec2 =",99(1x,f7.3))') al_vec2
        write(*, '("# al_vec1 =",99(1x,f7.3))') al_vec1
        write(*, '("# nstate =",99(1x,i7))') nstate
        write(*, '("# nelec =",99(1x,i7))') nelec
    end subroutine read_input
end module input_parameter

