! This file is automatically created by input_parameter.py
module input_parameter
    use salmon_file, only: open_filehandle
    implicit none

    character(64) :: directory
    character(64) :: sysname
    character(64) :: read_sbe_gs_bin
    real(8) :: dt
    integer :: nt
    real(8) :: epdir_re1(1:3)
    real(8) :: phi_cep1
    real(8) :: rlaser_int_wcm2_1
    real(8) :: pulse_tw1
    real(8) :: omega1
    real(8) :: epdir_im1(1:3)
    real(8) :: al_vec3(1:3)
    real(8) :: al_vec2(1:3)
    real(8) :: al_vec1(1:3)
    integer :: nstate
    integer :: nelec
    integer :: out_rt_step
    integer :: n_dielec
    character(1) :: out_dielec
    real(8) :: e_max_dielec
    real(8) :: gamma_dielec
    real(8) :: e_min_dielec
    integer :: num_kgrid_sbe(1:3)
    integer :: num_kgrid_gs(1:3)

contains

    subroutine read_input()
        implicit none
        integer :: ret, fh
        character(256) :: tmp

        namelist/control/ &
        & directory, &
        & sysname, &
        & read_sbe_gs_bin
        namelist/tgrid/ &
        & dt, &
        & nt
        namelist/emfield/ &
        & epdir_re1, &
        & phi_cep1, &
        & rlaser_int_wcm2_1, &
        & pulse_tw1, &
        & omega1, &
        & epdir_im1
        namelist/system/ &
        & al_vec3, &
        & al_vec2, &
        & al_vec1, &
        & nstate, &
        & nelec
        namelist/analysis/ &
        & out_rt_step, &
        & n_dielec, &
        & out_dielec, &
        & e_max_dielec, &
        & gamma_dielec, &
        & e_min_dielec
        namelist/kgrid/ &
        & num_kgrid_sbe, &
        & num_kgrid_gs

        directory = './'
        sysname = 'untitled'
        read_sbe_gs_bin = 'n'
        dt = 0.08
        nt = 5000
        epdir_re1 = (/0.0, 0.0, 1.0/)
        phi_cep1 = 0.0
        rlaser_int_wcm2_1 = 10000000000.0
        pulse_tw1 = 413.5
        omega1 = 0.05698
        epdir_im1 = (/0.0, 0.0, 0.0/)
        al_vec3 = 0.0
        al_vec2 = 0.0
        al_vec1 = 0.0
        nstate = 0
        nelec = 0
        out_rt_step = 10
        n_dielec = 1000
        out_dielec = 'y'
        e_max_dielec = 1.0
        gamma_dielec = 0.005
        e_min_dielec = 0.0
        num_kgrid_sbe = 0
        num_kgrid_gs = 0

        fh = open_filehandle('.namelist.tmp')
        do while (.true.)
            read(*, '(a)', iostat=ret) tmp
            if (ret < 0) exit ! End of file
            tmp = adjustl(tmp)
            if (tmp(1:1) .ne. '!') write(fh, '(a)') trim(tmp)
        end do
        rewind(fh); read(fh, nml=control, iostat=ret)
        rewind(fh); read(fh, nml=tgrid, iostat=ret)
        rewind(fh); read(fh, nml=emfield, iostat=ret)
        rewind(fh); read(fh, nml=system, iostat=ret)
        rewind(fh); read(fh, nml=analysis, iostat=ret)
        rewind(fh); read(fh, nml=kgrid, iostat=ret)

        close(fh)

        write(*, '("# directory =",99(1x,a))') directory
        write(*, '("# sysname =",99(1x,a))') sysname
        write(*, '("# read_sbe_gs_bin =",99(1x,a))') read_sbe_gs_bin
        write(*, '("# dt =",99(1x,es12.4e3))') dt
        write(*, '("# nt =",99(1x,i7))') nt
        write(*, '("# epdir_re1 =",99(1x,es12.4e3))') epdir_re1
        write(*, '("# phi_cep1 =",99(1x,es12.4e3))') phi_cep1
        write(*, '("# rlaser_int_wcm2_1 =",99(1x,es12.4e3))') rlaser_int_wcm2_1
        write(*, '("# pulse_tw1 =",99(1x,es12.4e3))') pulse_tw1
        write(*, '("# omega1 =",99(1x,es12.4e3))') omega1
        write(*, '("# epdir_im1 =",99(1x,es12.4e3))') epdir_im1
        write(*, '("# al_vec3 =",99(1x,es12.4e3))') al_vec3
        write(*, '("# al_vec2 =",99(1x,es12.4e3))') al_vec2
        write(*, '("# al_vec1 =",99(1x,es12.4e3))') al_vec1
        write(*, '("# nstate =",99(1x,i7))') nstate
        write(*, '("# nelec =",99(1x,i7))') nelec
        write(*, '("# out_rt_step =",99(1x,i7))') out_rt_step
        write(*, '("# n_dielec =",99(1x,i7))') n_dielec
        write(*, '("# out_dielec =",99(1x,a))') out_dielec
        write(*, '("# e_max_dielec =",99(1x,es12.4e3))') e_max_dielec
        write(*, '("# gamma_dielec =",99(1x,es12.4e3))') gamma_dielec
        write(*, '("# e_min_dielec =",99(1x,es12.4e3))') e_min_dielec
        write(*, '("# num_kgrid_sbe =",99(1x,i7))') num_kgrid_sbe
        write(*, '("# num_kgrid_gs =",99(1x,i7))') num_kgrid_gs
    end subroutine read_input
end module input_parameter
