
module sbe
    use salmon_file, only: open_filehandle
    implicit none


    type s_sbe
        ! k-point 
        integer :: nk, nb
        real(8), allocatable :: kvec(:, :), kweight(:)
    end type


    type s_sbe_gs
        ! Lattice information
        real(8) :: a_matrix(1:3, 1:3)
        real(8) :: b_matrix(1:3, 1:3)
        real(8) :: volume_cell, volume_bz

        ! Ground state (GS) electronic system information
        integer :: nk, nb, ne
        real(8), allocatable :: kvec(:, :), kweight(:)
        real(8), allocatable :: eigen(:, :)
        real(8), allocatable :: occup(:, :)
        real(8), allocatable :: omega(:, :, :)
        real(8), allocatable :: p_matrix(:, :, :, :)
        real(8), allocatable :: d_matrix(:, :, :, :)
        real(8), allocatable :: rv_matrix(:, :, :)

        ! k-space grid and geometry information
        ! NOTE: prepred for uniformally distributed k-grid....
        integer :: nkgrid(1:3)
        integer, allocatable :: iktbl_grid(:, :, :)
    end type






    contains

! Calculate ground state k-point index (ik_gs) corresponding to a given k-vector:
function get_ik_gs(gs, kvec) result(ik_gs) 
    implicit none
    type(s_sbe_gs), intent(in) :: gs
    real(8), intent(in) :: kvec(3)
    real(8) :: rka(3)    ! Reduced coordinate of k*a
    integer :: ikgrid(3) ! Discretized 3D grid coordinate of k*a
    integer :: ik_gs
    real(8), parameter :: one_over_2pi = 0.5 / 3.141592653589793
    rka = one_over_2pi * matmul(kvec, gs%a_matrix) + 0.5d0
    ikgrid = (rka - floor(rka)) * gs%nkgrid
    ik_gs = iktbl_grid(ik_grid(1), ik_grid(2), ik_grid(3))
    return
end function get_ik_gs

subroutine init_gs_info(gs, sysname, directory, nkgrid, nb, ne)
    implicit none
    type(s_sbe_gs), intent(inout) :: gs
    character(*), intent(in) :: sysname
    character(*), intent(in) :: directory
    integer, intent(in) :: nkgrid(1:3)
    integer, intent(in) :: nb
    integer, intent(in) :: ne
    integer :: nk

    nk = nkgrid(1) * nkgrid(2) * nkgrid(3)

    gs%nk = nk
    gs%nb = nb
    gs%ne = ne
    gs%nkgrid(1:3) = nkgrid(1:3)

    allocate(gs%kvec(1:3, 1:nk))
    allocate(gs%kweight(1:nk))
    allocate(gs%eigen(1:nb, 1:nk))
    allocate(gs%occup(1:nb, 1:nk))
    allocate(gs%omega(1:nb, 1:nb, 1:nk))
    allocate(gs%p_matrix(1:nb, 1:nb, 3, 1:nk))
    allocate(gs%d_matrix(1:nb, 1:nb, 3, 1:nk))
    allocate(gs%rv_matrix(1:nb, 1:nb, 1:nk))
    allocate(gs%iktbl_grid(1:nkgrid(1), 1:nkgrid(2), 1:nkgrid(3)))
    
    ! Retrieve eigenenergies from 'SYSNAME_eigen.data':
    call read_eigen_data()
    ! Retrieve k-points from 'SYSNAME_k.data':
    call read_k_data()
    ! Retrieve transition matrix from 'SYSNAME_tm.data':
    call read_tm_data()
    ! Calculate iktbl_grid for uniform (non-symmetric) k-grid:
    call create_uniform_iktbl_grid()
    ! Calculate omega and d_matrix (neglecting diagonal part):
    call create_omega_dmatrix()

    ! Initial Occupation Number
    gs%occup(1:(ne/2), :) = 2d0 !! Experimental !!

contains

    subroutine read_eigen_data()
        implicit none
        character(256) :: dummy
        integer :: fh, ik, ib, idummy

        fh = open_filehandle(trim(directory) // trim(sysname) // '_eigen.data')
        do i=1, 3
            read(fh, *) dummy ! Skip
        end do
        do ik=1, nk
            read(fh, *) dummy ! Skip
            do ib=1, nb
                read(fh) idummy, gs%eigen(ib, ik)
            end do
        end do
        close(fh)
    end subroutine read_eigen_data


    subroutine read_k_data()
        implicit none
        character(256) :: dummy
        integer :: fh, ik, ib, idummy

        fh = open_filehandle(trim(directory) // trim(sysname) // '_k.data')
        do i=1, 8
            read(fh, *) dummy ! Skip
        end do
        do ik=1, nk
            read(fh, *) dummy ! Skip
            do ib=1, nb
                read(fh) idummy, gs%kvec(1:3, ik), gs%kweight(ik)
            end do !ib
        end do !ik
        close(fh)
    end subroutine read_k_data


    subroutine read_tm_data()
        implicit none
        character(256) :: dummy
        integer :: fh, ik, ib, idummy
        real(8) :: tmp(1:6)

        fh = open_filehandle(trim(directory) // trim(sysname) // '_tm.data')
        do i=1, 3
            read(fh, *) dummy ! Skip
        end do
        do ik=1, nk
            do ib=1, nb
                do jb=1, nb
                    read(fh) idummy, idummy, idummy, tmp(1:6)
                    gs%p_matrix(ib, jb, 1, ik) = dcmplx(tmp(1), tmp(2))
                    gs%p_matrix(ib, jb, 2, ik) = dcmplx(tmp(3), tmp(4))
                    gs%p_matrix(ib, jb, 3, ik) = dcmplx(tmp(5), tmp(6))
                end do !jb
            end do !ib
        end do !ik
        close(fh)
    end subroutine read_tm_data


    subroutine create_omega_dmatrix()
        implicit none
        integer :: ik, ib, jb
        real(8), parameter :: epsilon = 1e-5
        complex(8), parameter :: zi = dcmplx(0d0, 1d0)
        !$omp parallel do default(shared) private(ik, ib, jb)
        do ik=1, nb
            do ib=1, nb
                do jb=1, nb
                    gs%omega(ib, jb, ik) = gs%eigen(ib) - gs%eigen(jb)
                    if (epsilon < abs(gs%omega(ib, jb, ik))) then
                        gs%d_matrix(ib, jb, ik) = &
                            & zi * gs%pmatrix(ib, jb, ik) / gs%omega(ib, jb, ik)
                    end if
                end do
            end do
        end do
        !$omp end parallel do
    end subroutine create_omega_dmatrix

    


    subroutine create_uniform_iktbl_grid()
        implicit none
        ! Based on GCEED's periodic system setup...
        integer :: ik1, ik2, ik3, ikcount

        ikcount = 1
        do ik3=1, nkgrid(3)
            do ik2=1, nkgrid(2)
                do ik1=1, nkgrid(1)
                    gs%iktbl_grid(ik1, ik2, ik3) = ikcount
                    ikcount = ikcount + 1
                end do ! ik1
            end do ! ik2
        end do ! ik3
    end subroutine create_uniform_iktbl_grid




! Calculate [H, rho] commutation:
subroutine calc_hrho(sbe, gs, E, Ac, rho, Hrho)
    implicit none
    type(s_sbe), intent(in)    :: sbe
    type(s_sbe_gs), intent(in) :: gs
    real(8), intent(in)        :: E(3), Ac(3)
    complex(8), intent(in)     :: rho(sbe%nb, sbe%nb, sbe%nk)
    complex(8), intent(out)    :: hrho(sbe%nb, sbe%nb, sbe%nk)

    !$omp parallel do default(shared) private(ik, ikAc_gs, idir)
    do ik=1, sbe%nk
        ikAc_gs = get_ik_gs(gs, sbe%kvec(1:3, ik) + Ac(1:3))
        ! hrho(k) = omega(k + A/c) * rho(k) 
        hrho(:, :, ik) = gs%omega(:, :, ikAc_gs) * rho(:, :, ik)
        ! hrho = hrho - E(t) * (d * rho - rho * d)
        do idir=1, 3 ! 1:x, 2:y, 3:z
            call ZHEMM('L', 'U', nb, nb, nb, &
                & dcmplx(-E(idir), 0d0), &
                & gs%d_matrix(:, :, idir, ikAc_gs), nb, &
                & rho(:, :, ik), nb, &
                & dcmplx(1d0, 0d0), hrho(:, :, ik))
            call ZHEMM('L', 'U', nb, nb, nb, &
                & dcmplx(+E(idir), 0d0), &
                & rho(:, :, ik), nb, &
                & gs%d_matrix(:, :, idir, ikAc_gs), nb, &
                & dcmplx(1d0, 0d0), hrho(:, :, ik))
        end do !idir
    end do !ik
    !$omp end parallel do
end subroutine calc_hrho


subroutine calc_current(sbe, gs, )
    implicit none
    integer :: ik
    real(8) :: jcur(1:3)

    ! Local diagonal component
    !$omp parallel do default(shared) private(ik, idir, ib, jb, kAc, ikAc_gs) reduction(+: jcur)
    do ik=1, sbe%nk
        ikAc_gs = get_ik_gs(gs, sbe%kvec(ik, 1:3) + Ac(1:3))
        do idir=1, 3
            do ib=1, nb
                do jb=1, nb
                    jcur(idir) = jcur(idir) + real(gs%pmat(ib, jb, idir, ikAc_gs) * rho(jb, ib, ik))
                end do
                jcur(idir) = jcur(idir) + gs%kvec(idir, ikAc_gs) * rho(ib, ib, ik)
            end do
        end do
    end do
    !$omp end parallel do

    return
end subroutine calc_current


subroutine dt_evolve(sbe)
    implicit none
    type(s_sbe_info), intent(in) :: sbe

    complex(8) :: hrho1(1:sbe%nb, 1:sbe%nb, 1:sbe%nk)
    complex(8) :: hrho2(1:sbe%nb, 1:sbe%nb, 1:sbe%nk)
    complex(8) :: hrho3(1:sbe%nb, 1:sbe%nb, 1:sbe%nk)
    complex(8) :: hrho4(1:sbe%nb, 1:sbe%nb, 1:sbe%nk)

    call calc_hrho(sbe, gs, E, Ac, rho, hrho1)
    call calc_hrho(sbe, gs, E, Ac, hrho1, hrho2)
    call calc_hrho(sbe, gs, E, Ac, hrho2, hrho3)
    call calc_hrho(sbe, gs, E, Ac, hrho3, hrho4)
    
end subroutine

    

program SSBE
    implicit none
    integer, parameter :: N = 5
    integer :: i, j
    complex(8) :: A(N, N), B(N, N), C(N, N)
    complex(8), parameter :: z1 = dcmplx(1d0, 0d0)
    complex(8), parameter :: zi = dcmplx(0d0, 1d0)
    complex(8), parameter :: z0 = dcmplx(0d0, 0d0)

    A = 0
    B = 0
    C = 0
    do i=1, N
        do j=1, N
            A(i, j) = abs(i - j) * 1d0
            if (i == j) then
                B(i, j) = 1.0d0
            end if
        end do !j
    end do !i


    do i=1, N
        do j=1, N
            write(*, *) i, j, c(i, j)
        end do
    end do

    stop
end program SSBE
