
module sbe_solver
    use salmon_math, only: pi
    implicit none

    integer, parameter :: ncycle = 10

    type s_sbe_gs
        !Lattice information
        real(8) :: a_matrix(1:3, 1:3)
        real(8) :: b_matrix(1:3, 1:3)
        real(8) :: volume

        !Ground state (GS) electronic system information
        integer :: nk, nb, ne
        real(8), allocatable :: kvec(:, :), kweight(:)
        real(8), allocatable :: eigen(:, :)
        real(8), allocatable :: occup(:, :)
        real(8), allocatable :: omega(:, :, :)
        complex(8), allocatable :: p_matrix(:, :, :, :)
        complex(8), allocatable :: d_matrix(:, :, :, :)
    end type


    type s_sbe
        !k-points for real-time SBE calculation
        integer :: nk, nb
        complex(8), allocatable :: rho(:, :, :)
    end type



contains


subroutine init_sbe_gs(gs, directory, ne, a1, a2, a3)
    use salmon_file, only: open_filehandle
    implicit none
    type(s_sbe_gs), intent(inout) :: gs
    character(*), intent(in) :: directory
    integer, intent(in) :: ne
    real(8), intent(in) :: a1(1:3), a2(1:3), a3(1:3)

    !Calculate b_matrix, volume_cell and volume_bz from a1..a3 vector.
    call calc_lattice_info()    
    !Retrieve k-points from 'kpoint.txt':
    call read_k_data()
    !Retrieve eigenenergies from 'eigen.txt':
    call read_eigen_data()
    !Retrieve transition matrix from 'pmatrix.txt':
    call read_tm_data()
    !Calculate omega and d_matrix (neglecting diagonal part):
    call create_omega_d()

    !Initial Occupation Number
    allocate(gs%occup(1:gs%nb, 1:gs%nk))
    gs%occup = 0
    gs%occup(1:(ne/2),:) = 2d0 !!Experimental!!

contains

    subroutine calc_lattice_info()
        implicit none
        real(8) :: a12(1:3), a23(1:3), a31(1:3), vol
        real(8) :: b1(1:3), b2(1:3), b3(1:3)

        a12(1) = a1(2) * a2(3) - a1(3) * a2(2)
        a12(2) = a1(3) * a2(1) - a1(1) * a2(3)
        a12(3) = a1(1) * a2(2) - a1(2) * a2(1)
        a23(1) = a2(2) * a3(3) - a2(3) * a3(2)
        a23(2) = a2(3) * a3(1) - a2(1) * a3(3)
        a23(3) = a2(1) * a3(2) - a2(2) * a3(1)
        a31(1) = a3(2) * a1(3) - a3(3) * a1(2)
        a31(2) = a3(3) * a1(1) - a3(1) * a1(3)
        a31(3) = a3(1) * a1(2) - a3(2) * a1(1)
        vol = dot_product(a12, a3)
        b1(1:3) = (2d0 * pi / vol) * a23(1:3)
        b2(1:3) = (2d0 * pi / vol) * a31(1:3)
        b3(1:3) = (2d0 * pi / vol) * a12(1:3)
        
        gs%a_matrix(1:3, 1) = a1(1:3)
        gs%a_matrix(1:3, 2) = a2(1:3)
        gs%a_matrix(1:3, 3) = a3(1:3)    
        gs%b_matrix(1, 1:3) = b1(1:3)
        gs%b_matrix(2, 1:3) = b2(1:3)
        gs%b_matrix(3, 1:3) = b3(1:3)
        gs%volume = vol
    end subroutine calc_lattice_info


    subroutine read_k_data()
        implicit none
        integer :: fh, i, ik, iik

        fh = open_filehandle(trim(directory) // 'kpoint.txt', 'old')
        read(fh, *) gs%nk ! Number of k-points in BZ
        allocate(gs%kvec(1:3, 1:gs%nk))
        do ik=1, gs%nk
            read(fh, *) iik, gs%kvec(1:3, ik)
            if (ik .ne. iik) &
                & stop "ERROR! Invalid kpoint.txt"
        end do !ik
        close(fh)
        ! k-point integration weights:
        allocate(gs%kweight(1:gs%nk))
        gs%kweight = 1d0 / gs%nk
    end subroutine read_k_data


    subroutine read_eigen_data()
        implicit none
        integer :: fh, i, ik, ib, iik, iib

        fh = open_filehandle(trim(directory)  // 'eigen.txt', 'old')
        read(fh, *) gs%nk, gs%nb
        allocate(gs%eigen(1:gs%nb, 1:gs%nk))
        do ik=1, gs%nk
            do ib=1, gs%nb
                read(*, *) iik, iib, gs%eigen(ib, ik)
                if ((ik .ne. iik) &
                    & .or. (ib .ne. iib)) &
                    & stop "ERROR! Invalid eigen.txt"
            end do
        end do
        close(fh)
    end subroutine read_eigen_data





    subroutine read_tm_data()
        implicit none
        character(256) :: dummy
        integer :: fh, i, ik, ib, jb, iik, iib, jjb
        real(8) :: tmp(1:6)

        fh = open_filehandle(trim(directory) // 'pmatrix.txt', 'old')
        read(fh, *) gs%nk, gs%nb
        allocate(gs%p_matrix(1:gs%nb, 1:gs%nb, 1:3, 1:gs%nk))

        do ik=1, gs%nk
            do ib=1, gs%nb
                do jb=1, gs%nb
                    read(fh, *) iik, iib, jjb, tmp(1:6)
                    if ((ik .ne. iik) &
                        & .or. (ib .ne. iib) &
                        & .or. (jb .ne. jjb)) &
                        & stop "ERROR! Invalid pmatrix.txt"
                    gs%p_matrix(ib, jb, 1, ik) = dcmplx(tmp(1), tmp(2))
                    gs%p_matrix(ib, jb, 2, ik) = dcmplx(tmp(3), tmp(4))
                    gs%p_matrix(ib, jb, 3, ik) = dcmplx(tmp(5), tmp(6))
                end do !jb
            end do !ib
        end do !ik
        close(fh)
    end subroutine read_tm_data


    subroutine create_omega_d()
        implicit none
        integer :: ik, ib, jb
        real(8), parameter :: epsilon = 0.04
        complex(8), parameter :: zi = dcmplx(0d0, 1d0)
        allocate(gs%omega(1:gs%nb, 1:gs%nb, 1:gs%nk))
        allocate(gs%d_matrix(1:gs%nb, 1:gs%nb, 1:3, 1:gs%nk))
        do ik=1, gs%nk
            do ib=1, gs%nb
                do jb=1, gs%nb
                    gs%omega(ib, jb, ik) = gs%eigen(ib, ik) - gs%eigen(jb, ik)
                    if (epsilon < abs(gs%omega(ib, jb, ik))) then
                        gs%d_matrix(ib, jb, 1:3, ik) = &
                            & zi * gs%p_matrix(ib, jb, 1:3, ik) / gs%omega(ib, jb, ik)
                    else
                        gs%d_matrix(ib, jb, 1:3, ik) = 0d0
                    end if
                end do
            end do
        end do
    end subroutine create_omega_d

    

end subroutine init_sbe_gs



subroutine init_sbe(sbe, gs)
    implicit none
    type(s_sbe), intent(inout) :: sbe
    type(s_sbe_gs), intent(in) :: gs

    sbe%nk = gs%nk
    sbe%nb = gs%nb

    allocate(sbe%rho(1:sbe%nb, 1:sbe%nb, 1:sbe%nk))
    
    call init_rho()

 contains

    subroutine init_rho()
        implicit none
        integer :: ik, ib
        ! Initial state of density matrix:
        sbe%rho = 0d0
        do ib = 1, sbe%nb
            sbe%rho(ib, ib, :) = gs%occup(ib, 1)
        end do
    end subroutine

end subroutine


subroutine calc_current_bloch(sbe, gs, Ac, jmat)
    implicit none
    type(s_sbe), intent(in) :: sbe
    type(s_sbe_gs), intent(in) :: gs
    real(8), intent(in) :: Ac(1:3)
    real(8), intent(out) :: jmat(1:3)
    complex(8) :: pk(sbe%nb, sbe%nb, 3)
    integer :: ik, idir, ib, jb
    real(8) :: jtot(1:3)

    jtot(1:3) = 0d0

    !$omp parallel do default(shared) private(ik,ib,jb,idir,pk) reduction(+:jtot)
    do ik=1, sbe%nk
        pk(:, :, :) = gs%p_matrix(:, :, :, ik)
        do idir=1, 3
            do jb=1, sbe%nb
                do ib=1, sbe%nb
                    jtot(idir) = jtot(idir) + gs%kweight(ik) &
                        & * real(pk(ib, jb, idir) * sbe%rho(jb, ib, ik)) 
                end do
            end do
        end do
    end do
    !$omp end parallel do
    jmat(:) =  (jtot(:) + gs%ne * ac(:)) / gs%volume
    

    return
end subroutine calc_current_bloch










subroutine dt_evolve(sbe, gs, E, Ac, dt)
    implicit none
    type(s_sbe), intent(inout) :: sbe
    type(s_sbe_gs), intent(in) :: gs
    real(8), intent(in) :: E(1:3)
    real(8), intent(in) :: Ac(1:3)
    real(8), intent(in) :: dt
    complex(8), parameter :: zi = dcmplx(0d0, 1d0)

    complex(8) :: hrho1(1:sbe%nb, 1:sbe%nb, 1:sbe%nk)
    complex(8) :: hrho2(1:sbe%nb, 1:sbe%nb, 1:sbe%nk)
    complex(8) :: hrho3(1:sbe%nb, 1:sbe%nb, 1:sbe%nk)
    complex(8) :: hrho4(1:sbe%nb, 1:sbe%nb, 1:sbe%nk)

    call calc_hrho_bloch(sbe%rho, hrho1)
    call calc_hrho_bloch(hrho1, hrho2)
    call calc_hrho_bloch(hrho2, hrho3)
    call calc_hrho_bloch(hrho3, hrho4)

    sbe%rho = sbe%rho + hrho1 * (- zi * dt)
    sbe%rho = sbe%rho + hrho2 * (- zi * dt) ** 2 * (1d0 / 2d0)
    sbe%rho = sbe%rho + hrho3 * (- zi * dt) ** 3 * (1d0 / 6d0)
    sbe%rho = sbe%rho + hrho4 * (- zi * dt) ** 4 * (1d0 / 24d0)
    return

contains


    !Calculate [H, rho] commutation:
    subroutine calc_hrho_bloch(rho, hrho)
        implicit none
        complex(8), intent(in)     :: rho(sbe%nb, sbe%nb, sbe%nk)
        complex(8), intent(out)    :: hrho(sbe%nb, sbe%nb, sbe%nk)
        integer :: ik, ib, jb, idir
        real(8) :: ek(sbe%nb)
        complex(8) :: pk(sbe%nb, sbe%nb, 3)
        !$omp parallel do default(shared) private(ik,ib,jb,idir,ek,pk)
        do ik=1, gs%nk
            ek(:) = gs%eigen(:, ik)
            pk(:, :, :) = gs%p_matrix(:, :, :, ik)
            !hrho(k) = omega(k + A/c) * rho(k) 
            do ib = 1, sbe%nb
                do jb = 1, sbe%nb
                    hrho(ib, jb, ik) = (ek(ib) - ek(jb)) * rho(ib, jb, ik)
                end do
            end do
            !hrho = hrho + Ac(t) * (p * rho - rho * p)
            do idir=1, 3 !1:x, 2:y, 3:z
                call ZHEMM('L', 'U', sbe%nb, sbe%nb, &
                    & dcmplx(+Ac(idir), 0d0), &
                    & pk, sbe%nb, &
                    & rho(:, :, ik), sbe%nb, &
                    & dcmplx(1d0, 0d0), hrho(:, :, ik), sbe%nb)
                call ZHEMM('L', 'U', sbe%nb, sbe%nb, &
                    & dcmplx(-Ac(idir), 0d0), &
                    & rho(:, :, ik), sbe%nb, &
                    & pk, sbe%nb, &
                    & dcmplx(1d0, 0d0), hrho(:, :, ik), sbe%nb)
            end do !idir
        end do !ik
        !$omp end parallel do
        return
    end subroutine calc_hrho_bloch


end subroutine

function calc_trace(sbe, gs, nb_max) result(tr)
    implicit none
    type(s_sbe), intent(in) :: sbe
    type(s_sbe_gs), intent(in) :: gs
    integer, intent(in) :: nb_max
    complex(8), parameter :: zi = dcmplx(0d0, 1d0)
    integer :: ik, ib
    real(8) :: tr
    tr = 0d0
    !$omp parallel do default(shared) private(ik, ib) reduction(+: tr) collapse(2) 
    do ik = 1, sbe%nk
        do ib = 1, nb_max
            tr = tr + real(sbe%rho(ib, ib, ik)) * gs%kweight(ik)
        end do
    end do
    !$omp end parallel do
    tr = tr / sum(gs%kweight)
    return 
end function calc_trace


end module



