
module sbe_solver
    use salmon_math, only: pi
    use sbe_gs
    implicit none



    type s_sbe_solver
        !k-points for real-time SBE calculation
        integer :: nk, nb
        real(8), allocatable :: kvec(:, :), kweight(:)
        complex(8), allocatable :: rho(:, :, :)
        ! integer :: ik_s, ik_e, icomm_k
    end type



contains



subroutine init_sbe(sbe, gs, nkgrid)
    implicit none
    type(s_sbe_solver), intent(inout) :: sbe
    type(s_sbe_gs), intent(in) :: gs
    integer, intent(in) :: nkgrid(1:3)

    integer :: nk, nb

    nk = nkgrid(1) * nkgrid(2) * nkgrid(3)
    nb = gs%nb

    sbe%nk = nk
    sbe%nb = nb

    allocate(sbe%rho(1:nb, 1:nb, 1:nk))
    allocate(sbe%kvec(1:3, 1:nk))
    allocate(sbe%kweight(1:nk))
    
    call init_rho()
    call create_uniform_kgrid()

 contains

    subroutine create_uniform_kgrid()
        implicit none 
        integer :: ik1, ik2, ik3, ik
        real(8) :: b1(1:3), b2(1:3), b3(1:3)
        real(8) :: h1, h2, h3
        real(8) :: dh(1:3)

        b1(1:3) = gs%b_matrix(1, 1:3)
        b2(1:3) = gs%b_matrix(2, 1:3)
        b3(1:3) = gs%b_matrix(3, 1:3)
        dh(1:3) = 1d0 / nkgrid(1:3)

        ik = 1
        do ik3=1, nkgrid(3)
            h3 = dh(3) * (ik3 - 0.5d0) - 0.5d0
            do ik2=1, nkgrid(2)
                h2 = dh(2) * (ik2 - 0.5d0) - 0.5d0
                do ik1=1, nkgrid(1)
                    h1 = dh(1) * (ik1 - 0.5d0) - 0.5d0
                    ! Uniformally sampled k-grid point:
                    sbe%kvec(1:3, ik) = h1 * b1(1:3) + h2 * b2(1:3) + h3 * b3(1:3)
                    sbe%kweight(ik) = (1d0 / nk)
                    ik = ik + 1
                end do
            end do
        end do
    end subroutine


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




subroutine interp_gs(gs, kvec, e_k, d_k, p_k)
    implicit none
    type(s_sbe_gs), intent(in) :: gs
    real(8), intent(in) :: kvec(3)
    real(8), intent(out), optional :: e_k(gs%nb)
    complex(8), intent(out), optional :: d_k(gs%nb, gs%nb, 1:3)
    complex(8), intent(out), optional :: p_k(gs%nb, gs%nb, 1:3)
    real(8) :: rkg(1:3), wkg(1:2,1:3), wj
    integer :: ikg(1:2, 1:3), jkg(1:3), j1, j2, j3, jk1, jk2, jk3, jk
    
    ! Calculate reduced coordinate of kvec:
    rkg(1:3) = (matmul(kvec, gs%a_matrix) / (2d0 * pi) + 0.5d0) * gs%nkgrid + 0.5d0

    if (present(e_k)) then
        ! Grid coordinates and weights for linear interpolation:
        ikg(1, 1:3) = int(floor(rkg(1:3)))      ! Lower index of data point
        ikg(2, 1:3) = ikg(1, 1:3) + 1           ! Upper index of data point
        wkg(2, 1:3) = rkg(1:3) - ikg(1, 1:3)
        wkg(1, 1:3) = 1d0 - wkg(2, 1:3)

        e_k = 0d0

        do j3 = 1, 2
            jk3 = gs%ikcycle_tbl3(ikg(j3, 3))
            do j2 = 1, 2
                jk2 = gs%ikcycle_tbl2(ikg(j2, 2))
                do j1 = 1, 2
                    jk1 = gs%ikcycle_tbl1(ikg(j1, 1))
                    ! write(*,*) j1, j2, j3, jk1, jk2, jk3, "j"
                    ! write(*,*) wkg(j1, 1), wkg(j2, 2) , wkg(j3, 3), wkg(j1, 1) * wkg(j2, 2) * wkg(j3, 3), "w"
                    ! Calculate interpolate coefficients:
                    jk = gs%iktbl_grid(jk1, jk2, jk3)
                    wj = wkg(j1, 1) * wkg(j2, 2) * wkg(j3, 3)
                    ! Interpolate e, d and p at given point:
                    e_k(:) = e_k(:) + wj * gs%eigen(:, jk)
                end do
            end do
        end do
    end if

    ! Nearlest grid coordinate
    jkg(1:3) = int(floor(rkg(1:3) + 0.5d0))    ! Rounded index of data point
    jk3 = gs%ikcycle_tbl3(jkg(3))
    jk2 = gs%ikcycle_tbl2(jkg(2))
    jk1 = gs%ikcycle_tbl1(jkg(1))
    jk = gs%iktbl_grid(jk1, jk2, jk3)
    if (present(d_k)) d_k(:, :, :) = gs%d_matrix(:, :, :, jk)
    if (present(p_k)) p_k(:, :, :) = gs%p_matrix(:, :, :, jk)

    return
end subroutine interp_gs


subroutine band_test(gs)
    implicit none
    type(s_sbe_gs), intent(in) :: gs
    integer :: i
    real(8) :: k(3), e(gs%nb)

    do i = -100, 100
        k(:) = (i * 0.01) * gs%b_matrix(1, :)
        call interp_gs(gs, k, e_k=e)
        write(*,*) i, e
    end do
end subroutine band_test






subroutine calc_current(sbe, gs, Ac, jmat)
    implicit none
    type(s_sbe_solver), intent(in) :: sbe
    type(s_sbe_gs), intent(in) :: gs
    real(8), intent(in) :: Ac(1:3)
    real(8), intent(out) :: jmat(1:3)
    complex(8) :: pk(sbe%nb, sbe%nb, 3)
    integer :: ik, idir, ib, jb
    real(8) :: jtot(1:3)

    jtot(1:3) = 0d0

    !$omp parallel do default(shared) private(ik,ib,jb,idir,pk) reduction(+:jtot)
    do ik=1, sbe%nk
        call interp_gs(gs, sbe%kvec(1:3, ik) + Ac(1:3), p_k=pk)
        do idir=1, 3
            do jb=1, sbe%nb
                do ib=1, sbe%nb
                    jtot(idir) = jtot(idir) + sbe%kweight(ik) &
                        & * real(pk(ib, jb, idir) * sbe%rho(jb, ib, ik)) 
                end do
            end do
        end do
    end do
    !$omp end parallel do
    jmat(:) =  (jtot(:) ) / gs%volume

    return
end subroutine calc_current








subroutine calc_current_bloch(sbe, gs, Ac, jmat)
    implicit none
    type(s_sbe_solver), intent(in) :: sbe
    type(s_sbe_gs), intent(in) :: gs
    real(8), intent(in) :: Ac(1:3)
    real(8), intent(out) :: jmat(1:3)
    complex(8) :: pk(sbe%nb, sbe%nb, 3)
    integer :: ik, idir, ib, jb
    real(8) :: jtot(1:3)

    jtot(1:3) = 0d0

    !$omp parallel do default(shared) private(ik,ib,jb,idir,pk) reduction(+:jtot)
    do ik=1, sbe%nk
        call interp_gs(gs, sbe%kvec(1:3, ik) + Ac(1:3), p_k=pk)
        do idir=1, 3
            do jb=1, sbe%nb
                do ib=1, sbe%nb
                    jtot(idir) = jtot(idir) + sbe%kweight(ik) &
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
    type(s_sbe_solver), intent(inout) :: sbe
    type(s_sbe_gs), intent(in) :: gs
    real(8), intent(in) :: E(1:3)
    real(8), intent(in) :: Ac(1:3)
    real(8), intent(in) :: dt
    complex(8), parameter :: zi = dcmplx(0d0, 1d0)

    complex(8) :: hrho1(1:sbe%nb, 1:sbe%nb, 1:sbe%nk)
    complex(8) :: hrho2(1:sbe%nb, 1:sbe%nb, 1:sbe%nk)
    complex(8) :: hrho3(1:sbe%nb, 1:sbe%nb, 1:sbe%nk)
    complex(8) :: hrho4(1:sbe%nb, 1:sbe%nb, 1:sbe%nk)

    call calc_hrho(sbe%rho, hrho1)
    call calc_hrho(hrho1, hrho2)
    call calc_hrho(hrho2, hrho3)
    call calc_hrho(hrho3, hrho4)

    sbe%rho = sbe%rho + hrho1 * (- zi * dt)
    sbe%rho = sbe%rho + hrho2 * (- zi * dt) ** 2 * (1d0 / 2d0)
    sbe%rho = sbe%rho + hrho3 * (- zi * dt) ** 3 * (1d0 / 6d0)
    sbe%rho = sbe%rho + hrho4 * (- zi * dt) ** 4 * (1d0 / 24d0)
    return

contains

    !Calculate [H, rho] commutation:
    subroutine calc_hrho(rho, hrho)
        implicit none
        complex(8), intent(in)     :: rho(sbe%nb, sbe%nb, sbe%nk)
        complex(8), intent(out)    :: hrho(sbe%nb, sbe%nb, sbe%nk)
        real(8) :: ek(sbe%nb)
        complex(8) :: dk(sbe%nb, sbe%nb, 3)
        integer :: ik, idir, ib, jb
        !$omp parallel do default(shared) private(ik,ib,jb,idir,ek,dk)
        do ik=1, sbe%nk
            call interp_gs(gs, sbe%kvec(1:3, ik) + Ac(1:3),  e_k=ek, d_k=dk)
            !hrho(k) = omega(k + A/c) * rho(k) 
            do ib = 1, sbe%nb
                do jb = 1, sbe%nb
                    hrho(ib, jb, ik) = (ek(ib) - ek(jb)) * rho(ib, jb, ik)
                end do
            end do
            !hrho = hrho - E(t) * (d * rho - rho * d)
            do idir=1, 3 !1:x, 2:y, 3:z
                call ZHEMM('L', 'U', sbe%nb, sbe%nb, &
                    & dcmplx(-E(idir), 0d0), &
                    & dk(:, :, idir), sbe%nb, &
                    & rho(:, :, ik), sbe%nb, &
                    & dcmplx(1d0, 0d0), hrho(:, :, ik), sbe%nb)
                call ZHEMM('L', 'U', sbe%nb, sbe%nb, &
                    & dcmplx(+E(idir), 0d0), &
                    & rho(:, :, ik), sbe%nb, &
                    & dk(:, :, idir), sbe%nb, &
                    & dcmplx(1d0, 0d0), hrho(:, :, ik), sbe%nb)
            end do !idir
        end do !ik
        !$omp end parallel do
        return
    end subroutine calc_hrho


    !Calculate [H, rho] commutation:
    subroutine calc_hrho_bloch(rho, hrho)
        implicit none
        complex(8), intent(in)     :: rho(sbe%nb, sbe%nb, sbe%nk)
        complex(8), intent(out)    :: hrho(sbe%nb, sbe%nb, sbe%nk)
        integer :: ik, ib, jb, idir
        real(8) :: ek(sbe%nb)
        complex(8) :: pk(sbe%nb, sbe%nb, 3)
        !$omp parallel do default(shared) private(ik,ib,jb,idir,ek,pk)
        do ik=1, sbe%nk
            call interp_gs(gs, sbe%kvec(1:3, ik),  e_k=ek, p_k=pk)
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

function calc_trace(sbe, nb_max) result(tr)
    implicit none
    type(s_sbe_solver), intent(in) :: sbe
    integer, intent(in) :: nb_max
    complex(8), parameter :: zi = dcmplx(0d0, 1d0)
    integer :: ik, ib
    real(8) :: tr
    tr = 0d0
    !$omp parallel do default(shared) private(ik, ib) reduction(+: tr) collapse(2) 
    do ik = 1, sbe%nk
        do ib = 1, nb_max
            tr = tr + real(sbe%rho(ib, ib, ik)) * sbe%kweight(ik)
        end do
    end do
    !$omp end parallel do
    tr = tr / sum(sbe%kweight)
    return 
end function calc_trace


end module



