subroutine test01(sbe, gs)
    use sbe_solver
    implicit none
    type(s_sbe), intent(in) :: sbe
    type(s_sbe_gs), intent(in) :: gs
    
    real(8), parameter :: e_min =  0.00
    real(8), parameter :: e_max = +1.00
    real(8), parameter :: gamma = 0.002
    integer, parameter :: ne = 1000
    real(8) :: de, e
    complex(8) :: eps

    integer :: ie, ik, ib, jb

    de = (e_max - e_min) / ne

    do ie = 1, ne
        e = e_min + de * ie
        eps = 1d0
        do ik = 1, gs%nk
            do ib = 1, gs%nb
                do jb = 1, gs%nb
                    ! write(*,*) ik, ib, jb, gs%omega(ib, jb, ik), "#OMG2"
                    if (gs%occup(ib, ik) < 1.0 .and. 1.0 < gs%occup(jb, ik)) then
                        !write(*,*) "accum", ib, jb
                        eps = eps + 2 * (4.0 * 3.14159265 / gs%volume) &
                        & * zabs(gs%d_matrix(ib, jb, 1, ik)) ** 2 &
                        & / (gs%omega(ib, jb, ik) - e - dcmplx(0d0, gamma)) / gs%nk
                    end if
                end do
            end do
        end do
        write(*,*) e, real(eps), aimag(eps), "#EPS"
    end do
end subroutine test01
