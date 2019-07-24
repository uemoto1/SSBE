module test
    implicit none
    contains

subroutine calc_dielec(file_dielec_data, gs, e_min, e_max, ne, gamma)
    use sbe_solver
    use salmon_math, only: pi
    use salmon_file, only: open_filehandle
    implicit none
    character(*), intent(in) :: file_dielec_data
    type(s_sbe_gs), intent(in) :: gs
    real(8), intent(in) :: e_min
    real(8), intent(in) :: e_max
    real(8), intent(in) :: gamma
    integer, intent(in) :: ne
    real(8) :: de, e
    complex(8) :: eps

    integer :: fh, ie, ik, ib, jb

    de = (e_max - e_min) / ne

    fh = open_filehandle(file_dielec_data, 'replace')
    do ie = 1, ne
        e = e_min + de * ie
        eps = 1.00
        do ik = 1, gs%nk
            do ib = 1, gs%nb
                do jb = 1, gs%nb
                    if (gs%occup(ib, ik) < 1.0 .and. 1.0 < gs%occup(jb, ik)) then
                        eps = eps + 2 * (4.0 * pi / gs%volume) &
                        & * zabs(gs%d_matrix(ib, jb, 1, ik)) ** 2 &
                        & / (gs%omega(ib, jb, ik) - e - dcmplx(0d0, gamma)) / gs%nk
                    end if
                end do
            end do
        end do
        write(fh, '(f9.4, 2(1x,e23.15e3))') e, real(eps), aimag(eps)
    end do
    close(fh)
end subroutine calc_dielec





end module test