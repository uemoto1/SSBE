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




subroutine test_interp(gs)
    use sbe_solver
    type(s_sbe_gs), intent(in) :: gs

    integer :: ik
    real(8) :: kv(3), ek(gs%nb)
    complex(8) :: dk(gs%nb, gs%nb, 3), pk(gs%nb, gs%nb, 3)

    write(99,*) "# Total size", gs%nb

    do ik = -100, 100
        kv(:) = gs%b_matrix(1,:) * ik * 0.01
        call interp_gs(gs, kv, e_k=ek, d_k=dk, p_k=pk)
        write(99, *) ik, ek
        write(98, *) ik, real(pk(1, 2, 1)), aimag(pk(1, 2, 1)), abs(pk(1, 2, 1))
        write(97, *) ik, real(pk(4, 5, 1)), aimag(pk(4, 5, 1)), abs(pk(4, 5, 1))
    end do

        
end subroutine test_interp




end module test