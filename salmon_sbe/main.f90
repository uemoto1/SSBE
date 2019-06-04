
module sbe

    implicit none


! Calculate [H, rho] commutation:
subroutine calc_hrho()
    implicit none
    complex(8), parameter :: z1 = dcmplx(1d0, 0d0)
    complex(8) :: zE

    !$omp parallel do default(shared) private(ik, idim, omega, dmat, pmat)
    do ik = ik_s, ik_e
        ! retrieve: omega, pmat(momentum), dmat(dipole)
        call get_matrix_at_k( &
            & kvec(1:3, ik) + ac(1:3), &
            & omega(:, :), dmat(:, :), pmat(:, :))

        ! hrho = de_{ij}(k + A/c) * rho_{ij}(k) 
        hrho(:, :, ik) = omega(:, :) * rho(:, :, ik)

        ! hrho = hrho - E(t) * (d * rho - rho * d)
        do idim = 1, 3
            call ZHEMM('L', 'U', nb, nb, nb, &
                & -E(idim), dmat(:, :, idim), nb, &
                & rho(:, :, ik), nb, &
                & z1, hrho(:, :, ik))
            call ZHEMM('L', 'U', nb, nb, nb, &
                & +E(dim), dmat(:, :, idim), nb, &
                & rho(:, :, ik), nb, &
                & z1, hrho(:, :, ik))
        end do !idim
    end do !ik
    !$omp end parallel do
end subroutine calc_hrho



subroutine dt_evolve()
    impclit none

    complex(8) :: hrho1(nb, nb, ik_s:ik_e)
    complex(8) :: hrho2(nb, nb, ik_s:ik_e)
    complex(8) :: hrho3(nb, nb, ik_s:ik_e)
    complex(8) :: hrho4(nb, nb, ik_s:ik_e)

    call calc_hrho(rho, hrho)
    rho = rho + zI ** 2 /
    call calc_hrho(hrho, hrho)
    call calc_hrho(hrho, hrho)
    call calc_hrho(hrho, hrho)
    

    

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
    do i = 1, N
        do j = 1, N
            A(i, j) = abs(i - j) * 1d0
            if (i == j) then
                B(i, j) = 1.0d0
            end if
        end do !j
    end do !i


    do i = 1, N
        do j = 1, N
            write(*, *) i, j, c(i, j)
        end do
    end do

    stop
end program SSBE

