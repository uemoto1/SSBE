module pulse
    use salmon_math, only: pi
    implicit none

    contains


    subroutine calc_pulse_ac(t, tw, iwcm2, omega, phi_cep, epdir_re, epdir_im, Ac)
        use salmon_math, only: pi
        implicit none
        real(8), intent(in) :: t, tw, iwcm2, omega, phi_cep
        real(8), intent(in) :: epdir_re(1:3), epdir_im(1:3)
        real(8), intent(out) :: Ac(1:3)
        real(8) :: f0, tt
        complex(8), parameter :: zI = dcmplx(0d0, 1d0)

        if (0d0 < t .and. t < tw) then
            f0 = (5.338d-9 * sqrt(iwcm2))
            tt = t - tw * 0.5d0
            Ac(1:3) = - f0 / omega * ( &
                & cos(pi * tt / tw) ** 2 * aimag( &
                & (epdir_re(:) + zI * epdir_im(:)) * exp( &
                & zI * (omega * tt + 2d0 * pi * phi_cep))))
        else
            Ac(1:3) = 0d0
        end if
        return
    end subroutine


end module pulse

