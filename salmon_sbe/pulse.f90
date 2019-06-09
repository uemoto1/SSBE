module pulse
    use salmon_math, only: pi
    implicit none

    contains


    subroutine calc_cos2_pulse(t, tw, iwcm2, omega, phi_cep, epdir_re, epdir_im, Ac, E)
        use salmon_math, only: pi
        implicit none
        real(8), intent(in) :: t, tw, iwcm2, omega, phi_cep
        real(8), intent(in) :: epdir_re(1:3), epdir_im(1:3)
        real(8), intent(out) :: Ac(1:3), E(1:3)
        real(8) :: e0, tt
        real(8), parameter :: factor = 5.338026811839179d-09
        complex(8), parameter :: zI = dcmplx(0d0, 1d0)

        if (0d0 < t .and. t < tw) then
            e0 = factor * sqrt(iwcm2)
            tt = t - tw * 0.5d0
            Ac(1:3) = (e0 / omega) * real( &
                - zI * cos(pi / tw * tt) ** 2 & 
                * (epdir_re + zI * epdir_im) &
                * exp(-zI * omega * tt - 2 * pi * zI * phi_cep))
            E(1:3) = e0 * real( & 
                & ( - 2.0 * pi * zI / (omega * tw) &
                & * sin(pi / tw * tt) * cos(pi / tw * tt) &
                & + cos(pi / tw * tt) ** 2) * (epdir_re + zI * epdir_im) &
                & * exp(-zI * omega * tt - 2 * pi * zI * phi_cep)) 
        else
            Ac(1:3) = 0d0
            E(1:3) = 0d0
        end if
        return
    end subroutine



end module pulse

