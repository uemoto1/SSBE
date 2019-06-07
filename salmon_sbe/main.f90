program main
    use sbe_solver
    use input_parameter
    use pulse
    implicit none

    type(s_sbe) :: sbe
    type(s_sbe_gs) :: gs
    real(8) :: t, E(3), Ac(3), Ac_dt(3)
    integer :: ib, it

    call read_input()

    ! Read ground state electronic system:
    call init_sbe_gs(gs, sysname, directory, &
        & num_kgrid, nstate, nelec, &
        & al_vec1, al_vec2, al_vec3)

    ! Export dielectric spectra to SYSNAME_dielec.data:
    call calc_dielec(trim(directory) // trim(sysname) // '_dielec.data', &
        & gs, e_min_dielec, e_max_dielec, n_dielec, gamma_dielec)

    ! Initialization of SBE solver:
    call init_sbe(sbe, gs, num_kgrid, nt, dt)
    ! Initialization of density matrix:
    sbe%rho = 0d0
    do ib = 1, sbe%nb
        sbe%rho(ib, ib, :) = gs%occup(ib, 1)
    end do

    ! Realtime calculation
    do it = 1, sbe%nt
        t = sbe%dt * it
        call calc_pulse_ac(t, pulse_tw1, rlaser_int_wcm2_1, &
            & omega1, phi_cep1, epdir_re1, epdir_im1, Ac)
        call calc_pulse_ac(t + dt, pulse_tw1, rlaser_int_wcm2_1, &
            & omega1, phi_cep1, epdir_re1, epdir_im1, Ac_dt)
        E = - (Ac_dt - Ac) / dt
        call dt_evolve(sbe, gs, E, Ac)

        if (mod(it, out_rt_step) == 0) then
            write(*, '(i6,1x,f7.3,99(1x,e23.15e3))') &
            & it, t, Ac, E, calc_total_elec(sbe, gs, sbe%nb/2)
        end if
    end do

    stop
end program 
