program main
    use sbe_solver
    use input_parameter
    use pulse
    use test
    implicit none

    type(s_sbe) :: sbe
    type(s_sbe_gs) :: gs
    real(8) :: t, E(3), Ac(3), jmat(3)
    integer :: it

    call read_input()

    ! Read ground state electronic system:
    call init_sbe_gs(gs, sysname, directory, &
        & num_kgrid_gs, nstate, nelec, &
        & al_vec1, al_vec2, al_vec3)
    
    call test_interp(gs); !stop
        
    ! Calculate dielectric spectra and save as SYSNAME_dielec.data:
    if (out_dielec == 'y') then
        call calc_dielec(trim(directory) // trim(sysname) // '_dielec.data', &
            & gs, e_min_dielec, e_max_dielec, n_dielec, gamma_dielec)
    end if

    ! Initialization of SBE solver and density matrix:
    call init_sbe(sbe, gs, num_kgrid_sbe)

    write(*, '("#",99(1x,a))') "1:Step", "2:Time[au]", "3:Ac_x", "4:Ac_y", "5:Ac_z", &
        & "6:E_x", "7:E_y", "8:E_z", "9:Jmat_x", "10:Jmat_y", "11:Jmat_z", "12:n_v", "13:n_all" 

    ! Realtime calculation
    do it = 1, nt
        t = dt * it
        call calc_cos2_pulse(t - dt * 0.5d0, pulse_tw1, &
            & rlaser_int_wcm2_1, omega1, phi_cep1, epdir_re1, epdir_im1, &
            & Ac, E)
        call dt_evolve(sbe, gs, E, Ac, dt)

        if (mod(it, out_rt_step) == 0) then
            call calc_cos2_pulse(t, pulse_tw1, &
                & rlaser_int_wcm2_1, omega1, phi_cep1, epdir_re1, epdir_im1, &
                & Ac, E)
            call calc_current_bloch(sbe, gs, Ac, jmat)
            write(*, '(i6,1x,f9.3,99(1x,es23.15e3))') &
            & it, t, Ac, E, jmat, calc_trace(sbe, sbe%nb/2), calc_trace(sbe, sbe%nb)
        end if
    end do

    stop
end program 
