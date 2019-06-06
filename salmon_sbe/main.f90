program main
    use sbe_solver
    use input_parameter
    implicit none

    type(s_sbe) :: sbe
    type(s_sbe_gs) :: gs
    real(8) :: E(3), Ac(3)

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
        E = 0d0
        Ac = 0d0
        call dt_evolve(sbe, gs, E, Ac)
        write(*,*) calc_total_elec(sbe, gs, sbe%nb/2)
    end do

    stop
end program 
