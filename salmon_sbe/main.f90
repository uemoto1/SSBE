
program main
    use sbe_solver
    use input_parameter
    implicit none

    type(s_sbe) :: sbe
    type(s_sbe_gs) :: gs
    
    call read_input()

    call init_sbe_gs(gs, sysname, directory, &
        & num_kgrid, nstate, nelec, &
        & al_vec1, al_vec2, al_vec3)

    call init_sbe(sbe, gs, num_kgrid)
        
    stop
end program 