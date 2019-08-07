module sbe_gs
    use salmon_math, only: pi
    implicit none

    integer, parameter :: ncycle = 2

    type s_sbe_gs
        !Lattice information
        real(8) :: a_matrix(1:3, 1:3)
        real(8) :: b_matrix(1:3, 1:3)
        real(8) :: volume

        !Ground state (GS) electronic system information
        integer :: nk, nb, ne
        real(8), allocatable :: kvec(:, :), kweight(:)
        real(8), allocatable :: eigen(:, :)
        real(8), allocatable :: occup(:, :)
        real(8), allocatable :: omega(:, :, :)
        complex(8), allocatable :: p_matrix(:, :, :, :)
        complex(8), allocatable :: d_matrix(:, :, :, :)
        complex(8), allocatable :: rv_matrix(:, :, :)
        complex(8), allocatable :: prod_dk(:, :, :, :, :, :)

        !k-space grid and geometry information
        !NOTE: prepred for uniformally distributed k-grid....
        integer :: nkgrid(1:3)
        integer, allocatable :: iktbl_grid(:, :, :)
        integer, allocatable :: ikcycle_tbl1(:)
        integer, allocatable :: ikcycle_tbl2(:)
        integer, allocatable :: ikcycle_tbl3(:)
    end type


contains


subroutine init_sbe_gs(gs, sysname, directory, nkgrid, nb, ne, a1, a2, a3, read_bin)
    use salmon_file, only: open_filehandle, get_filehandle
    implicit none
    type(s_sbe_gs), intent(inout) :: gs
    character(*), intent(in) :: sysname
    character(*), intent(in) :: directory
    integer, intent(in) :: nkgrid(1:3)
    integer, intent(in) :: nb
    integer, intent(in) :: ne
    real(8), intent(in) :: a1(1:3), a2(1:3), a3(1:3)
    logical, intent(in) :: read_bin
    integer :: nk

    nk = nkgrid(1) * nkgrid(2) * nkgrid(3)

    gs%nk = nk
    gs%nb = nb
    gs%ne = ne
    gs%nkgrid(1:3) = nkgrid(1:3)

    !Calculate b_matrix, volume_cell and volume_bz from a1..a3 vector.
    call calc_lattice_info()

    allocate(gs%kvec(1:3, 1:nk))
    allocate(gs%kweight(1:nk))
    allocate(gs%eigen(1:nb, 1:nk))
    allocate(gs%occup(1:nb, 1:nk))
    allocate(gs%omega(1:nb, 1:nb, 1:nk))
    allocate(gs%p_matrix(1:nb, 1:nb, 1:3, 1:nk))
    allocate(gs%d_matrix(1:nb, 1:nb, 1:3, 1:nk))
    allocate(gs%rv_matrix(1:nb, 1:nb, 1:nk))
    allocate(gs%prod_dk(1:nb, 1:nb, -1:1, -1:1, -1:1, 1:nk))
    allocate(gs%iktbl_grid(1:nkgrid(1),1:nkgrid(2),1:nkgrid(3)))
    allocate(gs%ikcycle_tbl1(1-nkgrid(1)*ncycle:nkgrid(1)*(ncycle+1)))
    allocate(gs%ikcycle_tbl2(1-nkgrid(2)*ncycle:nkgrid(2)*(ncycle+1)))
    allocate(gs%ikcycle_tbl3(1-nkgrid(3)*ncycle:nkgrid(3)*(ncycle+1)))

    if (read_bin) then
        !Retrieve all data from binray
        write(*,*) "read_sbe_gs_bin"
        call read_sbe_gs_bin()
    else
        !Retrieve eigenenergies from 'SYSNAME_eigen.data':
        write(*,*) "read_eigen_data"
        call read_eigen_data()
        !Retrieve k-points from 'SYSNAME_k.data':
        write(*,*) "read_k_data"
        call read_k_data()
        !Retrieve transition matrix from 'SYSNAME_tm.data':
        write(*,*) "read_tm_data"
        call read_tm_data()
        !Retrieve transition matrix from 'SYSNAME_tm.data':
        write(*,*) "read_prod_dk"
        call read_prod_dk()
        !Export all data from binray
        write(*,*) "save_sbe_gs_bin"
        call save_sbe_gs_bin()
    end if

    !Calculate iktbl_grid for uniform (non-symmetric) k-grid:
    write(*,*) "create_uniform_iktbl_grid"
    call create_uniform_iktbl_grid()
    stop
    !Calculate omega and d_matrix (neglecting diagonal part):
    write(*,*) "create_omega_d"
    call create_omega_d()

    !Initial Occupation Number
    gs%occup(1:(ne/2),:) = 2d0 !!Experimental!!

contains

    ! Calculate lattice and reciprocal vectors
    subroutine calc_lattice_info()
        implicit none
        real(8) :: a12(1:3), a23(1:3), a31(1:3), vol
        real(8) :: b1(1:3), b2(1:3), b3(1:3)

        a12(1) = a1(2) * a2(3) - a1(3) * a2(2)
        a12(2) = a1(3) * a2(1) - a1(1) * a2(3)
        a12(3) = a1(1) * a2(2) - a1(2) * a2(1)
        a23(1) = a2(2) * a3(3) - a2(3) * a3(2)
        a23(2) = a2(3) * a3(1) - a2(1) * a3(3)
        a23(3) = a2(1) * a3(2) - a2(2) * a3(1)
        a31(1) = a3(2) * a1(3) - a3(3) * a1(2)
        a31(2) = a3(3) * a1(1) - a3(1) * a1(3)
        a31(3) = a3(1) * a1(2) - a3(2) * a1(1)
        vol = dot_product(a12, a3)
        b1(1:3) = (2d0 * pi / vol) * a23(1:3)
        b2(1:3) = (2d0 * pi / vol) * a31(1:3)
        b3(1:3) = (2d0 * pi / vol) * a12(1:3)
        
        gs%a_matrix(1:3, 1) = a1(1:3)
        gs%a_matrix(1:3, 2) = a2(1:3)
        gs%a_matrix(1:3, 3) = a3(1:3)    
        gs%b_matrix(1, 1:3) = b1(1:3)
        gs%b_matrix(2, 1:3) = b2(1:3)
        gs%b_matrix(3, 1:3) = b3(1:3)
        gs%volume = vol
    end subroutine calc_lattice_info


    ! Read k-point coordinates from SALMON's output file
    subroutine read_k_data()
        implicit none
        character(256) :: dummy
        integer :: fh, i, ik, iik
        real(8) :: tmp(4)

        fh = open_filehandle(trim(directory) // trim(sysname) // '_k.data', 'old')
        do i=1, 8
            read(fh, '(a)') dummy !Skip
        end do
        do iik=1, nk
            read(fh, *) ik, tmp(1:4)
            if (ik .ne. iik) stop "ERROR! Invalid SYSNAME_k.data"
            gs%kvec(1:3, ik) = tmp(1:3)
            gs%kweight(ik) = tmp(4)
        end do !iik
        close(fh)
    end subroutine read_k_data


    ! Read eigenvalue data from SALMON's output file
    subroutine read_eigen_data()
        implicit none
        character(256) :: dummy
        integer :: fh, i, iik, iib, ib
        real(8) :: tmp

        fh = open_filehandle(trim(directory) // trim(sysname) // '_eigen.data', 'old')
        do i=1, 3
            read(fh, '(a)') dummy !Skip
        end do
        do iik=1, nk
            read(fh, '(a)') dummy !Skip
            do iib=1, nb
                read(fh, *) ib, tmp
                if (ib .ne. iib) stop "ERROR! Invalid SYSNAME_eigen.data"
                gs%eigen(ib, iik) = tmp
            end do
        end do
        close(fh)
    end subroutine read_eigen_data




    ! Read transition dipole moment from SALMON's output file
    subroutine read_tm_data()
        implicit none
        character(256) :: dummy
        integer :: fh, i, ik, ib, jb, iik, iib, jjb
        real(8) :: tmp(1:6)

        fh = open_filehandle(trim(directory) // trim(sysname) // '_tm.data', 'old')
        do i=1, 3
            read(fh, '(a)') dummy !Skip
        end do
        do iik=1, nk
            do iib=1, nb
                do jjb=1, nb
                    read(fh, *) ik, ib, jb, tmp(1:6)
                    if ((ik .ne. iik) .or. (ib .ne. iib) .or. (jb .ne. jjb)) &
                        stop "ERROR! Invalid SYSNAME_tm.data"
                    gs%p_matrix(ib, jb, 1, ik) = dcmplx(tmp(1), tmp(2))
                    gs%p_matrix(ib, jb, 2, ik) = dcmplx(tmp(3), tmp(4))
                    gs%p_matrix(ib, jb, 3, ik) = dcmplx(tmp(5), tmp(6))
                end do !jb
            end do !ib
        end do !ik
        close(fh)
    end subroutine read_tm_data


    ! Read band interconnection from SALMON's output file
    subroutine read_prod_dk()
        implicit none
        character(256) :: dummy
        integer :: fh, i, ik, ib, jb, iik, iib, jjb
        integer :: ik1, ik2, ik3, jk1, jk2, jk3, jjk
        real(8) :: tmp(1:6)

        fh = open_filehandle(trim(directory) // trim(sysname) // '_prod_dk.data', 'old')
        read(fh, '(a)') dummy !Skip
        do iik=1, nk
            do jjk=1, 3*3*3
                do jjb=1, nb
                    do iib=1, nb
                        read(fh, *) ik, ik1, ik2, ik3, & 
                            & jk1, jk2, jk3, ib, jb, tmp(1:2)
                        if ((ik .ne. iik) .or. (ib .ne. iib) .or. (jb .ne. jjb)) &
                            stop "ERROR! Invalid SYSNAME_prod_dk.data"
                        gs%prod_dk(ib, jb, jk1, jk2, jk3, ik) = &
                            & dcmplx(tmp(1), tmp(2))
                    end do
                end do
            end do
        end do
        close(fh)
    end subroutine read_prod_dk


    subroutine read_sbe_gs_bin()
        implicit none
        integer :: fh

        fh = get_filehandle()
        open(fh, file=trim(directory) // trim(sysname) // '_sbe_gs.bin', form='unformatted', status='old')
        read(fh) gs%kvec
        read(fh) gs%kweight
        read(fh) gs%eigen
        read(fh) gs%p_matrix
        read(fh) gs%prod_dk
        close(fh)
        return
    end subroutine read_sbe_gs_bin


    subroutine save_sbe_gs_bin()
        implicit none
        integer :: fh

        fh = get_filehandle()
        open(fh, file=trim(directory) // trim(sysname) // '_sbe_gs.bin', form='unformatted', status='replace')
        write(fh) gs%kvec
        write(fh) gs%kweight
        write(fh) gs%eigen
        write(fh) gs%p_matrix
        write(fh) gs%prod_dk
        close(fh)
        return
    end subroutine save_sbe_gs_bin




    subroutine create_omega_d()
        implicit none
        integer :: ik, ib, jb
        real(8), parameter :: epsilon = 0.04
        complex(8), parameter :: zi = dcmplx(0d0, 1d0)
        do ik=1, nk
            do ib=1, nb
                do jb=1, nb
                    gs%omega(ib, jb, ik) = gs%eigen(ib, ik) - gs%eigen(jb, ik)
                    if (epsilon < abs(gs%omega(ib, jb, ik))) then
                        gs%d_matrix(ib, jb, 1:3, ik) = &
                            & zi * gs%p_matrix(ib, jb, 1:3, ik) / gs%omega(ib, jb, ik)
                    else
                        gs%d_matrix(ib, jb, 1:3, ik) = 0d0
                    end if
                end do
            end do
        end do
    end subroutine create_omega_d

    
    subroutine create_uniform_iktbl_grid()
        implicit none
        integer :: ik1, ik2, ik3, ik, icycle
        ik = 1
        do ik3=1, nkgrid(3)
            do ik2=1, nkgrid(2)
                do ik1=1, nkgrid(1)
                    gs%iktbl_grid(ik1, ik2, ik3) = ik
                    ik = ik + 1
                end do !ik1
            end do !ik2
        end do !ik3

        do ik1=1, nkgrid(1)
            do icycle = -ncycle, ncycle
                gs%ikcycle_tbl1(ik1 + icycle * nkgrid(1)) = ik1
            end do
        end do
        
        do ik2=1, nkgrid(2)
            do icycle = -ncycle, ncycle
                gs%ikcycle_tbl2(ik2 + icycle * nkgrid(2)) = ik2
            end do
        end do
        
        do ik3=1, nkgrid(3)
            do icycle = -ncycle, ncycle
                gs%ikcycle_tbl3(ik3 + icycle * nkgrid(3)) = ik3
            end do
        end do
    end subroutine create_uniform_iktbl_grid


end subroutine init_sbe_gs





end module sbe_gs

