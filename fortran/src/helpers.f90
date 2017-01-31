module helpers

use parameters

real loc_lxd
real loc_lyd
real loc_lzd

contains

    !-----------------------------------------------------------
    !   Return the x coordinate of (the center of) cell i
    !     Note: based on the location of this MPI domain (loc_lxd)
    real function xc(i)
        integer i
        xc = loc_lxd + (i - 0.5)*dx
    end function xc

    !-----------------------------------------------------------

    real function yc(j)
        integer j
        yc = loc_lyd + (j - 0.5)*dy
    end function yc

    !-----------------------------------------------------------

    real function zc(k)
        integer k
        zc = loc_lzd + (k - 0.5)*dz
    end function zc

    !-----------------------------------------------------------

    real function rz(i,j)
        integer i,j
        rz = sqrt(yc(j)**2)
    end function rz

    !-----------------------------------------------------------

    real function r(i,j)
        integer i,j,k
        r = sqrt(xc(i)**2 + yc(j)**2)
    end function r

    !-----------------------------------------------------------

    real function theta(i,j)
        integer i,j
        theta = atan2(yc(j),xc(i))
    end function theta

    !-----------------------------------------------------------

    real function xvtk(i)
        integer i
        xvtk = loc_lxd + (i - 0.5)*dxvtk
    end function xvtk

    !-----------------------------------------------------------

    real function yvtk(j)
        integer j
        yvtk = loc_lyd + (j - 0.5)*dyvtk
    end function yvtk

    !-----------------------------------------------------------

    real function zvtk(k)
        integer k
        zvtk = loc_lzd + (k - 0.5)*dzvtk
    end function zvtk

    !-----------------------------------------------------------

    real function rvtk(i,j)
        integer i,j,k
        rvtk = sqrt(xvtk(i)**2 + yvtk(j)**2)
    end function rvtk

    !-----------------------------------------------------------

    real function thetavtk(i,j)
        integer i,j
        thetavtk = atan2(yvtk(j),xvtk(i))
    end function thetavtk

    !-----------------------------------------------------------


    !-----------------------------------------------------------
    ! Get current time from system_clock
    real function get_clock_time()
        integer :: ticks, count_rate, count_max
        call system_clock( ticks, count_rate, count_max )
        get_clock_time = 1.0 * ticks / count_rate
    end function get_clock_time

    !-----------------------------------------------------------
    ! Print a message from the MPI rank with ID mpi_id
    subroutine mpi_print(mpi_id, message)
        integer, intent(in) :: mpi_id
        character(*) :: message
        if (mpi_id .eq. print_mpi) then
            print *,message
            print *,''       ! print a new line after message by default
        endif
    end subroutine mpi_print


    !-----------------------------------------------------------
    ! Print a message from the MPI rank with ID mpi_id
    subroutine whoami_print(mpi_id, message)
        integer, intent(in) :: mpi_id
        character(*) :: message
        print *,'Rank',mpi_id, ':', message
        print *,''       ! print a new line after message by default
    end subroutine whoami_print


    !------------------------------------------
    subroutine init_random_seed(iseed)

        implicit none
        integer :: i, n, clock, iseed
        integer, dimension(:), allocatable :: seed

        call random_seed(size = n)
        allocate(seed(n))

        call system_clock(count=clock)

        if (iseed .eq. 0) then
            seed =  clock*(iam+1) + 37 * (/ (i - 1, i = 1, n) /)
        else
            seed =  iseed*(iam+1)
        endif
        call random_seed(put = seed)

        ! print *,seed(1)
        deallocate(seed)

    end subroutine

end module helpers
