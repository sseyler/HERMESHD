! Module hermeshd defined in file hermeshd.f90

subroutine f90wrap_main(comm)
    use hermeshd, only: main
    implicit none
    
    integer, intent(inout) :: comm
    call main(comm=comm)
end subroutine f90wrap_main

! End of module hermeshd defined in file hermeshd.f90

