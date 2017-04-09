! Module hermeshd defined in file hermeshd.f90

subroutine f90wrap_main(comm)
    use hermeshd, only: main
    implicit none
    
    integer, intent(inout) :: comm
    call main(comm=comm)
end subroutine f90wrap_main

subroutine f90wrap_temp(a, b, c, d)
    use hermeshd, only: temp
    implicit none
    
    integer, intent(in) :: a
    integer, intent(inout) :: b
    real :: c
    integer, intent(out) :: d
    call temp(a=a, b=b, c=c, d=d)
end subroutine f90wrap_temp

subroutine f90wrap_step(q_io, q_1, q_2, t, dt, n0, n1, n2, n3, n4, n5, n6, n7, &
    n8, n9, n10, n11, n12, n13, n14)
    use hermeshd, only: step
    implicit none
    
    real, intent(inout), dimension(n0,n1,n2,n3,n4) :: q_io
    real, intent(inout), dimension(n5,n6,n7,n8,n9) :: q_1
    real, intent(inout), dimension(n10,n11,n12,n13,n14) :: q_2
    real, intent(inout) :: t
    real, intent(inout) :: dt
    integer :: n0
    !f2py intent(hide), depend(q_io) :: n0 = shape(q_io,0)
    integer :: n1
    !f2py intent(hide), depend(q_io) :: n1 = shape(q_io,1)
    integer :: n2
    !f2py intent(hide), depend(q_io) :: n2 = shape(q_io,2)
    integer :: n3
    !f2py intent(hide), depend(q_io) :: n3 = shape(q_io,3)
    integer :: n4
    !f2py intent(hide), depend(q_io) :: n4 = shape(q_io,4)
    integer :: n5
    !f2py intent(hide), depend(q_1) :: n5 = shape(q_1,0)
    integer :: n6
    !f2py intent(hide), depend(q_1) :: n6 = shape(q_1,1)
    integer :: n7
    !f2py intent(hide), depend(q_1) :: n7 = shape(q_1,2)
    integer :: n8
    !f2py intent(hide), depend(q_1) :: n8 = shape(q_1,3)
    integer :: n9
    !f2py intent(hide), depend(q_1) :: n9 = shape(q_1,4)
    integer :: n10
    !f2py intent(hide), depend(q_2) :: n10 = shape(q_2,0)
    integer :: n11
    !f2py intent(hide), depend(q_2) :: n11 = shape(q_2,1)
    integer :: n12
    !f2py intent(hide), depend(q_2) :: n12 = shape(q_2,2)
    integer :: n13
    !f2py intent(hide), depend(q_2) :: n13 = shape(q_2,3)
    integer :: n14
    !f2py intent(hide), depend(q_2) :: n14 = shape(q_2,4)
    call step(Q_io=q_io, Q_1=q_1, Q_2=q_2, t=t, dt=dt)
end subroutine f90wrap_step

subroutine f90wrap_setup(q_io, t, dt, t1, t_start, dtout, nout, comm, n0, n1, &
    n2, n3, n4)
    use hermeshd, only: setup
    implicit none
    
    real, intent(inout), dimension(n0,n1,n2,n3,n4) :: q_io
    real, intent(inout) :: t
    real, intent(inout) :: dt
    real, intent(inout) :: t1
    real, intent(inout) :: t_start
    real, intent(inout) :: dtout
    integer, intent(inout) :: nout
    integer, intent(inout) :: comm
    integer :: n0
    !f2py intent(hide), depend(q_io) :: n0 = shape(q_io,0)
    integer :: n1
    !f2py intent(hide), depend(q_io) :: n1 = shape(q_io,1)
    integer :: n2
    !f2py intent(hide), depend(q_io) :: n2 = shape(q_io,2)
    integer :: n3
    !f2py intent(hide), depend(q_io) :: n3 = shape(q_io,3)
    integer :: n4
    !f2py intent(hide), depend(q_io) :: n4 = shape(q_io,4)
    call setup(Q_io=q_io, t=t, dt=dt, t1=t1, t_start=t_start, dtout=dtout, &
        nout=nout, comm=comm)
end subroutine f90wrap_setup

subroutine f90wrap_cleanup(t_start)
    use hermeshd, only: cleanup
    implicit none
    
    real, intent(in) :: t_start
    call cleanup(t_start=t_start)
end subroutine f90wrap_cleanup

subroutine f90wrap_generate_output(q_r, t, dt, t1, nout, n0, n1, n2, n3, n4)
    use hermeshd, only: generate_output
    implicit none
    
    real, intent(in), dimension(n0,n1,n2,n3,n4) :: q_r
    real, intent(inout) :: t
    real, intent(inout) :: dt
    real, intent(inout) :: t1
    integer, intent(inout) :: nout
    integer :: n0
    !f2py intent(hide), depend(q_r) :: n0 = shape(q_r,0)
    integer :: n1
    !f2py intent(hide), depend(q_r) :: n1 = shape(q_r,1)
    integer :: n2
    !f2py intent(hide), depend(q_r) :: n2 = shape(q_r,2)
    integer :: n3
    !f2py intent(hide), depend(q_r) :: n3 = shape(q_r,3)
    integer :: n4
    !f2py intent(hide), depend(q_r) :: n4 = shape(q_r,4)
    call generate_output(Q_r=q_r, t=t, dt=dt, t1=t1, nout=nout)
end subroutine f90wrap_generate_output

! End of module hermeshd defined in file hermeshd.f90

