from mpi4py import MPI
import hermeshd

me = MPI.COMM_WORLD.Get_rank()
nprocs = MPI.COMM_WORLD.Get_size()
print "Proc %d out of %d procs" % (me,nprocs)
hermeshd.hermeshd.main()
