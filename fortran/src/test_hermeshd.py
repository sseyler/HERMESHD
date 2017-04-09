from mpi4py import MPI
import hermeshd

comm = MPI.COMM_WORLD

rank = comm.Get_rank()
size = comm.Get_size()
print "Proc %d out of %d procs" % (rank,size)

# Get the communicator
fcomm = MPI.COMM_WORLD.py2f()


hermeshd.hermeshd.main(fcomm)
