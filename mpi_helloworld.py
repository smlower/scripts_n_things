from mpi4py import MPI
import sys

#=========================================================#
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()
name = MPI.Get_processor_name()
n_tasks = size
this_task = rank
sys.stdout.write(
   "Hello, World! I am process %d of %d on %s.\n"
   % (rank, size, name))
#=========================================================#
