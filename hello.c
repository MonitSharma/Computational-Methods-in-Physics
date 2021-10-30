#include <mpi.h>
#include <stdio.h>

int main (int argc, char *argv[])   {
    int rank, size;

    MPI Init (&argc, &argv);
    // Initialising the MPI Library
    MPI_Comms_size(MPI_COMM_WORLD, &size);
    // It gets the number of processors 
    MPI_Comms_rank(MPI_COMM_WORLD, &rank);
    // It gets the process id and determine the process rank within
    
    // Enter your code here
    printf(" Hello World from the rank %d\n", rank);
    if (rank == 0) printf("MPI World Size = %d processes\n", size);
    MPI_Finalize();
    // MPI cleanup
    return 0;
    
}