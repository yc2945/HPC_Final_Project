#include <iostream>
#include <vector>
#include <mpi.h>
#include <math.h>

using namespace std;

const int seed = 2019;
const int iterationCount = 3;
const int gridSize = 8;

int *grid;

// We define 0 as dead, 1 as alive

void runTick(int piece) {
    int *newGrid = new int[piece * piece];

    for (int i = 0; i < piece * piece; i++) {
        int liveCount = 0;
        int row = i / piece;
        int col = i % piece;

        for (int m = -1; m <= 1; m++) {
            for (int n = -1; n <= 1; n++) {
                if (m == 0 && n == 0) continue;
                int neighborRow = row + m;
                int neighborCol = col + n;

                if (neighborRow < 0 || neighborCol >= piece || neighborCol < 0 || neighborCol >= piece) continue;

                if (grid[neighborRow * piece + neighborCol] == 1) liveCount++;
            }
        }

        newGrid[i] = grid[i];
        if (grid[i] == 1 && (liveCount < 2 || liveCount > 3)) {
            newGrid[i] = 0;
        }
        if (grid[i] == 0 && liveCount == 3) {
            newGrid[i] = 1;
        }
    }

    grid = newGrid;
}

void printGrid(int iteration, int piece, int rank) {
    printf("Process %d: \n", rank);
    printf("Current Iteration: %d \n", iteration);
    for (int i = 0; i < piece; i++) {
        for (int j = 0; j < piece; j++) {
            printf("%d ", grid[i * piece + j]);
        }
        printf("\n");
    }
    printf("\n");
}

// void fillcube(subcube, piece)

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);
    
    int rank;
    MPI_Comm comm = MPI_COMM_WORLD;
    MPI_Comm_rank(comm, &rank);
    int world_size;
    MPI_Comm_size(comm, &world_size);
    MPI_Status status;

    // calculate the length of each subcube
    
    if (fabs((int)sqrt(world_size) - sqrt(world_size)) > 0.00001){
        printf("Please use nodes number whose square root is an integer. \n \n");
        abort();
    }
    int piece;
    piece = (int)(gridSize / (int)sqrt(world_size));


    //initiate each subcube
    srand(seed+rank);
    grid = new int[piece * piece];
    for (int i = 0; i < piece; i++) {
        for (int j = 0; j < piece; j++) {
            grid[i * piece + j] = rand() % 2;
        }
    }

    //update each subcube
    for (int i = 0; i < iterationCount; i++) {
        for (int j = 0; j < world_size; j++){
            if (rank == j)
                printGrid(i, piece, rank);
            MPI_Barrier(comm);
        }
        runTick(piece);
    }

    // master node gather subcubes 
    // if rank != 0 {
    //     MPI_Send(grid, piece * piece, MPI_INT, 0, rank, comm);
    // }
    // else{

    //     for (int j = 0; j < world_size; j++){
    //         MPI_Recv(grid, piece * piece, MPI_INT, j, j, comm, &status);
    //     }
    // }

//     MPI Gather(grid, piece * piece, MPI_INT, grid, piece * piece, 0, comm);
//     (void* sendbuff, int sendcount, MPI Datatype
// sendtype, void* recvbuf, int recvcount, MPI Datatype
// recvtype, int root, MPI Comm communicator)    
//     printGrid(iterationCount - 1, piece, rank);

    delete[] grid;
    
    MPI_Finalize();
    return 0;
}