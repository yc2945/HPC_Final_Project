#include <iostream>
#include <vector>
#include <mpi.h>
#include <math.h>

using namespace std;

const int seed = 2019;
const int iterationCount = 3;
const int gridSize = 9;

MPI_Status status;

MPI_Request request_out1, request_in1;
MPI_Request request_out2, request_in2;
MPI_Request request_out3, request_in3;
MPI_Request request_out4, request_in4;

// We define 0 as dead, 1 as alive

void runTick(int *grid, int piece) {

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

        if (grid[i] == 1 && (liveCount < 2 || liveCount > 3)) {
            grid[i] = 0;
        }
        if (grid[i] == 0 && liveCount == 3) {
            grid[i] = 1;
        }
    }

}

void printGrid(int *grid, int iteration, int piece, int rank) {
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

void printAllGrid(int *allgrid,int iteration) {
    for (int i = 0; i < gridSize; i++) {
        for (int j = 0; j < gridSize; j++) {
            printf("%d ", allgrid[i * gridSize + j]);
        }
        printf("\n");
    }
    printf("\n");
}

void fillcube(int *grid, int *allgrid, int rank, int rp, int piece){
    // e.g. rank = 4,rp = 3, piece = 2, then row = 2, col_start = 2
    int row_start = (int)(rank / rp * piece);
    int col_start = (int)(rank % rp) * piece;
    for (int i = 0; i < piece ; i++) {
        for (int j = 0; j < piece; j++) {
            allgrid[(row_start + i) * gridSize + col_start + j] = grid[i * piece + j];
        }
    }
}

void sendmargin(int *grid, int* top, int* bottom, int* left, int* right, int rank, int rp,\
 int piece, MPI_Comm comm)
{
    // e.g. rank = 4,rp = 3, then row_ind = 1, col_ind = 1.
    int row_ind = (int)(rank / rp);
    int col_ind = (int)(rank % rp);

    //send
    // not at the top
    if (row_ind != 0){
        for (int i=0;i<piece;i++) top[i] = grid[i];
        MPI_Isend(top, piece, MPI_INT, rank - rp, rank, comm, &request_out1);  
    }
    // not at the bottom
    if (row_ind != rp - 1){
        for (int i=0;i<piece;i++) bottom[i] = grid[piece * (rp - 1) + i];
        MPI_Isend(bottom, piece, MPI_INT, rank + rp, rank, comm, &request_out2);  
    }
    // //not at the left side
    // if (col_ind != 0){
    //     for (int i=0;i<piece;i++) left[i] = grid[piece * i];
    //     MPI_Isend(left, piece, MPI_INT, rank - 1, rank, comm, &request_out3);  
    // }
    // //not at the right side
    // if (col_ind != rp - 1){
    //     for (int i=0;i<piece;i++) right[i] = grid[piece * i + piece - 1];
    //     MPI_Isend(right, piece, MPI_INT, rank + 1, rank, comm, &request_out4);  
    // }


    //receive
    // not at the top, receive info from the grid above, bottom here is the part above the grid
    if (row_ind != 0){
        MPI_Irecv(top, piece, MPI_INT, rank + rp, rank + rp, comm, &request_in1);
    }
    // not at the bottom, receive info from the grid below, top here is the part below the grid
    if (row_ind != rp - 1){
        MPI_Irecv(bottom, piece, MPI_INT, rank - rp, rank - rp, comm, &request_in2);  
    }
    // // not at the right side, receive info from the grid right, left here is the part right of the grid
    // if (col_ind != rp - 1){
    //     MPI_Irecv(left, piece, MPI_INT, rank + 1, rank + 1, comm, &request_in3);  
    // }
    // // not at the left side, receive info from the grid left, right here is the part left of the grid
    // if (col_ind != 0){
    //     MPI_Irecv(right, piece, MPI_INT, rank - 1, rank - 1, comm, &request_in4);  
    // }
    printf("start waiting");
    if (row_ind != 0){
        MPI_Wait(&request_out1, &status);
        MPI_Wait(&request_in1, &status);
    }
    printf("top done");
    if (row_ind != rp - 1){
        MPI_Wait(&request_out2, &status);
        MPI_Wait(&request_in2, &status);
    }
    printf("bottom done");
    // if (col_ind != 0){
    //     MPI_Wait(&request_out3, &status);
    //     MPI_Wait(&request_in4, &status);
    //}
    // printf("left done");
    // if (col_ind != rp - 1){
    //     MPI_Wait(&request_out4, &status);
    //     MPI_Wait(&request_in3, &status);
    // }
    // printf("right done");
    
    // if (row_ind != rp - 1){
    //     for (int i=0;i<piece;i++) printf("rank = %d, top = %d\n", rank, top[i]);
    // }
    if (row_ind != rp - 1){
        for (int i=0;i<piece;i++) printf("rank = %d, bottom = %d\n", rank, bottom[i]);
    }
}



int main(int argc, char** argv) {

    MPI_Init(&argc, &argv);
    
    int rank;
    MPI_Comm comm = MPI_COMM_WORLD;
    MPI_Comm_rank(comm, &rank);
    int world_size;
    MPI_Comm_size(comm, &world_size);

    int name_len;
    char processor_name[MPI_MAX_PROCESSOR_NAME];
    MPI_Get_processor_name(processor_name, &name_len);
    printf("Rank %d/%d running on %s.\n", rank, world_size, processor_name);

    // calculate the rp of each subcube
    int rp = (int)sqrt(world_size);
    if (fabs(rp - sqrt(world_size)) > 0.00001){
        printf("Please use nodes number whose square root is an integer. \n \n");
        abort();
    }
    int piece;
    piece = (int)(gridSize / rp);


    //initiate each subcube
    srand(seed+rank);
    int* grid = (int*) malloc(piece * piece * sizeof(int));
    for (int i = 0; i < piece; i++) {
        for (int j = 0; j < piece; j++) {
            grid[i * piece + j] = rand() % 2;
        }
    }


    int* top = (int*) malloc(piece * sizeof(int));
    int* bottom = (int*) malloc(piece * sizeof(int));
    int* left = (int*) malloc(piece * sizeof(int));
    int* right = (int*) malloc(piece * sizeof(int));

    //update each subcube
    for (int i = 0; i < iterationCount; i++) {

        for (int j = 0; j < world_size; j++){
            if (rank == j)
                printGrid(grid, i, piece, rank);
            MPI_Barrier(comm);
        }

        sendmargin(grid, top, bottom, left, right, rank, rp, piece, comm);
        MPI_Barrier(comm);
        runTick(grid, piece);

    }
    MPI_Barrier(comm);
    // master node gather subcubes

    
    if (rank != 0) {
        MPI_Send(grid, piece * piece, MPI_INT, 0, rank, comm);
    }
    else{
        int* allgrid = (int*) malloc(gridSize * gridSize * sizeof(int));
        for (int i = 0; i < gridSize; i++) {
            for (int j = 0; j < gridSize; j++) {
            allgrid[i * gridSize + j] = 2;
            }
        }

        fillcube(grid, allgrid, rank, rp, piece);
        // printf("first cube filled");
        for (int j = 1; j < world_size; j++){
            free(grid);
            grid = (int*) malloc(piece * piece * sizeof(int));
            MPI_Recv(grid, piece * piece, MPI_INT, j, j, comm, &status);
            fillcube(grid, allgrid, j, rp, piece);
        }
        
        printAllGrid(allgrid, iterationCount);   
        free(allgrid); 
    }   
  
    MPI_Barrier(comm);
    free(grid);
    free(top);
    free(bottom);
    free(left);
    free(right);
    MPI_Finalize();
    return 0;
}