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

void transform(int* biggrid, int*smallgrid, int piece){
    for (int i = 1; i < piece + 1; i++) {
        for (int j = 1; j < piece + 1; j++) {
            smallgrid[(i - 1) * piece + j - 1] = biggrid[i * (piece + 2)  + j];

        }
    }
}
void runTick(int *grid, int piece, int rank) {
    
    int *newGrid = (int*) malloc((piece + 2) * (piece + 2) * sizeof(int));
    int liveCount;
    int row, col, val, ind;

    //for (int i = 1; i < (piece  + 1) * (piece  + 1); i++) 
    for (int i = 0; i < piece + 2; i++) {
        for (int j = 0; j < piece + 2; j++) {
            if (i == 0 || i == piece + 1 || j == 0 || j == piece + 1 )
                val = 0;

            else{
                ind = i * (piece + 2) + j;
                liveCount = 0;
                row = i;
                col = j;

                for (int m = -1; m <= 1; m++) {
                    for (int n = -1; n <= 1; n++) {
                        if (m == 0 && n == 0) continue;
                        int neighborRow = row + m;
                        int neighborCol = col + n;

                        if (neighborRow < 0 || neighborCol >= piece || neighborCol < 0 || neighborCol >= piece) continue;

                        if (grid[neighborRow * (piece + 2) + neighborCol] == 1) liveCount++;
                    }
                }

                val = grid[ind];
                if (grid[ind] == 1 && (liveCount < 2 || liveCount > 3)) {
                    val = 0;
                }
                if (grid[ind] == 0 && liveCount == 3) {
                    val = 1;
                }
            }
            newGrid[ind] = val;
        }
    }
    for (int i = 0; i < (piece+2) * (piece+2); i++) {grid[i] = newGrid[i];}
    // grid = newGrid;
}

void printGrid(int *grid, int iteration, int piece, int rank) {
    printf("Process %d: \n", rank);
    printf("Current Iteration: %d \n", iteration);
    for (int i = 0; i < piece + 2; i++) {
        for (int j = 0; j < piece + 2; j++) {
            printf("%d ", grid[i * (piece + 2) + j]);
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
    int* temp  = (int*) malloc(piece * piece * sizeof(int));
    transform(grid, temp, piece);
    grid = temp;
    int row_start = (int)(rank / rp * piece);
    int col_start = (int)(rank % rp) * piece;
    for (int i = 0; i < piece ; i++) {
        for (int j = 0; j < piece; j++) {
            allgrid[(row_start + i) * gridSize + col_start + j] = grid[i * piece + j];
        }
    }
    free(temp);
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
        for (int i=0;i<piece;i++) top[i] = grid[piece + 2 + 1 + i];
        MPI_Isend(top, piece, MPI_INT, rank - rp, rank, comm, &request_out1);  
    }
    // not at the bottom
    if (row_ind != rp - 1){
        for (int i=0;i<piece;i++) bottom[i] = grid[(piece + 2) * (rp - 1 + 1) + 1 + i];
        MPI_Isend(bottom, piece, MPI_INT, rank + rp, rank, comm, &request_out2);  
    }
    //not at the left side
    if (col_ind != 0){
        for (int i=0;i<piece;i++) left[i] = grid[(piece + 2) * (i + 1) + 1];
        MPI_Isend(left, piece, MPI_INT, rank - 1, rank, comm, &request_out3);  
    }
    //not at the right side
    if (col_ind != rp - 1){
        for (int i=0;i<piece;i++) right[i] = grid[(piece + 2) * (i + 1) + piece];
        MPI_Isend(right, piece, MPI_INT, rank + 1, rank, comm, &request_out4);  
    }


    //receive
    // not at the bottom, receive info from the grid below, top here is the part below the grid
    if (row_ind != rp - 1){
        MPI_Irecv(top, piece, MPI_INT, rank + rp, rank + rp, comm, &request_in1);
    }
    // not at the top, receive info from the grid above, bottom here is the part above the grid
    if (row_ind != 0){
        MPI_Irecv(bottom, piece, MPI_INT, rank - rp, rank - rp, comm, &request_in2);  
    }

    // not at the right side, receive info from the grid right, left here is the part right of the grid
    if (col_ind != rp - 1){
        MPI_Irecv(left, piece, MPI_INT, rank + 1, rank + 1, comm, &request_in3);  
    }
    // not at the left side, receive info from the grid left, right here is the part left of the grid
    if (col_ind != 0){
        MPI_Irecv(right, piece, MPI_INT, rank - 1, rank - 1, comm, &request_in4);  
    }
    if (row_ind != 0){
        MPI_Wait(&request_out1, &status);
        MPI_Wait(&request_in2, &status);
    }
    if (row_ind != rp - 1){
        MPI_Wait(&request_out2, &status);
        MPI_Wait(&request_in1, &status);
    }
    if (col_ind != rp - 1){
        MPI_Wait(&request_out4, &status);
        MPI_Wait(&request_in3, &status);
    }
    if (col_ind != 0){
        MPI_Wait(&request_out3, &status);
        MPI_Wait(&request_in4, &status);
    }
    if (row_ind != rp - 1){
        for (int i = 0; i < piece; i++) {grid[(piece + 2) * (rp + 1) + 1 + i] = top[i];}
    }
    if (row_ind != 0){
        for (int i = 0; i < piece; i++) {grid[1 + i] = bottom[i];}
    }
    if (col_ind != rp - 1){
        for (int i = 0; i < piece; i++) {grid[(piece + 2) * (1 + i) + piece + 1] = left[i];}
    }
    if (col_ind != 0){
        for (int i = 0; i < piece; i++) {grid[(piece + 2) * (1 + i)] = right[i];}
    }


    // if (row_ind != rp - 1){
    //     for (int i=0;i<piece;i++) printf("rank = %d, top = %d\n", rank, top[i]);
    // }
    // if (row_ind != 0){
    //     for (int i=0;i<piece;i++) printf("rank = %d, bottom = %d\n", rank, bottom[i]);
    // }
    // if (col_ind != rp - 1){
    //     for (int i=0;i<piece;i++) printf("rank = %d, left = %d\n", rank, left[i]);
    // }
    // if (col_ind != 0){
    //     for (int i=0;i<piece;i++) printf("rank = %d, right = %d\n", rank, right[i]);
    // }
}

// master node gather subcubes


void gather(int *allgrid, int *grid, int rank, int piece, int rp, int world_size, MPI_Comm comm){
    if (rank != 0) {
        MPI_Send(grid, (piece + 2) *(piece + 2), MPI_INT, 0, rank, comm);
    }
    else{
        fillcube(grid, allgrid, rank, rp, piece);
        int *other_grid = (int*) malloc((piece + 2) * (piece + 2) * sizeof(int));
        for (int j = 1; j < world_size; j++){
            
            MPI_Recv(other_grid, (piece + 2) * (piece + 2), MPI_INT, j, j, comm, &status);

            fillcube(other_grid, allgrid, j, rp, piece);
        }        
        printAllGrid(allgrid, iterationCount);
        free(other_grid);

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

    // the length of the subcube
    int piece;
    piece = (int)(gridSize / rp);

    int* top = (int*) malloc(piece * sizeof(int));
    int* bottom = (int*) malloc(piece * sizeof(int));
    int* left = (int*) malloc(piece * sizeof(int));
    int* right = (int*) malloc(piece * sizeof(int));
    int* allgrid = (int*) malloc(gridSize * gridSize * sizeof(int));





    //initiate each subcube
    srand(seed+rank);
    int* grid = (int*) malloc((piece + 2) * (piece + 2) * sizeof(int));
    int val;
    for (int i = 0; i < piece + 2; i++) {
        for (int j = 0; j < piece + 2; j++) {
            if (i == 0 || i == piece + 1 || j == 0 || j == piece + 1 )
                val = 0;
            else
                val = rand() % 2;
            grid[i * (piece + 2) + j] = val;
        }
    }


    //update each subcube
    for (int i = 0; i < iterationCount; i++) {


        gather(allgrid, grid, rank, piece, rp, world_size, comm);
        sendmargin(grid, top, bottom, left, right, rank, rp, piece, comm);
        MPI_Barrier(comm);
        for (int j = 0; j < world_size; j++){
            if (rank == j)
                printGrid(grid, i, piece, rank);
            MPI_Barrier(comm);
        }
        runTick(grid, piece, rank);

    }
    MPI_Barrier(comm);
    gather(allgrid, grid, rank, piece, rp, world_size, comm);

  
    MPI_Barrier(comm);
    free(allgrid);
    free(grid);
    free(top);
    free(bottom);
    free(left);
    free(right);
    MPI_Finalize();
    return 0;
}