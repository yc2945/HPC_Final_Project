#include <iostream>
#include <vector>

#include "Game.h"

using namespace std;

//const int seed = 2019;
const int iterationCount = 10;
const int gridSize = 20;

int *grid;

// We define 0 as dead, 1 as alive


int getStatus(int row, int col) {
    return grid[row * gridSize + col]; 
}

void runTick() {
    int *newGrid = new int[gridSize * gridSize];

    for (int i = 0; i < gridSize * gridSize; i++) {
        int liveCount = 0;
        int row = i / gridSize;
        int col = i % gridSize;

        for (int m = -1; m <= 1; m++) {
            for (int n = -1; n <= 1; n++) {
                if (m == 0 && n == 0) continue;
                int neighborRow = row + m;
                int neighborCol = col + n;

                if (neighborRow < 0 || neighborCol >= gridSize || neighborCol < 0 || neighborCol >= gridSize) continue;

                if (grid[neighborRow * gridSize + neighborCol] == 1) liveCount++;
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

void printGrid(int iteration) {
    printf("Current Iteration: %d \n", iteration);
    for (int i = 0; i < gridSize; i++) {
        for (int j = 0; j < gridSize; j++) {
            printf("%d ", grid[i * gridSize + j]);
        }
        printf("\n");
    }
    printf("\n");
}


int main() {
    Game *game = new Game(10);
    for (int itr = 0; itr < 10; itr++) {
        game->runTick();
        game->printGame();
    }

    return 0;


//    grid = new int[gridSize * gridSize];
//
//    srand(seed);
//
//    for (int i = 0; i < gridSize; i++) {
//        for (int j = 0; j < gridSize; j++) {
//            grid[i * gridSize + j] = rand() % 2;
//        }
//    }
//
//
//    for (int i = 0; i < iterationCount; i++) {
//        runTick();
//        printGrid(i);
//    }
//
//    delete[] grid;
//    return 0;
}