//
// Created by Yichang Chen on 2019-05-18.
//

#include "Game.h"

Game::Game(int gridSize) {
    this->gridSize = gridSize;
    this->iterationCount = 0;
    for (int i = 0; i < gridSize; i++) {
        for (int j = 0; j < gridSize; j++) {
            this->grid.push_back(rand() % 2);
        }
    }
}

void Game::runTick() {
    std::vector<int> newGrid(gridSize * gridSize);

    for (int row = 0; row < gridSize; row++) {
        for (int col = 0; col < gridSize; col++) {
            int liveCount = 0;
            int index = row * gridSize + col;

            for (int m = -1; m <= 1; m++) {
                for (int n = -1; n <= 1; n++) {
                    if (m == 0 && n == 0) continue;
                    int neighborRow = row + m;
                    int neighborCol = col + n;

                    if (neighborRow < 0 || neighborCol >= gridSize || neighborCol < 0 || neighborCol >= gridSize) {
                        continue;
                    }

                    if (getStatus(neighborRow, neighborCol) == 1) liveCount++;
                }
            }

            newGrid[index] = grid[index];
            if (grid[index] == 1 && (liveCount < 2 || liveCount > 3)) {
                newGrid[index] = 0;
            }
            if (grid[index] == 0 && liveCount == 3) {
                newGrid[index] = 1;
            }
        }
    }

    grid = newGrid;
    iterationCount++;
}

void Game::printGame() {
    std::cout << "Iteration: " << iterationCount << std::endl;
    for (int i = 0; i < gridSize; i++) {
        for (int j = 0; j < gridSize; j++) {
            std::cout << getStatus(i, j) << " ";
        }
        std::cout << std::endl;
    }
}

int Game::getStatus(int row, int col) {
    return grid[row * gridSize + col];
}