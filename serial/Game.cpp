//
// Created by Yichang Chen on 2019-05-18.
//

#include "Game.h"

using namespace std; 

Game::Game(int gridSize) {
    this->gridSize = gridSize;
    this->iterationCount = 0;
}

void Game::initialize() {
    for (int i = 0; i < gridSize; i++) {
        for (int j = 0; j < gridSize; j++) {
            int val; 
            if (rand() % 2 == 0) {
                val = 0; 
            } else {
                int c1 = rand() % 255; 
                int c2 = rand() % 255; 
                int c3 = rand() % 255; 
                val = c1 * 1000 * 1000 + c2 * 1000 + c3; 
            }
            this->grid.push_back(val);
        }
    }
}

void Game::runTick() {
    // printGame(); 

    std::vector<int> newGrid(gridSize * gridSize);

    for (int row = 0; row < gridSize; row++) {
        for (int col = 0; col < gridSize; col++) {
            int liveCount = 0;
            int index = row * gridSize + col;
            int crSum, cgSum, cbSum; 

            for (int m = -1; m <= 1; m++) {
                for (int n = -1; n <= 1; n++) {
                    if (m == 0 && n == 0) continue;
                    int neighborRow = row + m;
                    int neighborCol = col + n;

                    if (neighborRow < 0 || neighborRow >= gridSize || neighborCol < 0 || neighborCol >= gridSize) {
                        continue;
                    }

                    int neighborColor = getEncodedColor(neighborRow, neighborCol); 
                    int cr = neighborColor / (1000 * 1000);
                    int cg = (neighborColor - cr * (1000 * 1000)) / 1000;
                    int cb = neighborColor - cr * (1000 * 1000) - cg * 1000;
                    crSum += cr; 
                    cgSum += cg; 
                    cbSum += cb; 

                    if (neighborColor > 0) liveCount++; 
                }
            }


            newGrid[index] = grid[index];
            if (liveCount < 2 || liveCount > 3) {
                newGrid[index] = 0;
            }
            if (grid[index] == 0 && liveCount == 3) {
                // crSum = ceil((double) crSum / (double) liveCount); 
                // cgSum = ceil((double) cgSum / (double) liveCount); 
                // cbSum = ceil((double) cbSum / (double) liveCount); 
                crSum /= liveCount; 
                cgSum /= liveCount; 
                cbSum /= liveCount; 
                newGrid[index] = crSum * 1000 * 1000 + cgSum * 1000 + cbSum; 
                // newGrid[index] = 1;
            }
        }
    }

    grid = newGrid;
    iterationCount++;
    // if (iterationCount == 5) {
        // storeGrid(); 
    // }
}

void Game::printGame() {
    std::cout << "Iteration: " << iterationCount << std::endl;
    for (int i = 0; i < gridSize; i++) {
        for (int j = 0; j < gridSize; j++) {
            std::cout << getEncodedColor(i, j) << " ";
        }
        std::cout << std::endl;
    }
}

int Game::getEncodedColor(int row, int col) {
    return grid[row * gridSize + col];
}

void Game::storeGrid() {
    std::string fileName = "serial_" + std::to_string(gridSize) + "_" + std::to_string(iterationCount) + ".txt"; 
    std::ofstream outfile(fileName.c_str());
    for (int i = 0; i < gridSize; i++) {
        for (int j = 0; j < gridSize; j++) {
            outfile << i << " " << j << " " << getEncodedColor(i, j) << std::endl; 
        }
    }
    outfile.close(); 
    exit(0); 
}

std::vector<int> decodeColor(int color) {
    int c3 = color % 1000; 
    color = color / 1000; 
    int c2 = color % 1000; 
    color = color / 1000; 
    int c1 = color; 
    std::vector<int> result; 
    result.push_back(c1); 
    result.push_back(c2); 
    result.push_back(c3); 
    return result; 
}

int encodeColor(int c1, int c2, int c3) {
    return c1 * 1000 * 1000 + c2 * 1000 + c3; 
}