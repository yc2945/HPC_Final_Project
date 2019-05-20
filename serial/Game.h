//
// Created by Yichang Chen on 2019-05-18.
//

#define STATUS_DEAD 0;
#define STATUS_ALIVE 1;


#ifndef PROJECT_GAME_H
#define PROJECT_GAME_H

#include <vector>
#include <iostream>
#include <string> 
#include <fstream>
#include <math.h>

const int seed = 2019;


class Game {
private:
    int iterationCount;
    int gridSize;
    std::vector<int> grid;

public:
    Game(int gridSize);

    void initialize();
    void runTick();
    void printGame();
    void storeGrid(); 
    int getEncodedColor(int row, int col);
    std::vector<int> decodeColor(int color); 
    int encodeColor(int c1, int c2, int c3); 
};

#endif //PROJECT_GAME_H