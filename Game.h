//
// Created by Yichang Chen on 2019-05-18.
//

#define STATUS_DEAD 0;
#define STATUS_ALIVE 1;


#ifndef PROJECT_GAME_H
#define PROJECT_GAME_H

#include <vector>
#include <iostream>


const int seed = 2019;
//const int iterationCount = 10;
//const int gridSize = 20;


class Game {
private:
    int iterationCount;
    int gridSize;
    std::vector<int> grid;

public:
    Game(int gridSize);

    void runTick();
    void printGame();
    int getStatus(int row, int col);

};




#endif //PROJECT_GAME_H
