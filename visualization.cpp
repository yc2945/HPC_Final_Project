#include <OpenGL/gl.h>
#include <OpenGL/glu.h>
#include <GLUT/glut.h>
#include<math.h>
#include <iostream>
#include <algorithm>
#include "Game.h"

#define INITIAL_FPS 2;
#define ALIVE_CELL_COLOR glColor3f(0.0f, 0.0f, 0.0f);
#define DEAD_CELL_COLOR glColor3f(1.0f, 1.0f, 1.0f);
#define ROUND_ZERO(a) ((a) - 1 < 0 ? 0 : (a) - 1)
#define ROUND_MAX(a, max) ((a) + 1 > (max) ? (max) : (a) + 1)

Game *game;

const unsigned int WIDTH = 800;
const unsigned int HEIGHT = 800;
const unsigned int SIZE = 100;
const unsigned int ROWS = SIZE;
const unsigned int COLUMNS = SIZE;
const float CELL_WIDTH = (float) WIDTH / (float) ROWS;
const float CELL_HEIGHT = (float) HEIGHT / (float) COLUMNS;

unsigned int FPS = INITIAL_FPS;
bool continue_calc = true;


void init(void) {
    glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluOrtho2D(0.0, 1.0, 0.0, 1.0);
    glMatrixMode(GL_MODELVIEW);
    game->initialize();
}


void renderFunction(void) {
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    for (unsigned int i = 0; i < ROWS; i++) {
        for (unsigned int j = 0; j < COLUMNS; j++) {
            if (game->getStatus(i, j) == 1)
                ALIVE_CELL_COLOR
            else
                DEAD_CELL_COLOR

            glRectf(
                    (CELL_WIDTH * j) / WIDTH,
                    (CELL_HEIGHT * i) / HEIGHT,
                    (CELL_WIDTH * (j + 1)) / WIDTH,
                    (CELL_HEIGHT * (i + 1)) / HEIGHT
            );

        }
    }

    if (continue_calc) {
        game->runTick();
    }
    glutSwapBuffers();
}

void timer(int) {
    glutPostRedisplay();
    glutTimerFunc(1000 / FPS, timer, 0);
}


int main(int argc, char **argv) {
    game = new Game(SIZE);

    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_SINGLE);
    glutInitWindowSize(WIDTH, HEIGHT);
    glutCreateWindow("Game Of Life");
    glutDisplayFunc(renderFunction);

    init();

    glutTimerFunc(1000 / FPS, timer, 0);

    glutMainLoop();
    return 0;
}

