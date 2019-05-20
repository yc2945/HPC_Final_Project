#include <OpenGL/gl.h>
#include <OpenGL/glu.h>
#include <GLUT/glut.h>
#include<math.h>
#include <iostream>
#include <algorithm>
#include "Game.h"

#define INITIAL_FPS 2;

Game *game;

const unsigned int WIDTH = 800;
const unsigned int HEIGHT = 800;
const unsigned int SIZE = 500;
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
            int colors = game->getEncodedColor(i, j); 
            int cr, cg, cb; 
            cr = colors / (1000 * 1000);
            cg = (colors - cr * (1000 * 1000)) / 1000;
            cb = colors - cr * (1000 * 1000) - cg * 1000;
            glColor3f((float) cr / 255, (float) cg / 255, (float) cb / 255); 

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

