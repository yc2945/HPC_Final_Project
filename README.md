# HPC_Final_Project

To run viasualization of serial version on mac: 
$ vcpkg install opengl freeglut
$ g++ visualization.cpp Game.cpp Game.h -framework OpenGL -framework GLUT
$ ./a.out