# HPC_Final_Project


### Compile serial version visualization on mac

```sh
# install depedencies
$ vcpkg install opengl freeglut

# compile
$ g++ visualization.cpp Game.cpp Game.h -framework OpenGL -framework GLUT

# run
$ ./a.out
```