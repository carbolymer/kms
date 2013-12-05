#!/bin/bash
rm -f argon
g++ -O6 argon.cpp -o argon &
rm -f visualisation
g++ -O3 -lglut -lGL -lGLU visualisation.cpp -o visualisation &

