#!/bin/bash

qcc -O2 -Wall -disable-dimensions dropMove.c -o dropMove -lm
./dropMove