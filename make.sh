#!/bin/sh

cd src
gcc -Wall analysis.c memutils.c mathutils.c dsp.c -o ../analysis
