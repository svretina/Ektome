#!/usr/bin/zsh

cython cython_funcs.pyx
gcc -shared -pthread -fPIC -fwrapv -O2 -Wall -fno-strict-aliasing -I/work/stamatis.vretinaris/opt/python-3.8.7/include/python3.8/  -o cython_funcs.so cython_funcs.c
