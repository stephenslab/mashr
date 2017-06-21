#!/usr/bin/env sh

mkdir -p ~/.R

echo "CC=gcc-4.9"    >> ~/.R/Makevars
echo "CXX=g++-4.9"   >> ~/.R/Makevars
echo "CXX11=g++-4.9" >> ~/.R/Makevars
