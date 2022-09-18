#!/bin/bash

g++ BPP.cpp -o BPP


for f in ../test/EZ/*.BPP; do echo $f; ./BPP $f ../res/resEZ.csv -b; done
for f in ../test/EZ/*.BPP; do echo $f; ./BPP $f ../res/resEZ.csv; done
for f in ../test/N1/*.BPP; do echo $f; ./BPP $f ../res/resN1.csv; done
for f in ../test/N2/*.BPP; do echo $f; ./BPP $f ../res/resN2.csv; done
for f in ../test/H/*.BPP; do echo $f; ./BPP $f ../res/resH.csv; done

