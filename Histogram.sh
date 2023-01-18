#!/bin/bash


RATIO=$1
FILE="/Users/diggs/Desktop/TOPCon/OUT-FILES/out-long-1-17-23/long-${RATIO}.dump"
FILE2="/Users/diggs/Desktop/TOPCon/OUT-FILES/out-density-1-13-23/density-topcon-${RATIO}.dump"
FILE3="/Users/diggs/Desktop/TOPCon/OUT-FILES/out-nuc-1-16-23/nuc-${RATIO}.dump"

python Atom_Manip.py $FILE2