#!/bin/bash -eu

# Note: need to shift slightly as eigen value is zero.
./arpackmm --A issue215.mtx --mag LM --nbEV 1 --nbCV 4 --shiftReal 0.1
if [ "$?" -ne "0" ]; then exit 1; fi
./arpackmm --A issue215.mtx --mag LM --nbEV 1 --nbCV 4 --shiftReal 0.1 --restart
if [ "$?" -ne "0" ]; then exit 1; fi

echo "OK"
