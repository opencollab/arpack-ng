#!/bin/bash -eu

./arpackmm --A issue401.mtx --mag LA --nbEV 1 --nbCV 5
if [ "$?" -ne "0" ]; then exit 1; fi
./arpackmm --A issue401.mtx --mag LA --nbEV 1 --nbCV 5 --restart
if [ "$?" -ne "0" ]; then exit 1; fi
./arpackmm --A issue401.mtx --mag LA --nbEV 1 --nbCV 5 --restart --shift 0.9
if [ "$?" -ne "0" ]; then exit 1; fi

echo "OK"
