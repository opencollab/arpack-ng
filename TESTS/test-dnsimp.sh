#!/bin/sh
pwd=$(pwd)
cd "$srcdir" && exec "$pwd/dnsimp"
