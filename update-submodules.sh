#!/usr/bin/bash
cd JupiterMag/__data/libjupitermag/
make clean
git stash
git pull origin main
cd lib/libcon2020
git stash
git pull origin main
cd ../libspline
git stash
git pull origin main
cd ../libinternalfield
git stash
git pull origin main
cd ../../../../../