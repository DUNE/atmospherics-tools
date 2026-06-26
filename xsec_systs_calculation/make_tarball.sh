#!/bin/bash
# Script to build and package xsec_systs_calculation for grid submission

set -e

# 1. Compiling ndnusyst and copying custom libraries inside SL7 container
singularity exec -B /cvmfs,/exp /cvmfs/singularity.opensciencegrid.org/fermilab/fnal-wn-sl7:latest bash -c "./build_with_reweight.sh"

# 2. Packaging into tarball: ../xsec_systs_calculation_wmec.tar.gz
cd ..

echo "Packaging tarball, excluding temporary build directories..."
tar --exclude='xsec_systs_calculation/.git' \
    --exclude='xsec_systs_calculation/build/CMakeFiles' \
    --exclude='xsec_systs_calculation/build/_deps' \
    --exclude='xsec_systs_calculation/build/app/CMakeFiles' \
    --exclude='xsec_systs_calculation/build/app/cmake_install.cmake' \
    --exclude='xsec_systs_calculation/build/app/Makefile' \
    --exclude='xsec_systs_calculation/build/CMakeCache.txt' \
    --exclude='xsec_systs_calculation/build/Makefile' \
    --exclude='xsec_systs_calculation/build/cmake_install.cmake' \
    --exclude='xsec_systs_calculation/Reweight/.git' \
    --exclude='xsec_systs_calculation/Reweight/src' \
    -czf xsec_systs_calculation_wmec.tar.gz xsec_systs_calculation/

echo "=== Tarball successfully packaged at: $(readlink -f xsec_systs_calculation_wmec.tar.gz) ==="
