source /cvmfs/dune.opensciencegrid.org/products/dune/setup_dune.sh

## cmake
setup cmake v3_21_4

## for fhicl
setup boost v1_80_0 -q e26:prof

## GENIE

setup genie v3_04_00d -q e26:prof
setup genie_xsec v3_04_00 -q AR2320i00000:e1000:k250

## CAFs
setup duneanaobj v03_03_00 -q e26:prof
