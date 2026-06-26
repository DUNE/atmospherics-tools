source /cvmfs/dune.opensciencegrid.org/products/dune/setup_dune.sh

## cmake
setup cmake v3_21_4

## for fhicl
setup boost v1_80_0 -q e26:prof

## GENIE
setup genie v3_06_02e -q e26:prof
setup genie_xsec v3_06_00 -q AR2320i00000:e1000:k250

## CAFs
setup duneanaobj v03_15_00 -q e26:prof
setup sqlite v3_40_01_00

## SBN data
export PRODUCTS=$PRODUCTS:/cvmfs/sbn.opensciencegrid.org/products/sbn
setup sbndata v01_10
