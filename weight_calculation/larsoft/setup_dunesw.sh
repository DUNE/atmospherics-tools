# Environment (SL7 + UPS) to both run larsoft/extract_larsoft.C on LArSoft files and build/run the weightor.
# osclib v00.29 (prof:e26:n319:stanfree) is built against the same ROOT/eigen as dunesw v10_09_00d00.
source /cvmfs/dune.opensciencegrid.org/products/dune/setup_dune.sh
setup dunesw v10_09_00d00 -q e26:prof
setup osclib v00.29 -q e26:n319:prof:stanfree
setup cmake v3_27_4
