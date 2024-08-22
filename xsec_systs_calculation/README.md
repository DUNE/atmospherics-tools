# Xsec syst calculations

## Description

Slightly readapted from Jaesung Kim -> https://github.com/jedori0228/NDNuSyst/tree/master

Codes allowing to apply some xsec systs reweighting over CAF files.

## Setup

```bash
source setup.sh
mkdir build
cd build
cmake ..
make install
```

## Example

`./app/UpdateReweight -i /pnfs/dune/persistent/users/pgranger/atmospherics-data/atmospherics_prod_1M_events_cafs_hadd.root -N 5 -o out.root -c ../fcl/mach3_systs_updated.fcl` computes the new weights for several shifts in GENIE tweakable parameters on 5 files and prints the result.
