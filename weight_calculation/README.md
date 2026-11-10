# Atmospherics-tools

This repository contains basic scripts and tools allowing to produce analysis ntuples for the HEP atmospherics analysis.
The various tools are made available in the `src` directory while some executable sources are provided in `app`.
The supported input is the new hierarchical CAF format.

## Setup
On SL7:
```bash
source setup.sh
mkdir build
cd build
cmake ..
make install
```

## Apps
The `weightor.cxx` code produces a `weightor` binary that can be used to compute the various event weights necessary for the atmospherics oscillation analysis. Are provided in the final ntuple for each event:
- The xsec weight
- The nue(bar) flux weight
- The numu(bar) flux weight
- The osc probability from nue(bar) to detected flavour
- The osc probability from numu(bar) to detected flavour
- The final oscillated weight constructed from the previous ones

Multiple parameters can be tweaked to produce weight with different parameters or setups. The reference values are provided for the NuFIT 5.2 results and assumes the use of the 2023 HEP atmospherics sample.

The `split_channel.cxx` code produces a `split_channel` binary that is used to split the input sample into many subsamples, separated by initial flavour, final flavour, reco flavour that are used as inputs of MaCh3.

## LArSoft inputs
The weights can also be computed directly from LArSoft (art) files, without going through CAFs, in two steps:
1. Extract the needed truth information (neutrino energy, momentum and PDG from `MCTruth`, GENIE weight from `GTruth::fweight`, run/subrun/event) and the exposure (`POTSummary::totpot` of each subrun) into a flat ntuple. This needs a LArSoft environment for the dictionaries:
   ```bash
   source larsoft/setup_dunesw.sh   # dunesw + osclib, also usable to build the weightor
   ls /path/to/reco1/*.root > filelist.txt
   root -l -b -q 'larsoft/extract_larsoft.C("filelist.txt", "flat.root")'
   ```
2. Run the weightor on it with `--flat`. Remember to set the detector mass (`-m`) matching the simulated geometry. Only the `weights` tree is written. `sample_id` is not computed since there is no CVN information.
   ```bash
   ./app/weightor --flat -i flat.root -o weights.root -m <mass>
   ```
Events can be matched back to the LArSoft/hdf5 files through the run/subrun/eid branches.

## Example
An example input CAF file with ~1M events can be found at Fermilab at `/exp/dune/app/users/pgranger/weightor/caf_sum.root`
To compute the weights for it with the default parameters, one can simply run `./app/weightor -i /exp/dune/app/users/pgranger/weightor/caf_sum.root -o weighted_caf.root`
