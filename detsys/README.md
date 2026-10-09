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

## Example
An example input CAF file with ~1M events can be found at Fermilab at `/exp/dune/app/users/pgranger/weightor/caf_sum.root`
To compute the weights for it with the default parameters, one can simply run `./app/weightor -i /exp/dune/app/users/pgranger/weightor/caf_sum.root -o weighted_caf.root`

# Detsyst

This directory contains apps for comparing CAF files and build covariance matrix for detector systematics estimation.

## Comparison script
Execute with:
```./build/app/ComparisonScript Config.yaml```

This app plots various distribution of observables for all the samples given in Config.yaml.
It also computes the covariance matrix or the ratio of spectra in a root format to use as inputs for the fitters.

NB: If you use an old reco2 version (like for atm production in 2023), there might be a bug in the CVN (CVN numu and NC are inverted). In that case, use "FixCVN: true" in the YAML config for this sample.

## CompareEventByEvent
Execute with:
```./build/app/CompareEventByEvent ConfigCompareEventByEvent.yaml```

This app plots relative difference of observables between two samples. The samples must contain the same events: a matching event is done by run, subrun and event number between the two samples. This is useful to compare wiremod and detector variation for instance.

The file UnMatched.pdf plots the observables distributions for the first sample for events with particularly big discerepancies.

## root_covariance_to_yaml
Execute with:
```./build/app/root_covariance_to_yaml covariance.root covariance.yaml```

This app converts the covariance root file produced byt ComparisonScript to a YAML config file compatible with MaCh3.
