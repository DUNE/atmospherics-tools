# Description

The code hereafter is made to read in the xsec systematic variations extracted from `nusystematics`. It is composed of several tools in order to do so.

# Input

The input is produced using `nusystematics` by processing CAFs to which a `SystWeights` TTree has been added, recording the weight variations consequent from xsec syst dials variations.
Example of such a file can be found on the dunegpvms at `/exp/dune/data/users/pgranger/nusyst_new_sum.root` or on `lxplus` at `/afs/cern.ch/work/p/pigrange/public/nusyst_new_sum.root`

# extract_genie

`extract_genie.cc` is provide in order to extract and compute the relevant information from the GENIE NtpGHepRecord that is stored in the CAF files under an easy to read flat tree. With GENIE loaded, the following commands should allow to extract these information:

```bash
genie
.L extract_genie.cc
extract_genie("nusyst_new_sum.root", "nusyst_genie_extracted.root")
```
The output is a flat ROOT TTree in `nusyst_genie_extracted.root` that can be read in parallel to the CAFs information (`nusyst_new_sum.root`).
An example of such a file is provided on the dunegpvm at `/exp/dune/data/users/pgranger/nusyst_genie_extracted.root` (matching the `nusyst_new_sum.root` file). On lxplus, this file is located at `/afs/cern.ch/work/p/pigrange/public/nusyst_genie_extracted.root`

# read_systs.ipynb

This notebook shows how to make some example plots with the available data and save them.