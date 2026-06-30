#!/bin/bash

set -x

export X509_USER_PROXY=${CONDOR_DIR_INPUT}/jobscript-proxy.pem
chmod 0400 ${X509_USER_PROXY}
ls -l ${X509_USER_PROXY} #Check that the proxy is there

source /cvmfs/dune.opensciencegrid.org/products/dune/setup_dune.sh
setup ifdhc

CODEDIR=${INPUT_TAR_DIR_LOCAL}/xsec_systs_calculation
cd $CODEDIR
source $CODEDIR/setup.sh

if [ -f "$CODEDIR/build/Linux/bin/setup.fhicl_cpp_standalone.sh" ]; then
  source "$CODEDIR/build/Linux/bin/setup.fhicl_cpp_standalone.sh"
fi
if [ -f "$CODEDIR/build/Linux/bin/setup.systematicstools.sh" ]; then
  source "$CODEDIR/build/Linux/bin/setup.systematicstools.sh"
fi
if [ -f "$CODEDIR/build/Linux/bin/setup.nusystematics.sh" ]; then
  source "$CODEDIR/build/Linux/bin/setup.nusystematics.sh"
fi

# Override absolute paths from compilation environment with grid runtime paths
export fhicl_cpp_standalone_ROOT=$CODEDIR/build/Linux
export systematicstools_ROOT=$CODEDIR/build/Linux
export nusystematics_ROOT=$CODEDIR/build/Linux
export FHICL_FILE_PATH=".:${nusystematics_ROOT}/fcl:${CODEDIR}/fcl:${CODEDIR}/tools:${FHICL_FILE_PATH}"
export PATH="${CODEDIR}/build/Linux/bin:${PATH}"
export LD_LIBRARY_PATH="${CODEDIR}/build/Linux/lib:${LD_LIBRARY_PATH}"

# Capture original GENIE path from setup
ORIG_GENIE=$GENIE

export GINUKEHADRONDATA="${ORIG_GENIE}/data/evgen/intranuke"

#Creating the workdir
WORKDIR=${_CONDOR_SCRATCH_DIR}/work/
mkdir -p $WORKDIR && cd $WORKDIR

# Set up custom GENIE environment inside writeable WORKDIR by linking config and original data catalogues
mkdir -p genie_config
ln -s ${CODEDIR}/genie_config/config genie_config/config

# Build a MERGED data tree instead of a flat symlink to ${ORIG_GENIE}/data.
# The base ups GENIE data lacks the Martini MEC hadron tensors
# (data/evgen/hadron_tensors/martini), which the alternative-MEC reweight dials
# (XSecShape_CCMEC_Martini, ...) need. We mirror the base tree with symlinks and
# inject the custom Martini tensors shipped in the tarball.
mkdir -p genie_config/data/evgen/hadron_tensors
for p in ${ORIG_GENIE}/data/*; do
  bn=$(basename "$p"); [ "$bn" = "evgen" ] && continue
  ln -sfn "$p" genie_config/data/"$bn"
done
for p in ${ORIG_GENIE}/data/evgen/*; do
  bn=$(basename "$p"); [ "$bn" = "hadron_tensors" ] && continue
  ln -sfn "$p" genie_config/data/evgen/"$bn"
done
for p in ${ORIG_GENIE}/data/evgen/hadron_tensors/*; do
  ln -sfn "$p" genie_config/data/evgen/hadron_tensors/"$(basename "$p")"
done
# inject the custom Martini MEC hadron tensors (absent from the base ups release)
ln -sfn ${CODEDIR}/genie_config/hadron_tensors_custom/martini genie_config/data/evgen/hadron_tensors/martini

export GENIE="${WORKDIR}/genie_config"
export GENIE_REWEIGHT="${WORKDIR}/genie_config"

IFILE="root://fndca1.fnal.gov:1094/pnfs/fnal.gov/usr/dune/persistent/users/pgranger//caf_new_sum.2.6M_weighted.root"

cd $CODEDIR/build/
# export ROOT_INCLUDE_PATH=${CODEDIR}/build/ #Required to have StandardRecord.h found when loading CAFs

#Necessary libray path additions
echo "Listing the libs in ${CODEDIR}/build/Linux/lib/:"
ls -lht ${CODEDIR}/build/Linux/lib/

LD_LIBRARY_PATH=${CODEDIR}/build/Linux/lib/:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH

echo $LD_LIBRARY_PATH
ldd ./app/UpdateReweight

env

STEP_SIZE=3000
#STEP_SIZE=3000
# --- targeted resubmit: map condor PROCESS -> failed job id ---
FAILED_JOBS=(30 34 69 76 78 80 117 127 146 147 274 276 321 322 344 349 351 376 384 390 394 397 399 464 506 520 532 573 574 634 637 656 663 664 692 736 766 856 876 881 885 886 914 916 925 926 948)
JOB_ID=${FAILED_JOBS[$PROCESS]}
if [ -z "$JOB_ID" ]; then echo "No job id for PROCESS=$PROCESS, exiting"; exit 0; fi
START_EVENT=$(( JOB_ID * STEP_SIZE ))

LOCAL_OFILE="${WORKDIR}/xsec_systs_output_${JOB_ID}.root"
OFILE="/pnfs/dune/scratch/users/pgranger/xsec_systs_outputs_wsbn/xsec_systs_output_${JOB_ID}.root"

#Running code
CMD="./app/UpdateReweight -i $IFILE -c ../fcl/systs_atmospherics_v2.fcl -o $LOCAL_OFILE -s $START_EVENT -N $STEP_SIZE"
eval "$CMD"

#Copying output to dCache
ifdh cp $LOCAL_OFILE $OFILE
