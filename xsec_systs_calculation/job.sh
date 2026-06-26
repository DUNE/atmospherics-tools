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
ln -s ${ORIG_GENIE}/data genie_config/data
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
JOB_ID=$(( PROCESS ))
START_EVENT=$(( JOB_ID * STEP_SIZE ))

LOCAL_OFILE="${WORKDIR}/xsec_systs_output_${JOB_ID}.root"
OFILE="/pnfs/dune/scratch/users/pgranger/xsec_systs_outputs_wsbn/xsec_systs_output_${JOB_ID}.root"

#Running code
CMD="./app/UpdateReweight -i $IFILE -c ../fcl/systs_atmospherics_v2.fcl -o $LOCAL_OFILE -s $START_EVENT -N $STEP_SIZE"
eval "$CMD"

#Copying output to dCache
ifdh cp $LOCAL_OFILE $OFILE
