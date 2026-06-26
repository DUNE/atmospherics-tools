SOURCE="dropbox:///exp/dune/app/users/pgranger/systematics-lars/atmospherics-tools/xsec_systs_calculation_wmec.tar.gz"
GROUP="dune"
EXEC=$(readlink -f job.sh)
PROXYFILE=$PWD/jobscript-proxy.pem
# export X509_USER_PROXY=$PROXYFILE

# /exp/dune/app/users/pgranger/dune-justin/commands/justin --verbose get-token

jobsub_submit --expected-lifetime=4h \
  --group=${GROUP} \
  --resource-provides=usage_model=DEDICATED,OPPORTUNISTIC \
  --tar_file_name=${SOURCE} \
  --disk=4GB \
  --memory=4GB \
  --singularity-image="/cvmfs/singularity.opensciencegrid.org/fermilab/fnal-wn-sl7:latest" \
  -f dropbox://$PROXYFILE \
  --cpu=1 \
  -N 1000 \
  file://$EXEC
