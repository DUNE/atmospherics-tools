#!/bin/bash

# Source the environment setup
source ../../setup_update_reweight.sh

rm -rf build
mkdir build
cd build

export GENIE_REWEIGHT=/exp/dune/app/users/pgranger/systematics-lars/Reweight

cmake .. \
    -DCMAKE_EXE_LINKER_FLAGS="-L${GENIE_REWEIGHT}/lib -Wl,-rpath,'\$ORIGIN/../../../Reweight/lib:\$ORIGIN/../lib'" \
    -DCMAKE_INSTALL_RPATH_USE_LINK_PATH=FALSE \
    -DCMAKE_SKIP_BUILD_RPATH=FALSE \
    -DCMAKE_BUILD_WITH_INSTALL_RPATH=TRUE \
    -DSQLite3_LIBRARY="${SQLITE_LIB}/libsqlite3_ups.so" \
    -DSQLite3_INCLUDE_DIR="${SQLITE_INC}"
make -j$(nproc)
make install

echo ""
echo "Copying custom GENIE libraries to grid build directory..."
mkdir -p Linux/lib
# Copy all custom libraries (including generator, reweight, flux and geometry drivers)
cp -av /exp/dune/app/users/pgranger/systematics-lars/local_install/lib/libG*.so* Linux/lib/
# Copy the .pcm dictionary files (required for ROOT class decoding)
cp -av /exp/dune/app/users/pgranger/systematics-lars/local_install/lib/*_rdict.pcm Linux/lib/
# Copy rootmap files
cp -av /exp/dune/app/users/pgranger/systematics-lars/local_install/lib/*.rootmap Linux/lib/
