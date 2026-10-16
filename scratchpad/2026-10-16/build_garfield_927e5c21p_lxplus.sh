#!/bin/bash
# Rebuild Garfield++ 927e5c21 + the MX17 ComponentGrid patch for LCG_109 (lxplus + condor), 2026-10-16.
set -e
B=/tmp/dneff_gfbuild; I=$HOME/work/garfield_install/lcg109-927e5c21p
rm -rf $B; mkdir -p $B; cd $B
source /cvmfs/sft.cern.ch/lcg/views/LCG_109/x86_64-el9-gcc14-opt/setup.sh
git clone -q https://gitlab.cern.ch/garfield/garfieldpp.git src
cd src; git checkout -q 927e5c21; git apply $HOME/work/garfield_install/componentgrid_mx17.patch; git diff --stat
mkdir ../build; cd ../build
cmake ../src -DCMAKE_INSTALL_PREFIX=$I -DCMAKE_BUILD_TYPE=Release > cmake.log 2>&1
make -j8 > make.log 2>&1
make install > install.log 2>&1
cd $HOME/work/garfield_install && tar czf garfield-927e5c21p.tar.gz lcg109-927e5c21p
ls -la $HOME/work/garfield_install
echo BUILD DONE
