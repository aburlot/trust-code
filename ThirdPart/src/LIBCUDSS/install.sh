#!/bin/bash
build_dir=libcudss-linux-x86_64-0.6.0.5_cuda12-archive
mkdir -p $TRUST_ROOT/lib/src/LIBCUDSS
tar -xf $TRUST_ROOT/externalpackages/libcudss*
cp -r $build_dir/* $TRUST_ROOT/lib/src/LIBCUDSS
rm -rf $build_dir
echo "\$TRUST_ROOT/lib/src/LIBCUDSS/libcudss_static.a installed"
exit 0
