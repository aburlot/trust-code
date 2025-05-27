#!/bin/bash
# Only install for Nvidia GPU for the moment
# Installed by PETSc so disabled here:
exit 0
[ "$TRUST_USE_CUDA" != 1 ] && exit 0
archive=$TRUST_ROOT/externalpackages/kokkos/kokkos-kernels-4.6.01.tar.gz

build_dir=$TRUST_ROOT/build/kokkos-kernels
KOKKOS_ROOT_DIR=$TRUST_ROOT/lib/src/LIBKOKKOS
# Log file of the process:
log_file=$TRUST_ROOT/kokkos-kernels_compile.log

if [ ! -f $KOKKOS_ROOT_DIR/lib64/libkokkos.a ]
then
    echo "# Installing `basename $archive` ..."
    # Sub shell to get back to correct dir:
    (  
      mkdir -p $build_dir;cd $build_dir
      tar xzf $archive || exit -1
      src_dir=$build_dir/`ls $build_dir | grep kokkos-kernels`

      BUILD_TYPES="Release `[ "$build_debug" = "1" ] && echo Debug`"
      for CMAKE_BUILD_TYPE in $BUILD_TYPES
      do
        rm -rf BUILD;mkdir -p BUILD;cd BUILD
        [ "$TRUST_USE_CUDA" = 1 ] && CMAKE_OPT="-DCMAKE_CXX_COMPILER=$TRUST_CC_BASE_EXTP"
        [ "$TRUST_USE_ROCM" = 1 ] && CMAKE_OPT="-DCMAKE_CXX_COMPILER=hipcc" # $TRUST_CC_BASE pour profiter de ccache ?
        CMAKE_OPT="$CMAKE_OPT -DCMAKE_CXX_FLAGS=-fPIC"

        CMAKE_INSTALL_PREFIX=$KOKKOS_ROOT_DIR/$TRUST_ARCH`[ $CMAKE_BUILD_TYPE = Release ] && echo _opt`
        CMAKE_OPT=$CMAKE_OPT" -DCMAKE_BUILD_TYPE=$CMAKE_BUILD_TYPE -DCMAKE_INSTALL_PREFIX=$CMAKE_INSTALL_PREFIX -DCMAKE_INSTALL_LIBDIR=lib64"
        CMAKE_OPT=$CMAKE_OPT" -DKokkos_ROOT=$CMAKE_INSTALL_PREFIX"
        CMAKE_OPT=$CMAKE_OPT" -DKokkosKernels_INST_DOUBLE=ON"
        CMAKE_OPT=$CMAKE_OPT" -DKokkosKernels_INST_COMPLEX_DOUBLE=ON" 
        CMAKE_OPT=$CMAKE_OPT" -DKokkosKernels_INST_ORDINAL_INT=ON"
        CMAKE_OPT=$CMAKE_OPT" -DKokkosKernels_INST_ORDINAL_INT64_T=ON"
        CMAKE_OPT=$CMAKE_OPT" -DKokkosKernels_INST_OFFSET_INT=ON"
        CMAKE_OPT=$CMAKE_OPT" -DKokkosKernels_INST_OFFSET_SIZE_T=ON"
        CMAKE_OPT=$CMAKE_OPT" -DKokkosKernels_INST_LAYOUTLEFT=ON"
        CMAKE_OPT=$CMAKE_OPT" -DKokkosKernels_ADD_DEFAULT_ETI=ON"
        CMAKE_OPT=$CMAKE_OPT" -DKokkosKernels_ENABLE_TESTS=OFF"
        # Configure
        [ "$NVHPC_CUDA_HOME" != "" ] && export CUDA_ROOT=$NVHPC_CUDA_HOME
        cmake $src_dir $CMAKE_OPT 2>&1 | tee -a $log_file
        [ ${PIPESTATUS[0]} != 0 ] && echo "Error when configuring Kokkos (CMake) - look at $log_file" && exit -1

        # Build
        make -j$TRUST_NB_PHYSICAL_PROCS install 2>&1 | tee -a $log_file
        [ ${PIPESTATUS[0]} != 0 ] && echo "Error when compiling Kokkos-kernels - look at $log_file" && exit -1
        echo "Kokkos-kernels $CMAKE_BUILD_TYPE installed under $CMAKE_INSTALL_PREFIX"
        cd ..
      done
      # Clean build:
      rm -rf $build_dir $log_file
    )
else
    echo "# Kokkos-kernels: already installed. Doing nothing."
fi

