#!/bin/bash
# Create a PETSc matrix from test case and run ex10 on GPU in parallel
echo "Usage: $0 [-n NP] [petsc_options]"

ROOT=`pwd`
jdd=`basename $ROOT`"_BENCH"
if [ $jdd.data -nt petsc.data ]
then
   cp $jdd.data petsc.data
   sed -i "1,$ s?impr }?save_matrix_petsc_format impr }?" petsc.data
   sed -i "1,$ s?nb_pas_dt_max 10?nb_pas_dt_max 0?" petsc.data
   echo "Saving matrix from $jdd test case..."
   trust petsc
fi
matrix=`ls -rt Matrix*petsc | tail -1`
# See https://www.mcs.anl.gov/petsc/petsc-current/src/ksp/ksp/tutorials/ex10.c.html
exec=`pwd`/ex10
if [ ! -f $exec ]
then
   dir=$PETSC_ROOT/linux_opt/share/petsc/examples/src/ksp/ksp/tutorials
   #cp $dir/ex10.c $dir/makefile .
   [ "$TRUST_CC_BASE_EXTP" != "" ] && export MPICH_CC=$TRUST_cc_BASE_EXTP && export OMPI_CC=$TRUST_cc_BASE_EXTP
   make ex10 1>/dev/null 2>&1 || exit -1
fi
NP=2 && [ "$1" = -n ] && shift && NP=$1 && shift
# See https://www.mcs.anl.gov/petsc/documentation/faq.html#conditionnumber: -pc_type jacobi pour voir l'effet du preconditionnement
petsc_options="-ksp_monitor -ksp_type cg -pc_type hypre -pc_hypre_type boomeramg -mat_type mpiaijcusparse -pc_hypre_boomeramg_strong_threshold 0.7"
petsc_options="-ksp_monitor -ksp_type cg -pc_type amgx"
petsc_options=$petsc_options" "$*
cmd="$NP -f $matrix -ksp_view -log_view $petsc_options"
echo "mpirun -np $NP ./ex10 $cmd"
touch dummy.data 
trust -nsys dummy $cmd 2>&1 | tee ex10.log

# $exec -help 2>&1 | grep cuda
