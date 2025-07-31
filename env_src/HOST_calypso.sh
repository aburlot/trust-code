#!/bin/bash
##################################
# Variables for configure script #
##################################
define_modules_config()
{
   env=$TRUST_ROOT/env/machine.env
   # Initialisation de l environnement module $MODULE_PATH
   #echo "source /etc/profile.d/modules.sh " >> $env
   # Load modules
   module="python/3.12.10 tools/cmake/3.28.2_arm nvhpc/24.1"  # openmpi/4.1.7_gcc114_cuda124 gcc/12.3.0_arm
   [ "$TRUST_CUDA_CC" = "" ] && TRUST_CUDA_CC=90 # H100
   echo "# Module $module detected and loaded on $HOST." 
   echo "module purge 1>/dev/null" >> $env
   echo "module load $module" >> $env
   echo "[ \$? != 0 ] && echo \"Error: $module not found; we exit...\" && echo \"Contat TRUST support team or system administrator\" && exit -1" >> $env
   echo $source >> $env
   . $env
   # Creation wrapper qstat -> squeue
   echo "#!/bin/bash
squeue" > $TRUST_ROOT/bin/qstat
   chmod +x $TRUST_ROOT/bin/qstat
}

##############################
# Variables for trust script #
##############################
define_soumission_batch()
{
   soumission=2
   [ "$prod" = 1 ] && soumission=1
   [ "$gpu"  = 1 ] && soumission=1
   queue=grace && gpus_per_node=`echo $NB_PROCS | awk '{print $1<1?$1:1}'` && noeuds=`echo "1+($NB_PROCS-1)/1" | bc` # 1GPU/node
   if [ "$prod" = 1 ] || [ "$NB_PROCS" -gt 40 ]
   then
      qos=2jours && cpu=2880
      [ "$gpu" != 1 ] && node=1 # exclusif uniquement sur cpu
   else
      qos=test	&& cpu=60   && node=0 
   fi
   # Le binding ameliore fortement les performances sur AMD quelque soit MPI:
   if [ "$I_MPI_ROOT" != "" ] # IntelMPI
   then
      binding="-m block:block --cpu-bind=rank"
      mpirun="srun $binding -n \$SLURM_NTASKS"
   elif  [ "$HPCX_DIR" != "" ] # HPC-X
   then 
      binding="--map-by numa --bind-to core"
      mpirun="mpirun $binding -n \$SLURM_NTASKS"
   else
      mpirun="srun --mpi=pmix -n \$SLURM_NTASKS"
   fi
   sub=SLURM
}
