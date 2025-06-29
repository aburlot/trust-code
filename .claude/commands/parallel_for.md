Please analyze and use a Kokkos parallel_for instead of a loop in $ARGUMENTS for GPU port

Follow these steps:

1. Analyze the $ARGUMENTS method for potential loops
2. Search the codebase for similar port with parallel_for
3. Implement the parallel_for in the file by using similar pattern, and use timers: start_gpu_timer end_gpu_timer
4. Don't add additional #include
4. Ensure that Kokkos loop looks similar about names with previous loop
5. To keep same names inside parallel_for, rename TRUST arrays parameter of $ARGUMENTS by adding prefix tab_ if necessary
6. To keep same names inside parallel_for you may supress temporary references used before the loop:

    const DoubleTab& xv = domaine().xv();
    for(...)
    
    Becomes:
    
    CDoubleTabView xv = domaine().xv().view_ro();
    Kokkos::parallel_for(...)
    
7. Use static_cast<ArrOfDouble&> resp. static_cast<const ArrOfDouble&>  to create 1D view on DoubleTab resp const DoubleTab if necessary
8. Look for race conditions inside the parallel_for and fix with atomics if necessary
9. Build the code, fast way: cd MonoDir_mpi_opt/src && make -j20 
10. Run tests with: trust -check $ARGUMENTS 


