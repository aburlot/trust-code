#include <iostream>
#include <omp.h>
#include <chrono>
#include <cuda_runtime.h>

int main() {
    #pragma omp target
    {;}
    // Define the size of the array
    const int N = 1000000;

    // Allocate Unified Memory – accessible from CPU or GPU
    double *data;
    cudaMallocManaged(&data, N*sizeof(double));

    // Initialize data to 1 on the host
    for (int i = 0; i < N; i++) {
        data[i] = 1;
    }
    if (data[0]!=1) std::cerr << "Error: data[0]=" << data[0] << " should be 1." << std::endl;

    // Perform the computation on the device to 2
    // data is automatically copied to device
    #pragma omp target teams distribute parallel for
    for (int i = 0; i < N; i++)
    {
        data[i] = data[i] * 2;
    }
    // Wait for GPU to finish before accessing on host
    cudaDeviceSynchronize();
    
    // Normally thanks to UM, data is updated on the host:
    if (data[0]!=2) std::cerr << "Error: data[0]=" << data[0] << " should be 2." << std::endl;
    
    // Free memory
    cudaFree(data);
     
   /*
    // Exemple: 1 millions mailles 100^3, CL 6xn^2, donc 60000 
    const int N = 60000; // Gain interessant sur les petits envoi de donnees:
// [Pageable memory] Mean time taken to copy data H->D: 0.000295352 s 1.51356 GB/s
// [Pageable memory] Mean time taken to copy data D->H: 0.000196696 s 2.27272 GB/s
// [Pinned memory]   Mean time taken to copy data H->D: 0.000163316 s 2.73724 GB/s
// [Pinned memory]   Mean time taken to copy data D->H: 0.000150052 s 2.9792 GB/s

    const int samples = 10;
    // Initialize data on the host
    double *host_data = new double[N]();

    for (int i = 0; i < N; i++)
        host_data[i] = 0;
    #pragma omp target enter data map(alloc:host_data[0:N])
    
    double HtoD=0;
    double DtoH=0; 
    for (int sample=0;sample<samples;sample++)
    {
       auto start = std::chrono::high_resolution_clock::now();
       #pragma omp target update to(host_data[0:N])
       auto end = std::chrono::high_resolution_clock::now();
       std::chrono::duration<double> diff = end - start;
       HtoD += diff.count();       

       #pragma omp target teams distribute parallel for
       for (int i = 0; i < N; i++)
          host_data[i] = 1;

       // Copy data from device to host
       start = std::chrono::high_resolution_clock::now();
       #pragma omp target update from(host_data[0:N])
       end = std::chrono::high_resolution_clock::now();
       diff = end - start;
       DtoH += diff.count();
       if (host_data[10]!=1) abort();
    }
    #pragma omp target exit data map(release:host_data[0:N])
    HtoD/=samples;
    DtoH/=samples;
    std::cout << "[Pageable memory] Mean time taken to copy data H->D: " << HtoD << " s " << 8.0*N/(1024*1024*1024)/HtoD << " GB/s\n"; 
    std::cout << "[Pageable memory] Mean time taken to copy data D->H: " << DtoH << " s " << 8.0*N/(1024*1024*1024)/DtoH << " GB/s\n"; 

    // Pin memory on the host
    cudaError_t err = cudaHostRegister(host_data, N * sizeof(double), cudaHostRegisterDefault);
    double* pinned_mem = host_data;
    
    for (int i = 0; i < N; i++)
        pinned_mem[i] = 0;
    #pragma omp target enter data map(alloc:pinned_mem[0:N])

    HtoD=0;
    DtoH=0;
    for (int sample=0;sample<samples;sample++)
    {
       // Copy data from host to device using OpenMP target directives
       auto start = std::chrono::high_resolution_clock::now();
       #pragma omp target update to(pinned_mem[0:N])
       auto end = std::chrono::high_resolution_clock::now();
       std::chrono::duration<double> diff = end - start;
       HtoD += diff.count();
        
       // Compute data on device:
       #pragma omp target teams distribute parallel for
       for (int i = 0; i < N; i++)
          pinned_mem[i] = 1;
    
       // Copy data from device to host
       start = std::chrono::high_resolution_clock::now();
       #pragma omp target update from(pinned_mem[0:N])
       end = std::chrono::high_resolution_clock::now();
       diff = end - start;
       DtoH += diff.count();
       if (pinned_mem[10]!=1) abort();
    }
    #pragma omp target exit data map(release:pinned_mem[0:N])
    HtoD/=samples;
    DtoH/=samples;
    std::cout << "[Pinned memory]   Mean time taken to copy data H->D: " << HtoD << " s " << 8.0*N/(1024*1024*1024)/HtoD << " GB/s\n"; 
    std::cout << "[Pinned memory]   Mean time taken to copy data D->H: " << DtoH << " s " << 8.0*N/(1024*1024*1024)/DtoH << " GB/s\n"; 

    cudaHostUnregister(pinned_mem); */
    return 0;
}

