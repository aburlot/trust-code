# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

TRUST is a comprehensive computational fluid dynamics (CFD) platform developed by CEA for thermohydraulic simulations. Originally designed for incompressible single-phase and Low Mach Number flows, it now supports compressible multi-phase flows and is being progressively ported to support GPU acceleration (NVIDIA/AMD) through Kokkos.

**This repository (`review_gpu`) is specifically focused on GPU acceleration development and testing**, containing extensive GPU test suites and performance benchmarks for validating GPU implementations against CPU reference solutions.

The current development branch `ledac/next_gpu_vdf2` appears to be actively working on VDF (Finite Difference) GPU implementations, with multiple GPU test cases including:
- GAMELAN_AMG: GPU linear solver benchmarks
- JEL_bous: Boussinesq flow GPU tests  
- OpenMP integration with GPU acceleration
- Performance comparison scripts and solver benchmarking

## Core Build Commands

The project uses a combination of configure scripts, makefiles, and CMake for building:

### Initial Configuration
```bash
# Initialize TRUST environment
source ./env_TRUST.sh

# Configure with options (see ./configure -help for full list)
./configure [OPTIONS]

# Examples:
./configure -cuda -openmp                    # Enable CUDA and OpenMP
./configure -with-coolprop=/path/to/coolprop # Link with CoolProp
./configure -disable-optionals               # Minimal build
```

### Building TRUST
```bash
# Build optimized version only (for users)
make optim

# Build both optimized and debug versions (for developers)
make

# Build specific targets
make debug              # Debug version only
make semi_opt          # Semi-optimized with assertions
make tools             # Build tools only
make clean             # Clean build artifacts
```

### Testing
```bash
# Parallel testing (recommended)
make ctest_optim       # Run tests in parallel (optimized)
make ctest_debug       # Run tests in parallel (debug)

# Sequential testing
make check             # Run tests sequentially (optimized)
make check_debug       # Run tests sequentially (debug)

# Rerun only failed tests
make ctest_rerun_failed_optim  # Rerun failed tests (optimized)
make ctest_rerun_failed_debug  # Rerun failed tests (debug)

# Validation
make validation        # Run validation test suite

# GPU-specific testing
cd tests/GPU && ./check_perf.sh  # GPU performance regression tests
```

### Environment Management
```bash
# Initialize Python environment (MEDCoupling, swig, ICoCo)
source ./env_for_python.sh

# Go to TRUST root from anywhere
cd $TRUST_ROOT

# Go to TRUST tests directory
cd $TRUST_TESTS
```

## High-Level Architecture

TRUST follows a modular architecture with clear separation between physics, discretization methods, and computational kernels:

### Discretization Methods
- **VDF (Volumes Différentiels Finis)**: Finite difference on structured grids - `/src/VDF/`
- **VEF (Volumes Éléments Finis)**: Finite element volumes with P1NC elements - `/src/VEF/`
- **PolyMAC**: Polygonal/polyhedral MAC method with variants (P0, P0P1NC) - `/src/PolyMAC/`
- **EF**: Classical finite elements - `/src/EF/`
- **DG**: Discontinuous Galerkin - `/src/DG/`
- **IJK**: Optimized structured Cartesian grids - `/src/IJK/`

### Physics Modules
- **ThHyd**: Complete thermohydraulics framework including:
  - Incompressible flows (`Pb_Hydraulique`, `Pb_Thermohydraulique`)
  - Quasi-compressible flows (Low Mach)
  - Weakly compressible flows
  - Multiphase Euler-Euler with extensive correlation libraries
  - Turbulence models (RANS, LES) and wall laws
- **ThSol**: Solid heat conduction (`Pb_Conduction`)

### Core Framework (`/src/Kernel/`)
- **Framework**: Base classes (`Probleme_base`, `Equation_base`, `Champ_base`)
- **Math**: Array system (`TRUSTArray`, `TRUSTTab`), linear algebra
- **Geometry**: Domain management, meshing tools, parallel decomposition
- **Fields/BCs**: Generic fields, boundary conditions
- **GPU Support**: Kokkos integration, device management

### Class Naming Patterns
- `*_base`: Abstract base classes defining interfaces
- `*_VDF`, `*_VEF`, `*_PolyMAC`: Discretization-specific implementations
- `Op_*`: Discrete operators (convection, diffusion, divergence, gradient)
- `Pb_*`: Problem definitions
- `Terme_*`: Source terms

## GPU Acceleration (Kokkos Integration)

TRUST is being ported to GPU using Kokkos for performance portability:

### Configuration for GPU
```bash
./configure -cuda -openmp          # NVIDIA GPU support
./configure -rocm                  # AMD GPU support
./configure -kokkos_openmp         # OpenMP backend
./configure -kokkos_simd           # SIMD optimization
```

### Key GPU Components
- **Device Management**: `/src/Kernel/Utilitaires/Device.*`
- **Kokkos Wrappers**: `kokkos++.h`, view types
- **Performance Tools**: GPU timers, NVTX integration
- **Array System**: GPU-compatible `TRUSTArray` with view access

### GPU Development Conventions
Based on `/src/CLAUDE_UPDATE.md`, the established patterns are:

1. **Includes**: Kokkos headers (`kokkos++.h`, `TRUSTArray_kokkos.tpp`) are often already included transitively - check before adding
2. **Data Containers**: 
   - Prefer `DoubleTrav` over `DoubleTab` for GPU compatibility and performance
   - `DoubleTrav` uses memory pooling - when deleted, memory returns to pool instead of deallocation, reducing allocation overhead
   - Prefer `DoubleTrav` over `DoubleVect` for GPU compatibility
3. **View Creation**:
   ```cpp
   CDoubleArrView input_view = static_cast<const ArrOfDouble&>(array).view_ro();
   DoubleArrView output_view = static_cast<ArrOfDouble&>(array).view_wo();
   DoubleTabView3 output_3d = array.view_rw<3>();  // For 3D arrays
   ```
4. **Parallel Loops**:
   ```cpp
   Kokkos::parallel_for(start_gpu_timer(__KERNEL_NAME__), range_1D(start, end), KOKKOS_LAMBDA(const int i)
   {
     output_view(i) = Kokkos::sqrt(input_view(i));
   });
   end_gpu_timer(__KERNEL_NAME__);
   ```
5. **Atomic Operations**: Only use `Kokkos::atomic_add()` when actual race conditions exist (multiple threads writing to same memory location)
6. **Minimal Changes**: Keep variable naming and logic as close to original as possible in Kokkos regions

### GPU Testing and Performance
```bash
# Verify GPU/CPU consistency
trust -check ClassName::method_name

# Run GPU performance tests
cd tests/GPU && ./check_perf.sh          # All GPU performance tests
cd tests/GPU/[TestCase] && ./check_perf.sh  # Individual test case

# GPU profiling with NSight Systems
cd tests/GPU/[TestCase] && ./check_perf.sh -nsys

# Test specific solver performance
cd tests/GPU && ./check_solver.sh       # Check different GPU solver options
```

## Development Workflow

### Code Organization
- Each discretization method has complete parallel implementation
- Physics modules are discretization-independent where possible
- Template-based programming for performance
- Deep inheritance hierarchies for code reuse

### Testing and Validation
- Comprehensive test suite in `/tests/`
- Non-regression testing with reference solutions
- Validation forms for physical verification
- GPU/CPU consistency checking

### Third-Party Dependencies
- **PETSc**: Linear solvers, preconditioners
- **MED/HDF5**: Mesh and field I/O
- **MEDCoupling**: Advanced mesh operations
- **Metis/ParMetis**: Domain decomposition
- **CoolProp/EOS**: Thermophysical properties
- **Kokkos**: GPU portability
- **CGNS**: Mesh format support

### External Tools Integration
- **VisIt**: Visualization (install with `-download-visit` or `-without-visit`)
- **Gmsh**: Mesh generation
- **Valgrind**: Memory debugging
- **Doxygen**: Documentation generation

## Key Configuration Options

### Compilation Options
```bash
-c++=<compiler>           # Force C++ compiler
-cxxflags=<flags>         # Add C++ compiler flags
-with-32-bit-indices      # Use 32-bit integers (default: 64-bit)
-native                   # Architecture-specific optimization
-std=c++17               # C++ standard (default)
```

### Physics Libraries
```bash
-with-eos=<path>         # Link external EOS library
-with-coolprop=<path>    # Link external CoolProp library
-with-refprop=<path>     # REFPROP support via CoolProp
```

### Parallel/GPU Options
```bash
-openmp                  # Enable OpenMP
-cuda[=path/download]    # CUDA support
-rocm[=gfx...]          # AMD ROCm support
-force_system_mpi       # Use system MPI
-force_mpi_gpu_aware    # GPU-aware MPI
```

### Build Control
```bash
-disable-optionals      # Minimal build (no optional libraries)
-disable-tools         # Build binary only
-without-visit         # Skip VisIt installation
-without-conda         # Don't use conda for Python tools
```

## Performance and Debugging

### Performance Analysis
```bash
# GPU performance check
cd tests/GPU && ./check_perf.sh

# Individual GPU test performance
cd tests/GPU/[TestCase] && ./check_perf.sh

# GPU profiling with NSight Systems (for supported tests)
cd tests/GPU/[TestCase] && ./check_perf.sh -nsys

# GPU solver benchmarking
cd tests/GPU && ./check_solver.sh

# CPU profiling tools
trust -heaptrack datafile    # Memory profiling
trust -valgrind datafile     # Memory debugging
```

### Common Issues
- GPU builds require CUDA toolkit or ROCm installation
- MPI versions must be compatible (check `type mpicxx`)
- Large path names (>90 chars) can cause conda issues
- LaTeX packages needed for PDF documentation generation

This architecture represents a mature, research-grade CFD platform with comprehensive physics modeling, multiple discretization schemes, and modern GPU acceleration capabilities suitable for high-performance computing applications.