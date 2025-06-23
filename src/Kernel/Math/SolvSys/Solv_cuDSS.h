/****************************************************************************
* Copyright (c) 2025, CEA
* All rights reserved.
*
* Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
* 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
* 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
* 3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
*
* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
* IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
* OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*
*****************************************************************************/

#ifndef Solv_cuDSS_included
#define Solv_cuDSS_included

#include <Solv_Externe.h>
#include <Motcle.h>
#include <Device.h>
#include <cudss.h>
#include <cuda.h>
class EChaine;

#define cuDSS_

class Solv_cuDSS : public Solv_Externe
{
  Declare_instanciable_sans_constructeur_ni_destructeur(Solv_cuDSS);

public :
  Solv_cuDSS();
  Solv_cuDSS(const Solv_cuDSS&);
  ~Solv_cuDSS() override;

  inline int solveur_direct() const override { return 1; };
  inline int resoudre_systeme(const Matrice_Base& a, const DoubleVect& b, DoubleVect& x, int niter_max) override { return resoudre_systeme(a,b,x); };

  int resoudre_systeme(const Matrice_Base& a, const DoubleVect& b, DoubleVect& x) override;
  void create_solver(Entree& entree);

private :
  void initialize();
  void Create_objects(const Matrice_Morse&, const DoubleVect& bvect, DoubleVect& xvect);

  cudaError_t cuda_error = cudaSuccess;


  int *csr_offsets_d = nullptr;
  int *csr_columns_d = nullptr;
  double *csr_values_d = nullptr;
  double *x_values_d = nullptr, *b_values_d = nullptr;

  int nrhs=1; //For batched solve
  int n;
  int nnz;

#ifdef cuDSS_
  cudssStatus_t status = CUDSS_STATUS_SUCCESS;
  cudssAlgType_t reorder_alg;
  cudssMatrixType_t mtype;
  cudssMatrixViewType_t mview;
  cudssConfig_t solverConfig;
  cudssData_t solverData;
  cudssMatrix_t x, b, A;
  cudssIndexBase_t base = CUDSS_BASE_ZERO;
  cudssHandle_t handle;

  bool A_is_built = false;
  bool b_is_built = false;
  bool x is_built = false;
  bool handle_is_built = false;
  bool config_is_built=false;
  bool data_is_built=false;

#endif
};

#endif


