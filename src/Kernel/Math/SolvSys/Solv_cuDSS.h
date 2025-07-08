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
#include <Matrice_Morse.h>
class EChaine;

#define cuDSS_

class Solv_cuDSS : public Solv_Externe
{
  Declare_instanciable_sans_constructeur_ni_destructeur(Solv_cuDSS);

public :
  Solv_cuDSS() {};
  Solv_cuDSS(const Solv_cuDSS&);
  ~Solv_cuDSS() override;

  inline int solveur_direct() const override { return 1; };
  inline int resoudre_systeme(const Matrice_Base& a, const DoubleVect& b, DoubleVect& x, int niter_max) override { return resoudre_systeme(a,b,x); };

  int resoudre_systeme(const Matrice_Base& a, const DoubleVect& b, DoubleVect& x) override;

private :
  void Create_objects(const Matrice_Morse&);
  void set_pointers_A(const Matrice_Morse&);
  void set_pointers_xb(const DoubleVect& bvect, DoubleVect& xvect);


  int nrhs=1; //For batched solve
  int n=-1;
  int nnz=-1;

  double * b_values_d=nullptr;
  double * x_values_d=nullptr;
  double * x_values_h=nullptr;
  int * csr_offsets_d=nullptr;
  int * csr_columns_d =nullptr;
  double * csr_values_d =nullptr;

  bool Axb_are_built = false;
  bool solver_is_built = false;
  bool first_solve=true;
  Matrice_Morse csr_;

#ifdef cuDSS_

  cudssAlgType_t reorder_alg = CUDSS_ALG_DEFAULT; //Can be 0->5. Default / recommended is 0
  cudssConfig_t solverConfig;
  cudssHandle_t handle;
  cudssData_t solverData;
  cudssMatrixType_t mtype;
  cudssMatrixViewType_t mview;
  cudssMatrix_t x, b, A;
  cudssIndexBase_t base = CUDSS_BASE_ONE; //Fortran indexing
  cudssStatus_t status;
  cudaStream_t stream = nullptr;
#endif
};

#endif


