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

#include <Solv_cuDSS.h>
#include <Matrice_Morse_Sym.h>
#include <TRUSTVect.h>
#include <EChaine.h>
#include <Motcle.h>
#include <SFichier.h>
#include <Comm_Group_MPI.h>
#include <MD_Vector_std.h>
#include <MD_Vector_composite.h>
#include <Device.h>
#include <stat_counters.h>
#include <Array_tools.h>
#include <cuda.h>

Implemente_instanciable_sans_constructeur_ni_destructeur(Solv_cuDSS, "Solv_cuDSS", Solv_Externe);
// XD cuDSS Solv_Externe cuDSS 0 Solver via cuDSS API
// XD attr solveur chaine solveur 0 not_set
// XD attr option_solveur bloc_lecture option_solveur 0 not_set

// printOn
Sortie& Solv_cuDSS::printOn(Sortie& s ) const
{
  s << chaine_lue_;
  return s;
}

// readOn
Entree& Solv_cuDSS::readOn(Entree& is)
{
  create_solver(is);
  return is;
}

Solv_cuDSS::Solv_cuDSS()
{
  initialize();
}

Solv_cuDSS::Solv_cuDSS(const Solv_cuDSS& org)
{
  initialize();
  // on relance la lecture ....
  EChaine recup(org.get_chaine_lue());
  readOn(recup);
}

Solv_cuDSS::~Solv_cuDSS()
{
#ifdef cuDSS_

  if (A_is_built)
	  cudssMatrixDestroy(A);

  if (b_is_built)
	  cudssMatrixDestroy(b);

  if (x_is_built)
	  cudssMatrixDestroy(x);

  if (data_is_built)
	  cudssDataDestroy(handle, solverData);
  cudssConfigDestroy(solverConfig);
  cudssDestroy(handle);
#endif
}

void Solv_cuDSS::initialize()
{
#ifdef cuDSS_
#endif
}

// Lecture et creation du solveur
void Solv_cuDSS::create_solver(Entree& entree)
{
#ifdef cuDSS_

  cudssAlgType_t reorder_alg;
  reorder_alg = CUDSS_ALG_DEFAULT;

  lecture(entree);
  EChaine is(get_chaine_lue());
  Motcle accolade_ouverte("{"), accolade_fermee("}");
  Motcle solver, motlu;
  is >> solver;   // On lit le solveur en premier puis les options du solveur
  is >> motlu; // On lit l'accolade
  if (motlu != accolade_ouverte)
    {
      Cerr << "Error while reading the parameters of the solver " << solver << " { ... }" << finl;
      Cerr << "We expected " << accolade_ouverte << " instead of " << motlu << finl;
      exit();
    }
  // Lecture des parametres du solver (LU non symetric, cholesky symmetric)
  // LU|Cholesky { algo name [impr] }
  is >> motlu;
  while (motlu!=accolade_fermee)
    {
      if (motlu==(Motcle)"impr")
        {
          fixer_limpr(1);
        }
      else if (motlu==(Motcle)"algo")
        {
          is >> motlu; // Lecture de l'algo
        }
      else
        {
          Cerr << motlu << " keyword not recognized for cuDSS solver " << solver << finl;
          Process::exit();
        }
      is >> motlu;
    }



  /* Create handle */
  cudssCreate(&handle);

  /* Create config */
  cudssConfigCreate(&solverConfig);

  reorder_alg=CUDSS_ALG_DEFAULT;
  mtype=CUDSS_MTYPE_GENERAL;
  mview=CUDSS_MVIEW_FULL;

  cudssConfigSet(solverConfig, CUDSS_CONFIG_REORDERING_ALG,
                 &reorder_alg, sizeof(cudssAlgType_t));

  cudssDataCreate(handle, &solverData);

#else
  Process::exit("Sorry, cuDSS solvers not available with this build.");
#endif
}


int Solv_cuDSS::resoudre_systeme(const Matrice_Base& a, const DoubleVect& bvect, DoubleVect& xvect)
{
#ifdef cuDSS_
  if (nouvelle_matrice())
    {

      /* conversion matric base to csr */
      Matrice_Morse csr;
      construit_matrice_morse_intermediaire(a, csr);
      /* create cudss matrixfrom pointers */
      Create_objects(csr, bvect, xvect);

      /* Symbolic factorization */
      cudssExecute(handle, CUDSS_PHASE_ANALYSIS, solverConfig, solverData,
                   A, x, b);

      /* Factorization */
      cudssExecute(handle, CUDSS_PHASE_FACTORIZATION, solverConfig,
                   solverData, A, x, b);
    }

  /* Solving */
  cudssExecute(handle, CUDSS_PHASE_SOLVE, solverConfig, solverData,
               A, x, b);
  return 0;
#else
  Process::exit("Sorry, cuDSS solvers not available with this build.");
  return -1;
#endif
}

void Solv_cuDSS::Create_objects(const Matrice_Morse& csr, const DoubleVect& bvect, DoubleVect& xvect)
{
#ifdef cuDSS_
  if (Process::is_parallel())
    {
      Process::exit("Sorry, cuDSS is a sequential solver.");
    }

  /* get dimensions and nnz */
  n=csr.get_tab1().size();
  nnz=csr.get_coeff().size();

  /* get pointers */
  csr_offsets_d = const_cast<int*>(csr.get_tab1().view_ro<1>().data());
  csr_columns_d = const_cast<int*>(csr.get_tab2().view_ro<1>().data());
  csr_values_d  = const_cast<double*>(csr.get_coeff().view_ro<1>().data());
  x_values_d = const_cast<double*>(xvect.view_rw<1>().data());
  b_values_d = const_cast<double*>(bvect.view_ro<1>().data());

  /* create matrix */
  cudssMatrixCreateCsr(&A, n, n, nnz, csr_offsets_d, nullptr,
                       csr_columns_d, csr_values_d, CUDA_R_32I, CUDA_R_64F, mtype, mview,
                       base);

  /* create rhs and vector */
  cudssMatrixCreateDn(&b, n, nrhs, n, b_values_d, CUDA_R_64F,
                      CUDSS_LAYOUT_COL_MAJOR);
  cudssMatrixCreateDn(&x, n, nrhs, n, x_values_d, CUDA_R_64F,
                      CUDSS_LAYOUT_COL_MAJOR);
#endif
}

//Comment différencier la premiere fois ou on cree les objet et les fois dapres ou on change les valeurs ?
//recreer matrices à chaque fois ?
