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

#define CUDSS_CALL_AND_CHECK(call, status, msg) \
    do { \
        status = call; \
        std::cout<<call<<std::endl; \
        if (status != CUDSS_STATUS_SUCCESS) { \
            printf("Example FAILED: CUDSS call ended unsuccessfully with status = %d, details: " #msg "\n", status); \
            return -2; \
        } \
    } while(0);

#define CUDSS_CALL_AND_CHECK_VOID(call, status, msg) \
    do { \
        status = call; \
                std::cout<<call<<std::endl; \
        if (status != CUDSS_STATUS_SUCCESS) { \
            printf("Example FAILED: CUDSS call ended unsuccessfully with status = %d, details: " #msg "\n", status); \
        } \
    } while(0);


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

  if (Axb_are_built)
    {
      cudssMatrixDestroy(A);
      cudssMatrixDestroy(b);
      cudssMatrixDestroy(x);
    }

  if (solver_is_built)
    {
      cudssDataDestroy(handle, solverData);
      cudssConfigDestroy(solverConfig);
      cudssDestroy(handle);
    }
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

  Cerr <<"[cuDSS] reading jdd and creating solver "<< finl;

  lecture(entree);
  EChaine is(get_chaine_lue());
  Motcle accolade_ouverte("{"), accolade_fermee("}");
  Motcle solver, motlu;
  is >> solver;   // On lit le solveur en premier puis les options du solveur
  is >> motlu; // On lit l'accolade
  if (motlu != accolade_ouverte)
    {
      Cerr << "[cuDSS] Error while reading the parameters of the solver " << solver << " { ... }" << finl;
      Cerr << "We expected " << accolade_ouverte << " instead of " << motlu << finl;
      exit();
    }

  if (solver==(Motcle)"LU")
    {
      std::cout<<"[cuDSS] LU is read, matrix is assumed non symmetric"<<std::endl;
      mtype=CUDSS_MTYPE_GENERAL;
      mview=CUDSS_MVIEW_FULL;
    }
  else if (solver==(Motcle)"Cholesky")
    {
      std::cout<<"[cuDSS] Cholesky is read, matrix is assumed symmetric"<<std::endl;
      mtype=CUDSS_MTYPE_SYMMETRIC;
      mview=CUDSS_MVIEW_UPPER;
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

          if (motlu==(Motcle)"0") {  reorder_alg = CUDSS_ALG_DEFAULT;}
          else if (motlu==(Motcle)"1") {  reorder_alg = CUDSS_ALG_1;}
          else if (motlu==(Motcle)"2") {  reorder_alg = CUDSS_ALG_2;}
          else if (motlu==(Motcle)"3") {  reorder_alg = CUDSS_ALG_3;}
          else if (motlu==(Motcle)"4") {  reorder_alg = CUDSS_ALG_4;}
          else if (motlu==(Motcle)"5") {  reorder_alg = CUDSS_ALG_5;}
          else
            {
              Cerr << "[cuDSS] algo number "<<motlu<<" unrecognized (0-5) only"<<finl;
              Process::exit();
            }
          Cerr << "[cuDSS] algo number "<<motlu<<" read. "<<finl;
          Cerr <<"[cuDSS] Note from cuDSS doc: Different values represent different algorithms (for reordering, factorization, etc.) and can lead to significant differences in accuracy and performance. It is currently recommended to use CUDSS_ALG_DEFAULT (0) and only in case accuracy or performance are not sufficient, one can experiment with other values."<<finl;

        }
      else
        {
          Cerr << motlu << "[cuDSS] keyword not recognized for cuDSS solver " << solver << finl;
          Process::exit();
        }
      is >> motlu;
    }

#else
  Process::exit("Sorry, cuDSS solvers not available with this build (create solver).");
#endif
}


int Solv_cuDSS::resoudre_systeme(const Matrice_Base& a, const DoubleVect& bvect, DoubleVect& xvect)
{
#ifdef cuDSS_

  if ((nouvelle_matrice()||first_solve))
    {

      /* conversion matric base to csr */
      Matrice_Morse csr;
      /*build the csr matrix on host*/
      construit_matrice_morse_intermediaire(a, csr);
      /* create cudss matrix and size them */
      Create_objects(csr);
      /* give the device data / csr pointers to the cudss matrix */
      set_pointers_A(csr);
      /* Symbolic factorization x, and b are unused*/
      CUDSS_CALL_AND_CHECK(cudssExecute(handle, CUDSS_PHASE_ANALYSIS, solverConfig, solverData,
                                        A, x, b), status, "cudssExecute for analysis");

      /* Factorization x, and b are unused*/
      CUDSS_CALL_AND_CHECK(cudssExecute(handle, CUDSS_PHASE_FACTORIZATION, solverConfig,
                                        solverData, A, x, b), status, "cudssExecute for facto");
    }

  /*some checks*/
  /* sizes should never change */
  assert(a.nb_lignes()==n);
  assert(bvect.size_totale()==n);
  assert(xvect.size_totale()==n);

  /* give the device data / csr pointers to the cudss x and b */
  /* does nothing after the first solve*/
  set_pointers_xb(bvect, xvect);

  /* Solving */
  CUDSS_CALL_AND_CHECK(cudssExecute(handle, CUDSS_PHASE_SOLVE, solverConfig, solverData,
                                    A, x, b), status, "cudssExecute for solve");

#ifndef NDEBUG
  DoubleVect test(bvect);
  test*=-1;
  a.ajouter_multvect(xvect,test);
  double vrai_residu = mp_norme_vect(test);
  Cout << "||Ax-b||=" << vrai_residu << finl;
#endif

  first_solve = false;

  return 0;
#else
  Process::exit("Sorry, cuDSS solvers not available with this build (resoudre_systeme).");
  return -1;
#endif
}

void Solv_cuDSS::Create_objects(const Matrice_Morse& csr)
{
#ifdef cuDSS_

  if (Process::is_parallel())
    {
      Process::exit("Sorry, cuDSS is a sequential solver.");
    }

  /* get dimensions and nnz */
  /* check that n does not change between solves */

#ifndef NDEBUG
  int new_n = csr.get_tab1().size()-1;
  assert(((new_n==n)||(first_solve))); // n should never change between solves
#endif

  n = csr.get_tab1().size()-1;

  int new_nnz = csr.get_coeff().size();

  //Re-build A only if new nnz is different from old one, possible if matrix changes
  if( new_nnz != nnz)
    {
      nnz=new_nnz;

      /* destroy if A is built */
      if (Axb_are_built)
        cudssMatrixDestroy(A);

      /* create matrix, we give nullptrs*/
      CUDSS_CALL_AND_CHECK_VOID(cudssMatrixCreateCsr(&A, n, n, nnz, nullptr, nullptr,
                                                     nullptr, nullptr, CUDA_R_32I, CUDA_R_64F, mtype, mview,
                                                     base), status, "cudssMatrixCreateCsr");
    }

  if (first_solve)
    {
      /* create rhs and vector we give nullptrs*/
      // n should never change between solves, no need to resize

      CUDSS_CALL_AND_CHECK_VOID(cudssMatrixCreateDn(&b, n, nrhs, n, nullptr, CUDA_R_64F,
                                                    CUDSS_LAYOUT_COL_MAJOR), status, "cudssMatrixCreateDn for b");
      CUDSS_CALL_AND_CHECK_VOID(cudssMatrixCreateDn(&x, n, nrhs, n, nullptr, CUDA_R_64F,
                                                    CUDSS_LAYOUT_COL_MAJOR), status, "cudssMatrixCreateDn for x");
    }

  Axb_are_built=true;


  if (first_solve)
    {
      /* Create a CUDA stream */
      cudaStreamCreate(&stream);
      /* Create handle */
      CUDSS_CALL_AND_CHECK_VOID(cudssCreate(&handle), status, "cudssCreate");
      /* (optional) Setting the custom stream for the library handle */
      CUDSS_CALL_AND_CHECK_VOID(cudssSetStream(handle, stream), status, "cudssSetStream");
      /* Create config */
      CUDSS_CALL_AND_CHECK_VOID(cudssConfigCreate(&solverConfig), status, "cudssConfigCreate");
      /* set options to config */
      CUDSS_CALL_AND_CHECK_VOID(cudssConfigSet(solverConfig, CUDSS_CONFIG_REORDERING_ALG,
                                               &reorder_alg, sizeof(cudssAlgType_t)), status, "cudssConfigSet");
      /* create solver data */
      CUDSS_CALL_AND_CHECK_VOID(cudssDataCreate(handle, &solverData), status, "cudssDataCreate");
      /*log*/
      solver_is_built = true;
    }
#endif
}

void Solv_cuDSS::set_pointers_A(const Matrice_Morse& csr)
{

  /* get pointers */
  csr_offsets_d = const_cast<int*>(csr.get_tab1().view_ro<1>().data());
  csr_columns_d = const_cast<int*>(csr.get_tab2().view_ro<1>().data());
  csr_values_d = const_cast<double*>(csr.get_coeff().view_ro<1>().data());

  /* give pointers to A*/
  CUDSS_CALL_AND_CHECK_VOID(cudssMatrixSetCsrPointers(A,
                                                      csr_offsets_d, /*rowStart*/
                                                      nullptr,/*rowEnd*/
                                                      csr_columns_d,
                                                      csr_values_d), status, "cudssMatrixSetCsrPointers")
}

void Solv_cuDSS::set_pointers_xb(const DoubleVect& bvect, DoubleVect& xvect)
{
  //The pointers to x and b should never change between calls to resoudre
  //In debug we check that this is the case

  if (first_solve)
    {
      /* get pointers */
      x_values_d = const_cast<double*>(xvect.view_rw<1>().data()); //maybe wo ?
      b_values_d = const_cast<double*>(bvect.view_ro<1>().data());

      /* give pointers to vectors*/
      CUDSS_CALL_AND_CHECK_VOID(cudssMatrixSetValues(x, x_values_d), status, "cudssMatrixSetValues for x");
      CUDSS_CALL_AND_CHECK_VOID(cudssMatrixSetValues(b, b_values_d), status, "cudssMatrixSetValues for b");
    }
#ifndef NDEBUG
  else
    {
      assert(x_values_d==const_cast<double*>(xvect.view_rw<1>().data()));
      assert(b_values_d==const_cast<double*>(bvect.view_ro<1>().data()));
    }
#endif
}

