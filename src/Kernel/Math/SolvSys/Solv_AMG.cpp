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

#include <Solv_AMG.h>
#include <EChaine.h>
#include <Motcle.h>
#include <Solv_AMGX.h>
#include <Solv_Petsc_GPU.h>
#ifdef TRUST_USE_ROCM
#include <rocm-core/rocm_version.h>
#endif
#include <comm_incl.h> // Mandatory to have MPIX_CUDA_AWARE_SUPPORT defined or not
#include <MD_Vector_composite.h>

Implemente_instanciable(Solv_AMG,"Solv_AMG",SolveurSys_base);
// XD amg solveur_sys_base amg 0 Wrapper for AMG preconditioner-based solver which switch for the best one on CPU/GPU Nvidia/GPU AMD
// XD attr solveur chaine solveur 0 not_set
// XD attr option_solveur bloc_lecture option_solveur 0 not_set

// printOn
Sortie& Solv_AMG::printOn(Sortie& s ) const
{
  s << chaine_lue_;
  return s;
}

/**
 * @brief Reads the configuration for the AMG solver from the input stream.
 *
 * This function parses the input stream to configure the AMG solver parameters,
 * including the relative tolerance (RTOL), absolute tolerance (ATOL), and whether
 * to print additional information (IMPR). It supports different solver libraries
 * based on the available hardware (CPU or GPU).
 *
 * The expected input format is:
 *   amg solver { rtol value [impr] }
 *
 * @param is The input stream from which to read the solver configuration.
 * @return The input stream after reading the configuration.
 *
 * @throws Process::exit if the input syntax is incorrect or if an unsupported
 * library is specified.
 */

Nom boomeramg(double st)
{
  Nom chaine(" { precond boomeramg { }");
  if (st>=0)
    {
      chaine += " cli { -pc_hypre_boomeramg_strong_threshold";
      chaine += Nom(st, "%e");
      chaine += " }";
    }
  return chaine;
}

Nom pcfieldsplit_amgx(int n, double rtol, double atol)
{
  Nom chaine("cli { -ksp_type cg ");
  chaine+=rtol>0 ? "-ksp_rtol " : "-ksp_atol ";
  chaine+=rtol>0 ? Nom(rtol, "%e ") : Nom(atol, "%e ");
  chaine+="-ksp_norm_type UNPRECONDITIONED \
-pc_type fieldsplit \
-pc_fieldsplit_type additive \
-fieldsplit_P0_ksp_type preonly \
-fieldsplit_P0_pc_type amgx \
-fieldsplit_P1_ksp_type preonly \
-fieldsplit_P1_pc_type amgx ";
  if (n==3)
    {
      // ToDo: not efficient on P0P1Pa
      chaine+="-fieldsplit_P2_ksp_type preonly \
-fieldsplit_P2_pc_type amgx ";
    }
  chaine +="}";
  return chaine;
}

Nom pcfieldsplit_boomeramg(int n, double rtol, double atol)
{
  Nom chaine("cli { -ksp_type cg ");
  chaine+=rtol>0 ? "-ksp_rtol " : "-ksp_atol ";
  chaine+=rtol>0 ? Nom(rtol, "%e ") : Nom(atol, "%e ");
  chaine+="-ksp_norm_type UNPRECONDITIONED \
-pc_type fieldsplit \
-pc_fieldsplit_type additive \
-fieldsplit_P0_ksp_type preonly \
-fieldsplit_P0_pc_type hypre \
-fieldsplit_P0_pc_hypre_type boomeramg \
-fieldsplit_P0_pc_hypre_boomeramg_strong_threshold 0.1 \
-fieldsplit_P0_pc_hypre_boomeramg_print_statistics 1 \
-fieldsplit_P1_ksp_type preonly \
-fieldsplit_P1_pc_type hypre \
-fieldsplit_P1_pc_hypre_type boomeramg \
-fieldsplit_P1_pc_hypre_boomeramg_strong_threshold 0.1 \
-fieldsplit_P1_pc_hypre_boomeramg_print_statistics 1 ";
  if (n==3)
    {
      // ToDo: not efficient on P0P1Pa
      chaine+="-fieldsplit_P2_ksp_type preonly \
-fieldsplit_P2_pc_type hypre \
-fieldsplit_P2_pc_hypre_type boomeramg \
-fieldsplit_P2_pc_hypre_boomeramg_strong_threshold 0.1 \
-fieldsplit_P2_pc_hypre_boomeramg_print_statistics 1 ";
    }
  chaine +="}";
  return chaine;
}

Nom pcfieldsplit_gamg(int n, double rtol, double atol)
{
  Nom chaine("cli { -ksp_type cg ");
  chaine+=rtol>0 ? "-ksp_rtol " : "-ksp_atol ";
  chaine+=rtol>0 ? Nom(rtol, "%e ") : Nom(atol, "%e ");
  chaine+="-ksp_norm_type UNPRECONDITIONED \
-pc_type fieldsplit \
-pc_fieldsplit_type additive \
-fieldsplit_P0_ksp_type preonly \
-fieldsplit_P0_pc_type gamg \
-fieldsplit_P0_pc_gamg_threshold 0.01 \
-fieldsplit_P0_pc_gamg_square_graph 1 \
-fieldsplit_P1_ksp_type preonly \
-fieldsplit_P1_pc_type gamg \
-fieldsplit_P1_pc_gamg_threshold 0.01 \
-fieldsplit_P1_pc_gamg_square_graph 1 ";
  if (n==3)
    {
      // ToDo: not efficient on P0P1Pa
      chaine+="-fieldsplit_P2_ksp_type preonly \
-fieldsplit_P2_ksp_type preonly \
-fieldsplit_P2_pc_type gamg \
-fieldsplit_P2_pc_gamg_threshold 0.01 \
-fieldsplit_P2_pc_gamg_square_graph 1 ";
    }
  chaine +="}";
  return chaine;
}

Entree& Solv_AMG::readOn(Entree& is)
{
  // amg GCP|BISGTSTAB|GMRES { atol|rtol doublee [st double] [impr]  }
  Nom options_petsc("");
  Motcle motcle;
  Nom solver;
  is >> motcle;
  while (motcle != "}")
    {
      if (motcle=="GCP" || motcle=="BICGSTAB" || motcle=="GMRES") solver=motcle;
      else if (motcle=="{") {}
      else if (motcle=="RTOL") is >> rtol_;
      else if (motcle=="ATOL") is >> atol_;
      else if (motcle=="ST") is >> st_;
      else if (motcle=="IMPR") impr_ = true;
      else
        {
          options_petsc+=motcle;
          options_petsc+=" ";
        }
      is >> motcle;
    }
  // We select the more efficient/robust one:
  chaine_lue_ = solver;
#if defined(TRUST_USE_CUDA)
  library_ = "petsc_gpu";
  chaine_lue_ += boomeramg(st_); // Best GPU solver
  // KSP divergence with cg+boomeramg on multi-node with MPI Cuda Aware so we switch to AmgX:
  // Or switch to bcgs from cg ? Works !!! Strangely KSPSolve is 2x-3x slower on A100X vs MI250X... rocsparse better than cusparse ? And Kokkos-Kernels ?
  // Or use cg+gamg cg+amgx ?
#if defined(MPIX_CUDA_AWARE_SUPPORT)
  if (Process::nproc()>4)
    {
      library_ = "amgx";
      chaine_lue_ = solver;
      chaine_lue_ += " { precond c-amg {";    // Best GPU solver on Nvidia if MPI-GPU Aware OpenMPI 4.x
      if (st_>=0)
        {
          chaine_lue_ += " p:strength_threshold ";
          chaine_lue_ += Nom(st_, "%e");
        }
      chaine_lue_ += " }";
    }
#endif
#elif defined(TRUST_USE_ROCM)
  library_ = "petsc_gpu";
  const char* value = std::getenv("ROCM_ARCH");
  if (value != nullptr && std::string(value) == "gfx1100")
    {
      if (st_>=0) Process::exit("st option not supported yet in Solv_AMG");
      if (Process::is_parallel())
        chaine_lue_ += " { precond ua-amg { }";  // Converge mais plus lent que sa-amg
      else
        chaine_lue_ += " { precond sa-amg { }";  // Crash en parallele
    }
  else
    chaine_lue_ += boomeramg(st_); // Best GPU solver (// sa-amg is slow...)
#else
  library_ = "petsc";
  chaine_lue_ += boomeramg(st_); // Best CPU solver
#endif
  if (rtol_>0)
    {
      chaine_lue_ += " rtol ";
      chaine_lue_ += Nom(rtol_, "%e");
    }
  if (atol_>0)
    {
      chaine_lue_ += " atol ";
      chaine_lue_ += Nom(atol_, "%e");
    }
  if (impr_) chaine_lue_ += " impr";
  if (options_petsc!="")
    {
      chaine_lue_ += " ";
      chaine_lue_ += options_petsc;
    }
  chaine_lue_ += " }";
  return is;
}

int Solv_AMG::resoudre_systeme(const Matrice_Base& mat, const DoubleVect& b, DoubleVect& x)
{
  // We don't create solver during readOn as usual but just before solve to get more infos about matrix/vectors to fine tune
  if (solveur_.est_nul())
    {
      if (sub_type(MD_Vector_composite, b.get_md_vector().valeur()))
        {
          int nb_blocks = ref_cast(MD_Vector_composite, b.get_md_vector().valeur()).nb_parts();
          if (nb_blocks>1)
            {
              // Block matrix : we use PCFieldsplit (eg: VEF) for preconditioner
              // Much better convergence for P0P1 for instance
              Cerr << "Detecting " << nb_blocks << "x" << nb_blocks << " blocks into the matrix..." << finl;
              if (chaine_lue_.contient("boomeramg"))
                chaine_lue_ = pcfieldsplit_boomeramg(nb_blocks, rtol_, atol_);
              else if (chaine_lue_.contient("gamg"))
                chaine_lue_ = pcfieldsplit_gamg(nb_blocks, rtol_, atol_);
              else if (library_=="amgx")
                {
                  library_ = "petsc_gpu";
                  chaine_lue_ = pcfieldsplit_amgx(nb_blocks, rtol_, atol_);
                }
            }
        }
      Cerr << "====================================================================" << finl;
      Cerr << "Creating AMG solver: " << library_ << " " << chaine_lue_ << finl;
      Cerr << "====================================================================" << finl;
      EChaine entree(chaine_lue_);
      Nom nom_solveur("Solv_");
      nom_solveur+=library_;
      solveur_.typer(nom_solveur);
      solveur_.nommer("solveur_pression");
      if (library_=="amgx")
        ref_cast(Solv_AMGX, solveur_.valeur()).create_solver(entree);
      else if (library_=="petsc")
        ref_cast(Solv_Petsc, solveur_.valeur()).create_solver(entree);
      else if (library_=="petsc_gpu")
        ref_cast(Solv_Petsc_GPU, solveur_.valeur()).create_solver(entree);
      else
        Process::exit("Unsupported case in Solv_AMG::readOn");
    }
  return solveur_.resoudre_systeme(mat, b, x);
}
