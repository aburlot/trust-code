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
/////////////////////////////////////////////////////////////////////////////
//
// File      : Robin_VEF.cpp
// Directory : $TRUST_ROOT/src/Kernel/Cond_Lim
//
/////////////////////////////////////////////////////////////////////////////


#include <Robin_VEF.h>
#include <Param.h>
#include <Frontiere_dis_base.h>
#include <Front_VF.h>

Implemente_instanciable( Robin_VEF, "Robin_VEF", Cond_lim_base ) ;

Sortie& Robin_VEF::printOn( Sortie& os ) const
{
  Cond_lim_base::printOn( os );
  return os;
}

Entree& Robin_VEF::readOn( Entree& is )
{
  ;
  Param param(que_suis_je());
  // reading Robin border conditions parameters

  param.ajouter("alpha", 		&alpha_robin_cl_, 		Param::REQUIRED);
  param.ajouter("beta" , 		&beta_robin_cl_ , 		Param::REQUIRED);
  param.ajouter("champ_front_normal_et_tangentiel_robin",&le_champ_front,Param::REQUIRED); // Champ_front_tangentiel_robin::readOn(is);
  param.lire_avec_accolades_depuis(is);


  Cerr << "Reading Robin boundary condition, the value of the normal Robin coefficient is  "
       << alpha_robin_cl_ << " and the tangential Robin coefficient is  " << beta_robin_cl_ << finl;

  if (alpha_robin_cl_<=0 || beta_robin_cl_<=0)
    {
      Cerr << "Error of Robin boundary conditions : Alpha and Beta parameters should be positive  " << finl;
      exit() ;
    }


  return is; //Cond_lim_base::readOn(is);
}

int Robin_VEF::compatible_avec_eqn( const Equation_base& eqn) const
{
  // TODO : if we need to do some verification ?
  return 1;
}


double Robin_VEF::flux_robin_normal_et_trangentiel_imp(int i, int j) const
{
  // Si 2D : nb_dim = 3 : 1 col pour les indices des mailles, 1 col pour le champ normal, 1 col pour le champ tangent
  // Si 3D : nb_dim = 5 : 1 col pour les indices des mailles, 1 col pour le champ normal, 3 col pour le champ tangent



  if ((le_champ_front->valeurs().dimension(1)==2 && dimension == 2) || (le_champ_front->valeurs().dimension(1)==4 && dimension == 3))
    {
      if (le_champ_front->valeurs().dimension(0)==1) // a une valeur identique sur toutes les mailles du bord
        return le_champ_front->valeurs()(0,j);
      else
        {
          return le_champ_front->valeurs()(i,j);
        }
    }
  else
    Cerr << "Robin::flux_robin_normal_et_tangentiel error" << finl;
  Process::exit();
  return 0.;
}


double Robin_VEF::flux_normal_imp(int i) const {	return flux_robin_normal_et_trangentiel_imp(i, 0);}


double Robin_VEF::flux_tangentiel_imp(int i, int j ) const
{
  if ((dimension==2 && j<2) || (dimension ==3 && j<4))
    return flux_robin_normal_et_trangentiel_imp(i, j+1);
  else
    Cerr << " Robin::flux_tangentiel_imp error, dimension does not match" << finl;
  Process::exit();
  return 0.;
}





double Robin_VEF::flux_robin_imp_au_temps(double temps, int i) const
{
  return 0.;

}

double Robin_VEF::flux_robin_imp_au_temps(double temps, int i, int j ) const
{
  return 0.;

}



/*! @brief Retourne le tableau flux_robin_impose_ mis a jour
 *
 */
const DoubleTab& Robin_VEF::flux_robin_normal_imp() const
{
  const Front_VF& le_bord = ref_cast(Front_VF, frontiere_dis());
  int nb_faces_tot = le_bord.nb_faces_tot();
  if (nb_faces_tot>0)
    {
      if (flux_normal_impose_.dimension(0) != nb_faces_tot)
        flux_normal_impose_.resize(nb_faces_tot, 1);
      int size = flux_normal_impose_.dimension(0);
      //int nb_comp = flux_robin_impose_.dimension(1);
      for (int i = 0; i < size; i++)
        //for (int j = 0; j < nb_comp; j++)
        flux_normal_impose_(i, 0) = flux_normal_imp(i);
    }
  return flux_normal_impose_;
}


const DoubleTab& Robin_VEF::flux_robin_tangentiel_imp() const
{
  const Front_VF& le_bord = ref_cast(Front_VF, frontiere_dis());
  int nb_faces_tot = le_bord.nb_faces_tot();
  if (nb_faces_tot>0)
    {
      if (flux_tangentiel_impose_.dimension(0) != nb_faces_tot)
        flux_tangentiel_impose_.resize(nb_faces_tot, le_champ_front->valeurs().dimension(1));
      int size = flux_tangentiel_impose_.dimension(0);
      int nb_comp = flux_tangentiel_impose_.dimension(1);
      for (int i = 0; i < size; i++)
        for (int j = 0; j < nb_comp; j++)
          flux_tangentiel_impose_(i, 0) = flux_tangentiel_imp(i,j);
    }
  return flux_tangentiel_impose_;
}

