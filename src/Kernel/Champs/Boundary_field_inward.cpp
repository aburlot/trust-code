/****************************************************************************
* Copyright (c) 2015, CEA
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
//////////////////////////////////////////////////////////////////////////////
//
// File:        Boundary_field_inward.cpp
// Directory:   $TRUST_ROOT/src/Kernel/Champs
// Version:     /main/8
//
//////////////////////////////////////////////////////////////////////////////

#include <Boundary_field_inward.h>
#include <Motcle.h>
#include <Zone_VF.h>
#include <Linear_algebra_tools.h>


Implemente_instanciable(Boundary_field_inward,"Boundary_field_inward",Champ_front_normal);


Sortie& Boundary_field_inward::printOn(Sortie& os) const
{
  return Champ_front_normal::printOn(os);
}


Entree& Boundary_field_inward::readOn(Entree& is)
{
  return Champ_front_normal::readOn(is);
}


// Description:
// Initialize normal inward field given on boundaries

// Remark :
// This code is not in "Champ_front_normal::initialiser" function because in the case of using
// Champ_front_debit, the readOn function of Champ_front_normal is not called
// so vit_norm is not initialized and it can leads to critical values detected
// in debug mode.
// Moreover, when using Champ_front_debit we don't need to do all these calculations
// That is why Boundary_field_inward was created
int Boundary_field_inward::initialiser(double tps, const Champ_Inc_base& inco)
{
  if( !Champ_front_normal::initialiser(tps,inco) )
    return 0;

  const Zone_VF& zone_VF = ref_cast(Zone_VF,zone_dis());
  const Front_VF& le_bord= ref_cast(Front_VF,frontiere_dis());
  int first=le_bord.num_premiere_face();

  int nb_cases=les_valeurs->nb_cases();

  DoubleTab tab_normal_vector(le_bord.nb_faces(),dimension);

  for(int i=0; i<le_bord.nb_faces(); i++)
    {

      int signe = -1;//we want the opposite of the normal vector
      int global_boundary_face_number=i+first;
      //outward vector associated to the fac global_boundary_face_number
      DoubleTab normal_vector = zone_VF.normalized_boundaries_outward_vector(global_boundary_face_number,signe);

      for (int k=0; k<dimension; k++)
        {
          tab_normal_vector(i,k) = normal_vector(k);
        }

      // Fill all moments
      for (int t=0; t<nb_cases; t++)
        {
          DoubleTab& tab=les_valeurs[t].valeurs();
          for(int j=0; j<dimension; j++)
            {
              tab(i,j)=tab_normal_vector(i,j);
              tab(i,j)*=vit_norm;
            }
        }
    }
  return 1;
}


// Obsolete soon:
Implemente_instanciable(Champ_front_normal_VEF,"Champ_front_normal_VEF",Boundary_field_inward);

Sortie& Champ_front_normal_VEF::printOn(Sortie& os) const
{
  return Boundary_field_inward::printOn(os);
}

Entree& Champ_front_normal_VEF::readOn(Entree& is)
{
  Motcle motlu;
  Motcles les_mots(1);
  les_mots[0]="valeur_normale";
  is>>motlu;
  if (motlu != les_mots[0])
    {
      Cerr << "Error while reading a Champ_front_normal" << finl;
      Cerr << "We expected " << les_mots[0] << " instead of " << motlu << finl;
      exit();
    }

  is >> vit_norm;
  fixer_nb_comp(dimension);
  return is;
}
