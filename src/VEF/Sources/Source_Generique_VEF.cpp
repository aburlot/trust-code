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
// File:        Source_Generique_VEF.cpp
// Directory:   $TRUST_ROOT/src/VEF/Sources
// Version:     /main/5
//
//////////////////////////////////////////////////////////////////////////////

#include <Source_Generique_VEF.h>
#include <Zone_VEF.h>
#include <Zone_Cl_VEF.h>
#include <Zone_Cl_dis.h>
#include <Dirichlet.h>
#include <Dirichlet_homogene.h>

Implemente_instanciable(Source_Generique_VEF,"Source_Generique_VEF_P1NC",Source_Generique_base);


Sortie& Source_Generique_VEF::printOn(Sortie& os) const
{
  return os << que_suis_je() ;
}

Entree& Source_Generique_VEF::readOn(Entree& is)
{
  Source_Generique_base::readOn(is);
  return is;
}

DoubleTab& Source_Generique_VEF::ajouter(DoubleTab& resu) const
{
  Champ espace_stockage;
  const Champ_base& champ_calc = ch_source_->get_champ(espace_stockage);
  const DoubleTab& valeurs_calc = champ_calc.valeurs();

  int nb_faces = la_zone_VEF->nb_faces();
  const DoubleVect& vol_entrelaces = la_zone_VEF->volumes_entrelaces();
  const DoubleVect& vol_entrelaces_Cl = la_zcl_VEF->volumes_entrelaces_Cl();
  const DoubleVect& poro_face = la_zone_VEF->porosite_face();
  int nb_front_cl = la_zone_VEF->nb_front_Cl();
  int premiere_face_interne = la_zone_VEF->premiere_face_int();
  int num_face;

  int nb_comp = 1;
  int nb_dim = resu.nb_dim();

  if(nb_dim==2)
    nb_comp = resu.dimension(1);

  if (nb_dim==1)
    {
      for (int num_cl =0; num_cl<nb_front_cl; num_cl++)
        {
          const Cond_lim& la_cl = la_zcl_VEF->les_conditions_limites(num_cl);
          const Front_VF& le_bord = ref_cast(Front_VF,la_cl.frontiere_dis());
          int ndeb = le_bord.num_premiere_face();
          int nfin = ndeb + le_bord.nb_faces();

          if ((sub_type(Dirichlet,la_cl.valeur()))
              ||
              (sub_type(Dirichlet_homogene,la_cl.valeur()))
             )
            ;
          else
            for (num_face=ndeb; num_face<nfin; num_face++)
              resu(num_face) += valeurs_calc(num_face)*vol_entrelaces_Cl(num_face)*poro_face(num_face);
        }

      for (num_face=premiere_face_interne; num_face<nb_faces; num_face++)
        resu(num_face) += valeurs_calc(num_face)*vol_entrelaces(num_face)*poro_face(num_face);
    }
  else
    {
      for (int num_cl =0; num_cl<nb_front_cl; num_cl++)
        {
          const Cond_lim& la_cl = la_zcl_VEF->les_conditions_limites(num_cl);
          const Front_VF& le_bord = ref_cast(Front_VF,la_cl.frontiere_dis());
          int ndeb = le_bord.num_premiere_face();
          int nfin = ndeb + le_bord.nb_faces();

          if ((sub_type(Dirichlet,la_cl.valeur()))
              ||
              (sub_type(Dirichlet_homogene,la_cl.valeur()))
             )
            ;
          else
            for (num_face=ndeb; num_face<nfin; num_face++)
              {
                for (int nc=0; nc<nb_comp; nc++)
                  resu(num_face,nc) += valeurs_calc(num_face,nc)*vol_entrelaces_Cl(num_face)*poro_face(num_face);
              }
        }

      for (num_face=premiere_face_interne; num_face<nb_faces; num_face++)
        for (int nc=0; nc<nb_comp; nc++)
          resu(num_face,nc) += valeurs_calc(num_face,nc)*vol_entrelaces(num_face)*poro_face(num_face);
    }

  return resu;
}

void Source_Generique_VEF::associer_zones(const Zone_dis& zone_dis,
                                          const Zone_Cl_dis& zcl_dis)
{
  la_zone_VEF = ref_cast(Zone_VEF,zone_dis.valeur());
  la_zcl_VEF = ref_cast(Zone_Cl_VEF,zcl_dis.valeur());
}

Nom Source_Generique_VEF::localisation_source()
{
  return "faces";
}
