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

#include <Op_Conv_EF_Stab_PolyMAC_Face.h>
#include <Discretisation_base.h>
#include <Dirichlet_homogene.h>
#include <Champ_Face_PolyMAC.h>
#include <Domaine_Cl_PolyMAC.h>
#include <Schema_Temps_base.h>
#include <Domaine_PolyMAC.h>
#include <Probleme_base.h>
#include <Matrix_tools.h>
#include <Array_tools.h>
#include <TRUSTLists.h>
#include <Dirichlet.h>
#include <Param.h>
#include <cmath>

Implemente_instanciable( Op_Conv_EF_Stab_PolyMAC_Face, "Op_Conv_EF_Stab_PolyMAC_Face_PolyMAC", Op_Conv_PolyMAC_base );
Implemente_instanciable( Op_Conv_Amont_PolyMAC_Face, "Op_Conv_Amont_PolyMAC_Face_PolyMAC", Op_Conv_EF_Stab_PolyMAC_Face );
Implemente_instanciable( Op_Conv_Centre_PolyMAC_Face, "Op_Conv_Centre_PolyMAC_Face_PolyMAC", Op_Conv_EF_Stab_PolyMAC_Face );

// XD Op_Conv_EF_Stab_PolyMAC_Face interprete Op_Conv_EF_Stab_PolyMAC_Face 1 Class Op_Conv_EF_Stab_PolyMAC_Face_PolyMAC

Sortie& Op_Conv_EF_Stab_PolyMAC_Face::printOn(Sortie& os) const { return Op_Conv_PolyMAC_base::printOn(os); }
Sortie& Op_Conv_Amont_PolyMAC_Face::printOn(Sortie& os) const { return Op_Conv_PolyMAC_base::printOn(os); }
Sortie& Op_Conv_Centre_PolyMAC_Face::printOn(Sortie& os) const { return Op_Conv_PolyMAC_base::printOn(os); }

Entree& Op_Conv_EF_Stab_PolyMAC_Face::readOn(Entree& is)
{
  Op_Conv_PolyMAC_base::readOn(is);
  Param param(que_suis_je());
  param.ajouter("alpha", &alpha_);            // XD_ADD_P double parametre ajustant la stabilisation de 0 (schema centre) a 1 (schema amont)
  param.lire_avec_accolades_depuis(is);
  return is;
}

Entree& Op_Conv_Amont_PolyMAC_Face::readOn(Entree& is)
{
  alpha_ = 1.0;
  return Op_Conv_PolyMAC_base::readOn(is);
}

Entree& Op_Conv_Centre_PolyMAC_Face::readOn(Entree& is)
{
  alpha_ = 0.0;
  return Op_Conv_PolyMAC_base::readOn(is);
}

double Op_Conv_EF_Stab_PolyMAC_Face::calculer_dt_stab() const
{
  double dt = 1e10;
  const Domaine_Poly_base& domaine = le_dom_poly_.valeur();
  const DoubleVect& fs = domaine.face_surfaces(), &pf = equation().milieu().porosite_face(), &ve = domaine.volumes(), &pe = equation().milieu().porosite_elem();
  const DoubleTab& vit = vitesse_->valeurs();
  const IntTab& e_f = domaine.elem_faces(), &f_e = domaine.face_voisins();
  const int N = vit.line_size();
  DoubleTrav flux(N); //somme des flux pf * |f| * vf, volume minimal des mailles d'elements/faces affectes par ce flux

  for (int e = 0; e < domaine.nb_elem(); e++)
    {
      // Calcul du volume effectif de l'element
      const double vol = pe(e) * ve(e);
      flux = 0.;

      // Parcourt des faces associees a l'element
      for (int i = 0; i < e_f.dimension(1); i++)
        {
          int f = e_f(e, i);
          if (f < 0) continue; // face in-existante

          for (int n = 0; n < N; n++)
            {
              // Ajout du flux entrant pour la composante n : Seuls les flux entrants comptent
              double flux_f = pf(f) * fs(f) * std::max((e == f_e(f, 1) ? 1 : -1) * vit(f, n), 0.);
              flux(n) += flux_f;
            }
        }

      // Calcul du pas de temps pour chaque composante n
      for (int n = 0; n < N; n++)
        if (std::abs(flux(n)) > 1e-12)
          dt = std::min(dt, vol / flux(n));
    }

  return Process::mp_min(dt);
}

void Op_Conv_EF_Stab_PolyMAC_Face::completer()
{
  Op_Conv_PolyMAC_base::completer();

  /* au cas ou... */
  const Domaine_PolyMAC& domaine = le_dom_poly_.valeur();
  if (domaine.domaine().nb_joints() && domaine.domaine().joint(0).epaisseur() < 2)
    {
      Cerr << "Op_Conv_EF_Stab_PolyMAC_Face : largeur de joint insuffisante (minimum 2)!" << finl;
      Process::exit();
    }
  porosite_f.ref(mon_equation->milieu().porosite_face());
  porosite_e.ref(mon_equation->milieu().porosite_elem());
}

void Op_Conv_EF_Stab_PolyMAC_Face::dimensionner(Matrice_Morse& mat) const
{
  if (has_interface_blocs())
    {
      Operateur_base::dimensionner(mat);
      return;
    }

  const Domaine_PolyMAC& domaine = le_dom_poly_.valeur();
  const Champ_Face_PolyMAC& ch = ref_cast(Champ_Face_PolyMAC, equation().inconnue());
  const IntTab& e_f = domaine.elem_faces(), &f_e = domaine.face_voisins(), &equiv = domaine.equiv();
  const DoubleTab& xp = domaine.xp(), &xv = domaine.xv();
  const DoubleVect& fs = domaine.face_surfaces(), &vf = domaine.volumes_entrelaces();

  ch.fcl(), domaine.init_ve();

  IntTab stencil(0, 2);

  for (int e = 0; e < domaine.nb_elem_tot(); e++)
    {
      // Calcul de la divergence et contributions aux equations
      for (int i = 0; i < e_f.dimension(1); i++)
        {
          int f = e_f(e, i);
          if (f < 0) continue;

          // Contributions des faces voisines
          if (f < domaine.nb_faces() && ch.fcl()(f, 0) < 2)
            {
              for (int j = 0; j < e_f.dimension(1); j++)
                {
                  int fb = e_f(e, j);
                  if (fb < 0) continue;

                  for (int k = 0; k < 2; k++)
                    {
                      if (ch.fcl()(fb, 0) < 2 || ch.fcl()(fb, 0) == 3)
                        {
                          int eb = f_e(fb, k);
                          if (eb < 0) continue;

                          int fc = equiv(fb, e != f_e(fb, 0), i);
                          if (fc  >= 0 )
                            stencil.append_line(f, fc);
                          else
                            {
                              // Convection pour les faces internes sans equivalence
                              for (int l = domaine.vedeb(eb); l < domaine.vedeb(eb + 1); l++)
                                {
                                  fc = domaine.veji(l);
                                  if (ch.fcl()(fc, 0) < 2 && std::fabs(domaine.dot(&xv(f, 0), &domaine.veci(l, 0), &xp(e, 0))) > 1e-8 * vf(f) / fs(f))
                                    stencil.append_line(f, fc);
                                }
                            }
                        }
                    }
                }
            }
        }
    }

  // Tri et suppression des doublons
  tableau_trier_retirer_doublons(stencil);

  // Allocation de la matrice
  int taille = domaine.nb_faces_tot() + (dimension < 3 ? domaine.domaine().nb_som_tot() : domaine.domaine().nb_aretes_tot());
  Matrix_tools::allocate_morse_matrix(taille, taille, stencil, mat);
}

// ajoute la contribution de la convection au second membre resu
// renvoie resu
inline DoubleTab& Op_Conv_EF_Stab_PolyMAC_Face::ajouter(const DoubleTab& tab_inco, DoubleTab& tab_resu) const
{
  if (has_interface_blocs())
    return Operateur_base::ajouter(tab_inco, tab_resu);

  const Domaine_PolyMAC& domaine = le_dom_poly_.valeur();
  const Champ_Face_PolyMAC& ch = ref_cast(Champ_Face_PolyMAC, equation().inconnue());
  const Conds_lim& cls = la_zcl_poly_->les_conditions_limites();

  domaine.init_ve();

  const IntTab& f_e = domaine.face_voisins(), &e_f = domaine.elem_faces(), &equiv = domaine.equiv();
  const DoubleTab& xp = domaine.xp(), &xv = domaine.xv(), &vfd = domaine.volumes_entrelaces_dir(), &vit = vitesse_->valeurs();
  const DoubleVect& fs = domaine.face_surfaces(), &ve = domaine.volumes(),  &pf = porosite_f, &pe = porosite_e;
  const DoubleTab& nf = domaine.face_normales();

  for (int e = 0; e < domaine.nb_elem_tot(); e++)
    {
      for (int i = 0; i < e_f.dimension(1); i++)
        {
          int f = e_f(e, i);
          if (f < 0) continue;

          int dir_f = (e == f_e(f, 0) ? 1 : -1);

          // Traitement de la convection pour cette face
          if (f < domaine.nb_faces() && ch.fcl()(f, 0) < 2)
            for (int j = 0; j < e_f.dimension(1); j++)
              {
                int fb = e_f(e, j);
                if (fb < 0) continue;

                int dir_fb = (e == f_e(fb, 0) ? 1 : -1);

                for (int k = 0; k < 2; k++)
                  {
                    int eb = f_e(fb, k);
                    double fac = dir_f * tab_inco(fb) * dir_fb * fs(fb) * pe(eb >= 0 ? eb : f_e(fb, 0)) / ve(e) * (1. + (vit(fb) * (k ? -1 : 1) >= 0 ? 1. : vit[fb] ? -1. : 0.) * alpha_) / 2;

                    int fc = equiv(fb, e != f_e(fb, 0), i);
                    if (fc  >= 0 || f_e(fb, 0) < 0 || f_e(fb, 1) < 0)
                      {
                        if (eb >= 0)
                          {
                            int fd = (eb == e ? f : fc); //face source
                            double mult = (fd < 0 || domaine.dot(&nf(f, 0), &nf(fd, 0)) > 0) ? 1 : -1;
                            mult *= (fd >= 0) ? pf(fd) / pe(eb) : 1;

                            tab_resu(f) -= fac * mult * dir_f * vfd(f, e != f_e(f, 0)) * vit(fd);

                          }
                        else if (ch.fcl()(fb, 0) == 3 )
                          {
                            double masse = std::fabs(vit[fb]) > 1e-10 ? tab_inco(fb) / vit[fb] : 1.0;
                            for (int l = 0; l < dimension; l++)
                              tab_resu(f) -= fac * dir_f * vfd(f, e != f_e(f, 0))* nf(f, l) / fs(f) * ref_cast(Dirichlet, cls[ch.fcl()(fb, 1)].valeur()).val_imp(ch.fcl()(fb, 2), l) / masse;
                          }
                        if (!incompressible_)
                          tab_resu(f) += fac * dir_f * vfd(f, e != f_e(f, 0)) * vit(f);
                      }
                    else
                      {
                        for (int l = domaine.vedeb(eb); l < domaine.vedeb(eb + 1); l++)
                          {
                            fc = domaine.veji(l);
                            tab_resu(f) -= fac * fs(f) * domaine.dot(&xv(f, 0), &domaine.veci(l, 0), &xp(e, 0)) * vit(fc);
                          }
                        if (!incompressible_)
                          {
                            for (int l = domaine.vedeb(e); l < domaine.vedeb(e + 1); l++)
                              {
                                fc = domaine.veji(l);
                                tab_resu(f) += fac * fs(f) * domaine.dot(&xv(f, 0), &domaine.veci(l, 0), &xp(e, 0)) * vit(fc);
                              }
                          }
                      }
                  }
              }
        }
    }

  return tab_resu;
}

/*! @brief on assemble la matrice.
 *
 */
inline void Op_Conv_EF_Stab_PolyMAC_Face::contribuer_a_avec(const DoubleTab& tab_inco, Matrice_Morse& matrice) const
{
  if (has_interface_blocs())
    {
      Operateur_base::contribuer_a_avec(tab_inco, matrice);
      return;
    }

  const Domaine_PolyMAC& domaine = le_dom_poly_.valeur();
  const Champ_Face_PolyMAC& ch = ref_cast(Champ_Face_PolyMAC, equation().inconnue());
  const IntTab& f_e = domaine.face_voisins(), &e_f = domaine.elem_faces(), &equiv = domaine.equiv();
  const DoubleTab& xp = domaine.xp(), &xv = domaine.xv(), &vfd = domaine.volumes_entrelaces_dir(), &vit = vitesse_->valeurs();
  const DoubleVect& fs = domaine.face_surfaces(), &vf = domaine.volumes_entrelaces(), &ve = domaine.volumes(), &pf = porosite_f, &pe = porosite_e;
  const DoubleTab& nf = domaine.face_normales();

  for (int e = 0; e < domaine.nb_elem_tot(); e++)
    {
      // Calcul de la divergence et contributions aux equations
      for (int i = 0; i < e_f.dimension(1); i++)
        {
          int f = e_f(e, i);
          if (f < 0) continue;

          int dir_f = (e == f_e(f, 0) ? 1 : -1);

          // Contributions des faces voisines
          if (f < domaine.nb_faces() && ch.fcl()(f, 0) < 2)
            {
              for (int j = 0; j < e_f.dimension(1); j++)
                {
                  int fb = e_f(e, j);
                  if (fb < 0) continue;

                  int dir_fb = (e == f_e(fb, 0) ? 1 : -1);

                  for (int k = 0; k < 2; k++)
                    {
                      if (ch.fcl()(fb, 0) < 2 || ch.fcl()(fb, 0) == 3)
                        {
                          int eb = f_e(fb, k);
                          double fac = dir_f * tab_inco(fb) * dir_fb * fs(fb) * pe(eb >= 0 ? eb : f_e(fb, 0)) / ve(e) * (1. + (vit(fb) * (k ? -1 : 1) >= 0 ? 1. : vit[fb] ? -1. : 0.) * alpha_) / 2;

                          // Convection pour les equivalences ou bords


                          int fc = equiv(fb, e != f_e(fb, 0), i);
                          if (fc  >= 0 || f_e(fb, 0) < 0 || f_e(fb, 1) < 0)
                            {
                              if (eb >= 0)
                                {
                                  int fd = (eb == e ? f : fc); //face source
                                  double mult = (fd < 0 || domaine.dot(&nf(f, 0), &nf(fd, 0)) > 0) ? 1 : -1;
                                  mult *= (fd >= 0) ? pf(fd) / pe(eb) : 1;

                                  matrice(f, fd) += fac * mult * dir_f * vfd(f, e != f_e(f, 0));
                                }

                              if (!incompressible_)
                                matrice(f, f) -= fac * dir_f * vfd(f, e != f_e(f, 0));
                            }
                          else
                            {
                              // Convection pour les faces internes sans equivalence
                              for (int l = domaine.vedeb(eb); l < domaine.vedeb(eb + 1); l++)
                                {
                                  fc = domaine.veji(l);
                                  if (ch.fcl()(fc, 0) < 2 && std::fabs(domaine.dot(&xv(f, 0), &domaine.veci(l, 0), &xp(e, 0))) > 1e-8 * vf(f) / fs(f))
                                    matrice(f, fc) += fac * fs(f) * domaine.dot(&xv(f, 0), &domaine.veci(l, 0), &xp(e, 0));
                                }
                              if (!incompressible_)
                                // Convection pour les faces internes sans equivalence
                                for (int l = domaine.vedeb(e); l < domaine.vedeb(e + 1); l++)
                                  {
                                    fc = domaine.veji(l);
                                    if (ch.fcl()(fc, 0) < 2 && std::fabs(domaine.dot(&xv(f, 0), &domaine.veci(l, 0), &xp(e, 0))) > 1e-8 * vf(f) / fs(f))
                                      matrice(f, fc) += fac * fs(f) * domaine.dot(&xv(f, 0), &domaine.veci(l, 0), &xp(e, 0));
                                  }
                            }
                        }
                    }
                }
            }
        }
    }
}

void Op_Conv_EF_Stab_PolyMAC_Face::set_incompressible(const int flag)
{
  if (flag == 0)
    {
      Cerr << "Compressible form of operator \"" << que_suis_je() << "\" :" << finl;
      Cerr << "Discretization of \u2207(inco \u2297 v) - v \u2207.(inco)" << finl;
    }
  incompressible_ = flag;
}
