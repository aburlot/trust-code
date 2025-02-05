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
// File      : Robin_VEF.h
// Directory : $TRUST_ROOT/src/Kernel/Cond_Lim
//
/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
//
// .DESCRIPTION : class Robin_VEF
//
// <Description of class Robin_VEF>
//
// g = \alpha( \nabla(u)n.n - p ) + u.n
// xi = \beta(\nabla(u)n\times n) + u \times n
/////////////////////////////////////////////////////////////////////////////

#ifndef Robin_VEF_included
#define Robin_VEF_included

#include <Cond_lim_base.h>

/*! @brief Class Robin_VEF for Robin boundary condititions
 *
 *    Robin boundary conditions functions are set as it follows :
 *      g = \alpha( \nabla(u)n.n - p ) + u.n
 *      xi = \beta(\nabla(u)n\times n) + u \times n
 *
 *      The parameters \alpha, \beta, and the functions g and xi comes from the data file e.g.
 *
 *      Conditions_limites {
 *                          Robin_boudary_name {
 *                                                  alpha val
 *                                                  beta val
 *                                                  champ_front_normal_et_tangentiel_robin champ_front_fonc_xyz 4 (in 3D, 2 in 2D) normal_field1D, tangential_field_X, tangential_field_Y, tangential_field_Z (tangential_fiel_1D for 2D)
 *                                              }
 *                          }
 *
 * @sa Cond_lim_base Robin_VEF
 */

class Robin_VEF : public Cond_lim_base
{

  Declare_instanciable( Robin_VEF ) ;

public :
  int compatible_avec_eqn(const Equation_base&) const override ;
  inline  double get_alpha_cl() const { return alpha_robin_cl_; } ;
  inline double get_beta_cl() const { return beta_robin_cl_; };


  // recuperation du champ complet
  double flux_robin_normal_et_trangentiel_imp(int i, int j) const;

  // champ normal
  double flux_normal_imp(int i) const;

  // champ tangentiel
  double flux_tangentiel_imp(int i, int j ) const;


  double flux_robin_imp_au_temps(double temps, int i) const ;
  double flux_robin_imp_au_temps(double temps, int i, int j ) const;

  const DoubleTab& flux_robin_normal_imp() const;
  const DoubleTab& flux_robin_tangentiel_imp() const;

  inline double increment_pression_bord(int face) const { return 0.; } ;

protected :


  double alpha_robin_cl_  ;
  double beta_robin_cl_   ;
  // Champ_front_normal_robin champ_normal_robin_cl_;
  // Champ_front_tangentiel_robin champ_tangent_robin_cl_;
  mutable DoubleTab flux_normal_impose_; // Stocke toutes les valeurs du flux sur toutes les faces de la frontiere (pas d'hypothese sur un champ uniforme). Utile pour le GPU.
  mutable DoubleTab flux_tangentiel_impose_; // Stocke toutes les valeurs du flux sur toutes les faces de la frontiere (pas d'hypothese sur un champ uniforme). Utile pour le GPU.


};

#endif /* Robin_VEF_included */
