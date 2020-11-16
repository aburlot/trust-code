/****************************************************************************
* Copyright (c) 2020, CEA
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
// File:        Interpolation_IBM_mean_gradient.h
// Directory:   $TRUST_ROOT/src/Kernel/Geometrie/Interpolation_IBM
// Version:     1
//
//////////////////////////////////////////////////////////////////////////////

#ifndef Interpolation_IBM_mean_gradient_included
#define Interpolation_IBM_mean_gradient_included

#include <Interpolation_IBM_base.h>
#include <Champ_Don.h>
#include <IntList.h>
#include <Zone.h>

/////////////////////////////////////////////////////////////////////////////
//
// .DESCRIPTION : class Interpolation_IBM_mean_gradient
//
// <Description of class Interpolation_IBM_mean_gradient>
//
/////////////////////////////////////////////////////////////////////////////

class Interpolation_IBM_mean_gradient : public Interpolation_IBM_base
{

  Declare_instanciable( Interpolation_IBM_mean_gradient ) ;

public :
  void discretise(const Discretisation_base&, Zone_dis_base& la_zone_EF);
  inline IntList& getSommetsVoisinsOf(int i)
  {
    return sommets_voisins_[i];
  };
protected :
  void computeSommetsVoisins(Zone_dis_base&);
  Champ_Don solid_points_lu_;
  Champ_Don solid_points_;
  Champ_Don is_dirichlet_lu_;
  Champ_Don solid_elems_lu_;
  Champ_Don corresp_elems_lu_;
  Champ_Don is_dirichlet_;
  IntList* sommets_voisins_;
  Champ_Don corresp_elems_;
  Champ_Don solid_elems_;
  friend class Source_PDF_EF;
  friend class Interpolation_IBM_hybrid;
};

#endif /* Interpolation_IBM_mean_gradient_included */
