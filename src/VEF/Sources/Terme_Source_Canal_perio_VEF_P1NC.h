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
// File:        Terme_Source_Canal_perio_VEF_P1NC.h
// Directory:   $TRUST_ROOT/src/VEF/Sources
// Version:     /main/17
//
//////////////////////////////////////////////////////////////////////////////



#ifndef Terme_Source_Canal_perio_VEF_P1NC_included
#define Terme_Source_Canal_perio_VEF_P1NC_included



//
// .DESCRIPTION class Terme_Source_Canal_perio_VEF_P1NC
//  Cette classe permet de conserver le debit dans une simulation
//  temporelle de Canal
//
// .SECTION voir aussi
//  Terme_Source_Canal_perio

#include <Terme_Source_Canal_perio.h>
#include <Ref_Zone_VEF.h>
#include <Ref_Zone_Cl_VEF.h>
class Probleme_base;
class Navier_Stokes_std;
class DoubleTab;

// La classe derive de Source_base et peut etre d'un terme source
class Terme_Source_Canal_perio_VEF_P1NC : public Terme_Source_Canal_perio
{
  Declare_instanciable(Terme_Source_Canal_perio_VEF_P1NC);

public :
  DoubleTab& ajouter(DoubleTab& ) const;
  DoubleTab& calculer(DoubleTab& ) const;

protected :

  REF(Zone_VEF) la_zone_VEF;
  REF(Zone_Cl_VEF) la_zone_Cl_VEF;
  void associer_zones(const Zone_dis& ,const Zone_Cl_dis& );
  virtual void calculer_debit(double&) const;
  // les attributs ont ete mis dans la classe mere

};

class Terme_Source_Canal_perio_QC_VEF_P1NC : public Terme_Source_Canal_perio_VEF_P1NC
{
  Declare_instanciable(Terme_Source_Canal_perio_QC_VEF_P1NC);
};


#endif
