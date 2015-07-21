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
// File:        Loi_2couches_base.h
// Directory:   $TRUST_ROOT/src/ThHyd/Turbulence
// Version:     /main/8
//
//////////////////////////////////////////////////////////////////////////////

#ifndef Loi_2couches_base_included
#define Loi_2couches_base_included

#include <Objet_U.h>


//////////////////////////////////////////////////////////////////////////////
//
// .DESCRIPTION
//    classe Loi_2couches_base
//    Cette classe de base represente les modeles 1 equation pour le modele a deux couches.
// .SECTION voir aussi
//////////////////////////////////////////////////////////////////////////////
class Loi_2couches_base : public Objet_U
{

  Declare_base(Loi_2couches_base);

public:
  /*
   * Cette methode est charge de calculer Leps, Lmu et sqrt(vv). LEs parametres signifient :
   * k, l'energie cintetique.
   * nu, la viscosite
   * dist, la distnace a la paroi.
   * y_etoile, la valeur de y_etoile.
   * Leps, la variable qui contiendra LEps,
   * Lmu la variable qui contiendra Lmu.
   * vvSqRt, la variable qui contiendra sqrt(vv).
   */
  virtual void  LepsLmu(double k, double nu, double dist, double y_etoile, double& Leps, double& Lmu, double& vvSqRt) = 0;

};


#endif


