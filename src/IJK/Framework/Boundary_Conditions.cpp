/****************************************************************************
* Copyright (c) 2022, CEA
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

#include <Boundary_Conditions.h>
#include <Param.h>

Implemente_instanciable_sans_constructeur(Boundary_Conditions,"Boundary_Conditions",Objet_U);
Implemente_ref(Boundary_Conditions);

Boundary_Conditions::Boundary_Conditions()
{
  bctype_kmin_ = Paroi;
  bctype_kmax_ = Paroi;
}

// SYNTAXE:
//  {
//    BCTYPE_KMIN paroi|symetrie|perio
//    BCTYPE_KMAX paroi|symetrie|perio
//  }

Entree& Boundary_Conditions::readOn(Entree& is)
{
  Param param(que_suis_je());

  param.ajouter("bctype_kmin", &bctype_kmin_, Param::REQUIRED);
  param.dictionnaire("Paroi", Paroi);
  param.dictionnaire("Symetrie", Symetrie);
  param.dictionnaire("Perio", Perio);

  param.ajouter("bctype_kmax", &bctype_kmax_, Param::REQUIRED);
  param.dictionnaire("Paroi", Paroi);
  param.dictionnaire("Symetrie", Symetrie);
  param.dictionnaire("Perio", Perio);

  param.lire_avec_accolades(is);

  return is;
}

Sortie& Boundary_Conditions::printOn(Sortie& os) const
{
  return os;
}
