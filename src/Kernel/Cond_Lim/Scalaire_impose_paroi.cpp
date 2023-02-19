/****************************************************************************
* Copyright (c) 2023, CEA
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

#include <Scalaire_impose_paroi.h>

Implemente_instanciable(Scalaire_impose_paroi, "Scalaire_impose_paroi", Dirichlet);
// XD scalaire_impose_paroi dirichlet scalaire_impose_paroi 0 Imposed temperature condition at the wall called bord (edge).
// XD   attr ch front_field_base ch 0 Boundary field type.

Sortie& Scalaire_impose_paroi::printOn(Sortie& s) const { return s << que_suis_je() << finl; }

Entree& Scalaire_impose_paroi::readOn(Entree& s)
{
  if (app_domains.size() == 0) app_domains = { Motcle("Thermique"), Motcle("Concentration"), Motcle("Turbulence"), Motcle("indetermine") };
  if (supp_discs.size() == 0) supp_discs = { Nom("VEF"), Nom("EF"), Nom("EF_axi"), Nom("VEF_P1_P1"), Nom("VEFPreP1B"), Nom("PolyMAC"), Nom("PolyMAC_P0"), Nom("VDF") };

  return Dirichlet::readOn(s);
}
