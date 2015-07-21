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
// File:        Sous_zone_dis_base.cpp
// Directory:   $TRUST_ROOT/src/Kernel/Framework
// Version:     /main/5
//
//////////////////////////////////////////////////////////////////////////////

#include <Sous_zone_dis_base.h>

Implemente_base(Sous_zone_dis_base,"Sous_zone_dis_base",Objet_U);


Sortie& Sous_zone_dis_base::printOn(Sortie& os) const
{
  return os ;
}


Entree& Sous_zone_dis_base::readOn(Entree& is)
{
  return is ;
}

// Description:
//    Associe une Sous_Zone a l'objet.
void Sous_zone_dis_base::associer_sous_zone(const Sous_Zone& une_sous_zone)
{
  la_sous_zone=une_sous_zone;
}

// Description:
//    Associe une Zone_dis a l'objet.
void Sous_zone_dis_base::associer_zone_dis(const Zone_dis_base& une_zone_dis)
{
  la_zone_dis=une_zone_dis;
}

