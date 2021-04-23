/****************************************************************************
* Copyright (c) 2021, CEA
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
// File:        Perte_Charge_Singuliere_CoviMAC_Face.h
// Directory:   $TRUST_ROOT/src/CoviMAC/Sources
// Version:     /main/9
//
//////////////////////////////////////////////////////////////////////////////


#ifndef Perte_Charge_Singuliere_CoviMAC_Face_included
#define Perte_Charge_Singuliere_CoviMAC_Face_included

#include <Perte_Charge_CoviMAC_Face.h>
#include <Perte_Charge_Singuliere.h>

class Domaine;
class IntList;

//
// .DESCRIPTION class Perte_Charge_Singuliere_CoviMAC_Face
//

//
// .SECTION voir aussi Perte_Charge_CoviMAC_Face
//
//

class Perte_Charge_Singuliere_CoviMAC_Face : public Perte_Charge_CoviMAC_Face,
  public Perte_Charge_Singuliere
{

  Declare_instanciable(Perte_Charge_Singuliere_CoviMAC_Face);

public:
  int has_interface_blocs() const
  {
    return 1;
  };
  void check_multiphase_compatibility() const { };
  void dimensionner_blocs(matrices_t matrices, const tabs_t& semi_impl = {}) const { }; //rien
  void ajouter_blocs(matrices_t matrices, DoubleTab& secmem, const tabs_t& semi_impl = {}) const;
  void remplir_num_faces(Entree& );
  void mettre_a_jour(double temps);

protected:
  IntVect sgn;
};
#endif
