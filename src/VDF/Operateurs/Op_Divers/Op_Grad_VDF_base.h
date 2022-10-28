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

#ifndef Op_Grad_VDF_base_included
#define Op_Grad_VDF_base_included

#include <Operateur_Grad.h>
#include <Iterateur_VDF.h>

/*! @brief class Op_Grad_VDF_base Classe de base des operateurs de gradient VDF
 *
 */
class Op_Grad_VDF_base : public Operateur_Grad_base
{
  Declare_base(Op_Grad_VDF_base);
public:
  inline Op_Grad_VDF_base(const Iterateur_VDF_base& iter_base) : iter(iter_base) { }

  void check_multiphase_compatibility() const override { }
  void completer() override;
  int impr(Sortie& os) const override;

  inline DoubleTab& calculer(const DoubleTab& inco, DoubleTab& resu ) const override { return iter->calculer(inco,resu); } // calcule la contribution du gradient
  inline int has_interface_blocs() const override
  {
    return 1;
  };
  inline void dimensionner_blocs(matrices_t matrices, const tabs_t& semi_impl) const override {}
  inline void ajouter_blocs(matrices_t matrices, DoubleTab& secmem, const tabs_t& semi_impl) const override
  {
    const DoubleTab& inco = equation().inconnue().valeur().valeurs();
    iter->ajouter(inco,secmem);
  }


protected:
  Iterateur_VDF iter;
};

#endif /* Op_Grad_VDF_base_included */
