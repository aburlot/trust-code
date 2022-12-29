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

#ifndef Zone_dis_included
#define Zone_dis_included

#include <Zone_dis_base.h>
#include <TRUST_Deriv.h>

class  Zone_Cl_dis_base;

/*! @brief classe Zone_dis Classe generique de la hierarchie des zones discretisees.
 *
 * Un objet
 *      de type Zone_dis peut referencer n'importe quel objet derivant de
 *      Zone_dis_base.
 *      La plupart des methodes appellent les methodes de l'objet Zone
 *      sous-jacent via la methode valeur() declaree grace a la macro
 *
 * @sa Zone_dis_base
 */
class Zone_dis : public DERIV(Zone_dis_base)
{

  Declare_instanciable(Zone_dis);

public :

  inline void associer_zone(const Zone&);
  inline void associer_domaine_dis(const Domaine_dis&);
  inline const Zone& zone() const;
  inline Zone& zone();
  inline void discretiser();
  inline void creer_elements_fictifs(const Zone_Cl_dis_base&);
  inline const Frontiere_dis_base& frontiere_dis(int ) const;
};


/*! @brief Appel a l'objet sous-jacent.
 *
 * Associe une zone (non discretisee) a l'objet.
 *
 * @param (Zone& une_zone)
 */
inline void Zone_dis::associer_zone(const Zone& une_zone)
{
  valeur().associer_zone(une_zone);
}

inline void Zone_dis::associer_domaine_dis(const Domaine_dis& un_domaine_dis)
{
  valeur().associer_domaine_dis(un_domaine_dis);
}

/*! @brief Appel a l'objet sous-jacent.
 *
 * Renvoie la zone associee.
 *     (version const)
 *
 * @return (Zone&) la zone associee
 */
inline const Zone& Zone_dis::zone() const
{
  return valeur().zone();
}

/*! @brief Appel a l'objet sous-jacent.
 *
 * Renvoie la zone associee.
 *
 * @return (Zone&) la zone associee
 */
inline Zone& Zone_dis::zone()
{
  return valeur().zone();
}

/*! @brief Appel a l'objet sous-jacent.
 *
 * Se discretise cf Zone_dis_base
 *
 */
inline void Zone_dis::discretiser()
{
  valeur().discretiser();
}

inline void Zone_dis::creer_elements_fictifs(const Zone_Cl_dis_base& zcl)
{
  valeur().creer_elements_fictifs( zcl);
}
/*! @brief Appel a l'objet sous-jacent.
 *
 * Renvoie la ieme frontiere
 *
 * @param (int i) l'index de la frontiere a renvoyer
 * @return (Frontiere_dis_base&) la i-eme frontiere de la zone discretisee
 */
inline const Frontiere_dis_base& Zone_dis::frontiere_dis(int i) const
{
  return valeur().frontiere_dis(i);
}


#endif
