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
// File:        vect.cpp
// Directory:   $TRUST_ROOT/src/Kernel/Utilitaires
// Version:     /main/23
//
//////////////////////////////////////////////////////////////////////////////

#include <vect.h>

Implemente_instanciable_sans_constructeur_ni_destructeur(vect,"vect",Objet_U);
Sortie& vect::printOn(Sortie& s) const
{
  int i=sz;
  s<<sz<<space;
  for(i=0; i<sz; i++)
    s<<data[i].valeur()<<space;
  assert(sz>=0);
  if(sz) assert(data!=0);
  return s<<finl;
}
Entree& vect::readOn(Entree& s)
{
  int i;
  s >> i;
  dimensionner(i);
  for(i=0; i<sz; i++)
    s>>data[i].valeur();
  assert(sz>=0);
  if(sz) assert(data!=0);
  return s;
}
Entree& vect::lit(Entree& s)
{
  int i;
  s >> i;
  dimensionner_force(i);
  for(i=0; i<sz; i++)
    s>>data[i].valeur();
  assert(sz>=0);
  if(sz) assert(data!=0);
  return s;
}

// Description: methhode surcharge par les classes filles pour fournir la bonne
// instance aux constructeurs au readOn....
Objet_U* vect::cree_une_instance() const
{
  Cerr<<"cree_une_instance() must be redefined in "<<que_suis_je()<<finl;
  exit();
  return cree_une_instance() ;
}
/* VECT constructors */
vect::vect() : sz(0), data(0) {}
void vect::build_vect(int i)
{
  sz=i;
  if(i==0) data=0;
  else
    {
      data=0;
      dimensionner_force(i);
      if(!data)
        {
          Cerr << "Allocation is not possible : " << finl;
          exit();
        }
    }
  assert(sz>=0);
  if(sz) assert(data!=0);
}
/* Constructor by copy for VECT */
vect::vect(const vect& avect)
  : Objet_U(avect),sz(avect.sz), data(new DerObjU[avect.sz])
{
  if(!data)
    {
      Cerr << "Allocation is not possible : " << finl;
      exit();
    }
  int i=sz;
  if(i==0)
    {
      delete[] data;
      data=0;
    }
  else while(i--)
      {
        (data[i])=(avect[i]);
      }
  assert(sz>=0);
  if(sz) assert(data!=0);
}
/* VECT destructor */
vect::~vect()
{
  reset();
  assert (data==0) ;
}
/* VECT operator= */
vect& vect::operator=(const vect& avect)
{
  if (this==&(ref_cast(vect,avect)))
    return *this;
  if(data)
    {
      reset();
    }
  sz=avect.sz;
  if(sz)
    {
      data = new DerObjU[sz];
      if(!data)
        {
          Cerr << "Allocation is not possible : " << finl;
          exit();
        }
      int i=sz;
      while(i--)
        {
          data[i]=avect[i];
        }
    }
  else data=0;
  assert(sz>=0);
  if(sz) assert(data!=0);
  return *this;
}
/* Size a VECT to i elements */
void vect::dimensionner(int i)
{
  if(sz==i) return;
  if(sz!=0)
    {
      Cerr << "WARNING : dimensionner method of a vect" << finl;
      Cerr << "old size : " << sz << finl;
      Cerr << "new size : " << i << finl;
    }
  assert(data==0);
  Objet_U* ptr_inst=cree_une_instance();
  const Objet_U& une_instance=*ptr_inst;
  sz=i;
  data = new DerObjU[i];
  if(!data)
    {
      Cerr << "Allocation is not possible. " << finl;
      exit();
    }
  for (int ii=0; ii<sz; ii++)
    data[ii]=une_instance;
  assert(sz>=0);
  if(sz) assert(data!=0);
  delete ptr_inst;
}
/* Empty then size a VECT to i elements */
void vect::dimensionner_force(int i)
{
  if(data)
    {
      reset();
    }
  Objet_U* ptr_inst=cree_une_instance();
  const Objet_U& une_instance=*ptr_inst;
  sz=i;
  data = new DerObjU[i];
  if(!data)
    {
      Cerr << "Allocation is not possible. " << finl;
      exit();
    }
  for (int ii=0; ii<sz; ii++)
    data[ii]=une_instance;
  assert(sz>=0);
  if(sz) assert(data!=0);
  delete ptr_inst;

}
/* Empty a VECT (delete data and reset size to 0) */
void vect::reset()
{
  if (data)
    {
      delete [] data;
      data=0;
      sz=0;
    }
  assert(sz==0);
}
/* Add a new element to the VECT */
vect& vect::add(const Objet_U& data_to_add)
{
  DerObjU* tempo;
  int new_size=sz+1;
  tempo = new DerObjU[new_size];
  for(int i=0; i<sz; i++)
    (tempo[i]).deplace(data[i]);
  if (data) delete [] data;
  data=tempo;
  sz=new_size;
  data[sz-1]= data_to_add;
  assert(sz>=0);
  if(sz) assert(data!=0);
  return *this;
}
/* Append a VECT to a VECT */
vect& vect::add(const vect& data_to_add)
{
  int sz_data_to_add = data_to_add.size();
  DerObjU* tempo;
  int new_size=sz+sz_data_to_add;
  tempo = new DerObjU[new_size];
  for(int i=0; i<sz; i++)
    (tempo[i]).deplace( data[i]);
  if (data) delete [] data;
  data=tempo;
  sz=new_size;
  for(int j=0; j<sz_data_to_add; j++)
    {
      data[sz-sz_data_to_add+j]= data_to_add[j];
    }
  assert(sz>=0);
  if(sz) assert(data!=0);
  return *this;
}

