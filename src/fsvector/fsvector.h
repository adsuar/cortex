/***************************************************************************
 *   Copyright (C) 2008 by Antonio Jesús Adsuar Gómez                      *
 *   adsuar@lsi.upc.edu                                                    *
 *                                                                         *
 *   This file is part of CORTEX                                           *
 *   CORTEX is free software; you can redistribute it and/or modify        *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   CORTEX is distributed in the hope that it will be useful,             *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/

/****************************************************************************
 ** CORTEX
 ** Tool designed by:
 **    Antonio Jesús Adsuar Gómez
 **    adsuar@lsi.upc.edu
 ** Sevilla 2008
 **
 **
 ** File: fsvector.h
 ** Description: Definition of the Vector Class
 ***************************************************************************/
#ifndef _VECTOR_H_
#define _VECTOR_H_

#include <iostream>
#include <vector>
#include <cassert>
#include <algorithm>
#include <math.h>
#include "../fsdefinitions/fstypes.h"

using namespace std;

class fsVector
{
   // Atributos
   private:
      vector<fsReal> v;
   
   // Metodos
   public:
      fsVector();
      fsVector(fsUInt);
      fsVector(const fsVector &);
      fsVector(const vector<fsReal> &);
      ~fsVector();
      
      vector<fsReal>* getSTL();

      fsUInt sizeOf();

      fsVoid add(fsReal);

      fsVoid randomInitialization(fsSInt,fsSInt);

      fsVoid initialization(fsReal);

      friend ostream& operator<<(ostream &,const fsVector &);
      float& operator[](fsUInt i);
      fsVector& operator=(const fsVector &);
      friend fsBool operator==(const fsVector&,const fsVector&);

      //vector<float>::iterator cerca(float);                         // Busca un valor en el vector

      //bool isOutOfBounds(vector<float>::iterator i);                // Metodo que indica si una posicion esta fuera de los limites

      //void add(float,int);                                          // Anade un valor en el vector
      //float getValue(unsigned int);                                 // Obtiene el i-esimo elemento del vector
      //void eraseValue(float);                                       // Elimina el elemento i del vector
      //void setValue(unsigned int,float);                            // Instancia el i-esimo elemento del vector
      
      //int normalize();                                              // Normaliza el vector
      //int normalize(unsigned int);                                  // Normaliza el vector hasta una determinada dimension
      
      //float module();                                               // Calcula el modulo de un vector
      //float module(unsigned int);                                   // Calcula el modulo de un vector hasta una determinada dimension
      
      //float& operator()(unsigned int i);                            // Sobrecarga del operador de indexacion
};

#endif // _VECTOR_H_
