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

/***************************************************************************
 ** CORTEX
 ** Tool designed by:
 **    Antonio Jesús Adsuar Gómez
 **    adsuar@lsi.upc.edu
 ** Sevilla 2008
 **
 **
 ** File: fsdata.h
 ** Description: Definition of the Data class, that implements a new
 **              environment to manage data values and classes
 ***************************************************************************/

#ifndef _DATA_H_
#define _DATA_H_

using namespace std;

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <cassert>
#include <vector>

#include "../fsvector/fsvector.h"
#include "../fsmatrix/fsmatrix.h"
#include "../fsdefinitions/fstypes.h"
#include "../fsdefinitions/fsconstants.h"
#include "../fsdefinitions/fsmacros.h"

class fsData
{
   // Propiedades
   private:
      fsVector data_features;
      fsVector data_classes;
      fsMatrix features;
      fsMatrix classes;
      fsDataTypes features_domain;
      fsDataTypes classes_domain;

      fsBool dataLoaded;
      
   // Métodos
   public:
      fsData();                               // Constructor
      fsData(fsChar *);                         // Constructor a partir de un fichero
      ~fsData();                              // Destructor
      
      fsVoid normalizeFeatures();                  // Normalizar el conjunto de datos
                                             // por columnas
      
      fsVoid showFeatures();                       // Metodo que muestra los datos
      fsVoid showClasses();                 // Metodo que muestra las categorias
      fsVoid showTypesData();                  // Metodo que muestra el tipo de
                                             // valor de cada columna
      
      fsDataTypes getTypeFeatures();                     // Muestra el tipo de los datos
      
      fsDataTypes getTypeClasses();               // Muestra el tipo de las categorias
      
      fsMatrix* getClasses();              // Metodo que devuelve un puntero
                                             // a la matriz de las categorias
      fsMatrix* getFeatures();                    // Metodo que devuelve un puntero
                                             // a la matriz de datos
      
      fsUInt getInstances();                    // Metodo que indica el numero de
                                             // instancias que tiene el problema
      
      fsData& operator=(const fsData &);       // Sobrecarga del operador asignacion
      
      fsBool successfullyLoaded();
      
      fsUInt getNumberOfFeatures();
      
      fsUInt getNumberOfClasses();
};

#endif // _DATA_H_
