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
 ** File: fsvotes.h
 ** Description: Definition of the Votes Class, that gets the votes for each
 **              concept voted
 **
 ***************************************************************************/
#ifndef _VOTES_H_
#define _VOTES_H_

#include <iostream>
#include <map>
#include <cassert>
#include <vector>
#include "../fsdefinitions/fstypes.h"
#include "../fsmatrix/fsmatrix.h"

class fsVotes
{
   // Atributos
   private:
      // Empleare el tipo MAP para almacenar los elementos de la votacion
      // La declaracion contendra dos elementos, un TVector como indice (es la
      // clase) y un entero, que guarda el numero de votos.
      map<fsReal, fsInt> votation;
      fsReal winner;
      fsInt maximum;

      map<fsReal,fsInt>::iterator position;                                  // Iterador interno que me ayuda
                                                                          // a realizar un recorrido
                                                                          // por el mapa sin necesidad
                                                                          // de realizar una busqueda de una
                                                                          // posicion en un momento determinado
                                                                          // Empleare metodos del tipo SECUENCIA
      fsBool opened;                                                        // Booleano que indica si una secuencia
                                                                          // esta abierta o cerrada
   
   // Métodos
   public:
      fsVotes();
      ~fsVotes();
      
      map<fsReal,fsInt>::iterator search(fsReal);       // Busca un valor en el
                                                       // conjunto
      fsBool isOutOfBounds(map<fsReal,fsInt>::iterator);  // Metodo que indica si
                                                       // una posicion esta
                                                       // fuera de los limites
      fsVoid add(fsReal);                        // Vota a un vector STL
                                                       // determinado
      fsUInt sizeOf();                                    // Tamaño del universo de
                                                       // votación
      fsUInt getVotesOf(fsReal);                  // Obtiene los votos de i
      fsReal getTheBest();               // Obtiene el ganador
      
      fsUInt getClasses();                                // Devuelve el numero de
                                                       // clases
      fsVector getSet();                                                   // Devuelve un vector con las clases
      fsVector getVotes();                                                   // Devuelve un vector con las clases
      
      fsVoid emptyVotes();         // Método que elimina las votaciones, sin eliminar
                                 // las clases
      
      friend ostream& operator<<(ostream &,const fsVotes &);
      fsVotes& operator=(const fsVotes &);       // Sobrecarga del operador asignacion

      fsVoid open();                                                        // Metodo que abre el mapa
      fsVoid close();                                                       // Metodo que cierra el mapa
      fsReal read();                                              // Metodo que lee un elemento y pasa al siguiente
      fsBool atTheEnd();                                                    // Metodo que indica si estoy al final
};

#endif // _VOTES_H_

