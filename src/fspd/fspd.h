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
 ** File: fspd.h
 ** Description: Definition of the Problem Data Class, that controls all
 **              the variables involved in any problem
 ***************************************************************************/

#ifndef _PD_H_
#define _PD_H_

using namespace std;

#include "../fsvector/fsvector.h"
#include "../fsmatrix/fsmatrix.h"
#include "../fsevaluation/fsevaluation.h"
#include "../fsdefinitions/fstypes.h"
#include "../fsdefinitions/fsconstants.h"
#include "../fsdefinitions/fsmacros.h"

class fsPD
{
   public:
      fsVector *X;
      
      fsVector XD;
      
      fsVector vJ;
      
      fsUInt k;
      
      fsReal add_est,sub_est;
      
      fsInt add, sub;
      
      fsInt last;
      
      fsUInt features;
      
      fsReal best_of_all_J;                     // Valor de la mejor J

      fsInt position_best_J;                     // Momento en que se ha obtenido la mejor J
      
      fsPD();
      fsPD(fsUInt);
      ~fsPD();

      friend ostream& operator<<(ostream &,const fsPD &);
      fsPD& operator=(const fsPD &);       // Sobrecarga del operador asignacion
};

#endif // _PD_H
