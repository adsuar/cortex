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
 ** File: fsmatrix.h
 ** Description: Definition of the Matrix Class
 ***************************************************************************/
#ifndef _MATRIX_H_
#define _MATRIX_H_

#include <iostream>
#include <vector>
#include <cassert>
#include <algorithm>
#include <math.h>
#include "../fsvector/fsvector.h"
#include "../fsdefinitions/fstypes.h"

using namespace std;

class fsMatrix
{
   private:
      vector< vector<fsReal> > matrix;

   public:
      fsMatrix();
      fsMatrix(fsUInt,fsUInt);
      ~fsMatrix();
      
      fsVoid resize(fsUInt,fsUInt);
      fsVoid addRow(fsUInt,fsVector*);
      fsVoid addCol(fsUInt,fsVector*);
      
      fsVoid normalizeByCols(fsUInt,fsUInt);
      fsVoid normalizeByRows(fsUInt,fsUInt);
      
      fsVoid normalizeCol(fsUInt);
      fsVoid normalizeRow(fsUInt);
      
      fsReal moduleCol(fsUInt);
      fsReal moduleRow(fsUInt);
      
      fsReal averageCol(fsUInt);
      fsReal averageRow(fsUInt);
      
      fsReal varianceCol(fsUInt);
      fsReal varianceRow(fsUInt);
      
      fsVector extractRow(fsUInt);
      fsVector extractCol(fsUInt);
      
      vector<fsReal> *extractRowSTL(fsUInt);
      
      friend ostream& operator<<(ostream &,const fsMatrix &);
      
      fsUInt numRows();
      fsUInt numCols();

      fsMatrix& operator=(const fsMatrix &);
};

#endif // _MATRIX_H_

