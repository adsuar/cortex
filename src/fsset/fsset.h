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
 ** File: fsset.h
 ** Description: Definition of the Set Class, that manages sets
 **
 ***************************************************************************/
#ifndef _SET_H_
#define _SET_H_

#include <iostream>
#include <cassert>
#include <vector>
#include "../fsvector/fsvector.h"
#include "../fsdefinitions/fstypes.h"

using namespace std;

class fsSet
{
   // Attribs
   private:
      vector<fsReal> s_in;
      fsUInt size;

   // Methods
   public:
      fsSet();
      ~fsSet();
      fsSet(fsVector,fsReal,fsUInt k=0);
      fsVector getBinaryVector();

      fsUInt sizeOfDominion();
      fsUInt sizeOfSet();

      friend fsSet operator-(const fsSet &,const fsSet &);
      friend fsSet operator-(const fsSet &,const fsReal &);
      friend fsSet operator+(const fsSet &,const fsSet &);
      friend fsSet operator+(const fsSet &,const fsReal &);
      friend fsBool operator==(const fsSet &,const fsSet &);
      friend fsBool operator!=(const fsSet &,const fsSet &);
      fsSet& operator=(const fsSet &);
      fsReal& operator[](fsUInt);

      friend ostream& operator<<(ostream &,const fsSet &);

      friend fsSet getRandomSubselection(fsSet *,fsDouble);
};

#endif // _SET_H_

