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
 ** File: fsdistance.h
 ** Description: Definition of the Distance Class, that estimates the
 **              distance (difference) between two items (vectors)
 ***************************************************************************/

#ifndef _DISTANCE_H_
#define _DISTANCE_H_

#include "../fsvector/fsvector.h"
#include "../fsset/fsset.h"
#include "../fsdefinitions/fstypes.h"
#include "../fsdefinitions/fsconstants.h"

using namespace std;

class fsDistance 
{
   private:
      fsReal (fsDistance::*distance_vector)(vector<fsReal> *,vector<fsReal> *,vector<fsReal> *);
      fsReal (fsDistance::*distance_set)(vector<fsReal> *,vector<fsReal> *,fsSet *);
      fsReal euclideanDistance(vector<fsReal> *,vector<fsReal> *,vector<fsReal> *);
      fsReal squareEuclideanDistance(vector<fsReal> *,vector<fsReal> *,vector<fsReal> *);
      fsReal euclideanDistance(vector<fsReal> *,vector<fsReal> *,fsSet *);
      fsReal squareEuclideanDistance(vector<fsReal> *,vector<fsReal> *,fsSet *);
   public:
      fsDistance(fsDistanceAlgorithm d=SQUAREEUCLIDEANDISTANCE);
      ~fsDistance();
      fsReal get(fsVector*,fsVector*,fsVector*);
      fsReal get(vector<fsReal> *,vector<fsReal> *,vector<fsReal> *);
      fsReal get(fsVector*,fsVector*,fsSet*);
      fsReal get(vector<fsReal> *,vector<fsReal> *,fsSet *);
};

#endif // _DISTANCE_H_
