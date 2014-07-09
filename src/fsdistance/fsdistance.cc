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
 ** File: fsdistance.cc
 ** Description: Source code for the fsDistance class
 ***************************************************************************/

#include "fsdistance.h"

/****************************************************************************
 ** Name: fsDistance
 ** Class: fsDistance
 ** Created:           27/04/2008
 ** Last Modification: 27/04/2008
 ** Args:
 **    Input:
 **       - int d: It indicates the distance algorithm that'll be implemented
 **    Output:
 **       - None
 ** Return:
 **    Nothing
 ** 
 ** Description: Constructor of fsDistance class
 **
 ***************************************************************************/
fsDistance::fsDistance(fsDistanceAlgorithm d)
{
   switch(d)
   {
      case EUCLIDEANDISTANCE:
         distance_vector = &fsDistance::euclideanDistance;
         distance_set = &fsDistance::euclideanDistance;
         break;
      case SQUAREEUCLIDEANDISTANCE:
         distance_vector = &fsDistance::squareEuclideanDistance;
         distance_set = &fsDistance::squareEuclideanDistance;
         break;
      default:
         distance_vector = &fsDistance::squareEuclideanDistance;
         distance_set = &fsDistance::squareEuclideanDistance;
         break;
   }
}

/****************************************************************************
 ** Name: ~fsDistance
 ** Class: fsDistance
 ** Created:           27/04/2008
 ** Last Modification: 27/04/2008
 ** Args:
 **    Input:
 **       - None
 **    Output:
 **       - None
 ** Return:
 **    Nothing
 ** 
 ** Description: Destructor of fsDistance class
 **
 ***************************************************************************/
fsDistance::~fsDistance()
{
   distance_vector = NULL;
   distance_set = NULL;
}

/****************************************************************************
 ** Name: squareEuclideanDistance
 ** Class: fsDistance
 ** Created:           27/04/2008
 ** Last Modification: 03/05/2008
 ** Args:
 **    Input:
 **       - vector<fsReal> *v1: First vector
 **       - vector<fsReal> *v2: Second vector
 **       - vector<fsReal> *X: Vector of computable positions
 **    Output:
 **       - None
 ** Return:
 **    fsReal: Distance between vectors
 ** 
 ** Description: Method that calculates the square of the euclidean distance
 **              between the computable positions (given by vector X) of two
 **              given vectors.
 **              It's assumed that the three vectors have the same size
 **
 ***************************************************************************/
fsReal fsDistance::squareEuclideanDistance(
       vector<fsReal> *v1,
       vector<fsReal> *v2,
       vector<fsReal> *X
       )
{
#ifdef DEBUG_fsDistance
   cout << "Executing squareEuclideanDistance Algorithm, with:" << endl;
   cout << "\t- " << (*v1).size() << " features at origin" << endl;
   cout << "\t- " << (*v2).size() << " features at destination" << endl;
   cout << "\t- " << (*X).size() << " features mask max size" << endl;
#endif

   fsUInt sizeOfVectors = (*v1).size();
   
   fsReal sed=0.0;
   fsReal factor; 
   
   for(fsUInt i=0;i<sizeOfVectors;i++)
   {
      factor  = (*v1)[i]-(*v2)[i];
      sed    += (*X)[i]*factor*factor;
   }

   //cout << sed << endl;
   
   //exit(-1);
   return sed;
}

/****************************************************************************
 ** Name: euclideanDistance
 ** Class: fsDistance
 ** Created:           27/04/2008
 ** Last Modification: 03/05/2008
 ** Args:
 **    Input:
 **       - vector<fsReal> *v1: First vector
 **       - vector<fsReal> *v2: Second vector
 **       - vector<fsReal> *X: Vector of computable positions
 **    Output:
 **       - None
 ** Return:
 **    fsReal: Distance between vectors
 ** 
 ** Description: Method that calculates the euclidean distance between the
 **              computable positions (given by vector X) of two given
 **              vectors.
 **              It's assumed that the three vectors have the same size
 **
 ***************************************************************************/
fsReal fsDistance::euclideanDistance(
       vector<fsReal> *v1,
       vector<fsReal> *v2,
       vector<fsReal> *X
       )
{
   return sqrt(squareEuclideanDistance(v1,v2,X));
}

/****************************************************************************
 ** Name: squareEuclideanDistance
 ** Class: fsDistance
 ** Created:           05/06/2008
 ** Last Modification: 05/06/2008
 ** Args:
 **    Input:
 **       - vector<fsReal> *v1: First vector
 **       - vector<fsReal> *v2: Second vector
 **       - fsSet           *X: Set of taken features
 **    Output:
 **       - None
 ** Return:
 **    fsReal: Distance between vectors
 ** 
 ** Description: Method that calculates the square of the euclidean distance
 **              between the computable positions (given by vector X) of two
 **              given vectors.
 **              It's assumed that the three vectors have the same size
 **
 ***************************************************************************/
fsReal fsDistance::squareEuclideanDistance(
       vector<fsReal> *v1,
       vector<fsReal> *v2,
       fsSet           *X
       )
{
#ifdef DEBUG_fsDistance
   cout << "Executing squareEuclideanDistance Algorithm, with:" << endl;
   cout << "\t- " << (*v1).size() << " features at origin" << endl;
   cout << "\t- " << (*v2).size() << " features at destination" << endl;
   cout << "\t- " << (*X).sizeOfSet() << " features selected" << endl;
#endif

   fsUInt sizeOfFeatures = (*X).sizeOfSet();
   
   fsReal sed=0.0;
   fsReal factor; 
   
   for(fsUInt i=0;i<sizeOfFeatures;i++)
   {
      factor  = (*v1)[(fsUInt)(*X)[i]]-(*v2)[(fsUInt)(*X)[i]];
      sed    += factor*factor;
   }
   
   //exit(-1);
   return sed;
}

/****************************************************************************
 ** Name: euclideanDistance
 ** Class: fsDistance
 ** Created:           05/06/2008
 ** Last Modification: 05/06/2008
 ** Args:
 **    Input:
 **       - vector<fsReal> *v1: First vector
 **       - vector<fsReal> *v2: Second vector
 **       - fsSet           *X: Set of taken features
 **    Output:
 **       - None
 ** Return:
 **    fsReal: Distance between vectors
 ** 
 ** Description: Method that calculates the euclidean distance between the
 **              computable positions (given by vector X) of two given
 **              vectors.
 **              It's assumed that the three vectors have the same size
 **
 ***************************************************************************/
fsReal fsDistance::euclideanDistance(
       vector<fsReal> *v1,
       vector<fsReal> *v2,
       fsSet *X
       )
{
   return sqrt(squareEuclideanDistance(v1,v2,X));
}

/****************************************************************************
 ** Name: get
 ** Class: fsDistance
 ** Created:           27/04/2008
 ** Last Modification: 05/06/2008
 ** Args:
 **    Input:
 **       - fsVector *v1: First vector
 **       - fsVector *v2: Second vector
 **       - fsVector  *X:  Vector of computable positions
 **    Output:
 **       - None
 ** Return:
 **    fsReal: Distance between vectors
 ** 
 ** Description: Method that executes the distance function parametrized
 **              before at the Constructor Method
 **
 ***************************************************************************/
fsReal fsDistance::get(
       fsVector *v1,
       fsVector *v2,
       fsVector *X
       )
{
   return (this->*distance_vector)((*v1).getSTL(),(*v2).getSTL(),(*X).getSTL());
}

/****************************************************************************
 ** Name: get
 ** Class: fsDistance
 ** Created:           27/04/2008
 ** Last Modification: 05/06/2008
 ** Args:
 **    Input:
 **       - vector<fsReal> *v1: First vector
 **       - vector<fsReal> *v2: Second vector
 **       - vector<fsReal> *X: Vector of computable positions
 **    Output:
 **       - None
 ** Return:
 **    fsReal: Distance between vectors
 ** 
 ** Description: Method that executes the distance function parametrized
 **              before at the Constructor Method
 **
 ***************************************************************************/
fsReal fsDistance::get(
       vector<fsReal> *v1,
       vector<fsReal> *v2,
       vector<fsReal> *X
       )
{
   return (this->*distance_vector)(v1,v2,X);
}

/****************************************************************************
 ** Name: get
 ** Class: fsDistance
 ** Created:           05/06/2008
 ** Last Modification: 05/06/2008
 ** Args:
 **    Input:
 **       - fsVector *v1: First vector
 **       - fsVector *v2: Second vector
 **       - fsSet     *X: Set of features taken
 **    Output:
 **       - None
 ** Return:
 **    fsReal: Distance between vectors
 ** 
 ** Description: Method that executes the distance function parametrized
 **              before at the Constructor Method
 **
 ***************************************************************************/
fsReal fsDistance::get(
       fsVector *v1,
       fsVector *v2,
       fsSet *X
       )
{
   return (this->*distance_set)((*v1).getSTL(),(*v2).getSTL(),X);
}

/****************************************************************************
 ** Name: get
 ** Class: fsDistance
 ** Created:           05/06/2008
 ** Last Modification: 05/06/2008
 ** Args:
 **    Input:
 **       - vector<fsReal> *v1: First vector
 **       - vector<fsReal> *v2: Second vector
 **       - fsSet           *X: Set of features taken
 **    Output:
 **       - None
 ** Return:
 **    fsReal: Distance between vectors
 ** 
 ** Description: Method that executes the distance function parametrized
 **              before at the Constructor Method
 **
 ***************************************************************************/
fsReal fsDistance::get(
       vector<fsReal> *v1,
       vector<fsReal> *v2,
       fsSet *X
       )
{
   return (this->*distance_set)(v1,v2,X);
}

