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
 ** File: fsvector.cc
 ** Description: Source code for the fsVector class
 ***************************************************************************/

// fsvector.cc -- Instantiation of the methods of fsVector
#include "fsvector.h"

/****************************************************************************
 ** Name: fsVector
 ** Class: fsVector
 ** Created:           28/04/2008
 ** Last Modification: 28/04/2008
 ** Args:
 **    Input:
 **       - None
 **    Output:
 **       - None
 ** Return:
 **    Nothing
 ** 
 ** Description: Constructor of fsVector class
 **
 ***************************************************************************/

fsVector::fsVector()
{
   // Initially, the STL vector is already void
   // so we do not need to do anything more
}

/****************************************************************************
 ** Name: ~fsVector
 ** Class: fsVector
 ** Created:           28/04/2008
 ** Last Modification: 28/04/2008
 ** Args:
 **    Input:
 **       - None
 **    Output:
 **       - None
 ** Return:
 **    Nothing
 ** 
 ** Description: Destructor of fsVector class
 **
 ***************************************************************************/

fsVector::~fsVector()
{
   // We empty the vector
   v.clear();
}
/****************************************************************************
 ** Name: fsVector
 ** Class: fsVector
 ** Created:           28/04/2008
 ** Last Modification: 28/04/2008
 ** Args:
 **    Input:
 **       - fsUInt x: size of the vector
 **    Output:
 **       - None
 ** Return:
 **    Nothing
 ** 
 ** Description: Constructor of fsVector class, that creates an object where
 **              its size is previously determined
 **
 ***************************************************************************/

fsVector::fsVector(fsUInt x)
{
   // I make sure that x has to be greater than 0
   assert(x>=0);

   // Vector will be initialized with zeros
   v.insert(v.begin(), x, 0.0);
}

/****************************************************************************
 ** Name: fsVector
 ** Class: fsVector
 ** Created:           28/04/2008
 ** Last Modification: 28/04/2008
 ** Args:
 **    Input:
 **       - const fsVector &x: vector to be copied
 **    Output:
 **       - None
 ** Return:
 **    Nothing
 ** 
 ** Description: Copy constructor of fsVector class parting from another one
 **
 ***************************************************************************/

fsVector::fsVector(const fsVector &x)
{
   v=x.v;
}

/****************************************************************************
 ** Name: fsVector
 ** Class: fsVector
 ** Created:           28/04/2008
 ** Last Modification: 28/04/2008
 ** Args:
 **    Input:
 **       - const vector<fsReal> &x: STL vector to be copied
 **    Output:
 **       - None
 ** Return:
 **    Nothing
 ** 
 ** Description: Copy constructor of fsVector class parting from a STL vector
 **
 ***************************************************************************/

fsVector::fsVector(const vector<fsReal> &x)
{
   v=x;
}

/****************************************************************************
 ** Name: sizeOf
 ** Class: fsVector
 ** Created:           28/04/2008
 ** Last Modification: 28/04/2008
 ** Args:
 **    Input:
 **       - None
 **    Output:
 **       - None
 ** Return:
 **    - fsUInt
 ** 
 ** Description: Method that gets the size of a vector
 **
 ***************************************************************************/

fsUInt fsVector::sizeOf()
{
   return v.size();
}

/****************************************************************************
 ** Name: add
 ** Class: fsVector
 ** Created:           28/04/2008
 ** Last Modification: 28/04/2008
 ** Args:
 **    Input:
 **       - fsReal value: value to be added
 **    Output:
 **       - None
 ** Return:
 **    - None
 ** 
 ** Description: Method that adds a value at the end of the vector, at a
 **              position
 **
 ***************************************************************************/

fsVoid fsVector::add(fsReal value)
{
   v.push_back(value);
   
   return;
}
                

/****************************************************************************
 ** Name: randomInitialization
 ** Class: fsVector
 ** Created:           22/05/2008
 ** Last Modification: 22/05/2008
 ** Args:
 **    Input:
 **       - fsSInt min: minimum value of the initialization
 **       - fsSInt max: maximum value of the initialization
 **    Output:
 **       - None
 ** Return:
 **    - None
 ** 
 ** Description: Method that initializes the vector with random values
 **              between a minimum and a maximum value
 **
 ***************************************************************************/

fsVoid fsVector::randomInitialization(fsSInt min,fsSInt max)
{
   // Minimum value has to be smaller than the maximum one
   assert(min<max);

   // Establish a new seed for the generation of random integers
   srand(time(NULL));

   for(fsUInt i=0;i<v.size();i++)
   {
      v[i] = (fsReal)(min+(rand()%(max+1)));
   }

   return;
}

/****************************************************************************
 ** Name: initialization
 ** Class: fsVector
 ** Created:           22/05/2008
 ** Last Modification: 22/05/2008
 ** Args:
 **    Input:
 **       - fsReal value: value that'll be assigned to each position of the
 **                       vector
 **    Output:
 **       - None
 ** Return:
 **    - None
 ** 
 ** Description: Method that initializes the vector with a given value
 **
 ***************************************************************************/

fsVoid fsVector::initialization(fsReal value)
{
   for(fsUInt i=0;i<v.size();i++)
      v[i] = value;

   return;
}

/****************************************************************************
 ** Name: operator[]
 ** Class: fsVector
 ** Created:           28/04/2008
 ** Last Modification: 28/04/2008
 ** Args:
 **    Input:
 **       - fsUInt i: position to be indexed
 **    Output:
 **       - None
 ** Return:
 **    - fsReal&: memory address of the position indexed
 ** 
 ** Description: Method that overloads the indexation operator
 **
 ***************************************************************************/

float& fsVector::operator[](fsUInt i)
{
   // We access only to positions in the limits of the vector
   assert(i>=0 && i<v.size());

   return v[i];
}

/****************************************************************************
 ** Name: operator<<
 ** Class: fsVector
 ** Created:           28/04/2008
 ** Last Modification: 28/04/2008
 ** Args:
 **    Input:
 **       - ostream &os: output channel
 **       - const fsVector &c: vector to be printed
 **    Output:
 **       - None
 ** Return:
 **    ostream&
 ** 
 ** Description: Overload of the << operator
 **
 ***************************************************************************/

ostream& operator<<(ostream &os,const fsVector &c)
{
   if(c.v.size() == 0)
      os << "There's not any element; The vector is empty" << endl;
   else
   {
      for(fsUInt i=0;i<c.v.size();i++)
         os << c.v[i] << " ";
      os << endl;
   }

   return os;
}

/****************************************************************************
 ** Name: operator=
 ** Class: fsVector
 ** Created:           28/04/2008
 ** Last Modification: 03/05/2008
 ** Args:
 **    Input:
 **       - const fsVector &vect: vector to be assigned
 **    Output:
 **       - None
 ** Return:
 **    - fsVector&: memory address of the result of the assignation
 **
 ** Description: Method that assigns a vector to the propietary one
 **
 ***************************************************************************/

fsVector& fsVector::operator=(const fsVector &vect)
{
   v = vect.v;

   return *this;
}

/****************************************************************************
 ** Name: operator==
 ** Class: fsVector
 ** Created:           28/04/2008
 ** Last Modification: 28/04/2008
 ** Args:
 **    Input:
 **       - const fsVector &v1: first vector to be compared
 **       - const fsVector &v2: second vector to be compared
 **    Output:
 **       - None
 ** Return:
 **    - fsBool: result of the comparison
 **
 ** Description: Method that overloads the comparison operator between
 **              two vectors
 **
 ***************************************************************************/

fsBool operator==(const fsVector &v1,const fsVector &v2)
{
   return v1.v==v2.v;
}

/****************************************************************************
 ** Name: getSTL
 ** Class: fsVector
 ** Created:           28/04/2008
 ** Last Modification: 28/04/2008
 ** Args:
 **    Input:
 **       - None
 **    Output:
 **       - None
 ** Return:
 **    - vector<fsReal>*: memory address of stl vector
 **
 ** Description: Method that returns the pointer of the stl vector
 **
 ***************************************************************************/

vector<fsReal>* fsVector::getSTL()
{
   return &v;
}
