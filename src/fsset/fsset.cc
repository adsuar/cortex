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
 ** File: fsset.cc
 ** Description: Source code for the fsSet class
 ***************************************************************************/

#include "fsset.h"

/****************************************************************************
 ** Name: fsSet
 ** Class: fsSet
 ** Created:           03/06/2008
 ** Last Modification: 03/06/2008
 ** Args:
 **    Input:
 **       - None
 **    Output:
 **       - None
 ** Return:
 **    Nothing
 ** 
 ** Description: Constructor of fsSet class
 **
 ***************************************************************************/

fsSet::fsSet()
{
   size = -1;
}

/****************************************************************************
 ** Name: ~fsSet
 ** Class: fsSet
 ** Created:           03/06/2008
 ** Last Modification: 03/06/2008
 ** Args:
 **    Input:
 **       - None
 **    Output:
 **       - None
 ** Return:
 **    Nothing
 ** 
 ** Description: Destructor of fsSet class
 **
 ***************************************************************************/

fsSet::~fsSet()
{
   ;
}

/****************************************************************************
 ** Name: fsSet
 ** Class: fsSet
 ** Created:           03/06/2008
 ** Last Modification: 03/06/2008
 ** Args:
 **    Input:
 **       - fsVector v: Vector to be copied to set
 **       - fsReal val: Value that identifies the positions to be copied
 **       - fsUInt k: Number of features
 **    Output:
 **       - None
 ** Return:
 **    Nothing
 ** 
 ** Description: Copy constructor of fsSet class
 **
 ***************************************************************************/
fsSet::fsSet(fsVector v,fsReal val, fsUInt k)
{
   if(k<=0) k=v.sizeOf();
   // The set it's initialized with those positions of the vector
   // where its value it's like val.
   // The purpose of k is to optimize the travelling across the vector
   // indicating the maximum number of val we'll count
   for(fsUInt i=0,nk=0;i<v.sizeOf() && nk<k;i++)
   {
      if(v[i] == val)
      {
         nk++;
         s_in.push_back(i);
      }
   }
   size = v.sizeOf();
}

/****************************************************************************
 ** Name: getBinaryVector
 ** Class: fsSet
 ** Created:           04/06/2008
 ** Last Modification: 04/06/2008
 ** Args:
 **    Input:
 **       - None
 **    Output:
 **       - None
 ** Return:
 **    fsVector: Binary vector
 ** 
 ** Description: Get the binary vector related with the set
 **
 ***************************************************************************/

fsVector fsSet::getBinaryVector()
{
   fsVector v(size);

   v.initialization(0.0);

   for(fsUInt i=0;i<s_in.size();i++)
      v[(fsUInt)s_in[i]] = 1;

   return v;
}

/****************************************************************************
 ** Name: operator-
 ** Class: fsSet
 ** Created:           04/06/2008
 ** Last Modification: 04/06/2008
 ** Args:
 **    Input:
 **       - fsSet &x: Origin set
 **       - fsSet &y: Set to be substracted
 **    Output:
 **       - None
 ** Return:
 **    fsSet: Result of the substraction
 ** 
 ** Description: Substraction between two sets
 **
 ***************************************************************************/
fsSet operator-(const fsSet &x, const fsSet &y)
{
   fsSet z(x);   // Set that'll save z = x-y

   for(fsUInt i=0;i<y.s_in.size();i++)
   {
      fsBool finish = false;
      // Delete every element from y at z
      for(vector<fsReal>::iterator j=z.s_in.begin();j!=z.s_in.end() && !finish;++j)
      {
         if(*j == y.s_in[i])
         {
            finish = true;
            z.s_in.erase(j);
         }
      }
   }

   return z;
}


/****************************************************************************
 ** Name: operator-
 ** Class: fsSet
 ** Created:           04/06/2008
 ** Last Modification: 04/06/2008
 ** Args:
 **    Input:
 **       - fsSet &x: Origin set
 **       - fsReal: Value to be substracted
 **    Output:
 **       - None
 ** Return:
 **    fsSet: Result of the substraction
 ** 
 ** Description: Substraction between two sets
 **
 ***************************************************************************/
fsSet operator-(const fsSet &x, const fsReal &y)
{
   fsSet z(x);   // Set that'll save z = x-y

   fsBool finish = false;

   // Delete the y element from z
   for(vector<fsReal>::iterator j=z.s_in.begin();j!=z.s_in.end() && !finish;++j)
   {
      if(*j == y)
      {
         finish = true;
         z.s_in.erase(j);
      }
   }

   return z;
}


/****************************************************************************
 ** Name: operator+
 ** Class: fsSet
 ** Created:           04/06/2008
 ** Last Modification: 04/06/2008
 ** Args:
 **    Input:
 **       - fsSet &x: Origin set
 **       - fsSet &y: Set to be added
 **    Output:
 **       - None
 ** Return:
 **    fsSet: Result of the union
 ** 
 ** Description: Union between two sets
 **
 ***************************************************************************/
fsSet operator+(const fsSet &x, const fsSet &y)
{
   fsSet z(x);   // Set that'll save z = x+y

   for(fsUInt i=0;i<y.s_in.size();i++)
   {
      // If a value of y cannot be finded at z, it has to be inserted
      if(find(z.s_in.begin(), z.s_in.end(), y.s_in[i]) == z.s_in.end())
         z.s_in.push_back(y.s_in[i]);
   }

   return z;
}


/****************************************************************************
 ** Name: operator+
 ** Class: fsSet
 ** Created:           04/06/2008
 ** Last Modification: 04/06/2008
 ** Args:
 **    Input:
 **       - fsSet &x: Origin set
 **       - fsReal: Value to be added
 **    Output:
 **       - None
 ** Return:
 **    fsSet: Result of the union
 ** 
 ** Description: Union between two sets
 **
 ***************************************************************************/
fsSet operator+(const fsSet &x, const fsReal &y)
{
   fsSet z(x);   // Set that'll save z = x-y

   // If a value of y cannot be finded at z, it has to be inserted
   if(find(z.s_in.begin(), z.s_in.end(), y) == z.s_in.end())
      z.s_in.push_back(y);

   return z;
}

/****************************************************************************
 ** Name: operator==
 ** Class: fsSet
 ** Created:           04/06/2008
 ** Last Modification: 04/06/2008
 ** Args:
 **    Input:
 **       - fsSet &s1: First set to be compared
 **       - fsSet &s2: Second set to be compared
 **    Output:
 **       - None
 ** Return:
 **    fsBool: Result of the comparison
 ** 
 ** Description: Comparison between two sets
 **
 ***************************************************************************/
fsBool operator==(const fsSet &s1,const fsSet &s2)
{
   return (s1.s_in)==(s2.s_in);
}

/****************************************************************************
 ** Name: operator!=
 ** Class: fsSet
 ** Created:           04/06/2008
 ** Last Modification: 04/06/2008
 ** Args:
 **    Input:
 **       - fsSet &s1: First set to be compared
 **       - fsSet &s2: Second set to be compared
 **    Output:
 **       - None
 ** Return:
 **    fsBool: Result of the comparison
 ** 
 ** Description: Comparison between two sets
 **
 ***************************************************************************/
fsBool operator!=(const fsSet &s1,const fsSet &s2)
{
   return s1.s_in != s2.s_in;
}

/****************************************************************************
 ** Name: operator=
 ** Class: fsSet
 ** Created:           04/06/2008
 ** Last Modification: 04/06/2008
 ** Args:
 **    Input:
 **       - fsSet &s1: Set to receive the set
 **       - fsSet &s2: Set to be assigned
 **    Output:
 **       - None
 ** Return:
 **    fsSet&: The expression will be evaluated to the value of the
 **            assignation
 ** 
 ** Description: Assignation of one set
 **
 ***************************************************************************/
fsSet& fsSet::operator=(const fsSet &s)
{
   s_in = s.s_in;

   size = s.size;
   
   return *this;
}

/****************************************************************************
 ** Name: operator[]
 ** Class: fsSet
 ** Created:           04/06/2008
 ** Last Modification: 04/06/2008
 ** Args:
 **    Input:
 **       - fsUInt i: index
 **    Output:
 **       - None
 ** Return:
 **    fsReal&: Value saved at the position indexed with i
 ** 
 ** Description: Overload of the indexation operator
 **
 ***************************************************************************/
fsReal& fsSet::operator[](fsUInt i)
{
   // Access only to position in the limits of the set
   assert(i>=0 && i<s_in.size());

   return s_in[i];
}

/****************************************************************************
 ** Name: sizeOfSet
 ** Class: fsSet
 ** Created:           05/06/2008
 ** Last Modification: 05/06/2008
 ** Args:
 **    Input:
 **       - None
 **    Output:
 **       - None
 ** Return:
 **    fsUInt: Number of elements at the set
 ** 
 ** Description: Method that returns the number of elements at the set
 **
 ***************************************************************************/
fsUInt fsSet::sizeOfSet()
{
   return s_in.size();
}

/****************************************************************************
 ** Name: sizeOfDominion
 ** Class: fsSet
 ** Created:           05/06/2008
 ** Last Modification: 05/06/2008
 ** Args:
 **    Input:
 **       - None
 **    Output:
 **       - None
 ** Return:
 **    fsUInt: Number of different elements that can be at the set
 ** 
 ** Description: Method that returns the number of different values that
 **              can be saved at the set
 **
 ***************************************************************************/
fsUInt fsSet::sizeOfDominion()
{
   return size;
}

/****************************************************************************
 ** Name: operator<<
 ** Class: fsSet
 ** Created:           05/06/2008
 ** Last Modification: 05/06/2008
 ** Args:
 **    Input:
 **       - ostream &os: output channel
 **       - const fsSet &c: set to be printed
 **    Output:
 **       - None
 ** Return:
 **    ostream&
 ** 
 ** Description: Overload of the << operator
 **
 ***************************************************************************/

ostream& operator<<(ostream &os,const fsSet &s)
{
   if(s.s_in.size() == 0)
      os << "There's not any element; The set is empty" << endl;
   else
   {
      for(fsUInt i=0;i<s.s_in.size();i++)
         os << s.s_in[i] << " ";
      os << endl;
   }
 
   return os;
}

/****************************************************************************
 ** Name: getRandomSubselection
 ** Class: fsSet
 ** Created:           08/06/2008
 ** Last Modification: 08/06/2008
 ** Args:
 **    Input:
 **       - fsSet *origin: Origin Set
 **       - fsDouble k: Max amount of features to be taken
 **    Output:
 **       - None
 ** Return:
 **    fsSet&: Set constructed
 ** 
 ** Description: Friend function that returns a subset of a given set
 **              with the size indicated
 **
 ***************************************************************************/

fsSet getRandomSubselection(fsSet *origin,fsDouble k)
{
   // Number of features to be taken must be greater than 0
   assert(k>0 && (*origin).sizeOfSet()>=k);

   fsSet result(*origin);

   if(result.sizeOfSet()>k)
   {
      fsUInt i=0;
      
      srand(time(NULL));

      while(i<k)
      {
         fsUInt position = rand()%(result.sizeOfSet());
         
         result = result - result[position];
      }
   }
   
   return result;
}
