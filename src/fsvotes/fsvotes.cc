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
 ** File: fsvotes.cc
 ** Description: Source code for the fsVotes class
 ***************************************************************************/

#include "fsvotes.h"

/****************************************************************************
 ** Name: fsVotes
 ** Class: fsVotes
 ** Created:           03/05/2008
 ** Last Modification: 03/05/2008
 ** Args:
 **    Input:
 **       - None
 **    Output:
 **       - None
 ** Return:
 **    Nothing
 ** 
 ** Description: Constructor of fsVotes class
 **
 ***************************************************************************/

fsVotes::fsVotes()
{
   // I create a new votation, empty at first, and that doesn't contain
   // any elements
   // There's no maximum value
   maximum = 0;
   // At first, map is closed
   opened = false;
}

/****************************************************************************
 ** Name: ~fsVotes
 ** Class: fsVotes
 ** Created:           03/05/2008
 ** Last Modification: 03/05/2008
 ** Args:
 **    Input:
 **       - None
 **    Output:
 **       - None
 ** Return:
 **    Nothing
 ** 
 ** Description: Destructor of fsVotes class
 **
 ***************************************************************************/

fsVotes::~fsVotes()
{
   ;
}

/****************************************************************************
 ** Name: search
 ** Class: fsVotes
 ** Created:           03/05/2008
 ** Last Modification: 03/05/2008
 ** Args:
 **    Input:
 **       - fsReal value: value to be search
 **    Output:
 **       - None
 ** Return:
 **    - map<fsReal,fsInt>::iterator: position of the value
 ** 
 ** Description: Method that returns an element if its inside the Universe
 **
 ***************************************************************************/

map<fsReal,fsInt>::iterator fsVotes::search(fsReal value)
{
   map<fsReal,fsInt>::iterator i;
   
   i=votation.find(value);
   
   return i;
}

/****************************************************************************
 ** Name: isOutOfBounds
 ** Class: fsVotes
 ** Created:           03/05/2008
 ** Last Modification: 03/05/2008
 ** Args:
 **    Input:
 **       - map<float,int>::iterator i: position to be searched
 **    Output:
 **       - None
 ** Return:
 **    - fsBool: boolean that indicates if the iterator is or not inside the
 **              bounds of the map
 ** 
 ** Description: Method that indicates if the iterator is or not inside the
 **              bounds of the map
 **
 ***************************************************************************/

fsBool fsVotes::isOutOfBounds(map<fsReal,fsInt>::iterator i)
{
   return (i==votation.end());
}

/****************************************************************************
 ** Name: add
 ** Class: fsVotes
 ** Created:           03/05/2008
 ** Last Modification: 03/05/2008
 ** Args:
 **    Input:
 **       - fsReal value: value to be added at the map
 **    Output:
 **       - None
 ** Return:
 **    - None
 ** 
 ** Description: Method that adds a value at the map. If it's already at the
 **              map, it gets one vote
 **
 ***************************************************************************/

fsVoid fsVotes::add(fsReal value)
{
   votation[value]++;
   
   if(votation[value] > maximum)
   {
      maximum = votation[value];
      winner = value;                   
   }                   

   // Map is closed
   opened = false;
}

/****************************************************************************
 ** Name: sizeOf
 ** Class: fsVotes
 ** Created:           04/05/2008
 ** Last Modification: 04/05/2008
 ** Args:
 **    Input:
 **       - None
 **    Output:
 **       - None
 ** Return:
 **    - fsUInt: Size of the map
 ** 
 ** Description: Method that returns the number of elements of the map
 **
 ***************************************************************************/

fsUInt fsVotes::sizeOf()
{
   return votation.size();
}

/****************************************************************************
 ** Name: getVotesOf
 ** Class: fsVotes
 ** Created:           04/05/2008
 ** Last Modification: 04/05/2008
 ** Args:
 **    Input:
 **       - fsReal i: Value from which we want to know its votes
 **    Output:
 **       - None
 ** Return:
 **    - fsUInt: Votes of a given value
 ** 
 ** Description: Method that returns the votes for a given value
 **
 ***************************************************************************/

fsUInt fsVotes::getVotesOf(fsReal i)
{
   return votation[i];
}

/****************************************************************************
 ** Name: getTheBest
 ** Class: fsVotes
 ** Created:           04/05/2008
 ** Last Modification: 04/05/2008
 ** Args:
 **    Input:
 **       - None
 **    Output:
 **       - None
 ** Return:
 **    - fsReal: Most voted item
 ** 
 ** Description: Method that returns the most voted item
 **
 ***************************************************************************/

fsReal fsVotes::getTheBest()
{
   // There has to be one vote at least        
   assert (maximum > 0);
           
   return winner;
}

/****************************************************************************
 ** Name: getClasses
 ** Class: fsVotes
 ** Created:           04/05/2008
 ** Last Modification: 04/05/2008
 ** Args:
 **    Input:
 **       - None
 **    Output:
 **       - None
 ** Return:
 **    - fsUInt: Number of different classes or items
 ** 
 ** Description: Method that returns the number of different classes or items
 **
 ***************************************************************************/

fsUInt fsVotes::getClasses()
{
   return votation.size();
}

/****************************************************************************
 ** Name: getSet
 ** Class: fsVotes
 ** Created:           04/05/2008
 ** Last Modification: 04/05/2008
 ** Args:
 **    Input:
 **       - None
 **    Output:
 **       - None
 ** Return:
 **    - fsVector: Set of classes or items
 ** 
 ** Description: Method that returns the set of classes or items as a vector
 **
 ***************************************************************************/

fsVector fsVotes::getSet()
{
   vector<fsReal> v;
   map<fsReal,fsInt>::iterator i = votation.begin();
   
   while(i!=votation.end())
   {
      fsReal a = i->first;
      v.push_back(a);
      i++;
   }
   
   fsVector vec(v);
   
   return vec;
}

/****************************************************************************
 ** Name: getVotes
 ** Class: fsVotes
 ** Created:           04/05/2008
 ** Last Modification: 04/05/2008
 ** Args:
 **    Input:
 **       - None
 **    Output:
 **       - None
 ** Return:
 **    - fsVector: Votes of the different classes or items
 ** 
 ** Description: Method that returns a vector with the votes of classes or
 **              items, and that's related with the one obtained with
 **              fsVotes::getClasses method
 **
 ***************************************************************************/

// Metodo que devuelve los votos transformados en vector
fsVector fsVotes::getVotes()
{
   vector<fsReal> v;
   map<fsReal,fsInt>::iterator i = votation.begin();
   
   while(i!=votation.end())
   {
      fsReal a = (fsReal)(i->second);
      v.push_back(a);
      i++;
   }
   
   fsVector vec(v);
   
   return vec;
}


/****************************************************************************
 ** Name: emptyVotes
 ** Class: fsVotes
 ** Created:           04/05/2008
 ** Last Modification: 04/05/2008
 ** Args:
 **    Input:
 **       - None
 **    Output:
 **       - None
 ** Return:
 **    - None
 ** 
 ** Description: Method that empties the votation, but not the items
 **
 ***************************************************************************/

fsVoid fsVotes::emptyVotes()
{
   map<fsReal,fsInt>::iterator x;
   
   for( x=votation.begin() ; x!=votation.end() ; ++x )
      x->second = 0;

   maximum=0;

   // At first, map is closed
   opened = false;
}

/****************************************************************************
 ** Name: operator<<
 ** Class: fsVotes
 ** Created:           04/05/2008
 ** Last Modification: 04/05/2008
 ** Args:
 **    Input:
 **       - ostream &os: output channel
 **       - const fsVotes &v: votation to be printed
 **    Output:
 **       - None
 ** Return:
 **    - ostream&: output
 ** 
 ** Description: Overload of the << operator
 **
 ***************************************************************************/

ostream& operator<<(ostream &os,const fsVotes &v)
{
   cout << "==========================================================" << endl;
   cout << "At the present votation, the results are: (class,votes)" << endl;
   cout << "==========================================================" << endl;

   if((v.votation).size() <= 0)
      os << "There's not any element; The votation is empty\n";
   else
   {
      map<fsReal,fsInt> mapa = v.votation;
      map<fsReal,fsInt>::iterator i = mapa.begin();
      
      for(;i!= mapa.end();++i)
      {
         os << (i->first) << ";" << (i->second) << endl;
      }
   }
   
   cout << "==========================================================" << endl << endl;

   return os;
}

/****************************************************************************
 ** Name: operator=
 ** Class: fsVotes
 ** Created:           04/05/2008
 ** Last Modification: 04/05/2008
 ** Args:
 **    Input:
 **       - const fsVotes &v: votation to be assigned
 **    Output:
 **       - None
 ** Return:
 **    - fsVotes&: Result of the assignation
 ** 
 ** Description: Overload of the = operator
 **
 ***************************************************************************/

fsVotes& fsVotes::operator=(const fsVotes &v)
{
   votation = v.votation;
   
   winner = v.winner;
   
   maximum = v.maximum;
   
   // It's important to take care of the fact that when we use the this
   // operator, it's like if we copy one object into another one. At this
   // point, it's like if we've not browsed the map, so we initialize the
   // object with opened property as false (closed)

   opened = false;
     
   return *this;
}

/****************************************************************************
 ** Name: open
 ** Class: fsVotes
 ** Created:           04/05/2008
 ** Last Modification: 04/05/2008
 ** Args:
 **    Input:
 **       - None
 **    Output:
 **       - None
 ** Return:
 **    - None
 ** 
 ** Description: Method that opens the map, and initializes the iterator
 **
 ***************************************************************************/

fsVoid fsVotes::open()
{
   opened = true;
   
   position = votation.begin();
}
/****************************************************************************
 ** Name: close
 ** Class: fsVotes
 ** Created:           04/05/2008
 ** Last Modification: 04/05/2008
 ** Args:
 **    Input:
 **       - None
 **    Output:
 **       - None
 ** Return:
 **    - None
 ** 
 ** Description: Method that closes the map
 **
 ***************************************************************************/

fsVoid fsVotes::close()
{
   opened = false;
}
/****************************************************************************
 ** Name: read
 ** Class: fsVotes
 ** Created:           04/05/2008
 ** Last Modification: 04/05/2008
 ** Args:
 **    Input:
 **       - None
 **    Output:
 **       - None
 ** Return:
 **    - fsReal: Item or class
 ** 
 ** Description: Method that returns the value of the item at the current
 **              position of the iterator
 **
 ***************************************************************************/

fsReal fsVotes::read()
{
   assert(position != votation.end());
   
   fsReal value = (position->first);
   
   ++position;
   
   return value;
}
/****************************************************************************
 ** Name: atTheEnd
 ** Class: fsVotes
 ** Created:           04/05/2008
 ** Last Modification: 04/05/2008
 ** Args:
 **    Input:
 **       - None
 **    Output:
 **       - None
 ** Return:
 **    - fsBool
 ** 
 ** Description: Method that indicates if iterator it's at the end, or map
 **              is closed
 **
 ***************************************************************************/

fsBool fsVotes::atTheEnd()
{
   return (position==votation.end()) || opened==false;
}
