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
 ** File: fspd.cc
 ** Description: Source code for the fspd class
 ***************************************************************************/

#include "fspd.h"

/****************************************************************************
 ** Name: fsPD
 ** Class: fsPD
 ** Created:           10/05/2008
 ** Last Modification: 10/05/2008
 ** Args:
 **    Input:
 **       - None
 **    Output:
 **       - None
 ** Return:
 **    Nothing
 ** 
 ** Description: Constructor of fsPD class
 **
 ***************************************************************************/

fsPD::fsPD() 
{
   // Clase without initialize
   
   X = NULL;
}

/****************************************************************************
 ** Name: fsPD
 ** Class: fsPD
 ** Created:           10/05/2008
 ** Last Modification: 31/05/2008
 ** Args:
 **    Input:
 **       - fsUInt numFeat: Number of features for the problem
 **    Output:
 **       - None
 ** Return:
 **    Nothing
 ** 
 ** Description: Constructor of fsPD class with the correct number of
 **              features
 **
 ***************************************************************************/

fsPD::fsPD(fsUInt numFeat) 
{
   // The number of features must be positive
   assert(numFeat>0);
   
   // We assign the features
   features = numFeat;
   
#ifdef DEBUG_fsPD
   cout << "Instantiated the number of features at " << features << endl;
#endif

   X        = new fsVector[features+1];

   for(fsUInt i=0;i<features+1;i++)
   {
      X[i] = fsVector(features);
   }
   
#ifdef DEBUG_fsPD
   cout << "Instantiated the vector X" << endl;
#endif
   
   XD       = fsVector(features);
   
#ifdef DEBUG_fsPD
   cout << "Instantiated the vector XD" << endl;
#endif

   vJ       = fsVector(features+1);
   
#ifdef DEBUG_fsPD
   cout << "Instantiated the vector vJ" << endl;
#endif

   k        = 0;

   add_est  = -1;
   
   sub_est  = -1;
      
   add      = -1;
   
   sub      = -1;
   
   last     = -1;

   best_of_all_J = -1;

   position_best_J = -1;

#ifdef DEBUG_fsPD
   cout << "Instantiated all control values" << endl;
#endif
}

/****************************************************************************
 ** Name: ~fsPD
 ** Class: fsPD
 ** Created:           10/05/2008
 ** Last Modification: 10/05/2008
 ** Args:
 **    Input:
 **       - None
 **    Output:
 **       - None
 ** Return:
 **    Nothing
 ** 
 ** Description: Destructor of fsPata class
 **
 ***************************************************************************/

fsPD::~fsPD()
{
#ifdef DEBUG_fsPD
   cout << "I go inside fsPD destructor" << endl;
#endif

   if(X != NULL)
   {
      // En teoria debiera hacer un delete sobre X, pero como da errores, solo igualo a NULL
      //delete X;
      
      X = NULL;
   }
#ifdef DEBUG_fsPD
   cout << "I go out fsPD destructor" << endl;
#endif
}

/****************************************************************************
 ** Name: operator<<
 ** Class: fsPD
 ** Created:           10/05/2008
 ** Last Modification: 10/05/2008
 ** Args:
 **    Input:
 **       - ostream &os: output channel
 **       - const fsPD &pd: problem data to be printed
 **    Output:
 **       - None
 ** Return:
 **    ostream&
 **
 ** Description: Overload of the << operator
 **
 ***************************************************************************/

ostream& operator<<(ostream &os,const fsPD &pd)
{
   // Si ya se ha iniciado el estudio, entonces si puedo mostrar la informacion
   // de todo
   if(pd.k>0)
   {
      cout << "============================================================" << endl;
      cout << "   k is " << pd.k << endl;
      cout << "   features are " << pd.features << endl;
      cout << "   XD is " << pd.XD;
      cout << "============================================================" << endl << endl << endl;

      fsVector tmp = pd.vJ;

      for(fsUInt i=0;i<=pd.features;i++)
      {
         cout << "J is " << tmp[i];
         cout << " and its X set is " << pd.X[i];
      }

      cout << "============================================================" << endl << endl << endl;

      cout << "add_est is " << pd.add_est << endl;
      cout << "sub_est is " << pd.sub_est << endl;
      
      cout << "add is " << pd.add << endl;
      cout << "sub is " << pd.sub << endl;
      
      cout << "last is " << pd.last << endl;

      cout << "the best J for this thread is " << pd.best_of_all_J << endl;

      cout << "the position were best J was got is " << pd.position_best_J << endl;
   }
   else
      cout << "No data can't be shown about this feature selection problem" << endl;

   return os;
}

/****************************************************************************
 ** Name: operator=
 ** Class: fsPD
 ** Created:           10/05/2008
 ** Last Modification: 10/05/2008
 ** Args:
 **    Input:
 **       - const fsPD &pd: problem data to be printed
 **    Output:
 **       - None
 ** Return:
 **    ostream&
 **
 ** Description: Overload of the << operator
 **
 ***************************************************************************/

// Sobecarga del operador de asignacion
fsPD& fsPD::operator=(const fsPD &pd)
{
#ifdef DEBUG_fsPD
   cout << "Assign between fsPD starts" << endl;
#endif   
   features = pd.features;
   
#ifdef DEBUG_fsPD
   cout << "I get " << pd.features << " features" << endl;
#endif   
   
   //if(X != NULL)
   //   delete X;
   
   X = NULL;
   
#ifdef DEBUG_fsPD
   cout << "Reserve memory" << endl;
#endif   
   X = new fsVector[features+1];

#ifdef DEBUG_fsPD
   cout << "Reserved the memory for " << features << " features" << endl;
#endif   

   for(fsUInt i=0;i<=features;i++)
   {
      cerr << i << endl;
      X[i]       = pd.X[i];
   }
   
   XD      = pd.XD;
      
   vJ      = pd.vJ;
      
   k       = pd.k;
      
   add_est = pd.add_est;
   sub_est = pd.sub_est;
     
   add     = pd.add;
   sub     = pd.sub;
   
   last    = pd.last;

   best_of_all_J = pd.best_of_all_J;

   position_best_J = pd.position_best_J;
   
#ifdef DEBUG_fsPD
   cout << "Assign between fsPD ended" << endl;
#endif   
   return  *this;
}
