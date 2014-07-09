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
 ** File: fsdata.cc
 ** Description: Source code for the fsData class
 ***************************************************************************/

#include "fsdata.h"

/****************************************************************************
 ** Name: fsData
 ** Class: fsData
 ** Created:           29/04/2008
 ** Last Modification: 29/04/2008
 ** Args:
 **    Input:
 **       - None
 **    Output:
 **       - None
 ** Return:
 **    Nothing
 ** 
 ** Description: Constructor of fsData class
 **
 ***************************************************************************/

fsData::fsData()
{
   dataLoaded = false;
}

/****************************************************************************
 ** Name: ~fsData
 ** Class: fsData
 ** Created:           29/04/2008
 ** Last Modification: 29/04/2008
 ** Args:
 **    Input:
 **       - None
 **    Output:
 **       - None
 ** Return:
 **    Nothing
 ** 
 ** Description: Destructor of fsData class
 **
 ***************************************************************************/

fsData::~fsData()
{
   dataLoaded = false;
}

/****************************************************************************
 ** Name: fsData
 ** Class: fsData
 ** Created:           01/05/2008
 ** Last Modification: 30/05/2008
 ** Args:
 **    Input:
 **       - fsChar*: File to be opened containing the items of the problem
 **    Output:
 **       - None
 ** Return:
 **    Nothing
 ** 
 ** Description: Constructor of fsData class
 **
 ***************************************************************************/

fsData::fsData(fsChar *fichero)
{
#ifdef DEBUG_fsData
   cout << "fsData with file: " << fichero << endl;
#endif

   // Descriptor of the file
   fstream alg_data;
   
   // d es el numero de columnas de datos, y c de clases
   fsUInt col_features=0,col_classes=0;
   
   // This vector will control which columns are features and which ones not
   fsVector column_type;
   
   // This vector will control the type of each column
   fsVector domains;
   
   // By default, there isn't defined any domain for classes or features
   (*this).features_domain = NOTDEFINED;
   (*this).classes_domain = NOTDEFINED;
   
   try
   {
      alg_data.open(fichero,ios::in);
   
      if(!alg_data)
      {
         throw "Data file can not be opened";
      }
   }
   catch(const fsChar *s)
   {
      cout << fichero << ": " << s << endl;
      exit(-1);
   }
   
   fsChar character;
   
   // Now, I have to read data, first of all, the type of columns, ranged
   // between '@' symbols
   // I read the first character
   character=alg_data.get();
   
   // The first character must be '@'
   assert(character == '@' );
   
   character='\0';
   
   while(character != '@')
   {
      fsChar c;
      alg_data >> c;
      
      switch(c)
      {
         case 'D':
            column_type.add(DATA);
            col_features++;
            break;
         case 'C':
            column_type.add(CLASS);
            col_classes++;
            break;
         case '@':
            character = c;
            break;
         default:
            // Error in the header
            assert(false);
            break;
      }
   }
   
   character=alg_data.get();

   // I must read a new line character
   assert(character == '\n' || ((fsInt)character == 13));
   
   //character=alg_data.get();
   alg_data >> character;

   // I must read the character '#'
   assert(character == '#');
   
   character='\0';
   
   while(character != '#')
   {
      fsChar c;
      alg_data >> c;
      
      switch(c)
      {
         case 'D':
            domains.add(DISCREET);
            break;
         case 'C':
            domains.add(CONTINUOUS);
            break;
         case 'B':
            domains.add(BINARY);
            break;
         case '#':
            character = c;
            break;
         default:
            assert(false);
            break;
      }
   }
   
   character=alg_data.get();
   
   // I must read a new line character
   assert(character == '\n' || ((fsInt)character == 13));
      
   // The two vectors have to have the same length
   assert (domains.sizeOf() == column_type.sizeOf());
   
#ifdef DEBUG_fsData
   cout << "I proceed to load data" << endl;
#endif
   
   //fsInt contador_filas = 0;
   
   // I proceed to load data
   while(!alg_data.eof())
   {
      // I create a vector to load features, and another to load classes
      fsVector row_features(col_features);
      fsVector row_classes(col_classes);
      
      fsInt num_datos=0,num_clases=0;
      
      // I read data
      for(fsUInt i=0;i<column_type.sizeOf() && !alg_data.eof();i++)
      {
         if(column_type[i] == DATA)
            alg_data >> row_features[num_datos++];
         else
            alg_data >> row_classes[num_clases++];
      }
      
      // Only if we've read an entire row without any problem, we'll be able
      // to insert data in matrices
      if(column_type.sizeOf() == (fsUInt)(num_datos + num_clases))
      {
         // I read the number of rows
         fsInt numRowsFeatures = features.numRows();
         fsInt numRowsClasses = classes.numRows();
         
         // Primero cambio el tamaño, redimensiono las matrices
         features.resize(numRowsFeatures+1,col_features);
         classes.resize(numRowsClasses+1,col_classes);
   
         features.addRow(numRowsFeatures,&row_features);
         classes.addRow(numRowsClasses,&row_classes);
         
         /*
         if(num_clases == 1)
            classes.add(row_classes[0]);
         else
            assert(false);
         */
      }
   }

#ifdef DEBUG_fsData
   cout << "Data loaded in matrices" << endl;
#endif
   
   // Once I've splitted data features and data classes, I've to know if they're
   // DISCREET or not
   for(fsUInt i=0;i<column_type.sizeOf();i++)
   {
      if(column_type[i] == DATA)
      {
         data_features.add(domains[i]);

         if(domains[i] == (*this).features_domain)
            ;
         else
         if((*this).features_domain == -1)
            (*this).features_domain = (fsDataTypes)domains[i];
         else
            (*this).features_domain = COMPOUND;         
      }
      else
      {
         data_classes.add(domains[i]);

         if(domains[i] == (*this).classes_domain)
            ;
         else
         if((*this).classes_domain == -1)
            (*this).classes_domain = (fsDataTypes)domains[i];
         else
            (*this).classes_domain = COMPOUND;         
      }
   }
   
   // I close the file
   alg_data.close();
   
   dataLoaded = true;
   
#ifdef DEBUG_fsData
   cout << "Structure loaded successfully" << endl;
   cout << "It has: " << endl;
   cout << " - " << data_features.sizeOf() << " feature column(s)" << endl;
   cout << " - " << data_classes.sizeOf() << " class column(s)" << endl;
   cout << " - " << features.numRows() << " instance(s)" << endl << endl;
#endif
}

/****************************************************************************
 ** Name: normalizeFeatures
 ** Class: fsData
 ** Created:           01/05/2008
 ** Last Modification: 01/05/2008
 ** Args:
 **    Input:
 **       - None
 **    Output:
 **       - None
 ** Return:
 **    - None
 ** 
 ** Description: Method that normalizes CONTINUOUS features, not DISCREET, nor
 **              BINARY
 **
 ***************************************************************************/

fsVoid fsData::normalizeFeatures()
{
   for(fsUInt i=0;i<data_features.sizeOf();i++)
   {
      if(data_features[i] == CONTINUOUS)
      {
         features.normalizeCol(i);
      }
   }
}

/****************************************************************************
 ** Name: showFeatures
 ** Class: fsData
 ** Created:           01/05/2008
 ** Last Modification: 01/05/2008
 ** Args:
 **    Input:
 **       - None
 **    Output:
 **       - None
 ** Return:
 **    - None
 ** 
 ** Description: Method that shows the features of each item
 **
 ***************************************************************************/

fsVoid fsData::showFeatures()
{
   cout << "=====================================================" << endl;
   cout << " Features values by items are:" << endl;
   cout << "=====================================================" << endl;
   cout << endl << endl;
   cout << features << endl;
   cout << "=====================================================" << endl;
}

/****************************************************************************
 ** Name: showClasses
 ** Class: fsData
 ** Created:           01/05/2008
 ** Last Modification: 01/05/2008
 ** Args:
 **    Input:
 **       - None
 **    Output:
 **       - None
 ** Return:
 **    - None
 ** 
 ** Description: Method that shows the classes of each item
 **
 ***************************************************************************/

fsVoid fsData::showClasses()
{
   cout << "=====================================================" << endl;
   cout << " Classes values by items are:" << endl;
   cout << "=====================================================" << endl;
   cout << endl << endl;
   cout << classes << endl;
   cout << "=====================================================" << endl;
}


/****************************************************************************
 ** Name: showTypesData
 ** Class: fsData
 ** Created:           01/05/2008
 ** Last Modification: 01/05/2008
 ** Args:
 **    Input:
 **       - None
 **    Output:
 **       - None
 ** Return:
 **    - None
 ** 
 ** Description: Method that shows the data type of every feature, column by
 **              column
 **
 ***************************************************************************/

fsVoid fsData::showTypesData()
{
   cout << "=====================================================" << endl;
   cout << " The domain of every feature, column by column, is:" << endl;
   cout << "=====================================================" << endl;
   cout << endl << endl;
   cout << data_features << endl;
   cout << "=====================================================" << endl;
}

/****************************************************************************
 ** Name: getClasses
 ** Class: fsData
 ** Created:           01/05/2008
 ** Last Modification: 01/05/2008
 ** Args:
 **    Input:
 **       - None
 **    Output:
 **       - None
 ** Return:
 **    - fsMatrix*: Pointer to classes data matrix
 ** 
 ** Description: Method that returns the pointer to the classes data matrix
 **
 ***************************************************************************/

fsMatrix* fsData::getClasses()
{
   return &classes;
}

/****************************************************************************
 ** Name: getFeatures
 ** Class: fsData
 ** Created:           01/05/2008
 ** Last Modification: 01/05/2008
 ** Args:
 **    Input:
 **       - None
 **    Output:
 **       - None
 ** Return:
 **    - fsMatrix*: Pointer to features data matrix
 ** 
 ** Description: Method that returns the pointer to the features data matrix
 **
 ***************************************************************************/

fsMatrix* fsData::getFeatures()
{
   return &features;
}


/****************************************************************************
 ** Name: getInstances
 ** Class: fsData
 ** Created:           29/04/2008
 ** Last Modification: 29/04/2008
 ** Args:
 **    Input:
 **       - None
 **    Output:
 **       - None
 ** Return:
 **    - fsInt: Number of instances (rows)
 ** 
 ** Description: Method that returns the number of instances (or rows) of
 **              the data
 **
 ***************************************************************************/

fsUInt fsData::getInstances()
{
   return features.numRows();
}

/****************************************************************************
 ** Name: getTypeFeatures
 ** Class: fsData
 ** Created:           29/04/2008
 ** Last Modification: 29/04/2008
 ** Args:
 **    Input:
 **       - None
 **    Output:
 **       - None
 ** Return:
 **    - fsDataTypes: Type of data features
 ** 
 ** Description: Method that returns the type of the features
 **
 ***************************************************************************/

fsDataTypes fsData::getTypeFeatures()
{
   return (*this).features_domain;
}

/****************************************************************************
 ** Name: getTypeClasses
 ** Class: fsData
 ** Created:           29/04/2008
 ** Last Modification: 29/04/2008
 ** Args:
 **    Input:
 **       - None
 **    Output:
 **       - None
 ** Return:
 **    - fsDataTypes: Type of data classes
 ** 
 ** Description: Method that returns the type of the classes
 **
 ***************************************************************************/

fsDataTypes fsData::getTypeClasses()
{
   return (*this).classes_domain;
}

/****************************************************************************
 ** Name: operator=
 ** Class: fsData
 ** Created:           29/04/2008
 ** Last Modification: 29/04/2008
 ** Args:
 **    Input:
 **       - const fsData &d: data to be assigned
 **    Output:
 **       - None
 ** Return:
 **    - fsData&: Memory address of the data assigned
 ** 
 ** Description: Overload of the assign operator
 **
 ***************************************************************************/

fsData& fsData::operator=(const fsData &d)
{
   data_features   = d.data_features;
   data_classes    = d.data_classes;
   features        = d.features;
   classes         = d.classes;
   features_domain = d.features_domain;
   classes_domain  = d.classes_domain;
  
   return *this;
}

/****************************************************************************
 ** Name: successfullyLoaded
 ** Class: fsData
 ** Created:           01/05/2008
 ** Last Modification: 01/05/2008
 ** Args:
 **    Input:
 **       - None
 **    Output:
 **       - None
 ** Return:
 **    - fsBool: Boolean that indicates if data has been or not loaded
 ** 
 ** Description: Method that indicates if data has been or not loaded
 **
 ***************************************************************************/

fsBool fsData::successfullyLoaded()
{
   return dataLoaded;
}

/****************************************************************************
 ** Name: getNumberOfClasses
 ** Class: fsData
 ** Created:           04/05/2008
 ** Last Modification: 04/05/2008
 ** Args:
 **    Input:
 **       - None
 **    Output:
 **       - None
 ** Return:
 **    - fsUInt: Number of classes
 ** 
 ** Description: Method that returns the number of classes of the problem
 **
 ***************************************************************************/

fsUInt fsData::getNumberOfClasses()
{
   assert(dataLoaded);
   
   return classes.numCols();
}

/****************************************************************************
 ** Name: getNumberOfFeatures
 ** Class: fsData
 ** Created:           04/05/2008
 ** Last Modification: 04/05/2008
 ** Args:
 **    Input:
 **       - None
 **    Output:
 **       - None
 ** Return:
 **    - fsUInt: Number of features
 ** 
 ** Description: Method that returns the number of features of the problem
 **
 ***************************************************************************/

fsUInt fsData::getNumberOfFeatures()
{
   assert(dataLoaded);
   
   return features.numCols();
}
