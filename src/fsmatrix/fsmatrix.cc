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
 ** File: fsmatrix.cc
 ** Description: Source code for the fsMatrix class
 ***************************************************************************/

#include "fsmatrix.h"

/****************************************************************************
 ** Name: fsMatrix
 ** Class: fsMatrix
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
 ** Description: Constructor of fsMatrix class
 **
 ***************************************************************************/

fsMatrix::fsMatrix()
{
   // There's no need to initialize anything
   // Initially, the STL vector is already void
   // so we do not need to do anything more
}

/****************************************************************************
 ** Name: ~fsMatrix
 ** Class: fsMatrix
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
 ** Description: Destructor of fsMatrix class
 **
 ***************************************************************************/

fsMatrix::~fsMatrix()
{
   // If matrix isn't empty, we empty every row
   if(matrix.size() != 0)
   {
      for(fsUInt i=0;i<matrix.size();i++)
         matrix[i].clear();
      matrix.clear();
   }

   // Destructor of matrix is called automatically at the exit of this method
}

/****************************************************************************
 ** Name: fsMatrix
 ** Class: fsMatrix
 ** Created:           28/04/2008
 ** Last Modification: 28/04/2008
 ** Args:
 **    Input:
 **       - fsUInt x: Number of rows of the matrix
 **       - fsUInt y: Number of columns of the matrix
 **    Output:
 **       - None
 ** Return:
 **    Nothing
 ** 
 ** Description: Constructor of fsMatrix class, given the dimensions
 **
 ***************************************************************************/

fsMatrix::fsMatrix(fsUInt x, fsUInt y)
{
   // Rows and Columns must be greater than 0
   assert(x>0 && y>0);
   
   // We create a Real vector, having as many positions as columns
   vector<fsReal> nuevo(y);
   
   // At every row, we push a Real vector with y columns
   while((fsUInt)matrix.size() != x)
   {
      matrix.push_back(nuevo);
   }
}

/****************************************************************************
 ** Name: resize
 ** Class: fsMatrix
 ** Created:           28/04/2008
 ** Last Modification: 28/04/2008
 ** Args:
 **    Input:
 **       - fsUInt x: Number of rows of the matrix
 **       - fsUInt y: Number of columns of the matrix
 **    Output:
 **       - None
 ** Return:
 **    Nothing
 ** 
 ** Description: Method that changes the dimensions of the matrix
 **
 ***************************************************************************/

fsVoid fsMatrix::resize(fsUInt x,fsUInt y)
{
   // Rows and Columns must be greater than 0
   assert(x>0 && y>0);
   
   // First of all, I adjust rows
   if(x<(fsUInt)matrix.size())                 // If there's more rows than
                                               // needed
      while((fsUInt)matrix.size() != x)
         matrix.pop_back();
   else
   if(x>(fsUInt)matrix.size())                 // If we need more rows
   {
      fsInt num_cols=0;                        // We insert needed rows with
                                               // previous number of columns
      
      if(matrix.size() > 0)
         num_cols = matrix[0].size();
      
      vector<fsReal> nuevo(num_cols);
      
      
      while((fsUInt)matrix.size() != x)
         matrix.push_back(nuevo);
   }
   
   // We adjust now the number of columns
   if((fsUInt)matrix[0].size() > y)            // If there's more columns
                                               // than needed
   {
      for(fsUInt i=0;i<(fsUInt)matrix.size();i++)
         while((fsUInt)matrix[i].size() != y)
            matrix[i].pop_back();
   }
   else
   if((fsUInt)matrix[0].size() < y)            // If we need more columns
   {
      for(fsUInt i=0;i<(fsUInt)matrix.size();i++)
         while((fsUInt)matrix[i].size() != y)
            matrix[i].push_back(0);
   }
}


/****************************************************************************
 ** Name: addRow
 ** Class: fsMatrix
 ** Created:           28/04/2008
 ** Last Modification: 28/04/2008
 ** Args:
 **    Input:
 **       - fsUInt pos: Position where vector will be inserted
 **       - fsVector *v: Vector to be inserted
 **    Output:
 **       - None
 ** Return:
 **    Nothing
 **
 ** Description: Method that inserts a vector in a given row
 **
 ***************************************************************************/

fsVoid fsMatrix::addRow(fsUInt pos, fsVector *v)
{
   // Position has to be between the ranges
   assert(pos>=0 && pos<(fsUInt)matrix.size());
   
   // Matrix columns and vector must have the same lentgh
   assert((fsUInt)matrix[pos].size()==(*v).sizeOf());
   
   for(fsUInt i=0;i<(fsUInt)matrix[pos].size();i++)
      matrix[pos][i] = (*v)[i];
}

/****************************************************************************
 ** Name: addCol
 ** Class: fsMatrix
 ** Created:           28/04/2008
 ** Last Modification: 28/04/2008
 ** Args:
 **    Input:
 **       - fsUInt pos: Position where vector will be inserted
 **       - fsVector *v: Vector to be inserted
 **    Output:
 **       - None
 ** Return:
 **    Nothing
 **
 ** Description: Method that inserts a vector in a given column
 **
 ***************************************************************************/

fsVoid fsMatrix::addCol(fsUInt pos, fsVector *v)
{
   // Matrix rows and vector must have the same lentgh
   assert((fsUInt)matrix.size()==(*v).sizeOf());
   
   // Position has to be between the ranges
   assert(pos>=0 && pos<(fsUInt)matrix[0].size());
   
   for(fsInt i=0;i<(fsSInt)matrix.size();i++)
      matrix[i][pos] = (*v)[i];
}

/****************************************************************************
 ** Name: extractRow
 ** Class: fsMatrix
 ** Created:           28/04/2008
 ** Last Modification: 28/04/2008
 ** Args:
 **    Input:
 **       - fsUInt row: Position to be extracted
 **    Output:
 **       - None
 ** Return:
 **    fsVector
 **
 ** Description: Method that gets one row of the matrix
 **
 ***************************************************************************/

fsVector fsMatrix::extractRow(fsUInt row)
{
   // Row must be between the ranges
   assert(row>=0 && row<(fsUInt)matrix.size());

   fsVector v(matrix[row]);
   /*
   fsVector v;
   for(fsInt i=0;i<(fsSInt)matrix[row].size();i++)
      v.add(matrix[row][i]);
   */
   return v;
}

/****************************************************************************
 ** Name: extractCol
 ** Class: fsMatrix
 ** Created:           28/04/2008
 ** Last Modification: 28/04/2008
 ** Args:
 **    Input:
 **       - fsUInt col: Column to be extracted
 **    Output:
 **       - None
 ** Return:
 **    fsVector
 **
 ** Description: Method that gets one column of the matrix
 **
 ***************************************************************************/

fsVector fsMatrix::extractCol(fsUInt col)
{
   // There must be one row at least
   assert((fsSInt)matrix.size()>0);
   
   // Column must be between the ranges
   assert(col>=0 && col<(fsUInt)matrix[0].size());
   
   fsVector v;

   for(fsInt i=0;i<(fsSInt)matrix.size();i++)
      v.add(matrix[i][col]);
   
   return v;
}

/****************************************************************************
 ** Name: extractRowSTL
 ** Class: fsMatrix
 ** Created:           28/04/2008
 ** Last Modification: 28/04/2008
 ** Args:
 **    Input:
 **       - fsUInt row: Row to be extracted
 **    Output:
 **       - None
 ** Return:
 **    vector<fsReal>*
 **
 ** Description: Method that gets the pointer to one row of the matrix
 **              avoiding to make a copy to be returned
 **
 ***************************************************************************/

vector<fsReal> *fsMatrix::extractRowSTL(fsUInt row)
{
   // Row must be between the ranges
   assert(row>=0 && row<(fsUInt)matrix.size());

   return &matrix[row];
}

/****************************************************************************
 ** Name: operator<<
 ** Class: fsMatrix
 ** Created:           28/04/2008
 ** Last Modification: 28/04/2008
 ** Args:
 **    Input:
 **       - ostream &os: output channel
 **       - const fsMatrix &c: matrix to be printed
 **    Output:
 **       - None
 ** Return:
 **    ostream&
 **
 ** Description: Overload of the << operator
 **
 ***************************************************************************/

ostream& operator<<(ostream &os,const fsMatrix &c)
{
   if(c.matrix.size() <= 0)
      os << "Empty matrix\n";
   else
   {
      for(fsUInt i=0;i<c.matrix.size();i++)
      {
         for(fsUInt j=0;j<c.matrix[i].size();j++)
            os << c.matrix[i][j] << " ";
         os << endl;
      }
      
      os << "\n";
   }
   
   return os;
}

/****************************************************************************
 ** Name: averageCol
 ** Class: fsMatrix
 ** Created:           28/04/2008
 ** Last Modification: 28/04/2008
 ** Args:
 **    Input:
 **       - fsUInt col: column to get the average value
 **    Output:
 **       - None
 ** Return:
 **    fsReal
 **
 ** Description: Method that get the average value of a given column
 **
 ***************************************************************************/

fsReal fsMatrix::averageCol(fsUInt col)
{
   // Matrix must be initialized
   assert((fsUInt)matrix.size()>0 && (fsUInt)matrix[0].size()>0);
   
   // Column must be between the ranges
   assert(col>=0 && col<(fsUInt)matrix[0].size());
      
   fsReal total=0,tamx=(fsUInt)matrix.size();
   
   for(fsInt i=0;i<tamx;i++)
      total+=matrix[i][col];
   
   return (total/tamx);
}

/****************************************************************************
 ** Name: averageRow
 ** Class: fsMatrix
 ** Created:           28/04/2008
 ** Last Modification: 28/04/2008
 ** Args:
 **    Input:
 **       - fsUInt row: row to get the average value
 **    Output:
 **       - None
 ** Return:
 **    fsReal
 **
 ** Description: Method that get the average value of a given row
 **
 ***************************************************************************/

fsReal fsMatrix::averageRow(fsUInt row)
{
   // Matrix must be initialized
   assert((fsUInt)matrix.size()>0 && (fsUInt)matrix[0].size()>0);
   
   // Row must be between the ranges
   assert(row>=0 && row<(fsUInt)matrix.size());

   fsReal total=0,tamy=(fsUInt)matrix[row].size();
   
   for(fsInt i=0;i<tamy;i++)
      total+=matrix[row][i];
   
   return (total/tamy);
}

/****************************************************************************
 ** Name: moduleCol
 ** Class: fsMatrix
 ** Created:           28/04/2008
 ** Last Modification: 28/04/2008
 ** Args:
 **    Input:
 **       - fsUInt col: col to get the module value
 **    Output:
 **       - None
 ** Return:
 **    fsReal
 **
 ** Description: Method that get the module value of a given column
 **
 ***************************************************************************/

fsReal fsMatrix::moduleCol(fsUInt col)
{
   // Matrix must be initialized
   assert((fsUInt)matrix.size()>0 && (fsUInt)matrix[0].size()>0);
   
   // Column must be between the ranges
   assert(col>=0 && col<(fsUInt)matrix[0].size());

   fsReal total=0,tamx=(fsUInt)matrix.size();
   
   for(fsInt i=0;i<tamx;i++)
      total+=matrix[i][col]*matrix[i][col];
   
   total = sqrt(total);
   
   return total;
}

/****************************************************************************
 ** Name: moduleRow
 ** Class: fsMatrix
 ** Created:           28/04/2008
 ** Last Modification: 28/04/2008
 ** Args:
 **    Input:
 **       - fsUInt row: row to get the module value
 **    Output:
 **       - None
 ** Return:
 **    fsReal
 **
 ** Description: Method that get the module value of a given row
 **
 ***************************************************************************/

fsReal fsMatrix::moduleRow(fsUInt row)
{
   // Matrix must be initialized
   assert((fsUInt)matrix.size()>0 && (fsUInt)matrix[0].size()>0);
   
   // Row must be between the ranges
   assert(row>=0 && row<(fsUInt)matrix.size());

   fsReal total=0,tamy=(fsUInt)matrix[0].size();
   
   for(fsInt i=0;i<tamy;i++)
      total+=matrix[row][i]*matrix[row][i];
   
   total = sqrt(total);
   
   return total;
}

/****************************************************************************
 ** Name: varianceCol
 ** Class: fsMatrix
 ** Created:           28/04/2008
 ** Last Modification: 28/04/2008
 ** Args:
 **    Input:
 **       - fsUInt col: col to get the module value
 **    Output:
 **       - None
 ** Return:
 **    fsReal
 **
 ** Description: Method that get the variance value of a given column
 **
 ***************************************************************************/

fsReal fsMatrix::varianceCol(fsUInt col)
{
   // Matrix must be initialized
   assert((fsUInt)matrix.size()>0 && (fsUInt)matrix[0].size()>0);
   
   // Column must be between the ranges
   assert(col>=0 && col<(fsUInt)matrix[0].size());
   
   fsReal varianza=0,tamx=(fsUInt)matrix.size();
   
   fsReal avg=this->averageCol(col);
   
   for(fsInt i=0;i<tamx;i++)
   {
      fsReal value = matrix[i][col]-avg;
      varianza+=value*value;
   }
   
   // Aqui aparece un error ya que la varianza la divide entre tamx-1
   // y debe ser entre tamx!!!
   varianza = varianza/tamx;
   
   return varianza;
}

/****************************************************************************
 ** Name: varianceRow
 ** Class: fsMatrix
 ** Created:           28/04/2008
 ** Last Modification: 28/04/2008
 ** Args:
 **    Input:
 **       - fsUInt row: row to get the module value
 **    Output:
 **       - None
 ** Return:
 **    fsReal
 **
 ** Description: Method that get the variance value of a given row
 **
 ***************************************************************************/

fsReal fsMatrix::varianceRow(fsUInt row)
{
   // Matrix must be initialized
   assert((fsUInt)matrix.size()>0 && (fsUInt)matrix[0].size()>0);
   
   // Row must be between the ranges
   assert(row>=0 && row<(fsUInt)matrix.size());

   fsReal varianza=0,tamy=(fsSInt)matrix[0].size();
   
   fsReal avg=this->averageRow(row);
   
   for(fsInt i=0;i<tamy;i++)
   {
      fsReal value = matrix[row][i]-avg;
      varianza+=value*value;
   }
   
   // Aqui aparece un error ya que la varianza la divide entre tamy-1
   // y debe ser entre tamy!!!
   varianza = varianza/tamy;
   
   return varianza;
}

/****************************************************************************
 ** Name: normalizeCol
 ** Class: fsMatrix
 ** Created:           28/04/2008
 ** Last Modification: 28/04/2008
 ** Args:
 **    Input:
 **       - fsUInt col: col to be normalized
 **    Output:
 **       - None
 ** Return:
 **    - None
 **
 ** Description: Method that normalizes one column
 **
 ***************************************************************************/

fsVoid fsMatrix::normalizeCol(fsUInt col)
{
   // Matrix must be initialized
   assert((fsUInt)matrix.size()>0 && (fsUInt)matrix[0].size()>0);
   
   // Column must be between the ranges
   assert(col>=0 && col<(fsUInt)matrix[0].size());

   // First of all, we get the average value
   fsReal avg = averageCol(col);
   
   // We get the typical deviation
   fsReal desvTip = sqrt(this->varianceCol(col));
   
   // We get the number of rows
   fsReal tamx=matrix.size();
   
   // We calculate the module but, if typical deviation is 0, the module is 0
   for(fsInt i=0;i<tamx;i++)
   {
      if(desvTip != 0)
         matrix[i][col] = (matrix[i][col]-avg)/desvTip;
      else
         matrix[i][col] = 0;
   }
}

/****************************************************************************
 ** Name: normalizeRow
 ** Class: fsMatrix
 ** Created:           28/04/2008
 ** Last Modification: 28/04/2008
 ** Args:
 **    Input:
 **       - fsUInt row: row to be normalized
 **    Output:
 **       - None
 ** Return:
 **    - None
 **
 ** Description: Method that normalizes one row
 **
 ***************************************************************************/

fsVoid fsMatrix::normalizeRow(fsUInt row)
{
   // Matrix must be initialized
   assert((fsSInt)matrix.size()>0 && (fsSInt)matrix[0].size()>0);
   
   // Row must be between the ranges
   assert(row>=0 && row<(fsUInt)matrix.size());

   // First of all, we get the average value
   fsReal avg = this->averageRow(row);
   
   // We get the typical deviation
   fsReal desvTip = sqrt(this->varianceRow(row));
   
   // We get the number of columns
   fsReal tamy=matrix[0].size();
   
   // We calculate the module but, if typical deviation is 0, the module is 0
   for(fsInt i=0;i<tamy;i++)
   {
      if(desvTip != 0)
         matrix[row][i] = (matrix[row][i]-avg)/desvTip;
      else
         matrix[row][i] = 0;
   }
}

/****************************************************************************
 ** Name: normalizeByCols
 ** Class: fsMatrix
 ** Created:           28/04/2008
 ** Last Modification: 28/04/2008
 ** Args:
 **    Input:
 **       - fsUInt start: first column to be normalized
 **       - fsUInt end: last column to be normalized
 **    Output:
 **       - None
 ** Return:
 **    - None
 **
 ** Description: Method that normalizes a range of columns
 **
 ***************************************************************************/

fsVoid fsMatrix::normalizeByCols(fsUInt start,fsUInt end)
{
   // Matrix must be initialized
   assert((fsSInt)matrix.size()>0 && (fsUInt)matrix[0].size()>0);

   // start and end indices must be correct, kept on the right range
   assert(start>=0 && start<=end && end<(fsUInt)matrix[0].size());
    
   for(fsUInt i=start;i<=end;i++)
      this->normalizeCol(i);
}

/****************************************************************************
 ** Name: normalizeByRows
 ** Class: fsMatrix
 ** Created:           28/04/2008
 ** Last Modification: 28/04/2008
 ** Args:
 **    Input:
 **       - fsUInt start: first row to be normalized
 **       - fsUInt end: last row to be normalized
 **    Output:
 **       - None
 ** Return:
 **    - None
 **
 ** Description: Method that normalizes a range of rows
 **
 ***************************************************************************/

fsVoid fsMatrix::normalizeByRows(fsUInt start,fsUInt end)
{
   // Matrix must be initialized
   assert((fsSInt)matrix.size()>0 && (fsUInt)matrix[0].size()>0);

   // start and end indices must be correct, kept on the right range
   assert(start>=0 && start<=end && end<(fsUInt)matrix.size());

   for(fsUInt i=start;i<=end;i++)
      this->normalizeRow(i);
}

/****************************************************************************
 ** Name: numRows
 ** Class: fsMatrix
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
 ** Description: Method that gets the number of rows of the matrix
 **
 ***************************************************************************/

fsUInt fsMatrix::numRows()
{
   return matrix.size();
}

/****************************************************************************
 ** Name: numCols
 ** Class: fsMatrix
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
 ** Description: Method that gets the number of columns of the matrix
 **
 ***************************************************************************/

fsUInt fsMatrix::numCols()
{
   if(matrix.size() > 0)
      return matrix[0].size();
   
   return 0;
}

/****************************************************************************
 ** Name: operator=
 ** Class: fsMatrix
 ** Created:           28/04/2008
 ** Last Modification: 28/04/2008
 ** Args:
 **    Input:
 **       - const fsMatrix &m: matrix to be assigned
 **    Output:
 **       - None
 ** Return:
 **    - fsMatrix&: memory address of the result of the assignation
 **
 ** Description: Method that assigns a matrix to the propietary one
 **
 ***************************************************************************/

fsMatrix& fsMatrix::operator=(const fsMatrix &m)
{
   matrix = m.matrix;
     
   return *this;
}
