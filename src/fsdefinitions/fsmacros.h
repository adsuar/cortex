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
 ** File: fsmacros.h
 ** Description: Definition of global macros
 ***************************************************************************/

#ifndef _MACROS_H_
#define _MACROS_H_

#define PAUSE                             system("pause")

#define SWAP_fsInt(a,b)                   {fsInt tmp;tmp=a;a=b;b=tmp;}
#define SWAP_fsUInt(a,b)                  {fsUInt tmp;tmp=a;a=b;b=tmp;}
#define SWAP_fsReal(a,b)                  {fsReal tmp;tmp=a;a=b;b=tmp;}
#define Boundary_K(k)                     (k<0 || k>((*evaluator).numberOfFeatures()))?true:false

// Macro that increases the value of the current main_k
#define incK()                            main_k++

// Macro that decreases the value of the current main_k
#define decK()                            main_k--

// Macro that counts a new J evaluation
#define next_J()                          executedJs++

// Macro that updates the best J
#define J_best(value,opc) if(value>threads[opc].best_of_all_J) {threads[opc].best_of_all_J = value; threads[opc].position_best_J = executedJs;}; addOfAllJ += value;

// Macro posts information in a file
#define J_post(value,sizeOfFeat)        {alg_plot_J << executedJs << " " << value << " " << addOfAllJ/executedJs << " " << sizeOfFeat << endl;}


#endif //_MACROS_H_
