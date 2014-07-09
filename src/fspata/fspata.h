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
 ** File: fspata.h
 ** Description: Definition of the PATA Class, that has the Plus And Take 
 **              Away algorithms
 ***************************************************************************/

#ifndef _PATA_H_
#define _PATA_H_

using namespace std;

#include "../fsset/fsset.h"
#include "../fsvector/fsvector.h"
#include "../fsmatrix/fsmatrix.h"
#include "../fsevaluation/fsevaluation.h"
#include "../fspd/fspd.h"
#include "../fsdefinitions/fstypes.h"
#include "../fsdefinitions/fsconstants.h"
#include "../fsdefinitions/fsmacros.h"

class fsPata
{
   private:
      fsEvaluation *evaluator;
      fsVector X;
      fsVoid (fsPata::*algorithm)();
      fsVoid bestFeatureToAdd(fsUInt);
      fsVoid bestFeatureToSubstract(fsUInt);
      
      fsPD *threads;

      fsUInt goForward(fsUInt,fsUInt);
      fsUInt goBackward(fsUInt,fsUInt);

      fsVoid sfs();

      fsVoid sbs();
      
      fsVoid sbfs();
      
      fsVoid sffs();
      
      fsVoid pata();
      
      fsInt stepsForward,stepsBackward;
      
      fsUInt executedJs;
      
      fsReal addOfAllJ;
      
      fstream alg_plot_J;                      // Fichero que almacenara las Js
      
      fsUInt main_k;

   public:
      fsPata(
             fsAlgorithm fsalg=SFS,
             fsEvaluationAlgorithm fsea=NN,
             fsDistanceAlgorithm fsda=SQUAREEUCLIDEANDISTANCE,
             fsUInt kFCV=0,
             fsUInt kFCVTest=0,
             fsBool r=false,
             fsChar *fsFileData=NULL,
             fsInt sf=1,
             fsInt sb=0);
      ~fsPata();

      fsVoid process();
};

#endif // _PATA_H
