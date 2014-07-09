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
 ** File: fsaprs.h
 ** Description: Definition of the APRS Class, that has the Adaptive 
 **              Partitioned Random Search algorithms
 ***************************************************************************/

#ifndef _APRS_H_
#define _APRS_H_

using namespace std;

#include "../fsset/fsset.h"
#include "../fsvector/fsvector.h"
#include "../fsmatrix/fsmatrix.h"
#include "../fsevaluation/fsevaluation.h"
#include "../fsdefinitions/fstypes.h"
#include "../fsdefinitions/fsconstants.h"
#include "../fsdefinitions/fsmacros.h"

class fsAprs
{
   private:
      fsEvaluation *evaluator;

      fsSweepSearchConf sweep;
      
      fsVoid (fsAprs::*algorithm)();

      fsVoid granularitySearch();

      fsVoid sweepSearch();
      
      fsVoid codichotomousSearch();

      fsVoid mixedSearch();

      // I've to improve the apply of the BinarySelector
      fsVoid applyBinarySelector(fsVector*,fsVector*,fsVector*);
      fsVoid applyBinarySelector(fsUInt,fsUInt,fsVector*);
      
      fsVoid subSFS(fsUInt, fsUInt, fsVector*);
      fsVoid subSBS(fsUInt, fsUInt, fsVector*);
      fsVoid subExhaustive(fsUInt, fsUInt, fsVector*);

      fsVoid F(fsSet *,fsUInt);
      fsVoid B(fsSet *,fsUInt);

      fsSet setOfFeatures;

      fsVector Xdone;

      fsUInt Js;
      
   public:
      fsAprs(
             fsAlgorithm fsalg=GRASEA,
             fsEvaluationAlgorithm fsea=NN,
             fsDistanceAlgorithm fsda=SQUAREEUCLIDEANDISTANCE,
             fsUInt kFCV=0,
             fsUInt kFCVTest=0,
             fsBool r=false,
             fsChar *fsFileData=NULL,
             fsAlgorithm fsssalg=SFS,
             fsSubPartitionsTreatment fsssspt=TOONE,
             fsUInt k=5
      );
      ~fsAprs();

      fsVoid process();
};

#endif // _APRS_H
