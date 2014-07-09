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
 ** File: fstypes.h
 ** Description: Definition of the new data types
 ***************************************************************************/

#ifndef _DATA_TYPES_H_
#define _DATA_TYPES_H_

typedef int   fsInt;                      // Integer
typedef unsigned int fsUInt;              // Unsigned Integer
typedef signed int fsSInt;                // Signed Integer
typedef float fsReal;                     // Real
typedef double fsDouble;                  // Double
typedef bool fsBool;                      // Boolean
typedef char fsChar;                      // Character
typedef void fsVoid;                      // Void

typedef enum {
           EUCLIDEANDISTANCE =1,         // Euclidean Distance
           SQUAREEUCLIDEANDISTANCE=2,    // Square Euclidean Distance
           UNKNOWN=3                     // Unknown
        } fsDistanceAlgorithm;           // Distance Algorithm

typedef enum {
           NN  = 0,                      // Nearest Neighbours
           SC  = 1,                      // Scorer
           NB  = 2,                      // Naive Bayes
           LDA = 3,                      // Naive Bayes
           QDA = 4,                      // Naive Bayes
           UK  = 5                       // Unknown
        } fsEvaluationAlgorithm;

typedef enum {
           NOALGORITHM = -1,             // Algorithm not defined
           GRASEA      =  0,             // Granularity Search
           SWESEA      =  1,             // Sweep Search
           CODSEA      =  2,             // Co-dichotomous Search
           MIXSEA      =  3,             // Mixed Search
           SFS         =  4,             // Sequential Forward Search
           SBS         =  5,             // Sequential Backward Search
           SFFS        =  6,             // Sequential Forward Floating Search
           SBFS        =  7,             // Sequential Backward Floating Search
           PATA        =  8,             // Plus and Take Away Generic Algorithm
           EXHAUSTIVE  =  9,             // Exhaustive Algorithm
        } fsAlgorithm;

typedef enum {
           RANDOM      = -1,             // Data in partitions is randomized for
                                         // each partition
           TOONE       =  0,             // Data in each partition is 0
           TOZERO      =  1,             // Data in each partition is 1
           LASTDONE    =  2              // Data stays as is was at last 
                                         // optimization
        } fsSubPartitionsTreatment;

typedef enum {
           TRAINING    = 0,              // Data Training
           TEST        = 1               // Data Test
        } fsDataType;

typedef struct {
           fsAlgorithm algorithm;        // Algorithm to be used at each
                                         // partition
           fsSubPartitionsTreatment
                       treatment;        // How subpartitions will be treated
                                         // at each iteration where won't be
                                         // used
           fsUInt k;                     // Granularity of the partitions
        } fsSweepSearchConf;

// Definitions for fsData class
typedef enum {
           NOTDEFINED = -1,              // Unknown Data
           DISCREET   = 0,               // Discreet Data
           CONTINUOUS = 1,               // Continuous Data
           BINARY     = 2,               // Binary Data
           COMPOUND   = 3                // Compound Data
        } fsDataTypes;

typedef enum {
           DATA  = 0,                    // Column with Data
           CLASS = 1                     // Column with Class
        } fsDataColumns;

#endif // _DATA_TYPES_H_
