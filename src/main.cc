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
 ** File: main.cc
 ** Description: Execution manager that allows the launch of different
 **              algorithms of feature selection.
 ***************************************************************************/
#include "main.h"

/****************************************************************************
 ** Name: echoAlgorithmName
 ** Class: None
 ** Created:           31/05/2008
 ** Last Modification: 10/06/2008
 ** Args:
 **    Input:
 **       - fsAlgorithm alg: Code of the algorithm
 **    Output:
 **       - None
 ** Return:
 **    fsChar*: Name of the algorithm
 ** 
 ** Description: Function that returns the name of an algorithm, given its
 **              code
 **
 ***************************************************************************/
fsChar *echoAlgorithmName(fsAlgorithm alg)
{
   fsChar *algorithmName;
   algorithmName = (fsChar*)malloc(64*sizeof(fsChar));

   switch(alg)
   {
      case GRASEA:
         sprintf(algorithmName,"Granularity Search");
         break;
      case SWESEA:
         sprintf(algorithmName,"Sweep Search");
         break;
      case CODSEA:
         sprintf(algorithmName,"Co-dichotomous Search");
         break;
      case MIXSEA:
         sprintf(algorithmName,"Mixed Search");
         break;
      case SFS:
         sprintf(algorithmName,"Sequential Forward Search");
         break;
      case SBS:
         sprintf(algorithmName,"Sequential Backward Search");
         break;
      default:
         algorithmName[0] = '\0';
         break;
   }

   return algorithmName;
}

/****************************************************************************
 ** Name: showInputPattern
 ** Class: None
 ** Created:           31/05/2008
 ** Last Modification: 10/06/2008
 ** Args:
 **    Input:
 **       - None
 **    Output:
 **       - None
 ** Return:
 **    Nothing
 ** 
 ** Description: Function that show the correct pattern of the user input
 **
 ***************************************************************************/
fsVoid showInputPattern()
{
   cout << "NAME:" << endl;
   cout << "\tcortex - executes feature selection algorithms" << endl;

   cout << "SYNOPSIS:" << endl;
   cout << "\tcortex [OPTIONS] FILE" << endl;

   cout << "DESCRIPTION:" << endl;
   cout << "\tCORTEX is a system that executes a feature selection algorithm" << endl;
   cout << "\treturning for each problem four values:" << endl;
   cout << "\t\t#J    - Number of executions of the evaluator" << endl;
   cout << "\t\tJtest - Evaluation for the test data set, calculated" << endl;
   cout << "\t\t        from the best solution got with training data" << endl;
   cout << "\t\t        set." << endl;
   cout << "\t\ttime  - Time of the execution" << endl;
   cout << "\t\tXbest - Best subset of features" << endl << endl;
   cout << "\tEach value will be saved in an independent file:" << endl;
   cout << "\t\t#J    - FILE.J" << endl;
   cout << "\t\tJtest - FILE.Jtest" << endl;
   cout << "\t\ttime  - FILE.time" << endl;
   cout << "\t\tXbest - FILE.Xbest" << endl;
   cout << "\tThe mandatory arguments are:" << endl;
   cout << "\t\t--algorithm algID" << endl;
   cout << "\t\t\talgID identifies the algorithm that has to be used at" << endl;
   cout << "\t\t\tstudy of the file" << endl;
   cout << "\t\t\tThe Algorithms that can be used are:" << endl;
   cout << "\t\t\t\tSFS - Sequential Forward Search" << endl;
   cout << "\t\t\t\tSBS - Sequential Backward Search" << endl;
   cout << "\t\t\t\tGS  - Granularity Search" << endl;
   cout << "\t\t\t\tSS  - Sweep Search" << endl;
   cout << "\t\t\t\tCS  - Co-dichotomous Search" << endl;
   cout << "\t\t\t\tMS  - Mixed Search" << endl;
   cout << "\t\t-NN" << endl;
   cout << "\t\t\tNearest Neighbours evaluation" << endl;
   cout << "\t\t-NB" << endl;
   cout << "\t\t\tNaive Bayes evaluation" << endl;
   cout << "\t\t-SC" << endl;
   cout << "\t\t\tScorer evaluation" << endl;
   cout << "\t\t-EUCLIDEAN" << endl;
   cout << "\t\t\tEUCLIDEAN distance" << endl;
   cout << "\t\t-SQUAREEUCLIDEAN" << endl;
   cout << "\t\t\tSQUAREEUCLIDEAN distance" << endl;
   cout << "\t\t--RANDOMSPLIT=ON" << endl;
   cout << "\t\t\tData will be splitted randomly" << endl;
   cout << "\t\t--RANDOMSPLIT=OFF" << endl;
   cout << "\t\t\tData will not be splitted randomly" << endl;
   return;
}

/****************************************************************************
 ** Name: main
 ** Class: None
 ** Created:           31/05/2008
 ** Last Modification: 31/05/2008
 ** Args:
 **    Input:
 **       - int argc:    Number of arguments passed to the program
 **       - char **argv: String of strings that saves at each main
 **                      position, the value of the ith argument 
 **    Output:
 **       - None
 ** Return:
 **    Code of the end of the function
 ** 
 ** Description: Main function that calls the FS algorithms
 **
 ***************************************************************************/

int main(int argc, char**argv)
{
   fsAlgorithm cortexAlgorithm = NOALGORITHM;
   fsUInt fileDescriptor=0;
   fsEvaluationAlgorithm fsea=NN;
   fsDistanceAlgorithm fsda=SQUAREEUCLIDEANDISTANCE;
   fsBool rs=false;

   if(argc == 1)
   {
      showInputPattern();
      return -1;
   }

   for(int i=1;i<argc;i++)
   {
      if(strcmp(argv[i],"--algorithm") == 0)
      {
         i++;
         if(strcmp(argv[i],"SFS") == 0 || strcmp(argv[i],"sfs") == 0)
         {
            cortexAlgorithm = SFS;
         }
         else
         if(strcmp(argv[i],"SBS") == 0 || strcmp(argv[i],"sbs") == 0)
         {
            cortexAlgorithm = SBS;
         }
         else
         if(strcmp(argv[i],"GS") == 0 || strcmp(argv[i],"gs") == 0)
         {
            cortexAlgorithm = GRASEA;
         }
         else
         if(strcmp(argv[i],"SS") == 0 || strcmp(argv[i],"ss") == 0)
         {
            cortexAlgorithm = SWESEA;
         }
         else
         if(strcmp(argv[i],"CS") == 0 || strcmp(argv[i],"cs") == 0)
         {
            cortexAlgorithm = CODSEA;
         }
         else
         if(strcmp(argv[i],"MS") == 0 || strcmp(argv[i],"ms") == 0)
         {
            cortexAlgorithm = MIXSEA;
         }
      }
      else 
      if(strcmp(argv[i],"--file") == 0)
      {
         i++;
         fileDescriptor = i;
      }
      else 
      if(strcmp(argv[i],"-NN") == 0)
      {
         fsea = NN;
      }
      else 
      if(strcmp(argv[i],"-NB") == 0)
      {
         fsea = NB;
      }
      else 
      if(strcmp(argv[i],"-SC") == 0)
      {
         fsea = SC;
      }
      else 
      if(strcmp(argv[i],"-SQUAREEUCLIDEAN") == 0)
      {
         fsda = SQUAREEUCLIDEANDISTANCE;
      }
      else 
      if(strcmp(argv[i],"-EUCLIDEAN") == 0)
      {
         fsda = EUCLIDEANDISTANCE;
      }
      else 
      if(strcmp(argv[i],"--RANDOMSPLIT=ON") == 0)
      {
         rs = true;
      }
      else 
      if(strcmp(argv[i],"--RANDOMSPLIT=OFF") == 0)
      {
         rs = false;
      }
   }

   if(fileDescriptor == 0)
   {
      cout << "There's no data file defined" << endl;
      return -1;
   }

   fsPata *ss;
   fsAprs *as;

   switch(cortexAlgorithm)
   {
      case SFS:
         cout << echoAlgorithmName(cortexAlgorithm) << " in use." << endl;
         ss = new fsPata(cortexAlgorithm,fsea,fsda,0,0,rs,argv[fileDescriptor]);
         (*ss).process();
         break;
      case SBS:
         cout << echoAlgorithmName(cortexAlgorithm) << " in use." << endl;
         ss = new fsPata(cortexAlgorithm,fsea,fsda,0,0,rs,argv[fileDescriptor]);
         (*ss).process();
         break;
      case GRASEA:
         cout << echoAlgorithmName(cortexAlgorithm) << " in use." << endl;
         as = new fsAprs(cortexAlgorithm,fsea,fsda,0,0,rs,argv[fileDescriptor]);
         (*as).process();
         break;
      case SWESEA:
         cout << echoAlgorithmName(cortexAlgorithm) << " in use." << endl;
         as = new fsAprs(cortexAlgorithm,fsea,fsda,0,0,rs,argv[fileDescriptor],SFS,LASTDONE,10);
         (*as).process();
         break;
      case MIXSEA:
         cout << echoAlgorithmName(cortexAlgorithm) << " in use." << endl;
         as = new fsAprs(cortexAlgorithm,fsea,fsda,0,0,rs,argv[fileDescriptor],SFS,LASTDONE,20);
         (*as).process();
         break;
      case CODSEA:
         cout << echoAlgorithmName(cortexAlgorithm) << " in use." << endl;
         as = new fsAprs(cortexAlgorithm,fsea,fsda,0,0,rs,argv[fileDescriptor]);
         (*as).process();
         break;
      default:
         showInputPattern();
         return -1;
         break;
   }

   return EXIT_SUCCESS;
}
