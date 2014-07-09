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
 ** File: fsevaluation.cc
 ** Description: Source code for the fsEvaluation class
 ***************************************************************************/

#include "fsevaluation.h"

/****************************************************************************
 ** Name: fsEvaluation
 ** Class: fsEvaluation
 ** Created:           27/04/2008
 ** Last Modification: 05/06/2008
 ** Args:
 **    Input:
 **       - fsEvaluationAlgorithm fsea: Evaluation Algorithm for Data
 **       - fsDistanceAlgorithm fsda: Distance Algorithm to calculate de
 **                                   distance between two vectors
 **       - fsUInt kFCV: K-Fold Cross Validation value
 **       - fsBool r: Boolean that indicates whether Data is splitted randomly
 **                   or not
 **       - fsChar *fsFileData:  Data file
 **    Output:
 **       - None
 ** Return:
 **    Nothing
 ** 
 ** Description: Constructor of fsEvaluation class
 **
 ***************************************************************************/

fsEvaluation::fsEvaluation(
                 fsEvaluationAlgorithm fsea,
                 fsDistanceAlgorithm fsda,
                 fsUInt kFCV,
                 fsUInt kFCVTest,
                 fsBool r,
                 fsChar *fsFileData
                 )
{
#ifdef DEBUG_fsEvaluation
   cout << "Processing to create fsEvaluation" << endl;
#endif
   // K-Fold must be greater than 1 or equal to 0
   assert(kFCV==0 || kFCV>1);
   
   // We must assure that, if distance is initialized to a value different
   // from UNKNOWN, that's because the evaluation algorithm is NN (Nearest
   // Neighbours)
   assert((fsda!=UNKNOWN && fsea==NN) || fsda==UNKNOWN);
   
   // We assign the evaluation and distance algorithms
   fsEA = fsea;
   fsDA = fsda;

   // We initialize distance function if and only if it's different from
   // UNKNOWN
   if(fsDA != UNKNOWN)
      fsDist = new fsDistance(fsDA);
   else
      fsDist = NULL;

   fsChar fileName[1024];

#ifdef DEBUG_fsEvaluation
   cout << "\tProcessing Training Data for Evaluation" << endl;
#endif
   // We create the Training Data Structures
   // We must now get the training data
   sprintf(fileName,"%s.trn.ctx",fsFileData);
   study = new fsData(fileName);
   
   // Data must be loaded
   assert((*study).successfullyLoaded());

   // Normalize Study Data
   (*study).normalizeFeatures();
   
   // When kFCV it's 0, it means LOOCV
   if(kFCV == 0)
      kFCV = (*((*study).getFeatures())).numRows();

   // The number of features must be greater or equal than kFoldCrossValidation
   assert((*((*study).getFeatures())).numRows() >= kFCV);
   
   // If it's greater or equal, then kFCV is a valid value
   kFoldCrossValidation = kFCV;
   
   // We calculate the size of the blocks
   kFoldSteps = (fsUInt)
                (kFoldCrossValidation/((*((*study).getFeatures())).numRows()));
   
   // It's possible that the remainder of the division won't be 0. In this case,
   // Last block will be smaller than the rest
   
   
   // We get the randomized flag
   randomized = r;
   
   // If randomized is true, blocks will start at any position
   // On the other hand, blocks will start at 0 position
   if(randomized)
   {
      srand(time(NULL));
      firstBlock = rand()%((*((*study).getFeatures())).numRows());
   }
   else
      firstBlock = 0;

#ifdef DEBUG_fsEvaluation
   cout << "\tProcessing Test Data for Evaluation" << endl;
#endif
   // We create the Test Data Structures
   // We must now get the training data
   sprintf(fileName,"%s.tst.ctx",fsFileData);
   studyTest = new fsData(fileName);
   
   // Data must be loaded
   assert((*studyTest).successfullyLoaded());

   // Normalize Study Test Data
   (*studyTest).normalizeFeatures();
   
   // When kFCV it's 0, it means LOOCV
   if(kFCVTest == 0)
      kFCVTest = (*((*studyTest).getFeatures())).numRows();

   // The number of features must be greater or equal than kFoldCrossValidation
   assert((*((*studyTest).getFeatures())).numRows() >= kFCVTest);
   
   // If it's greater or equal, then kFCV is a valid value
   kFoldCrossValidationTest = kFCVTest;
   
   // We calculate the size of the blocks
   kFoldStepsTest = (fsUInt)
                (kFoldCrossValidationTest/((*((*studyTest).getFeatures())).numRows()));
   
   // It's possible that the remainder of the division won't be 0. In this case,
   // Last block will be smaller than the rest
   
   
   // If randomized is true, blocks will start at any position
   // On the other hand, blocks will start at 0 position
   if(randomized)
   {
      srand(time(NULL));
      firstBlockTest = rand()%((*((*studyTest).getFeatures())).numRows());
   }
   else
      firstBlockTest = 0;

   // We initialize the evaluation function
   switch(fsEA)
   {
      case NN:
         eval = &fsEvaluation::NearestNeighbours;
         eval_set = &fsEvaluation::NearestNeighbours;
         break;
      case SC:
         eval = &fsEvaluation::Scorer;
         break;
      case NB:
         eval = &fsEvaluation::NaiveBayes;
#ifdef DEBUG_fsEvaluation
   cout << "Developing Training Data" << endl;
#endif
         setActive(TRAINING);
         initializeNaiveBayes();
#ifdef DEBUG_fsEvaluation
   cout << "Developing Test Data" << endl;
#endif
         setActive(TEST);
         initializeNaiveBayes();
         break;
      case UK:
         eval = &fsEvaluation::NearestNeighbours;
         eval_set = &fsEvaluation::NearestNeighbours;
         break;
   }
   
   // We set by default data training as the data to be used at the studies
   setActive(TRAINING);
   
#ifdef DEBUG_fsEvaluation
   cout << "fsEvaluation created" << endl;
#endif
}

/****************************************************************************
 ** Name: ~fsEvaluation
 ** Class: fsEvaluation
 ** Created:           27/04/2008
 ** Last Modification: 05/06/2008
 ** Args:
 **    Input:
 **       - None
 **    Output:
 **       - None
 ** Return:
 **    Nothing
 ** 
 ** Description: Destructor of fsEvaluation class
 **
 ***************************************************************************/
fsEvaluation::~fsEvaluation()
{
   eval = NULL;
   eval_set = NULL;
}

/****************************************************************************
 ** Name: initializeNaiveBayes
 ** Class: fsEvaluation
 ** Created:           27/04/2008
 ** Last Modification: 27/04/2008
 ** Args:
 **    Input:
 **       - None
 **    Output:
 **       - None
 ** Return:
 **    - None
 ** 
 ** Description: Method that initializes that Naive Bayes weight values. With
 **              this matriz, when we need to make a calculation, we access
 **              directly to this matrix to get data
 **
 ***************************************************************************/

fsVoid fsEvaluation::initializeNaiveBayes()
{
#ifdef DEBUG_fsEvaluation
   cout << endl << "BAYES Method" << endl;
   cout << "Creation and instantiation of bayes preproccess matrix" << endl;
#endif
       
    try
    {
       if(!((*studyActive).successfullyLoaded()))
       {
          throw "initializeNaiveBayes: Data matriz hasn't been instantiated yet";
       }
       
       if((*studyActive).getTypeFeatures() != DISCREET ||
          (*studyActive).getTypeFeatures() != BINARY )
       {
#ifdef DEBUG_fsEvaluation
   cout << "Data type is: " << (*studyActive).getTypeFeatures() << endl;
   cout << "Remember that types are: " << endl;
   cout << "        NOTDEFINED = -1" << endl;
   cout << "        DISCREET   = 0" << endl;
   cout << "        CONTINUOUS = 1" << endl;
   cout << "        BINARY     = 2" << endl;
   cout << "        COMPOUND   = 3" << endl;
#endif
          throw "initializeNaiveBayes:We can't use Naive Bayes 'cause features values are continuous";
       }
       

       if((*studyActive).getTypeClasses() != DISCREET ||
          (*studyActive).getTypeClasses() != BINARY )
       {
#ifdef DEBUG_fsEvaluation
   cout << "Data type is: " << (*studyActive).getTypeClasses() << endl;
   cout << "Remember that types are: " << endl;
   cout << "        NOTDEFINED = -1" << endl;
   cout << "        DISCREET   = 0" << endl;
   cout << "        CONTINUOUS = 1" << endl;
   cout << "        BINARY     = 2" << endl;
   cout << "        COMPOUND   = 3" << endl;
#endif
          throw "initializeNaiveBayes:We can't use Naive Bayes 'cause classes values are continuous";
       }


       if((*(*studyActive).getClasses()).numCols() != 1)
       {
#ifdef DEBUG_fsEvaluation
   cout << "Classes has: " << (*(*studyActive).getClasses()).numCols();
   cout << " dimensions." << endl;
#endif
          throw "initializeNaiveBayes:Classes have more than one dimension";
       }
       
       cout << "Begining of the Bayes preproccess matrix instantiation " << endl;
       
       // I initialize the items votation
       // With this block, I discover the number of different classes, so
       // I can optimize the use of memory
       for(fsUInt i=0;i<(*studyActive).getInstances();i++)
       {
          //(*(*(*study).getClasses()).extractRowSTL(i))[0];
          (*itemsActive).add((*(*(*studyActive).getClasses()).extractRowSTL(i))[0]);
       }
       
       fsUInt itemsClasses=(*itemsActive).getClasses();
       
       cout << "Creation of the weight matrix for Naive Bayes" << endl;
       
       (*matrixNaiveBayesActive) = vector< vector < vector <float> > > (itemsClasses);
       
       fsUInt features = (*studyActive).getNumberOfFeatures();
       
       for(fsUInt i=0;i<itemsClasses;i++)
       {
          (*matrixNaiveBayesActive)[i] = vector < vector < float > > (features);
          
          for(fsUInt j=0;j<features;j++)
             (*matrixNaiveBayesActive)[i][j] = vector<float> (kMAXVALDISC,0);
       }
       
       fsReal count_iterations = 0.0;
       
       for(fsUInt i=0;i<(*studyActive).getInstances();i++)
       {
//          for(int j=0;j<((*(*data).getData()).numCols());j++)
          for(fsUInt j=0;j<(*studyActive).getNumberOfFeatures();j++)
          {
             count_iterations++;
             
             fsUInt val = (fsUInt)((*(*(*studyActive).getClasses()).extractRowSTL(i))[0]);
               
             (*matrixNaiveBayesActive)[val][j][(fsUInt)((*(*(*studyActive).getFeatures()).extractRowSTL(i))[j])]++;
          }
       }

       cout << "End of the Bayes preproccess matrix instantiation " << endl;
       
#ifdef DEBUG_fsEvaluation
       cout << "The number of iterations has been " << count_iterations << endl;
       cout << "It should've been " << (*studyActive).getInstances() * (*studyActive).getNumberOfFeatures() << endl;
#endif
    }
    catch(const char *s)
    {
       cout << s << endl;
       exit(-1);
    }
}

/****************************************************************************
 ** Name: J
 ** Class: fsEvaluation
 ** Created:           28/04/2008
 ** Last Modification: 02/05/2008
 ** Args:
 **    Input:
 **       - fsVector *X: Pointer to the Vector with features selected
 **    Output:
 **       - None
 ** Return:
 **    - fsReal: evaluation value
 ** 
 ** Description: Method that implements the call to the generic evaluation
 **              function. This method will made as many calls to the eval
 **              function as defined by kFoldCrossValidation
 **
 ***************************************************************************/
fsReal fsEvaluation::J(fsVector *X)
{
   // If data hasn't been loaded successfully, we return the minimum value
   if(!(*studyActive).successfullyLoaded())
      return 0.0;
   
   // If data has been loaded successfully, we process to execute the
   // evaluation funcion as many times as blocks made 'cause of 
   // kFoldCrossValidation
   
   //fsUInt nextBlock = firstBlock;
   
   fsUInt fold      = 0;
   
   fsReal value,average=0.0;

   
   for(fold = 0; fold < kFoldCrossValidationActive; fold++)
   {
      value = (this->*eval)(X,fold);
      
      //cout << " " << fold << "(" << value << ") ";
      average += value;
   }

   //cout << endl;   

   return (average/kFoldCrossValidationActive);
}

/****************************************************************************
 ** Name: J
 ** Class: fsEvaluation
 ** Created:           05/06/2008
 ** Last Modification: 05/06/2008
 ** Args:
 **    Input:
 **       - fsSet *X: Pointer to the set  with features selected
 **    Output:
 **       - None
 ** Return:
 **    - fsReal: evaluation value
 ** 
 ** Description: Method that implements the call to the generic evaluation
 **              function. This method will made as many calls to the eval
 **              function as defined by kFoldCrossValidation
 **
 ***************************************************************************/
fsReal fsEvaluation::J(fsSet *X)
{
   // Pointer to evaluation function has to be set
   assert(eval_set != NULL);

   // If data hasn't been loaded successfully, we return the minimum value
   if(!(*studyActive).successfullyLoaded())
      return 0.0;
   
   // If data has been loaded successfully, we process to execute the
   // evaluation funcion as many times as blocks made 'cause of 
   // kFoldCrossValidation
   
   //fsUInt nextBlock = firstBlock;
   
   fsUInt fold      = 0;
   
   fsReal value,average=0.0;

   
   for(fold = 0; fold < kFoldCrossValidationActive; fold++)
   {
      value = (this->*eval_set)(X,fold);
      
      //cout << " " << fold << "(" << value << ") ";
      average += value;
   }

   //cout << endl;   

   return (average/kFoldCrossValidationActive);
}

/****************************************************************************
 ** Name: NearestNeighbours
 ** Class: fsEvaluation
 ** Created:           28/04/2008
 ** Last Modification: 01/05/2008
 ** Args:
 **    Input:
 **       - fsVector *X: Pointer to the Vector with features selected
 **       - fsUInt fold: Number of training block
 **    Output:
 **       - None
 ** Return:
 **    - fsReal: evaluation value
 ** 
 ** Description: Method that implements the call to the evaluation function 
 **              with the Nearest Neighbours algorithm
 **
 ***************************************************************************/
fsReal fsEvaluation::NearestNeighbours(fsVector *X, fsUInt fold)
{
   // Number of Observations
   fsUInt maxObservations = (*((*studyActive).getFeatures())).numRows();
   
   // Start Position of Training
   fsUInt startPositionTest = 
          (firstBlockActive + fold*kFoldStepsActive) % maxObservations;
   
   // Start Position of Test
   fsUInt startPositionTraining = 
          (startPositionTest + kFoldStepsActive) % maxObservations;

#ifdef DEBUG_fsEvaluation
   cout << "Execution of NearestNeighbours distance algorithm" << endl;
   cout << "This is the " << fold << "th fold" << endl;
   cout << "The Training Position is " << startPositionTraining << endl;
   cout << "The Test Position is " << startPositionTest << endl;
   cout << "Size Of Training " << maxObservations-1 << endl;
   cout << "Size Of Test " << kFoldStepsActive << endl;
   cout << "X=" << *X;
#endif

   // minimum distance
   fsReal min=MAX_FSREAL;
   // position of the mininmum distance
   fsUInt min_pos;
   
   // Number of correct associations according to minimum distance
   fsReal correct=0.0;
   
   // Iterators for test and training
   fsUInt i,j;
   fsUInt position_test,position_training;
   
   // Distance between two vectors
   fsReal distance;
   
   // Matrices of features and classes
   fsMatrix *features,*classes;
   
   features = (*studyActive).getFeatures();
   classes  = (*studyActive).getClasses();
   
   // We create the loop of study
   // Outside Loop is the one that browses data test
   for(i=0,position_test=startPositionTest;
       i<kFoldStepsActive;
       i++,position_test=(position_test+1)%maxObservations)
   {
      // Inside Loop is the one that browses data training
      for(j=0,position_training=startPositionTraining;
          j<(maxObservations-1);
          j++,position_training=(position_training+1)%maxObservations)
      {
         //cout << *features;
         //cout << (*X);
         // We get the distance between two vectors
         distance = (*fsDist).get((*features).extractRowSTL(position_test),
                                  (*features).extractRowSTL(position_training),
                                  (*X).getSTL());
#ifdef DEBUG_fsEvaluation
         cout << "distance(" << position_test << "," << position_training;
         cout << ")= " << distance << endl;
#endif
         // If it's lower than minimum distance, this last becomes the new
         // minimum distance
         if(distance < min )
         {
            min = distance;
            min_pos = position_training;
         }
      }
      
      //cout << (*((*classes).extractRowSTL(min_pos)))[0] << endl;
      //cout << (*((*classes).extractRowSTL(position_test)))[0] << endl;
      //exit(-1);

      // If both vectors (minimum and test) belong to the same class, then we
      // count this iteration as correct
      if((*((*classes).extractRowSTL(min_pos))) == 
         (*((*classes).extractRowSTL(position_test))))
         correct++;
   }

#ifdef DEBUG_fsEvaluation
   cout << (correct/kFoldStepsActive)*100 << "% accuracy" << endl;
#endif
   // We return the average of accuracy of the algorithm
   return correct/kFoldStepsActive;
}

/****************************************************************************
 ** Name: NearestNeighbours
 ** Class: fsEvaluation
 ** Created:           05/06/2008
 ** Last Modification: 05/06/2008
 ** Args:
 **    Input:
 **       - fsSet    *X: Pointer to the Set with features selected
 **       - fsUInt fold: Number of training block
 **    Output:
 **       - None
 ** Return:
 **    - fsReal: evaluation value
 ** 
 ** Description: Method that implements the call to the evaluation function 
 **              with the Nearest Neighbours algorithm
 **
 ***************************************************************************/
fsReal fsEvaluation::NearestNeighbours(fsSet *X, fsUInt fold)
{
   // Number of Observations
   fsUInt maxObservations = (*((*studyActive).getFeatures())).numRows();
   
   // Start Position of Training
   fsUInt startPositionTest = 
          (firstBlockActive + fold*kFoldStepsActive) % maxObservations;
   
   // Start Position of Test
   fsUInt startPositionTraining = 
          (startPositionTest + kFoldStepsActive) % maxObservations;

#ifdef DEBUG_fsEvaluation
   cout << "Execution of NearestNeighbours distance algorithm" << endl;
   cout << "This is the " << fold << "th fold" << endl;
   cout << "The Training Position is " << startPositionTraining << endl;
   cout << "The Test Position is " << startPositionTest << endl;
   cout << "Size Of Training " << maxObservations-1 << endl;
   cout << "Size Of Test " << kFoldStepsActive << endl;
   cout << "X=" << *X;
#endif

   // minimum distance
   fsReal min=MAX_FSREAL;
   // position of the mininmum distance
   fsUInt min_pos;
   
   // Number of correct associations according to minimum distance
   fsReal correct=0.0;
   
   // Iterators for test and training
   fsUInt i,j;
   fsUInt position_test,position_training;
   
   // Distance between two vectors
   fsReal distance;
   
   // Matrices of features and classes
   fsMatrix *features,*classes;
   
   features = (*studyActive).getFeatures();
   classes  = (*studyActive).getClasses();
   
   // We create the loop of study
   // Outside Loop is the one that browses data test
   for(i=0,position_test=startPositionTest;
       i<kFoldStepsActive;
       i++,position_test=(position_test+1)%maxObservations)
   {
      // Inside Loop is the one that browses data training
      for(j=0,position_training=startPositionTraining;
          j<(maxObservations-1);
          j++,position_training=(position_training+1)%maxObservations)
      {
         //cout << *features;
         //cout << (*X);
         // We get the distance between two vectors
         distance = (*fsDist).get((*features).extractRowSTL(position_test),
                                  (*features).extractRowSTL(position_training),
                                  X);
#ifdef DEBUG_fsEvaluation
         cout << "distance(" << position_test << "," << position_training;
         cout << ")= " << distance << endl;
#endif
         // If it's lower than minimum distance, this last becomes the new
         // minimum distance
         if(distance < min )
         {
            min = distance;
            min_pos = position_training;
         }
      }
      
      //cout << (*((*classes).extractRowSTL(min_pos)))[0] << endl;
      //cout << (*((*classes).extractRowSTL(position_test)))[0] << endl;
      //exit(-1);

      // If both vectors (minimum and test) belong to the same class, then we
      // count this iteration as correct
      if((*((*classes).extractRowSTL(min_pos))) == 
         (*((*classes).extractRowSTL(position_test))))
         correct++;
   }

#ifdef DEBUG_fsEvaluation
   cout << (correct/kFoldStepsActive)*100 << "% accuracy" << endl;
#endif
   // We return the average of accuracy of the algorithm
   return correct/kFoldStepsActive;
}

/****************************************************************************
 ** Name: Scorer
 ** Class: fsEvaluation
 ** Created:           28/04/2008
 ** Last Modification: 31/05/2008
 ** Args:
 **    Input:
 **       - fsVector *X: Pointer to the Vector with features selected
 **       - fsUInt fold: Number of training block
 **    Output:
 **       - None
 ** Return:
 **    - fsReal: evaluation value
 ** 
 ** Description: Method that implements the call to the evaluation function 
 **              with the Scorer algorithm
 **
 ***************************************************************************/
fsReal fsEvaluation::Scorer(fsVector *X, fsUInt fold)
{
#ifdef DEBUG_fsEvaluation
   cout << "Execution of Scorer distance algorithm" << endl;
#endif

   return 0.0;
}

/****************************************************************************
 ** Name: NaiveBayes
 ** Class: fsEvaluation
 ** Created:           28/04/2008
 ** Last Modification: 31/05/2008
 ** Args:
 **    Input:
 **       - fsVector *X: Pointer to the Vector with features selected
 **       - fsUInt fold: Number of training block
 **    Output:
 **       - None
 ** Return:
 **    - fsReal: evaluation value
 ** 
 ** Description: Method that implements the call to the evaluation function 
 **              with the Naive Bayes algorithm
 **
 ***************************************************************************/
fsReal fsEvaluation::NaiveBayes(fsVector *X, fsUInt fold)
{
#ifdef DEBUG_fsEvaluation
   cout << "Execution of NaiveBayes distance algorithm" << endl;
#endif
   //cdebug_1("Ejecutamos NaiveBayes en modo normal");
   fsReal correct=0;
   
   //cdebug_2("El numero de filas de datos es ",dataRows);
   fsUInt ndata = (*studyActive).getInstances();
   fsUInt x_size,z_size;
   
   fsVector x,y,z;
   
   x = (*itemsActive).getSet();
   x_size = x.sizeOf();
   
   y = (*itemsActive).getVotes();
   
   z = (*X);
   z_size = (*X).sizeOf();

   // Variable that saves the best class
   fsReal best_class;
   
   for(fsUInt j=0;j<ndata;j++)
   {
      float best_probability=0;

      vector<float> v = (*(*(*studyActive).getFeatures()).extractRowSTL(j));
         
      for(fsUInt i=0;i<x_size;i++)
      {
         // We read the class
         fsReal itemClass = x[i];
         
         // We get the votes of a class
         fsReal v_itemClass = y[i];

         fsReal prod_prob = (fsReal)(v_itemClass-1)/(fsReal)(ndata-1);
         
         for(fsUInt k=0;k<z_size;k++)
         {
            if((*X)[k] == 1)
            {
               fsUInt data = (fsUInt)v[k];

               prod_prob *= ((fsReal)(*matrixNaiveBayesActive)[(fsUInt)itemClass][k][data]-1)/(v_itemClass-1);
            }
         }

         if(prod_prob>best_probability)
         {
            best_class = itemClass;
            best_probability = prod_prob;
         }
      }
      
      // If the best class fits in (matches) with the real one
      if(best_class == ((*(*(*studyActive).getClasses()).extractRowSTL(j))[0]))
         correct++;
   }

   return correct/ndata;
}


/****************************************************************************
 ** Name: numberOfFeatures
 ** Class: fsEvaluation
 ** Created:           07/05/2008
 ** Last Modification: 07/05/2008
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
fsUInt fsEvaluation::numberOfFeatures()
{
   // Data has to be successfully loaded
   assert((*study).successfullyLoaded());

   return (*study).getNumberOfFeatures();
}

/****************************************************************************
 ** Name: numberOfObservations
 ** Class: fsEvaluation
 ** Created:           07/05/2008
 ** Last Modification: 07/05/2008
 ** Args:
 **    Input:
 **       - None
 **    Output:
 **       - None
 ** Return:
 **    - fsUInt: Number of Observations
 ** 
 ** Description: Method that returns the number of observations of the
 **              problem
 **
 ***************************************************************************/
fsUInt fsEvaluation::numberOfObservations()
{
   // Data has to be successfully loaded
   assert((*studyActive).successfullyLoaded());

   return (*((*studyActive).getFeatures())).numRows();
}

/****************************************************************************
 ** Name: setActive
 ** Class: fsEvaluation
 ** Created:           07/05/2008
 ** Last Modification: 07/05/2008
 ** Args:
 **    Input:
 **       - fsDataType fsda: Data to be set as active
 **    Output:
 **       - None
 ** Return:
 **    - None
 ** 
 ** Description: Method that sets as the object of calculations TRAINING or
 **              TEST DATA
 **
 ***************************************************************************/
fsVoid fsEvaluation::setActive(fsDataType fsda)
{
   switch(fsda)
   {
      case TEST:
         studyActive                = studyTest;
         itemsActive                = &itemsTest;
         matrixNaiveBayesActive     = &matrixNaiveBayesTest;
         kFoldCrossValidationActive = kFoldCrossValidationTest;
         kFoldStepsActive           = kFoldStepsTest;
         firstBlockActive           = firstBlockTest;

         break;

      case TRAINING:
      default:
         studyActive                = study;
         itemsActive                = &items;
         matrixNaiveBayesActive     = &matrixNaiveBayes;
         kFoldCrossValidationActive = kFoldCrossValidation;
         kFoldStepsActive           = kFoldSteps;
         firstBlockActive           = firstBlock;
 
         break;
   }

   fsDataActive = fsda;
   
   return;
}
