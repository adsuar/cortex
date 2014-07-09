#Load MASS library for algorithms functions
library(MASS)

# DEBUG MODE
debugging <- FALSE

back <- function()
{
   setwd("../src/R")
}

printDebug <- function(text)
{
   if(debugging) print(text)
}

loadDataTraining <- function(data_name)
{
   printDebug("LOADING COLON TRAINING DATA")
   data_training_file_name <- paste("ma_",data_name,"_trn.csv",sep="")
   print(data_training_file_name)
   file_loaded <- read.csv(data_training_file_name,header=FALSE)
   printDebug("COLON TRAINING DATA LOADED")
   file_loaded
}

loadDataTest <- function(data_name)
{
   printDebug("LOADING COLON TEST DATA")
   data_test_file_name <- paste("ma_",data_name,"_tst.csv",sep="")
   print(data_test_file_name)
   file_loaded <- read.csv(data_test_file_name,header=FALSE)
   printDebug("COLON TEST DATA LOADED")
   file_loaded
}

J <- function(X,data.trn)
{
   ###############################################################
   # Read data columns
   data.trn.data <- NULL

   # Number of features
   k <- 0
   
   # Solution of the calculation
   solution <- vector("numeric",2)
   
   # There's no estimation 
   solution[1] <- -1

   for(i in 1:length(X))
   {
      if(as.matrix(X)[i] != 0)
      {
         data.trn.data <- cbind(data.trn.data,as.matrix(data.trn[,i]))
         k <- k + 1
      }
   }

   solution[2] <- k

   if(k > 0)
   {
      # Read the class column
      data.trn.class <- as.matrix(data.trn[,length(data.trn)])

      data.qda <- qda(data.trn.data,data.trn.class)
      upredict <- predict(data.qda,data.trn.data,CV=TRUE)

      result.boolean <- (upredict$class == data.trn.class)

      accuracy <- 0

      for(i in 1:length(result.boolean))
      {
         if(result.boolean[i] == TRUE)
            accuracy <- accuracy + 1
      }

      solution[1] <- accuracy/nrow(data.trn.data)
   }
  
   # Return the solution
   solution
}

Jtest <- function(X,data.trn,data.tst)
{
   ###############################################################
   # Read data columns
   data.trn.data <- NULL
   data.tst.data <- NULL

   # Number of features
   k <- 0
   
   # Solution of the calculation
   solution <- vector("numeric",2)
   
   # There's no estimation 
   solution[1] <- -1

   for(i in 1:length(X))
   {
      if(as.matrix(X)[i] != 0)
      {
         data.trn.data <- cbind(data.trn.data,as.matrix(data.trn[,i]))
         data.tst.data <- cbind(data.tst.data,as.matrix(data.tst[,i]))
         k <- k + 1
      }
   }

   solution[2] <- k

   if(k > 0)
   {
      # Read the class column
      data.trn.class <- as.matrix(data.trn[,length(data.trn)])
      data.tst.class <- as.matrix(data.tst[,length(data.trn)])

      data.qda <- qda(data.trn.data,data.trn.class)
      upredict <- predict(data.qda,data.tst.data,CV=TRUE)

      result.boolean <- (upredict$class == data.tst.class)

      accuracy <- 0

      for(i in 1:length(result.boolean))
      {
         if(result.boolean[i] == TRUE)
            accuracy <- accuracy + 1
      }

      solution[1] <- accuracy/nrow(data.tst.data)
   }
  
   # Return the solution
   solution
}

# Function that inverts the values of a given partition
inverse <- function(partition,sizeOfZ,X)
{
   features <- length(X)

   sizeOfPartitions <- ceiling(features / sizeOfZ)
   realSizeOfPartitions <- as.double(features) / as.double(sizeOfZ)

   offset <- floor(partition*realSizeOfPartitions) + 1
   #printDebug(paste("Offset antes de cuadrar",offset))

   if((offset+sizeOfPartitions) > features)
   {
      offset <- (features - sizeOfPartitions) + 1
   }

   #printDebug(paste("El offset para la particion",partition,"es",offset))

   for(j in 1:sizeOfPartitions)
   {
      X[offset] <- 1 - X[offset]
      offset <- offset + 1
   }

   X
}

granularitySearch <- function(data_name)
{
   actualwd <- getwd()
   #Change the directory to work
   setwd("/home/adsuar/Trabajo/Proyectos/Doctorado/CORTEX/data")

   # Read training data without header
   data.trn <- loadDataTraining(data_name)
   # Read test data without header
   data.tst <- loadDataTest(data_name)

   # Number of J
   Js <- 0

   # Data loaded successfully so, right now starts the initilization
   # of the algorithm
   fsAprsFeatures = ncol(as.matrix(data.trn))-1
   printDebug(fsAprsFeatures)

   # Declare the vector that'll save the features
   # X = {X1, ..., Xn), with Xi belonging to {0,1}
   # It's initialized to 0
   X <- rep(0.0,fsAprsFeatures)

   # Estimation value
   estimation <- -1

   # Best Estimation values
   bestEstimation <- -1
   bestEstimationOld <- -1

   # Best size of solution
   bestK <- -1
   bestKOld <- -1

   # We create a new vector, Z. Its function is to save the
   # information of the segments for each iteration of the main
   # loop, where X data is splitted in ki different segments.
   # Each position of Z indicates if the correspondent partition
   # of geatures is taken as it is (1) or its negation (0)
   Z <- vector("numeric",fsAprsFeatures)

   # Create a new vector that'll save the best Z configuration
   Zbest <- Z

   # The granularity value is defined at k
   # Initially there's onlu one set
   # The use of k will be in the way k=n div (i+1), being n the
   # number of features and i the present iteration
   k <- -1

   for(i in fsAprsFeatures:2)
   {
      # Only integer divisions will be studied by now
      if(fsAprsFeatures%%i == 0)
      {
         
         k <- ceiling(fsAprsFeatures/i)

         message <- paste("Size Of Each SubSet of Z is",k,"for",i,"partitions")

         print(message)

         #print(i)
         #print(k)

         # Get a copy of features' vector
         Xaux <- X

         solution <- J(X,data.trn)

         # Number of J is incremented
         Js <- Js + 1
   
         bestEstimation <- solution[1]
         bestK <- solution[2]

         repeat
         {
            # Construct the vector Z, that'll save the result,
            # and stays as it is
            Z <- rep(1,i)

            # Zbest it's Z
            Zbest <- Z

            # We save the current solution
            bestEstimationOld <- bestEstimation
            bestKOld <- bestK

            # We do not have to force a linear relation between the
            # features so the starting point will be always randomly
            # selected

            endPosition <- ceiling(runif(1,0,100000))%%i
            startPosition <- endPosition

            # For each position of Z, we select if we need its current
            # value or the negated one
            j <- startPosition

            repeat
            {
               Z[j] <- 0

               Xaux <- inverse(j,length(Z),Xaux)

               printDebug(paste("estoy con",i,"en",j,"start:",startPosition,"end:",endPosition))

               estimation <- J(Xaux,data.trn)
               Js <- Js + 1

               if(estimation > bestEstimation)
               {
                  bestEstimation <- estimation
                  Zbest <- Z
               }
               else
               {
                  # Undo the Xaux change
                  Xaux <- inverse(j,length(Z),Xaux)

                  # Undo the Z change
                  Z[j] <- 1
               }

               # Increase the index without going out of bounds
               j <- (j+1)%%i

               if(j == endPosition) break
            } 

            if(bestEstimation == bestEstimationOld)
            {
               break
            }
         }
         
         X <- Xaux
      }
   }

   # Here we'd have the solution

   print(paste("Best solution at TRAINING:",bestEstimation))
   solution <- Jtest(X,data.trn,data.tst)

   # Number of J is incremented
   Js <- Js + 1
   
   print(paste("Best solution at TEST:",solution[1]))
   print(paste("Size of the Best solution at TEST:",solution[2]))
   print(paste("The number of executions of J has been:",Js))

   # And here we'd execute Jtest

   ###############################################################
   # Read one data column
   #colon.trn.data <- as.matrix(colon.trn[,1])
   # Read the class column
   #colon.trn.class <- as.matrix(colon.trn[,2001])

   # Read one data column
   #colon.tst.data <- as.matrix(colon.tst[,1])
   # Read the class column
   #colon.tst.class <- as.matrix(colon.tst[,2001])

   #colon.lda <- lda(colon.trn.data,colon.trn.class)
   #upredict <- predict(colon.lda,colon.tst.data,CV=TRUE)
   #upredict$class == colon.tst.class

   # Training Data Rows
   #nrow(as.matrix(colon.trn))
   # Test Data Rows
   #nrow(as.matrix(colon.tst))

   #print(upredict)

   setwd(actualwd)
}
