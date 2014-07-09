#!/bin/bash
# $1 : number of features
# $2 : data type
# $3 : data file name

HEADER_FILE_NAME=.header.ctx.tmp

# Delete the previous header
rm -rf $HEADER_FILE_NAME

# Create a new header file
touch $HEADER_FILE_NAME

# Add the values 
echo -n "@" >> $HEADER_FILE_NAME
for ((k=0;k<$1;k++))
do
   echo -n " D" >> $HEADER_FILE_NAME;
done

echo " C @" >> $HEADER_FILE_NAME

echo -n "#">>$HEADER_FILE_NAME
for ((k=0;k<$1;k++))
do
   echo -n " $2" >> $HEADER_FILE_NAME
done

echo " D #" >> $HEADER_FILE_NAME

# Modify the data file to add the result
cp $3 ".$3"
cp $HEADER_FILE_NAME $3
cat ".$3" >> $3
rm -rf $HEADER_FILE_NAME ".$3"
