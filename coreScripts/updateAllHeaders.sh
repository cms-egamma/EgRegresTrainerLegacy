#!/usr/bin/env bash

packageWithHeaders=$1
filesToUpdate=$2

for fileAndPath in `ls ${packageWithHeaders}/include/*`
do
  fileAndPath=`echo $fileAndPath | sed 's|//|/|g'`
  file=`echo $fileAndPath | awk -F "/" '{print $NF}'`
  package=`echo $fileAndPath | awk -F "/" '{print $(NF-2)}'`
  echo $fileAndPath
  echo "updating $file to $package/$file"

  ./updateHeaders.sh $filesToUpdate $file $package/$file
#  ./updateHeaders.sh $packageToUpdate/src/\* $file $package/$file
 
done
