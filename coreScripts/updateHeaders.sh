#!/usr/bin/env bash

filesToUpdate=$1
oldHeader=$2
newHeader=$3

for fileAndPath in `cat $filesToUpdate`
do
 # echo "processing $fileAndPath old: $oldHeader new: $newHeader"
  cat $fileAndPath  | sed 's|#include "'$oldHeader'"|#include "'$newHeader'"|g' >$fileAndPath.tmp
  
  mv $fileAndPath.tmp $fileAndPath

#  if [[ `md5sum $fileAndPath.tmp | awk '{print $1}'` != `md5sum $fileAndPath | awk '{print $1}'` ]];then
  #    echo "md5sums do not match"
 #     mv $fileAndPath.tmp $fileAndPath
#  else
 #     echo "md5sums do match"
 #     rm $fileAndPath.tmp
#  fi
done
