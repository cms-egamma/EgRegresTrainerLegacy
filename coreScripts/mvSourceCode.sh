#!/usr/bin/env bash

pattern=$1
outputDir=$2

for fileAndPath in `ls $pattern`
do
 # echo $fileAndPath

  package=`echo $outputDir | sed 's|/$||g' | awk -F "/" '{print toupper($(NF-1))}'`

  file=`echo $fileAndPath | awk -F "/" '{print $NF}'`

 # echo "package "$package "file "$file
  classNameCaps=`echo $file | awk -F ".hh" '{print toupper($1)}'`

  cat $fileAndPath |  sed 's|#include "\(SH.*\).hh"|#include "SHEvent/\1.hh"|g' | \
      sed 's|#ifndef '$classNameCaps'|#ifndef '$package'_'$classNameCaps'_HH|' | \
      sed 's|#define '$classNameCaps'|#define '$package'_'$classNameCaps'_HH|' > $outputDir/$file > $outputDir/$file
      
done
