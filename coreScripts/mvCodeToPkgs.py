#!/usr/bin/env python

import shutil
def mvClass(className,package,oldBaseDir):

    #we want to copy the header, source and dict file, anything other than the header may not exist
    import os
    headerFile = oldBaseDir+"/include/"+className+".hh"
    dictFile = oldBaseDir+"/dict/"+className+"_LinkDef.h"
    srcFile = oldBaseDir+"/src/"+className+".cc"

    headDst = "packages/"+package+"/include"
    dictDst = "packages/"+package+"/dict"
    srcDst = "packages/"+package+"/src"
    
    if os.path.isfile(headerFile): shutil.copy(headerFile,headDst)
    if os.path.isfile(dictFile): shutil.copy(dictFile,dictDst)
    if os.path.isfile(srcFile): shutil.copy(srcFile,srcDst)   
#    print "would move",className,package,oldBaseDir

import argparse

parser = argparse.ArgumentParser(description='creates the package directory structure')
parser.add_argument('--oldCodeBaseDir',help='base directory of old code to migrate',required=True)
parser.add_argument('--fileToPkgIndex',help='text file with list of file to move to each package',required=True)
parser.add_argument('--mkPkgs',help='if true we remake the packages',default=False)

args = parser.parse_args()


validPackages="Utility SHUtility Dead Unknown Obsolete SHEvent AnaTrees Analysis Stats"

if args.mkPkgs:
    for package in validPackages.split():
        import subprocess
        subprocess.Popen(["./mkPkg.py","--pkgName",package],stdout=subprocess.PIPE)
    exit

indexFile = open(args.fileToPkgIndex,'r')

for line in indexFile:
    line = line.rstrip()
    if len(line.split())<2:
        print "error, line "+line+" doesnt have a package"
    else:
    
        package=line.split()[1]
        headerFile=line.split()[0]

        tempSplit=headerFile.split("/")
        className=tempSplit[len(tempSplit)-1].rstrip(".hh")
        
     #   print className,package

        if package in validPackages:
            mvClass(className,package,args.oldCodeBaseDir)
        else:
            print "package "+package+" for class "+className+" is not valid"


