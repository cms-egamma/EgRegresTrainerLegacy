#!/usr/bin/env python


import sys
import os

import argparse

parser = argparse.ArgumentParser(description='creates the package directory structure')
parser.add_argument('--pkgName',help='name of package',required=True)
parser.add_argument('--pkgDir',help='directory containing packages',default="packages")

args = parser.parse_args()

fullPkgDir=args.pkgDir+"/"+args.pkgName

subDirs=["src","include","dict","main"]

for subDir in subDirs:
    if not os.path.exists(fullPkgDir+"/"+subDir):
        os.makedirs(fullPkgDir+"/"+subDir)

if not os.path.isfile(fullPkgDir+"/package.mk"):

    f = open(fullPkgDir+"/package.mk","w")
    f.write("PKG_OBJ_DIR=$(OBJ_DIR)/"+args.pkgName+"/src\n\n")
    f.write(args.pkgName.upper()+"_LIBFILES \t= $(PKG_OBJ_DIR)/Dummy.o\n")
    f.write(args.pkgName.upper()+"_LIBNAME \t= $(LIB_DIR)/lib"+args.pkgName+".so\n\n")
    
    f.write("$("+args.pkgName.upper()+"_LIBNAME):\t$("+args.pkgName.upper()+"_LIBFILES)\n")
    f.write("\t\t@$(LD) $(SOFLAGS) $(LDFLAGS) $(LIBS) $^ -o $@\n")
    f.write("\t\t@echo \"$@ done\"\n\n")
    f.write(args.pkgName+"/src/%Dict.cc: PACKAGENAME:= "+args.pkgName+"\n")
    f.write(args.pkgName+"/src/%_LinkDef.h: packages/"+args.pkgName+"/dict/%_LinkDef.h packages/"+args.pkgName+"/include/%.hh\n")
    f.write("\t@sleep 0.001s #dummy command, okay its been a long battle and this hack works... \n\n")
    f.write("STD_LIBS\t+= $("+args.pkgName.upper()+"_LIBNAME)\n")
    f.close()
else:
    print "package makefile "+fullPkgDir+"/package.mk exists"

    
