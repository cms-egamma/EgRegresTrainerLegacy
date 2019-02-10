PKG_OBJ_DIR=$(OBJ_DIR)/Utility/src

UTILITY_LIBFILES        =$(PKG_OBJ_DIR)/AnaFuncs.o $(PKG_OBJ_DIR)/CmdLineInt.o $(PKG_OBJ_DIR)/ComCodes.o $(PKG_OBJ_DIR)/ComCodesBase.o $(PKG_OBJ_DIR)/DetIdTools.o $(PKG_OBJ_DIR)/EventListComp.o $(PKG_OBJ_DIR)/EvtIndexLUT.o $(PKG_OBJ_DIR)/EvtList.o $(PKG_OBJ_DIR)/GoodLumiChecker.o $(PKG_OBJ_DIR)/JsonMaker.o $(PKG_OBJ_DIR)/Lumi3DReWeighting.o $(PKG_OBJ_DIR)/LumiReWeightingRootWrapper.o $(PKG_OBJ_DIR)/MathFuncs.o $(PKG_OBJ_DIR)/AnaFuncsDict.o $(PKG_OBJ_DIR)/DetIdToolsDict.o $(PKG_OBJ_DIR)/EventListCompDict.o $(PKG_OBJ_DIR)/EvtIndexLUTDict.o $(PKG_OBJ_DIR)/MathFuncsDict.o $(PKG_OBJ_DIR)/RootBoost.o $(PKG_OBJ_DIR)/RootBoostDict.o $(PKG_OBJ_DIR)/CaloTools.o $(PKG_OBJ_DIR)/CaloToolsDict.o $(PKG_OBJ_DIR)/TempFuncsDict.o $(PKG_OBJ_DIR)/HistFuncs.o $(PKG_OBJ_DIR)/HistFuncsDict.o   $(PKG_OBJ_DIR)/EvtLUTDict.o $(PKG_OBJ_DIR)/EcalCodes.o $(PKG_OBJ_DIR)/BitPackerDict.o $(PKG_OBJ_DIR)/BXPUInfo.o $(PKG_OBJ_DIR)/BXPUInfoDict.o $(PKG_OBJ_DIR)/MiniFloatConverter.o

UTILITY_LIBNAME 	= $(LIB_DIR)/libUtility.so

$(UTILITY_LIBNAME):	$(UTILITY_LIBFILES)
		@$(LD) $(SOFLAGS) $(LDFLAGS) $(LIBS) $^ -o $@
		@echo "$@ done"

Utility/src/%Dict.cc: PACKAGENAME:= Utility
Utility/src/%_LinkDef.h:   packages/Utility/dict/%_LinkDef.h packages/Utility/include/%.hh
#	@echo "joy and happyness $*"
	@sleep 0.001s #dummy command, okay its been a long battle and this hack works... 

STD_LIBS	+= $(UTILITY_LIBNAME)
