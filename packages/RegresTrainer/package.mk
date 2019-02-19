PKG_OBJ_DIR=$(OBJ_DIR)/RegresTrainer/src

REGRESTRAINER_LIBFILES 	= $(PKG_OBJ_DIR)/ErrorCorrection.o $(PKG_OBJ_DIR)/GBRApply.o $(PKG_OBJ_DIR)/GBREvent.o $(PKG_OBJ_DIR)/GBRForest.o $(PKG_OBJ_DIR)/GBRMaker.o $(PKG_OBJ_DIR)/GBRTrainer.o $(PKG_OBJ_DIR)/GBRTree.o $(PKG_OBJ_DIR)/HybridGBRMaker.o $(PKG_OBJ_DIR)/ParReader.o $(PKG_OBJ_DIR)/RegressionManager.o $(PKG_OBJ_DIR)/SmearingCorrection.o $(PKG_OBJ_DIR)/TMVAMaker.o $(PKG_OBJ_DIR)/TrackMomentumCorrection.o $(PKG_OBJ_DIR)/Utilities.o $(PKG_OBJ_DIR)/VariableCorrectionApply.o $(PKG_OBJ_DIR)/CruijffPdf.o $(PKG_OBJ_DIR)/CruijffPdfDict.o 
REGRESTRAINER_LIBNAME 	= $(LIB_DIR)/libRegresTrainer.so

$(REGRESTRAINER_LIBNAME):	$(REGRESTRAINER_LIBFILES)
		@$(LD) $(SOFLAGS) $(LDFLAGS) $(LIBS) $^ -o $@
		@echo "$@ done"

RegresTrainer/src/%Dict.cc: PACKAGENAME:= RegresTrainer
RegresTrainer/src/%_LinkDef.h: packages/RegresTrainer/dict/%_LinkDef.h packages/RegresTrainer/include/%.h
	@sleep 0.001s #dummy command, okay its been a long battle and this hack works... 

STD_LIBS	+= $(REGRESTRAINER_LIBNAME)
