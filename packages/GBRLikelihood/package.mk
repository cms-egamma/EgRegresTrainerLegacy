PKG_OBJ_DIR=$(OBJ_DIR)/GBRLikelihood/src

GBRLIKELIHOOD_LIBFILES 	= $(PKG_OBJ_DIR)/GBRArrayUtils.o $(PKG_OBJ_DIR)/GBRMath.o $(PKG_OBJ_DIR)/HybridGBREvent.o $(PKG_OBJ_DIR)/HybridGBRForestD.o $(PKG_OBJ_DIR)/HybridGBRForestDDict.o $(PKG_OBJ_DIR)/HybridGBRForestFlex.o $(PKG_OBJ_DIR)/HybridGBRForestFlexDict.o $(PKG_OBJ_DIR)/HybridGBRTreeD.o $(PKG_OBJ_DIR)/HybridGBRTreeDDict.o $(PKG_OBJ_DIR)/RooDoubleCBFast.o $(PKG_OBJ_DIR)/RooDoubleCBFastDict.o  $(PKG_OBJ_DIR)/RooHybridBDTAutoPdf.o $(PKG_OBJ_DIR)/RooHybridBDTAutoPdfDict.o 
GBRLIKELIHOOD_LIBNAME 	= $(LIB_DIR)/libGBRLikelihood.so

$(GBRLIKELIHOOD_LIBNAME):	$(GBRLIKELIHOOD_LIBFILES)
		@$(LD) $(SOFLAGS) $(LDFLAGS) $(LIBS) $^ -o $@
		@echo "$@ done"

GBRLikelihood/src/%Dict.cc: PACKAGENAME:= GBRLikelihood
GBRLikelihood/src/%_LinkDef.h: packages/GBRLikelihood/dict/%_LinkDef.h packages/GBRLikelihood/include/%.h
	@sleep 0.001s #dummy command, okay its been a long battle and this hack works... 

STD_LIBS	+= $(GBRLIKELIHOOD_LIBNAME)
