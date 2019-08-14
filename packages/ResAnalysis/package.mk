PKG_OBJ_DIR=$(OBJ_DIR)/ResAnalysis/src

RESANALYSIS_LIBFILES 	= $(PKG_OBJ_DIR)/ResFitter.o $(PKG_OBJ_DIR)/ResFitterDict.o $(PKG_OBJ_DIR)/ResPlotter.o $(PKG_OBJ_DIR)/ResPlotterDict.o
RESANALYSIS_LIBNAME 	= $(LIB_DIR)/libResAnalysis.so

$(RESANALYSIS_LIBNAME):	$(RESANALYSIS_LIBFILES)
		@$(LD) $(SOFLAGS) $(LDFLAGS) $(LIBS) $^ -o $@
		@echo "$@ done"

ResAnalysis/src/%Dict.cc: PACKAGENAME:= ResAnalysis
ResAnalysis/src/%_LinkDef.h: packages/ResAnalysis/dict/%_LinkDef.h packages/ResAnalysis/include/%.hh
	@sleep 0.001s #dummy command, okay its been a long battle and this hack works... 

STD_LIBS	+= $(RESANALYSIS_LIBNAME)
