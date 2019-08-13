PKG_OBJ_DIR=$(OBJ_DIR)/RegresTester/src

REGRESTESTER_LIBFILES 	= $(PKG_OBJ_DIR)/ResFitter.o $(PKG_OBJ_DIR)/RegresValidator.o
REGRESTESTER_LIBNAME 	= $(LIB_DIR)/libRegresTester.so

$(REGRESTESTER_LIBNAME):	$(REGRESTESTER_LIBFILES)
		@$(LD) $(SOFLAGS) $(LDFLAGS) $(LIBS) $^ -o $@
		@echo "$@ done"

RegresTester/src/%Dict.cc: PACKAGENAME:= RegresTester
RegresTester/src/%_LinkDef.h: packages/RegresTester/dict/%_LinkDef.h packages/RegresTester/include/%.hh
	@sleep 0.001s #dummy command, okay its been a long battle and this hack works... 

STD_LIBS	+= $(REGRESTESTER_LIBNAME)
