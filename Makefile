#Sam's top level Makefile
#libaries are defined in seperate locations

ARCH=linux


#test to see if SCRAM_ARCH exists and if so use that as arch otherwise default to linux
ifeq ($(strip $(SCRAM_ARCH)),)
BFARCH		= Linux
else
BFARCH		= $(SCRAM_ARCH)
endif

#the input directories containing the source
SRC_DIR = main
LIBSRC_DIR = src
INCLUDE_DIR =include
DICT_DIR = dict

SOURCES := $(wildcard $(SRC_DIR)/*.cc $(LIBSRC_DIR)/*.cc )
#SOURCES := $(wildcard $(LIBSRC_DIR)/*.cc) 

#the output directories
LIB_DIR = libs/$(BFARCH)
OBJ_DIR = tmp/$(BFARCH)
BIN_DIR = bin/$(BFARCH)

#letting my dict make file know where to find source
export SRC_DIR
export LIBSRC_DIR
export INCLUDE_DIR

#making sure my output directories exist
$(shell test -d $(INCLUDE_DIR) || mkdir -p $(INCLUDE_DIR))
$(shell test -d $(BIN_DIR) || mkdir -p $(BIN_DIR))
$(shell test -d $(LIB_DIR) || mkdir -p $(LIB_DIR))
$(shell test -d $(OBJ_DIR) || mkdir -p $(OBJ_DIR))


PKG_DIR=packages
PKGSTOBUILD= $(shell ls $(PKG_DIR)/* -d  | sed 's|packages/||g' | sed 's|\(^.*\)|$(LIB_DIR)/lib\1.so|'  )
PKG_SRC_DIRS=$(shell ls $(PKG_DIR)/* -d  | sed 's|packages/||g' | sed 's|\(^.*\)|$(OBJ_DIR)/\1/src|' )
$(shell mkdir -p $(PKG_SRC_DIRS))
PKG_MAIN_DIRS=$(shell ls $(PKG_DIR)/* -d  | sed 's|packages/||g' | sed 's|\(^.*\)|$(OBJ_DIR)/\1/main|' )
$(shell mkdir -p $(PKG_MAIN_DIRS))
PACKAGES=$(shell ls $(PKG_DIR)/* -d  | sed 's|packages/||g')

#$(shell rm include/*)
$(foreach package,$(PACKAGES),$(shell test -h $(INCLUDE_DIR)/$(package) || ln -s $(PWD)/$(PKG_DIR)/$(package)/include $(INCLUDE_DIR)/$(package) );)

DEPEND_DIR 	= $(OBJ_DIR)
#DEPENDFILES 	:= $(patsubst $(SRC_DIR)/%.cc,$(DEPEND_DIR)/%.P,$(wildcard $(SRC_DIR)/*.cc ))
#DEPENDFILES 	+= $(patsubst $(LIBSRC_DIR)/%.cc,$(DEPEND_DIR)/%.P,$(wildcard $(LIBSRC_DIR)/*.cc ))
#DEPENDFILES 	+= $(patsubst $(DICT_DIR)/%.cc,$(DEPEND_DIR)/%.P,$(wildcard $(DICT_DIR)/*.cc ))

#DEPENDFILES     += $(foreach package,$(PACKAGES),$(shell ls $(OBJ_DIR)/$(package)/src/*.P))
DEPENDFILES     := $(shell find $(OBJ_DIR) -name "*.P")
#$(foreach package,$(PACKAGES),DEPIN$(shell echo "$package" );)

#telling make to look in these directories for what it needs
vpath  %Dict.cc $(DICT_DIR)
vpath  %Dict.h $(DICT_DIR)
vpath  %LinkDef.h $(PKG_DIR)
vpath  %.so $(LIB_DIR)
vpath  %.o  $(OBJ_DIR)
vpath  %.cc $(SRC_DIR):$(LIBSRC_DIR)
vpath  %.hh $(INCLUDE_DIR)


ifeq ($(ARCH),linux)
# Linux with egcs, gcc 2.9x, gcc 3.x (>= RedHat 5.2)
# -02 option heavily optimises the code but makes it difficult to debug
# -g prints out extra debuging info
# -W does extra but not all warnings
CXX           = g++
#CXXFLAGS      = -O2  -fPIC -Wall
ifeq ($(DEBUG),1)
CXXFLAGS      = -fPIC -Wall -g -O  -Werror  -Wfatal-errors
else
CXXFLAGS      = -fPIC -Wall -O3 -fopenmp -Werror -Wfatal-errors

endif
LD            = g++
#LDFLAGS       = -O2 
LDFLAGS       = ${CXXFLAGS}
SOFLAGS       = -shared
endif


ROOTCFLAGS   := $(shell root-config --cflags) -I$(ROOFITSYS)/include
ROOTCFLAGS   += -I$(ROOFITSYS)/include
#ROOTLIBS     := $(shell root-config --libs) -lDCache -L/raid/expt-sw/SL4/cms/slc4_ia32_gcc345/external/dcap/1.2.35-CMS3/lib/ -ldcap
ROOTLIBS     := $(shell root-config --libs) -lDCache 
ROOTLIBS     += -L$(ROOFITSYS)/lib -lRooFitCore -lRooFit -lRooStats -lFoam -lMinuit -lTMVA


#we dont always need to link this library in, only if we are going to run
#on files in the dcache pool
#so if the lib location isnt specified or doesnt exist, then its okay
ifeq ($(shell test -d $(DCAP_LIB_DIR) && echo true),true)
ROOTLIBS     += -L$(DCAP_LIB_DIR) -ldcap 
endif



ROOTGLIBS    := $(shell root-config --glibs)

CMSSWFLAGS    = -I$(CMSSW_BASE)/src -I$(CMSSW_RELEASE_BASE)/src
CMSSWLIBS     = -L${CMSSW_BASE}/lib/${SCRAM_ARCH} -L${CMSSW_RELEASE_BASE}/lib/${SCRAM_ARCH} -lCondFormatsEgammaObjects

#we dont know exactly where the boost directory is for CMSSW, we just pick the 
#the first one in the directory, dont think this will matter
BOOST_DIR     = $(shell ls $$CMSSW_DATA_PATH/../external/boost/* -d  | head -n 1)

CXXFLAGS     += $(ROOTCFLAGS) -I$(INCLUDE_DIR) $(CMSSWFLAGS)  -fexceptions  -I$(BOOST_DIR)/include -I/cvmfs/cms.cern.ch/slc6_amd64_gcc700/cms/vdt/0.4.0/include/



LIBS          = $(ROOTLIBS) $(SYSLIBS) $(USERLIBS) $(CMSSWLIBS)
GLIBS         = $(ROOTGLIBS) $(SYSLIBS)
#LDFLAGS      += $(shell root-config --nonew --ldflags)
LDFLAGS      += $(shell root-config --nonew ) -fno-exceptions

#if the file doesnt exist, we will ignore it, and often the file may not exist

#naughty, remove this asap
#CXXFLAGS	+= -I/home/sharper/lhapdfInstall/include 
#LIBS 		+= -L/home/sharper/lhapdfInstall/lib -lLHAPDF



-include $(shell ls $(PKG_DIR)/* -d | sed 's|\(.*$ \)|\1/package.mk|')


#default: 	libs

it.so:		all

libs:           $(STD_LIBS)


clean:		clean.tmp clean.dict
		rm -rf $(BIN_DIR) $(LIB_DIR) $(OBJ_DIR)


clean.tmp: 
		rm -f *~ $(SRC_DIR)/*~ $(LIBSRC_DIR)/*~ $(INCLUDE_DIR)/*~ $(DICT_DIR)/*~

clean.dict:	
		rm -rf $(DICT_DIR)/*

#$(LIB_DIR)/lib%.so:
#
#cd $(PKG_DIR)/$* && make
#	@echo $*

$(OBJ_DIR)/%Dict.o:%Dict.cc
#right I want make to rerun the dictionary defination if it needs to
#this next bit solely sorts out the dependences
#yes I know it looks like somebody mashed the keyboard 

#after much puzzling over this, I figured out what Past Sam was trying to do
#make the rule as normal for the Dict.o and then remake the rule with Dict.cc so Make sees that it doesnt need to be updated

#warning, it wont regenerate after a link def change and after a manual change of Dict.h (which shouldnt be done), can fix but for another time

#	@echo "dict comp started $@"
	@$(CXX) $(CXXFLAGS) -MM -o $(DEPEND_DIR)/$*Dict.d -c $(DICT_DIR)/$*Dict.cc; \
#	sed 's,$*Dict.o,$(OBJ_DIR)/$*Dict.o,' < $(DEPEND_DIR)/$*Dict.d  > $(DEPEND_DIR)/$*Dict.P; \
	sed -e '1s,^\(.*\):\(.*\),$(OBJ_DIR)/$*Dict.o:\2,' -e 's,dict/,,g' -e 's,$*Dict.h,,g' < $(DEPEND_DIR)/$*Dict.d > $(DEPEND_DIR)/$*Dict.P; \
	grep -v "dict/$*Dict" < $(DEPEND_DIR)/$*Dict.d | sed -e 's/#.*//' -e 's/^[^:]*: *//' -e 's/ *\\$$//' -e '/^$$/ d' -e 's/$$/ :/'  -e 's,dict/,,g' >> $(DEPEND_DIR)/$*Dict.P; \
	sed -e 's,dict/$*Dict.cc,,' -e '1s,^\(.*\):\(.*\),dict/$*Dict.cc:\2,' -e 's,dict/,,g' -e 's,$*Dict.h,,g'  < $(DEPEND_DIR)/$*Dict.d >> $(DEPEND_DIR)/$*Dict.P; \
	rm -f $(DEPEND_DIR)/$*Dict.d
#yay, dependences done
	@$(CXX) $(CXXFLAGS) -I$(DICT_DIR) -I./ -c $(DICT_DIR)/$*Dict.cc -o $@
#now clean up dict files (otherwise they wont be regened properly)
#	@rm $(DICT_DIR)/$*Dict.*

$(OBJ_DIR)/%.o: packages/%.cc
#this bit sorts out all the dependences
#	@echo "packages $*"
#so we had some fun when moving to packages, I ended up just deleting the first part of the gcc output and replacing it with the rule name it should have
	@$(CXX) $(CXXFLAGS) -MM -o $(DEPEND_DIR)/$*.d -c $<; \
	sed '1s,^\(.*\):\(.*\),$(OBJ_DIR)/$*.o:\2,' < $(DEPEND_DIR)/$*.d > $(DEPEND_DIR)/$*.P; \
	sed -e 's/#.*//' -e 's/^[^:]*: *//' -e 's/ *\\$$//' -e '/^$$/ d' -e 's/$$/ :/'  < $(DEPEND_DIR)/$*.d >> $(DEPEND_DIR)/$*.P; \
#	rm -f $(DEPEND_DIR)/$*.d
#	@echo "depend "$(DEPENDFILES)
#dependences sorted, build it
	@echo "Compiling $<..."
	@$(CXX) $(CXXFLAGS)  -c $< -o $@

$(OBJ_DIR)/%.o: %.cc
#this bit sorts out all the dependences
#	@echo "not packages $*"
	@$(CXX) $(CXXFLAGS) -MM -o $(DEPEND_DIR)/$*.d -c $<; \
	sed 's,$*.o,$(OBJ_DIR)/$*.o,' < $(DEPEND_DIR)/$*.d > $(DEPEND_DIR)/$*.P; \
	sed -e 's/#.*//' -e 's/^[^:]*: *//' -e 's/ *\\$$//' -e '/^$$/ d' -e 's/$$/ :/'  < $(DEPEND_DIR)/$*.d >> $(DEPEND_DIR)/$*.P; \
	rm -f $(DEPEND_DIR)/$*.d

#dependences sorted, build it
	@echo "Compiling $<..."
	$(CXX) $(CXXFLAGS)  -c $< -o $@


%Dict.cc:  # %_LinkDef.h
#	@echo "first dict $@" 
#so I needed to put this in dict/PACKAGENAME/src, seemed path of least resistance
	@test -d dict/$(PACKAGENAME)/src || mkdir -p dict/$(PACKAGENAME)/src
#now I allow ".h" or ".hh" hence the wild card
	@rootcint  -f dict/$@ -c -I$(ROOFITSYS)/include -Iinclude $(wildcard packages/$(PACKAGENAME)/include/$(subst $(PACKAGENAME)/src/,,$*).*h) packages/$(PACKAGENAME)/dict/$(subst $(PACKAGENAME)/src/,,$*)_LinkDef.h

#need to fix the include path in the generated header
	@sed 's|/include/|/|g' < dict/$@ >$ dict/$@.tmp
	@mv dict/$@.tmp dict/$@
	@sed 's|packages/||g' < dict/$@ >$ dict/$@.tmp
	@mv dict/$@.tmp dict/$@
	@echo "Finished dictionary Generation $@..."

#now need to move the .pcm files over to the libaray location
	@mv dict/$*Dict_rdict.pcm $(LIB_DIR)




%Exe:$(OBJ_DIR)/%.o $(STD_LIBS) 
	$(LD) $(LDFLAGS) $< $(LIBS) $(STD_LIBS) -o $(BIN_DIR)/$@	

-include $(DEPENDFILES)
# DO NOT DELETE
