# ======================================================================
# NOTE: make all
# -----
#
# ======================================================================

ifdef ROOTSYS
  ROOTCINT      = $(ROOTSYS)/bin/rootcint
  ROOTCFLAGS    = $(shell $(ROOTSYS)/bin/root-config --cflags)
  ROOTLIBS      = $(shell $(ROOTSYS)/bin/root-config --libs)
  ROOTGLIBS     = $(shell $(ROOTSYS)/bin/root-config --glibs)
else 
  ROOTCINT      = rootcint
  ROOTCFLAGS    = $(shell root-config --cflags)
  ROOTLIBS      = $(shell root-config --libs)
  ROOTGLIBS     = $(shell root-config --glibs)
endif

ifdef SCRAM_ARCH
  CXX         := $(shell scramv1 tool info cxxcompiler | grep CXX= | sed s/CXX=//)
else 
  ifdef CXXCOMPILER
    CXX         := $(CXXCOMPILER)
  else 
    CXX         := c++
  endif
endif

CXXFLAGS      = -g -Wall -fPIC
LD            = $(CXX)
LDFLAGS       = -g
SOFLAGS       = -shared

CXXFLAGS     += $(ROOTCFLAGS)
LIBS          = $(ROOTLIBS)
GLIBS         = $(filter-out -lz, $(ROOTGLIBS))
RELLIBPATH    = ../../../AnalysisDataFormats/HeavyFlavorObjects/lib


# ======================================================================
# -- Default rules
#$(addprefix obj/,%.o) : $(addprefix rootio/,%.cc )
#	$(CXX) $(CXXFLAGS) -c $< -o $@

$(addprefix obj/,%.o) : %.cc
	$(CXX) $(CXXFLAGS) -c $< -o $@

#$(addprefix obj/,%Dict.o) : %Dict.cc
#	echo "bla"

#%Dict.cc : $subst(Dict.cc,.hh,$1)
#	echo "$1"

# ======================================================================


# ======================================================================
READERCLASSES = treeReader01.o treeReader01Dict.o \
                bmmReader.o bmmReaderDict.o 


# ================================================================================
all: 
# --
#	@$(MAKE) readerclasses
	@$(MAKE) runTreeReaders



# ================================================================================
readerclasses: $(addprefix obj/,$(READERCLASSES))
# -----------------------------------------------
	$(CXX) $(SOFLAGS) $(addprefix obj/,$(READERCLASSES)) -o lib/libAnaClasses.so

bmmReaderDict.cc: bmmReader.hh 
	$(ROOTSYS)/bin/rootcint -f bmmReaderDict.cc -c bmmReader.hh 
treeReader01Dict.cc: treeReader01.hh 
	$(ROOTSYS)/bin/rootcint -f treeReader01Dict.cc -c treeReader01.hh 

# ================================================================================
runTreeReaders: readerclasses obj/runTreeReaders.o
	$(LD) $(LDFLAGS) -o bin/runTreeReaders $(GLIBS) obj/runTreeReaders.o obj/bmmReader.o obj/treeReader01.o $(RELLIBPATH)/libAna00.so


# ================================================================================
clean:
	rm -f obj/*.o *Dict.* *Dict.* lib/libAna00.so

# --------------------------------------------------------------------------------
cleanall:
	@$(MAKE) clean
	rm -f bin/runTreeReaders
	rm -f bin/writeA01Event bin/readA01Event
	rm -f lib/lib*.so
	rm -f ../../../lib/$(SCRAM_ARCH)/libAna00.so
	rm -f ../../../lib/$(SCRAM_ARCH)/libUtil.so
	rm -f ../../../lib/$(SCRAM_ARCH)/libAnaClasses.so

