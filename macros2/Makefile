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

CXXFLAGS      = -g -O3 -Wall -fPIC -pipe
LD            = $(CXX)
LDFLAGS       = -g
SOFLAGS       = -shared

CXXFLAGS     += $(ROOTCFLAGS)
LIBS          = $(ROOTLIBS)
GLIBS         = $(filter-out -lz, $(ROOTGLIBS))
RELLIBPATH    = ../../../../AnalysisDataFormats/HeavyFlavorObjects/lib

#########
# FILES #
#########
READER_HEADERS = bmm2Reader.hh candAna.hh candAnaMuMu.hh candAnaBu2JpsiKp.hh
READER = bmm2Reader.o candAna.o candAnaMuMu.o candAnaBu2JpsiKp.o

READERDICT = ${READER:.o=Dict.o}

#################
# PATTERN RULES #
#################

obj/%.o : %.cc $(READER_HEADERS)
	$(CXX) $(CXXFLAGS) -c $< -o $@

%Dict.cc : %.hh %LinkDef.h
	$(ROOTCINT) -f $@ -c $^

%Dict.cc : %.hh
	$(ROOTCINT) -f $@ -c $<


obj/%.o : %.c $(CSTUFF)
	$(CXX) $(CXXFLAGS) -c $< -o $@


###########
# TARGETS #
###########

all:
	@$(MAKE) bin/runBmm2

lib/libBmm2Classes.so: $(addprefix obj/,$(CSTUFF) $(READER) $(SUPPORTCLASSES) $(READERDICT)) lib/libAna00.so lib/libUtil.so lib/libAnaClasses.so
	$(CXX) $(SOFLAGS) $(addprefix obj/,$(CSTUFF) $(READER) $(READERDICT) $(SUPPORTCLASSES)) -o lib/libBmm2Classes.so $(GLIBS) -lMinuit lib/libAna00.so lib/libUtil.so lib/libAnaClasses.so

lib/libAna00.so:
	cd lib && rm -f libAna00.so && ln -s $(RELLIBPATH)/libAna00.so libAna00.so && cd -

lib/libUtil.so:
	cd lib && rm -f libUtil.so && ln -s $(RELLIBPATH)/libUtil.so libUtil.so && cd -

lib/libAnaClasses.so:
	cd lib && rm -f libAnaClasses.so && ln -s ../../macros/lib/libAnaClasses.so libAnaClasses.so && cd -

bin/runBmm2: lib/libBmm2Classes.so obj/runBmm2.o lib/libBmm2Classes.so
	$(LD) $(LDFLAGS) -o bin/runBmm2 $(GLIBS) obj/runBmm2.o lib/libBmm2Classes.so lib/libAna00.so lib/libUtil.so lib/libAnaClasses.so

clean:
	rm -f obj/*.o *Dict.* lib/lib*.so

cleanall: clean
	rm -f bin/runBmm2
	rm -f lib/lib*.so
	rm -f ../../../lib/$(SCRAM_ARCH)/libAna00.so
	rm -f ../../../lib/$(SCRAM_ARCH)/libUtil.so
	rm -f ../../../lib/$(SCRAM_ARCH)/libAnaClasses.so