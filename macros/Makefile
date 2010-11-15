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

READER_HEADERS = bmmReader.hh bmmSignalReader.hh bmmNormalizationReader.hh d0Reader.hh dumpReader.hh muCharmReader.hh copyReader.hh dpReader.hh massReader.hh treeReader01.hh triggerValidation.hh ksReader.hh lambdaReader.hh genLevel.hh kpReader.hh phiReader.hh decayCounter.hh impactReader.hh kpIsolationReader.hh phiIsolationReader.hh
READER = treeReader01.o bmmReader.o bmmSignalReader.o bmmNormalizationReader.o massReader.o dumpReader.o copyReader.o muCharmReader.o dpReader.o d0Reader.o triggerValidation.o ksReader.o lambdaReader.o genLevel.o kpReader.o phiReader.o decayCounter.o impactReader.o kpIsolationReader.o phiIsolationReader.o
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

###########
# TARGETS #
###########

all: 
	@$(MAKE) bin/runTreeReaders
	@$(MAKE) bin/runMergeTrees

lib/libAnaClasses.so: $(addprefix obj/,$(READER) $(READERDICT)) lib/libAna00.so lib/libUtil.so
	$(CXX) $(SOFLAGS) $(addprefix obj/,$(READER) $(READERDICT)) -o lib/libAnaClasses.so $(GLIBS) lib/libAna00.so lib/libUtil.so

lib/libAna00.so:
	cd lib && ln -s $(RELLIBPATH)/libAna00.so libAna00.so && cd -

lib/libUtil.so:
	cd lib && ln -s $(RELLIBPATH)/libUtil.so libUtil.so && cd -

bin/runTreeReaders: lib/libAnaClasses.so obj/runTreeReaders.o lib/libAna00.so
	$(LD) $(LDFLAGS) -o bin/runTreeReaders $(GLIBS) obj/runTreeReaders.o lib/libAnaClasses.so lib/libAna00.so lib/libUtil.so

bin/runMergeTrees: lib/libAnaClasses.so obj/runMergeTrees.o lib/libAna00.so
	$(LD) $(LDFLAGS) -o bin/runMergeTrees $(GLIBS) obj/runMergeTrees.o lib/libAnaClasses.so lib/libAna00.so lib/libUtil.so

clean:
	rm -f obj/*.o *Dict.* lib/libAna00.so lib/libUtil.so

cleanall: clean
	rm -f bin/runTreeReaders bin/runMergeTrees
	rm -f bin/writeA01Event bin/readA01Event
	rm -f lib/lib*.so
	rm -f ../../../lib/$(SCRAM_ARCH)/libAna00.so
	rm -f ../../../lib/$(SCRAM_ARCH)/libUtil.so
	rm -f ../../../lib/$(SCRAM_ARCH)/libAnaClasses.so
