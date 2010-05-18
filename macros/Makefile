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

READER_HEADERS = bmmReader.hh d0Reader.hh dumpReader.hh muCharmReader.hh copyReader.hh dpReader.hh massReader.hh treeReader01.hh triggerValidation.hh ksReader.hh
READER = treeReader01.o bmmReader.o massReader.o dumpReader.o copyReader.o muCharmReader.o dpReader.o d0Reader.o triggerValidation.o ksReader.o
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

lib/libAnaClasses.so: $(addprefix obj/,$(READER) $(READERDICT)) lib/libAna00.so
	$(CXX) $(SOFLAGS) $(addprefix obj/,$(READER) $(READERDICT)) -o lib/libAnaClasses.so $(GLIBS) lib/libAna00.so

lib/libAna00.so:
	cd lib && ln -s $(RELLIBPATH)/libAna00.so libAna00.so && cd -

bin/runTreeReaders: lib/libAnaClasses.so obj/runTreeReaders.o lib/libAna00.so
	$(LD) $(LDFLAGS) -o bin/runTreeReaders $(GLIBS) obj/runTreeReaders.o lib/libAnaClasses.so lib/libAna00.so

clean:
	rm -f obj/*.o *Dict.* lib/libAna00.so

cleanall: clean
	rm -f bin/runTreeReaders
	rm -f bin/writeA01Event bin/readA01Event
	rm -f lib/lib*.so
	rm -f ../../../lib/$(SCRAM_ARCH)/libAna00.so
	rm -f ../../../lib/$(SCRAM_ARCH)/libUtil.so
	rm -f ../../../lib/$(SCRAM_ARCH)/libAnaClasses.so
