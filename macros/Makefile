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
#CXXFLAGS      = -g -O0 -Wall -fPIC -pipe
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
READER = treeReader01.o massReader.o copyReader.o  genLevel.o impactReader.o AnalysisCuts.o AnalysisDistribution.o initFunc.o HistCutEfficiency.o pixelReader.o
READERDICT = ${READER:.o=Dict.o}

SUPPORTCLASSES = mclimit_csm.o
CSTUFF = bayesianlimit.o

#################
# PATTERN RULES #
#################

obj/%.o : %.cc
	$(CXX) $(CXXFLAGS) -c $< -o $@

%Dict.cc : %.hh %LinkDef.h
	$(ROOTCINT) -f $@ -c $^

%Dict.cc : %.hh
	$(ROOTCINT) -f $@ -c $<

obj/%.o : %.c $(CSTUFF)
	$(CXX) $(CXXFLAGS) -c $< -o $@

obj/%.d: %.cc
	@echo "Creating $@"
	$(shell $(CXX) -MM $(CXXFLAGS) $< | sed "s/.*\.o[ :]*/obj\/$(notdir $@) obj\/&/g" > $@)

###########
# TARGETS #
###########

all: bin/runTreeReaders bin/runMergeTrees
	@true

lib/libAnaClasses.so: $(addprefix obj/,$(CSTUFF) $(READER) $(SUPPORTCLASSES) $(READERDICT)) lib/libAna00.so lib/libUtil.so 
	$(CXX) $(SOFLAGS) $(addprefix obj/,$(CSTUFF) $(READER) $(READERDICT) $(SUPPORTCLASSES)) -o lib/libAnaClasses.so $(GLIBS) -lMinuit lib/libAna00.so lib/libUtil.so

lib/libAna00.so:
	cd lib && ln -s $(RELLIBPATH)/libAna00.so libAna00.so && cd -

lib/libUtil.so:
	cd lib && ln -s $(RELLIBPATH)/libUtil.so libUtil.so && cd -

bin/runTreeReaders: lib/libAnaClasses.so obj/runTreeReaders.o lib/libAna00.so
	$(LD) $(LDFLAGS) -o bin/runTreeReaders $(GLIBS) obj/runTreeReaders.o lib/libAnaClasses.so lib/libAna00.so lib/libUtil.so

bin/runMergeTrees: lib/libAnaClasses.so obj/runMergeTrees.o lib/libAna00.so
	$(LD) $(LDFLAGS) -o bin/runMergeTrees $(GLIBS) obj/runMergeTrees.o lib/libAnaClasses.so lib/libAna00.so lib/libUtil.so

bin/cls: lib/libAnaClasses.so obj/cls.o lib/libAna00.so
	$(LD) $(LDFLAGS) -o bin/cls $(GLIBS) obj/cls.o lib/libAnaClasses.so lib/libAna00.so lib/libUtil.so

clean:
	rm -f obj/*.o obj/*.d *Dict.* lib/libAna00.so lib/libUtil.so lib/libNCMVA.so lib/libAnaClasses.so

cleanall: clean
	rm -f bin/runTreeReaders bin/runMergeTrees 
	rm -f bin/writeA01Event bin/readA01Event
	rm -f bin/ncBatch
	rm -f lib/lib*.so
	rm -f ../../../lib/$(SCRAM_ARCH)/libAna00.so
	rm -f ../../../lib/$(SCRAM_ARCH)/libUtil.so
	rm -f ../../../lib/$(SCRAM_ARCH)/libAnaClasses.so

-include $(addprefix obj/,${READER:.o=.d})
