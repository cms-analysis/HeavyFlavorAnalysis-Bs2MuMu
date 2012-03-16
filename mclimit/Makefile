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
#CXXFLAGS      = -g -O3 -Wall -fPIC -pipe -fpermissive
LD            = $(CXX)
LDFLAGS       = -g 
SOFLAGS       = -shared

CXXFLAGS     += $(ROOTCFLAGS)
LIBS          = $(ROOTLIBS)
GLIBS         = $(filter-out -lz, $(ROOTGLIBS))
RELLIBPATH    = ../../../../AnalysisDataFormats/HeavyFlavorObjects/lib
ANALIBPATH    = ../../macros/lib
#########
# FILES #
#########
SUPPORTCLASSES = mclimit_csm.o


#################
# PATTERN RULES #
#################
obj/%.o : %.cc
	$(CXX) $(CXXFLAGS) -c $< -o $@


###########
# TARGETS #
###########

all:	
lib/libAnaClasses: $(addprefix obj/,$(SUPPORTCLASSES)) 
	$(CXX) $(SOFLAGS) $(addprefix obj/,$(SUPPORTCLASSES)) -o lib/libAnaClasses.so $(GLIBS) -lMinuit

bin/multibin_cls: lib/libAnaClasses.so obj/multibin_cls.o
	$(LD) $(LDFLAGS) -o bin/multibin_cls $(GLIBS) obj/multibin_cls.o lib/libAnaClasses.so

bin/cms_lhcb_comb_cls: lib/libAnaClasses.so obj/cms_lhcb_comb_cls.o 
	$(LD) $(LDFLAGS) -o bin/cms_lhcb_comb_cls $(GLIBS) obj/cms_lhcb_comb_cls.o lib/libAnaClasses.so

bin/upperlimitCalc: lib/libAnaClasses.so obj/upperlimitCalc.o 
	$(LD) $(LDFLAGS) -o bin/upperlimitCalc $(GLIBS) obj/upperlimitCalc.o lib/libAnaClasses.so

clean:
	rm -f obj/*.o 

cleanall: clean
	rm -f lib/lib*.so
#	rm -f ../../../lib/$(SCRAM_ARCH)/libAnaClasses.so
