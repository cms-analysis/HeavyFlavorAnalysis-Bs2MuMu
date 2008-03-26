ROOTCFLAGS    = $(shell $(ROOTSYS)/bin/root-config --cflags)
ROOTLIBS      = $(shell $(ROOTSYS)/bin/root-config --libs)
ROOTGLIBS     = $(shell $(ROOTSYS)/bin/root-config --glibs)

CXX           = c++
CXXFLAGS      = -O2 -Wall -fPIC -Wno-deprecated
CXXFLAGS      = -g -Wall -fPIC -Wno-deprecated
LD            = c++
LDFLAGS       = -O2
LDFLAGS       = -g
SOFLAGS       = -shared

CXXFLAGS     += $(ROOTCFLAGS)
LIBS          = $(ROOTLIBS)
MLIBS         = $(ROOTGLIBS) 
MLIBS        += -lMinuit
GLIBS         = $(filter-out -lz, $(MLIBS))


$(addprefix obj/,%.o) : %.cc 
	$(CXX) $(CXXFLAGS) -c $< -o $@

ANA00   = TAna00Event.o TAna00EventDict.o \
          TGenCand.o TGenCandDict.o \
          TAnaTrack.o TAnaTrackDict.o \
          TAnaCand.o TAnaCandDict.o \
          TAnaVertex.o TAnaVertexDict.o

ANACLASSES = anaBmm.o anaBmmDict.o 

#             \
#             blimit.o incompletebeta.o incompletegamma.o posterior.o postint.o

# ================================================================================
bmm: bmm.cc treeBmm.cc
# --
	$(CXX) $(CXXFLAGS) -I../ -c treeBmm.cc -o treeBmm.o
	$(ROOTSYS)/bin/rootcint -f treeBmmDict.cc -c treeBmm.hh LinkDef.h
	$(CXX) $(CXXFLAGS) -o bmm.o -c bmm.cc
	$(LD) $(LDFLAGS)  -o bmm $(GLIBS) ../../../AnalysisDataFormats/HeavyFlavorObjects/lib/libAna00.so bmm.o treeBmm.o


ana: $(ANACLASSES)
	$(CXX) $(SOFLAGS) $(ANACLASSES) -o ../../../AnalysisDataFormats/HeavyFlavorObjects/lib/libAnaClasses.so

anaBmmDict.cc: anaBmm.hh 
	$(ROOTSYS)/bin/rootcint  -f anaBmmDict.cc -c anaBmm.hh


# ======================================================================
clean:
	rm -f *.o *Dict.* *~
	rm -f bmm

cleanall:
	@$(MAKE) clean
	rm -f *.o *Dict.* *~
	rm -f bmm
	rm -f ../../../AnalysisDataFormats/HeavyFlavorObjects/lib/libAnaClasses.so
