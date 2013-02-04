# Makefile for ana-producer

#########
# SETUP #
#########
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

XMLCFLAGS     = $(shell pkg-config --cflags libxml-2.0)
XMLLIBS       = $(shell pkg-config --libs libxml-2.0)

CXXFLAGS      = -g -O3 -Wall -fPIC -pipe $(ROOTCFLAGS) $(XMLCFLAGS)
LIBS          = $(ROOTLIBS) -lTreePlayer $(XMLLIBS)
GLIBS		  = $(ROOTGLIBS) -lTreePlayer $(XMLLIBS)
LDFLAGS       = -g
SOFLAGS       = -shared
LD            = $(CXX)

################
# OBJECT FILES #
################

OBJS = ana-produce.o
LIB_OBJS = ncMVA.o ncAna.o ncAna2.o ncCut.o ncFormula.o ncEvaluate.o ncVarReader.o ncConfig.o
LIB_OBJSDICT = ${LIB_OBJS:.o=Dict.o}

#################
# PATTERN RULES #
#################

obj/%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

%Dict.cpp: %.h
	$(ROOTCINT) -f $@ -c $<

###########
# TARGETS #
###########

all: bin obj lib bin/ana-produce lib/libNCAna.so bin/ncBatch bin/ncEvalTMVA
	@true

lib/libNCAna.so: $(addprefix obj/,$(LIB_OBJS) $(LIB_OBJSDICT)) lib/libNCRootUtils.so
	$(LD) $(SOFLAGS) -o $@ $(addprefix obj/,$(LIB_OBJS) $(LIB_OBJSDICT)) $(GLIBS) lib/libNCRootUtils.so -lTMVA -lRooFit -lRooFitCore

lib/libNCRootUtils.so:
	rm -f lib/libNCRootUtils.so
	cd lib && ln -s ../../rootutils/lib/libNCRootUtils.so && cd -

bin/ana-produce: $(addprefix obj/,$(OBJS))
	$(CXX) -o $@ $(addprefix obj/,$(OBJS)) $(LIBS)

bin/ncBatch: obj/ncBatch.o lib/libNCAna.so
	$(LD) $(LDFLAGS) -o $@ obj/ncBatch.o lib/libNCAna.so lib/libNCRootUtils.so $(GLIBS)

bin/ncEvalTMVA: obj/ncEvalTMVA.o
	$(LD) $(LDFLAGS) -o $@ obj/ncEvalTMVA.o $(GLIBS)

bin:
	@mkdir -p bin

obj:
	@mkdir -p obj

lib:
	@mkdir -p lib

clean:
	rm -rf obj *Dict.h

distclean: clean
	rm -rf bin lib
